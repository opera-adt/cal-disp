"""Tests for UNR grid staging: URL/local name per grid type, download robustness.

The robustness tests run a local ``http.server`` in a thread that serves a
good file, a truncated body, HTML "error" pages and a slow endpoint; no
network access is needed.
"""

from __future__ import annotations

import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from unittest.mock import MagicMock

import pytest
import requests

from cal_disp.download import _stage_unr
from cal_disp.download._stage_unr import (
    DownloadError,
    create_session,
    download_file,
    download_grid_file,
    download_lookup_table,
)

GOOD_BODY = b"2022.0 0 0 0 1 1 1 0\n" * 50
HTML_BODY = (
    b"<!DOCTYPE html>\n<html><body><h1>503 Service Unavailable</h1></body></html>"
)


class _Handler(BaseHTTPRequestHandler):
    """Serve the fixed responses the robustness tests need."""

    protocol_version = "HTTP/1.0"

    def _send(self, body: bytes, content_type: str, length: int | None = None):
        self.send_response(200)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body) if length is None else length))
        self.end_headers()
        self.wfile.write(body)
        self.wfile.flush()

    def do_GET(self):  # noqa: N802
        route = self.path.split("/")[-1]
        if route == "good.tenv8" or route.endswith("_IGS20.tenv8"):
            self._send(GOOD_BODY, "text/plain")
        elif route == "truncated.tenv8":
            # Announce the full length but send only the first half, then close
            self._send(GOOD_BODY[: len(GOOD_BODY) // 2], "text/plain", len(GOOD_BODY))
        elif route == "html-typed.tenv8":
            self._send(HTML_BODY, "text/html; charset=utf-8")
        elif route == "html-untyped.tenv8":
            # A misconfigured server: HTML body with a binary content type
            self._send(HTML_BODY, "application/octet-stream")
        elif route == "slow.tenv8":
            time.sleep(3)
            self._send(GOOD_BODY, "text/plain")
        else:
            self.send_error(404)

    def log_message(self, *args):  # noqa: D102
        pass


class _Server(ThreadingHTTPServer):
    # Non-daemon handler threads so server_close() joins them (the slow one)
    daemon_threads = False

    def handle_error(self, request, client_address):  # noqa: D102
        pass  # the client hanging up on the slow endpoint is expected


@pytest.fixture(scope="module")
def http_url():
    server = _Server(("127.0.0.1", 0), _Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{server.server_address[1]}"
    server.shutdown()
    server.server_close()
    thread.join()


@pytest.fixture
def session():
    # No retries so failures surface immediately and the tests stay fast
    with create_session(retries=0) as s:
        yield s


def test_download_file(tmp_path: Path, http_url, session):
    out = tmp_path / "good.tenv8"
    assert download_file(f"{http_url}/good.tenv8", out, session=session) == out
    assert out.read_bytes() == GOOD_BODY
    assert not (tmp_path / "good.tenv8.part").exists()


def test_truncated_body_leaves_no_file(tmp_path: Path, http_url, session):
    out = tmp_path / "truncated.tenv8"
    with pytest.raises((DownloadError, requests.RequestException)):
        download_file(f"{http_url}/truncated.tenv8", out, session=session)
    assert not out.exists()
    assert not (tmp_path / "truncated.tenv8.part").exists()


@pytest.mark.parametrize("route", ["html-typed.tenv8", "html-untyped.tenv8"])
def test_html_error_page_rejected(tmp_path: Path, http_url, session, route):
    out = tmp_path / route
    with pytest.raises(DownloadError, match="HTML page"):
        download_file(f"{http_url}/{route}", out, session=session)
    assert not out.exists()
    assert not (tmp_path / f"{route}.part").exists()


def test_read_timeout_raises(tmp_path: Path, http_url, session):
    """A stalled server raises after the read timeout instead of hanging.

    Through the retry adapter urllib3 reports the exhausted read timeout as
    ``MaxRetryError``, which requests maps to ``ConnectionError``; without
    the adapter it is a ``ReadTimeout``. Either way the message names it.
    """
    out = tmp_path / "slow.tenv8"
    start = time.monotonic()
    with pytest.raises((requests.Timeout, requests.ConnectionError), match="timed out"):
        download_file(f"{http_url}/slow.tenv8", out, session=session, timeout=(5, 0.5))
    assert time.monotonic() - start < 2.5  # the endpoint sleeps for 3 s
    assert not out.exists()
    assert not (tmp_path / "slow.tenv8.part").exists()


def test_timeout_passed_to_every_request(tmp_path: Path):
    session = MagicMock()
    response = session.get.return_value.__enter__.return_value
    response.headers = {"Content-Type": "text/plain", "Content-Length": "4"}
    response.iter_content.return_value = [b"data"]
    download_file("https://example.invalid/x", tmp_path / "x", session=session)
    assert session.get.call_args.kwargs["timeout"] == _stage_unr.DEFAULT_TIMEOUT
    assert _stage_unr.DEFAULT_TIMEOUT == (10.0, 120.0)
    assert session.get.call_args.kwargs["stream"] is True


def test_retry_statuses():
    session = create_session(retries=3, backoff=2.0)
    retry = session.get_adapter("https://").max_retries
    assert {429, 500, 502, 503, 504} <= set(retry.status_forcelist)
    assert retry.total == 3
    assert retry.backoff_factor == 2.0
    # Plain http (e.g. a local mirror) gets the same retry policy
    assert session.get_adapter("http://").max_retries is retry


def test_leftover_part_is_redownloaded(tmp_path: Path, http_url, session, monkeypatch):
    monkeypatch.setitem(
        _stage_unr.GRID_BASE_URLS, "constant", f"{http_url}/v{{version}}"
    )
    monkeypatch.setattr(_stage_unr, "grid_file_name", lambda *_: "good.tenv8")
    part = tmp_path / "good.tenv8.part"
    part.write_bytes(b"half")

    path = download_grid_file(1, tmp_path, session=session)

    assert path == tmp_path / "good.tenv8"
    assert path.read_bytes() == GOOD_BODY
    assert not part.exists()


def test_complete_file_is_reused(tmp_path: Path):
    session = MagicMock()
    out = tmp_path / "grid_latlon_lookup_v0.3.txt"
    out.write_bytes(b"1 2 3\n")
    assert download_lookup_table(tmp_path, session=session) == out
    session.get.assert_not_called()


def test_empty_file_is_not_trusted(tmp_path: Path, http_url, session, monkeypatch):
    monkeypatch.setattr(_stage_unr, "LOOKUP_URL", f"{http_url}/good.tenv8")
    out = tmp_path / "grid_latlon_lookup_v0.3.txt"
    out.touch()
    assert download_lookup_table(tmp_path, session=session) == out
    assert out.read_bytes() == GOOD_BODY


class _FakeResponse:
    """Minimal streaming response for URL checks with a mock session."""

    def __init__(self, body: bytes, headers: dict[str, str] | None = None):
        self._body = body
        self.headers = {"Content-Type": "text/plain", **(headers or {})}

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def raise_for_status(self):
        pass

    def iter_content(self, chunk_size):  # noqa: ARG002
        yield self._body


@pytest.mark.parametrize(
    ("headers", "rejected"),
    [
        ({"Content-Length": "100"}, True),
        ({"Content-Length": "21"}, False),
        ({}, False),
        # Encoded bodies: the header counts encoded bytes, so it is not compared
        ({"Content-Length": "100", "Content-Encoding": "gzip"}, False),
    ],
)
def test_content_length_checked(tmp_path: Path, headers, rejected):
    """Our own length check, for clients that do not enforce Content-Length.

    urllib3 >= 2 raises on a short body by itself (the local-server test
    above), urllib3 1.x returned it silently.
    """
    session = MagicMock()
    session.get.return_value = _FakeResponse(b"2022.0 0 0 0 1 1 1 0\n", headers)
    out = tmp_path / "x.tenv8"

    if rejected:
        with pytest.raises(DownloadError, match="Content-Length is 100"):
            download_file("https://example.invalid/x", out, session=session)
        assert not out.exists()
    else:
        assert download_file("https://example.invalid/x", out, session=session) == out
        assert out.read_bytes() == b"2022.0 0 0 0 1 1 1 0\n"
    assert not (tmp_path / "x.tenv8.part").exists()


@pytest.mark.parametrize(
    ("grid_type", "remote_dir"),
    [("constant", "time_contsant_gridded"), ("variable", "time_variable_gridded")],
)
def test_download_grid_file(tmp_path: Path, grid_type, remote_dir):
    session = MagicMock()
    session.get.return_value = _FakeResponse(b"2022.0 0 0 0 1 1 1 0\n")

    path = download_grid_file(12, tmp_path, grid_type=grid_type, session=session)

    url = session.get.call_args.args[0]
    assert (
        url
        == f"https://geodesy.unr.edu/grid_timeseries/Version0.3/{remote_dir}"
        "/IGS20/000012_IGS20.tenv8"
    )
    # Same name Venti's download_station looks for, so staged files are reused
    assert path == tmp_path / f"000012_IGS20_{grid_type}.tenv8"
    assert path.read_bytes() == b"2022.0 0 0 0 1 1 1 0\n"


def test_constant_grid_unavailable(tmp_path: Path):
    with pytest.raises(ValueError, match="constant grid only"):
        download_grid_file(12, tmp_path, version="0.2", grid_type="constant")
