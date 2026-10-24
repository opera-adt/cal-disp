#!/usr/bin/env python
"""Check the operational install against the dependency budget (plan T26.4).

The cal-disp image must not carry heavy research packages (Venti PRD R-O6).
`budget.yaml` lists the forbidden distributions and the exceptions still
tolerated, each with a reason. The installed set is read either from the
running interpreter (default) or from a pip resolution report::

    # this environment
    python scripts/check_env_budget.py
    # what `pip install cal-disp` (no extras) pulls in, from installed metadata
    python scripts/check_env_budget.py --closure cal-disp
    # a pip resolution
    pip install --dry-run --ignore-installed --report r.json .
    python scripts/check_env_budget.py --report r.json

    python scripts/check_env_budget.py --peak-memory-log test_golden/validation.log

The last form checks the ``Maximum memory usage: X GB`` line of a golden run
against ``peak_memory_gb`` (plan T26.5).

Exit status 1 when a forbidden package is present without an exception, an
exception names a package that is no longer present (stale budget), or the
peak memory exceeds the budget.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import re
import sys
from pathlib import Path

from ruamel.yaml import YAML

ROOT = Path(__file__).resolve().parents[1]


def _norm(name: str) -> str:
    return name.lower().replace("-", "_").replace(".", "_")


def installed_from_environment() -> set[str]:
    """Return the distribution names installed in the running interpreter."""
    return {_norm(d.metadata["Name"]) for d in importlib.metadata.distributions()}


def installed_from_report(path: Path) -> set[str]:
    """Return the distribution names a ``pip install --report`` would install."""
    report = json.loads(Path(path).read_text())
    return {_norm(item["metadata"]["name"]) for item in report.get("install", [])}


def installed_closure(root: str) -> tuple[set[str], set[str]]:
    """Return the requirement closure of distribution `root` without extras.

    Walks ``Requires-Dist`` of the installed distributions, keeping only the
    requirements whose markers hold with no extra requested, so the result is
    what ``pip install <root>`` brings in. Returns ``(closure, missing)``.
    """
    from packaging.requirements import Requirement

    seen: set[str] = set()
    missing: set[str] = set()
    todo = [root]
    while todo:
        name = _norm(todo.pop())
        if name in seen:
            continue
        try:
            dist = importlib.metadata.distribution(name)
        except importlib.metadata.PackageNotFoundError:
            missing.add(name)
            continue
        seen.add(name)
        for spec in dist.requires or []:
            req = Requirement(spec)
            if req.marker is not None and not req.marker.evaluate({"extra": ""}):
                continue
            todo.append(req.name)
            # requested extras of a dependency (e.g. geepers[grid]) add their deps
            for extra in req.extras:
                for spec2 in dist_requires_for_extra(req.name, extra):
                    todo.append(spec2)
    return seen, missing


def dist_requires_for_extra(name: str, extra: str) -> list[str]:
    """Return the requirement names a dependency's `extra` adds."""
    from packaging.requirements import Requirement

    try:
        dist = importlib.metadata.distribution(name)
    except importlib.metadata.PackageNotFoundError:
        return []
    out = []
    for spec in dist.requires or []:
        req = Requirement(spec)
        if (
            req.marker is not None
            and req.marker.evaluate({"extra": extra})
            and not req.marker.evaluate({"extra": ""})
        ):
            out.append(req.name)
    return out


def peak_memory_gb(log_text: str) -> float:
    """Return the last ``Maximum memory usage: X GB`` value of a run log."""
    values = re.findall(r"Maximum memory usage:\s*([0-9.]+)\s*GB", log_text)
    if not values:
        msg = "no 'Maximum memory usage: X GB' line in the log"
        raise ValueError(msg)
    return float(values[-1])


def check(installed: set[str], budget: dict) -> tuple[list[str], list[str], list[str]]:
    """Return ``(violations, tolerated, stale)`` package names."""
    forbidden = [_norm(n) for n in budget.get("forbidden", [])]
    exceptions = {
        _norm(k): v for k, v in (budget.get("known_exceptions") or {}).items()
    }
    present = [n for n in forbidden if n in installed]
    violations = [n for n in present if n not in exceptions]
    tolerated = [n for n in present if n in exceptions]
    stale = [n for n in exceptions if n not in installed]
    return violations, tolerated, stale


def main(argv: list[str] | None = None) -> int:
    """Print the budget report; return the exit status."""
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--budget", type=Path, default=ROOT / "budget.yaml")
    src = p.add_mutually_exclusive_group()
    src.add_argument("--report", type=Path, help="pip install --dry-run --report JSON")
    src.add_argument(
        "--closure", metavar="DIST", help="requirement closure of DIST without extras"
    )
    src.add_argument(
        "--peak-memory-log", type=Path, help="log of a golden run (cal-disp run)"
    )
    args = p.parse_args(argv)
    budget = YAML(typ="safe").load(args.budget.read_text())
    if args.peak_memory_log:
        limit = float(budget["peak_memory_gb"])
        peak = peak_memory_gb(args.peak_memory_log.read_text())
        ok = peak <= limit
        print(
            f"peak memory {peak:.2f} GB, budget {limit:.2f} GB:"
            f" {'OK' if ok else 'FAILED'}"
        )
        return 0 if ok else 1
    if args.report:
        installed = installed_from_report(args.report)
    elif args.closure:
        installed, missing = installed_closure(args.closure)
        for name in sorted(missing):
            print(f"  (not installed here, not checked: {name})")
    else:
        installed = installed_from_environment()
    violations, tolerated, stale = check(installed, budget)
    exceptions = {
        _norm(k): v for k, v in (budget.get("known_exceptions") or {}).items()
    }
    source = (
        f"report {args.report}"
        if args.report
        else f"closure of {args.closure}" if args.closure else sys.executable
    )
    print(f"dependency budget: {len(installed)} distributions ({source})")
    for name in tolerated:
        print(f"  tolerated  {name}: {exceptions[name]}")
    for name in violations:
        print(f"  FORBIDDEN  {name}: not allowed in the operational install")
    for name in stale:
        print(
            f"  STALE      {name}: listed as an exception but not installed; remove it"
        )
    ok = not violations and not stale
    print("OK" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
