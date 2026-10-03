"""The core workflow must import without the optional ``download`` extra."""

import subprocess
import sys

CHECK = """
import sys
# Simulate the extra not being installed
for name in ("asf_search",):
    sys.modules[name] = None
import cal_disp.workflow
import cal_disp.main
import cal_disp.cli
"""


def test_core_imports_without_download_extra():
    result = subprocess.run(
        [sys.executable, "-c", CHECK], capture_output=True, text=True, check=False
    )
    assert result.returncode == 0, result.stderr
