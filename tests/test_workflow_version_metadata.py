"""The product must record the package's own version, not the dist-info's."""

import ast
from pathlib import Path

WORKFLOW = Path(__file__).resolve().parents[1] / "src" / "cal_disp" / "workflow.py"


def test_run_calibration_uses_package_version():
    tree = ast.parse(WORKFLOW.read_text())
    calls = [
        n
        for n in ast.walk(tree)
        if isinstance(n, ast.Call)
        and isinstance(n.func, ast.Name)
        and n.func.id == "_pkg_version"
    ]
    names = {ast.literal_eval(c.args[0]) for c in calls if c.args}
    # venti is a plain dependency; cal-disp itself must come from __version__
    assert names == {"venti"}
