"""Package the manuscript together with the code and data that produced it.

A reviewer who has the PDF and nothing else cannot check anything. This
assembles one self-contained folder (and optionally a zip) holding the
manuscript, its figures, the numbers behind every table, and the source that
generated them:

    python paper/make_submission.py            # build paper/submission/
    python paper/make_submission.py --zip      # and zip it
    python paper/make_submission.py --clean    # remove it again

ponytail: the bundle is *built*, not committed. Keeping a second copy of the
source inside paper/ would go stale the first time anything changed, and then
the paper would ship code that never produced its own figures. One copy in the
repository stays the truth; this copies it on demand.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from datetime import date
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = HERE / "submission"

# Only the results the manuscript actually cites travel with it.
DATA = [
    "exp06_shootout.csv",
    "exp07_formulations.csv",
    "exp08_optimizer_benchmark.csv",
    "exp08_certified_instance.csv",
    "exp08_summary.json",
    "exp09_science_cases.csv",
    "exp09_multifrequency.csv",
    "exp09_soft_vs_hard.csv",
    "exp09_summary.json",
    "exp10_multiepoch.csv",
    "exp10_window_sweep.csv",
    "exp10_summary.json",
]

FIGURES = [
    "fig16_formulation_comparison.png",
    "fig17_optimizer_comparison.png",
    "fig18_optimized_layout.png",
    "fig19_science_case_layouts.png",
    "fig21_multifrequency_gain.png",
    "fig22_multiepoch_tradeoff.png",
]

EXPERIMENTS = ["common.py", "06_layout_shootout.py", "07_formulation_comparison.py",
               "08_optimizer_benchmark.py", "09_science_cases.py", "10_multiepoch.py"]

TOP_LEVEL = ["requirements.txt", "pyproject.toml", "LICENSE"]


def copy(src: Path, dst: Path) -> bool:
    if not src.exists():
        print(f"  ! missing, skipped: {src.relative_to(ROOT)}")
        return False
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return True


def build(make_zip: bool) -> Path:
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)

    n = 0

    # the paper itself
    for name in ("manuscript.md", "README.md", "build_pdf.py"):
        n += copy(HERE / name, OUT / name)
    pdf = HERE / "manuscript.pdf"
    if not pdf.exists():
        print("  building manuscript.pdf first...")
        subprocess.run([sys.executable, str(HERE / "build_pdf.py")], check=True)
    n += copy(pdf, OUT / "manuscript.pdf")

    # figures, at the relative path the manuscript expects (../figures/...)
    for f in FIGURES:
        n += copy(ROOT / "figures" / f, OUT / "figures" / f)

    # the numbers behind every table
    for d in DATA:
        n += copy(ROOT / "data" / "results" / d, OUT / "data" / "results" / d)

    # the code: library, the experiments cited, the tests that check the maths
    shutil.copytree(ROOT / "src", OUT / "src",
                    ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    n += sum(1 for _ in (OUT / "src").rglob("*.py"))
    for e in EXPERIMENTS:
        n += copy(ROOT / "experiments" / e, OUT / "experiments" / e)
    shutil.copytree(ROOT / "tests", OUT / "tests",
                    ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    n += sum(1 for _ in (OUT / "tests").rglob("*.py"))
    shutil.copytree(ROOT / "configs", OUT / "configs")
    for t in TOP_LEVEL:
        n += copy(ROOT / t, OUT / t)

    (OUT / "HOW_TO_VERIFY.txt").write_text(_instructions(), encoding="utf-8")

    print(f"\n  built {OUT.relative_to(ROOT)} with {n} files")
    if make_zip:
        archive = shutil.make_archive(str(HERE / f"submission_{date.today():%Y%m%d}"),
                                      "zip", root_dir=OUT)
        size = Path(archive).stat().st_size / 1e6
        print(f"  zipped -> {Path(archive).name} ({size:.1f} MB)")
    return OUT


def _instructions() -> str:
    return f"""How to verify this manuscript
=============================
Built {date.today():%Y-%m-%d} from the project repository.

This bundle contains the paper, every figure in it, the numbers behind every
table, and the source code that produced all of them.

1. Install dependencies

       pip install -r requirements.txt

2. Check the mathematics (about one minute)

       python -m pytest -q

   144 checks. They verify, among other things, that the matrix QUBO equals a
   direct term-by-term implementation over all 2^8 assignments, that Parseval's
   identity holds for the gridded sampling function, that the Bonferroni
   expression is a genuine lower bound on unique-cell coverage, and that the
   Earth-rotation mapping reduces to the zenith snapshot.

3. Regenerate the figures and tables

       cd experiments
       python 07_formulation_comparison.py   # Figure 1,  Table 4.1
       python 08_optimizer_benchmark.py      # Figures 2-3, Tables 4.2 / 4.2.1
       python 09_science_cases.py            # Figures 4-5, Tables 4.3 / 4.4
       python 10_multiepoch.py               # Figure 6,  Tables 4.5
       python 06_layout_shootout.py          # Table 4.6

   Outputs overwrite figures/ and data/results/. Experiment 08 takes about ten
   minutes because it runs two mixed-integer solves.

4. Rebuild the PDF

       python paper/build_pdf.py

README.md maps every figure, table and equation in the manuscript to the exact
file that implements or produces it.

Note on scope: all candidate pad coordinates and atmospheric constants here are
synthetic controlled inputs, not a site survey. The manuscript states this.
"""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--zip", action="store_true", help="also produce a zip archive")
    ap.add_argument("--clean", action="store_true", help="delete the bundle and exit")
    a = ap.parse_args()

    if a.clean:
        if OUT.exists():
            shutil.rmtree(OUT)
            print(f"removed {OUT}")
        else:
            print("nothing to remove")
        return

    build(a.zip)


if __name__ == "__main__":
    main()
