"""Refresh the code, figures and data held alongside the manuscript.

The paper directory is a complete, self-contained artifact: a reader gets the
manuscript, the source that implements every equation in it, the scripts that
produced every figure, the numbers behind every table, and the test suite that
checks the mathematics -- without needing the rest of the repository.

Holding a second copy of the source risks it drifting from the code that
actually produced the results, so every synchronised file is checksummed into
``code/MANIFEST.md`` together with the commit it came from. ``--check`` compares
the two and reports any drift instead of hiding it.

    python paper/sync_code.py            # refresh, and rewrite the manifest
    python paper/sync_code.py --check    # verify against the repository, change nothing
    python paper/sync_code.py --archive  # refresh, then produce a single zip
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
import subprocess
from datetime import date
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
CODE = HERE / "code"

ARCHIVE_STEM = "thz-qubo-array-design-supplementary"

# Experiments cited by the manuscript. The rest of the repository's experiments
# are exploratory and are deliberately not shipped with the paper.
EXPERIMENTS = [
    "common.py",
    "06_layout_shootout.py",
    "07_formulation_comparison.py",
    "08_optimizer_benchmark.py",
    "09_science_cases.py",
    "10_multiepoch.py",
]

FIGURES = [
    "fig16_formulation_comparison.png",
    "fig17_optimizer_comparison.png",
    "fig18_optimized_layout.png",
    "fig19_science_case_layouts.png",
    "fig21_multifrequency_gain.png",
    "fig22_multiepoch_tradeoff.png",
]

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

PYTEST_INI = """[pytest]
testpaths = tests
pythonpath = src
addopts = -q
"""


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()[:16]


def commit() -> str:
    try:
        out = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
                             capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except Exception:
        return "unknown"


def planned() -> list:
    """Every (source, destination) pair the paper directory should hold."""
    pairs = []
    for p in sorted((ROOT / "src").rglob("*.py")):
        pairs.append((p, CODE / "src" / p.relative_to(ROOT / "src")))
    for p in sorted((ROOT / "tests").glob("*.py")):
        pairs.append((p, CODE / "tests" / p.name))
    for name in EXPERIMENTS:
        pairs.append((ROOT / "experiments" / name, CODE / "experiments" / name))
    for p in sorted((ROOT / "configs").glob("*.yaml")):
        pairs.append((p, CODE / "configs" / p.name))
    pairs.append((ROOT / "requirements.txt", CODE / "requirements.txt"))
    for name in FIGURES:
        pairs.append((ROOT / "figures" / name, HERE / "figures" / name))
    for name in DATA:
        pairs.append((ROOT / "data" / "results" / name, HERE / "data" / name))
    return pairs


def check() -> int:
    """Report files that differ from the repository. Returns the drift count."""
    drift, missing = [], []
    for src, dst in planned():
        if not src.exists():
            continue
        if not dst.exists():
            missing.append(dst)
        elif digest(src) != digest(dst):
            drift.append(dst)

    if not drift and not missing:
        print(f"paper/ is in sync with the repository at commit {commit()}")
        return 0
    for d in missing:
        print(f"  missing : {d.relative_to(HERE)}")
    for d in drift:
        print(f"  outdated: {d.relative_to(HERE)}")
    print(f"\n{len(missing)} missing, {len(drift)} outdated. "
          f"Run 'python paper/sync_code.py' to refresh.")
    return len(drift) + len(missing)


def sync() -> list:
    for folder in (CODE, HERE / "figures", HERE / "data"):
        if folder.exists():
            shutil.rmtree(folder)

    copied = []
    for src, dst in planned():
        if not src.exists():
            print(f"  ! missing in repository, skipped: {src.relative_to(ROOT)}")
            continue
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        copied.append(dst)

    (CODE / "pytest.ini").write_text(PYTEST_INI, encoding="utf-8")
    copied.append(CODE / "pytest.ini")
    write_manifest(copied)
    return copied


def write_manifest(copied: list) -> None:
    sha = commit()
    lines = [
        "# Code manifest",
        "",
        f"Synchronised from the project repository at commit `{sha}` "
        f"on {date.today():%Y-%m-%d}.",
        "",
        "Each entry is the SHA-256 prefix of the file as shipped. Regenerate and "
        "re-verify with:",
        "",
        "```bash",
        "python paper/sync_code.py --check",
        "```",
        "",
        "| File | SHA-256 (first 16) | Bytes |",
        "|---|---|---:|",
    ]
    for p in sorted(copied):
        if p.name == "MANIFEST.md":
            continue
        rel = p.relative_to(HERE).as_posix()
        lines.append(f"| `{rel}` | `{digest(p)}` | {p.stat().st_size} |")
    lines += ["", f"{len(copied)} files."]
    (CODE / "MANIFEST.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def archive() -> Path:
    staging = HERE.parent / f"_{ARCHIVE_STEM}"
    if staging.exists():
        shutil.rmtree(staging)
    shutil.copytree(HERE, staging,
                    ignore=shutil.ignore_patterns("__pycache__", "*.pyc",
                                                  ".pytest_cache", "_*"))
    out = shutil.make_archive(str(HERE / ARCHIVE_STEM), "zip", root_dir=staging)
    shutil.rmtree(staging)
    return Path(out)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true",
                    help="verify against the repository without changing anything")
    ap.add_argument("--archive", action="store_true",
                    help="also produce a single zip of the paper directory")
    a = ap.parse_args()

    if a.check:
        raise SystemExit(1 if check() else 0)

    copied = sync()
    print(f"synchronised {len(copied)} files from commit {commit()}")
    if a.archive:
        out = archive()
        print(f"archive: {out.name} ({out.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
