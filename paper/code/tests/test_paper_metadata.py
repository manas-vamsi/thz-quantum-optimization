"""The submission-metadata gate.

These tests exist so the manuscript cannot drift into claiming to be
submission-ready. Placeholder values are easy to leave behind and impossible to
spot in a 13-page PDF, and an invented affiliation or an uncredited coauthor on
a submitted paper is misconduct rather than an oversight.
"""

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "paper"))

import build_pdf  # noqa: E402

META = ROOT / "paper" / "metadata.yaml"
ZENODO = ROOT / ".zenodo.json"


def test_metadata_file_exists_and_parses():
    assert META.exists(), "paper/metadata.yaml is missing"
    meta = build_pdf.load_metadata(META)
    assert isinstance(meta, dict) and meta


def test_every_required_field_is_present_as_a_key():
    """Filling in a value must not require knowing which keys exist."""
    meta = build_pdf.load_metadata(META)
    for key in ("title", "authors", "funding", "acknowledgements",
                "repository_doi", "competing_interests", "target_venue"):
        assert key in meta, f"metadata.yaml has no {key!r} key"


def test_the_checker_detects_placeholders():
    """A checker that passes on TODO values would be worse than none."""
    stub = {"authors": [{"name": "A", "affiliation": "TODO", "email": "TODO"}],
            "funding": "TODO", "acknowledgements": "None",
            "repository_doi": "TODO", "competing_interests": "None",
            "target_venue": "TODO"}
    missing = build_pdf.missing_metadata(stub)
    assert "funding" in missing
    assert "repository_doi" in missing
    assert "A: affiliation" in missing
    assert "A: email" in missing
    assert "acknowledgements" not in missing      # legitimately answered "None"


def test_the_checker_passes_a_complete_record():
    full = {"authors": [{"name": "A", "affiliation": "Somewhere",
                         "email": "a@example.org", "corresponding": True}],
            "funding": "None", "acknowledgements": "None",
            "repository_doi": "10.5281/zenodo.0000000",
            "competing_interests": "None", "target_venue": "arXiv"}
    assert build_pdf.missing_metadata(full) == []


def test_an_empty_author_list_is_caught():
    assert "authors" in build_pdf.missing_metadata(
        {"authors": [], "funding": "None", "acknowledgements": "None",
         "repository_doi": "x", "competing_interests": "None",
         "target_venue": "x"})


def test_author_block_numbers_shared_affiliations_once():
    meta = {"authors": [
        {"name": "A", "affiliation": "Inst X", "email": "a@x", "corresponding": True},
        {"name": "B", "affiliation": "Inst Y", "email": "b@y"},
        {"name": "C", "affiliation": "Inst X", "email": "c@x"},
    ]}
    lines = build_pdf.author_block(meta)
    assert "1. Inst X" in lines and "2. Inst Y" in lines
    assert sum(l.startswith("1. ") for l in lines) == 1
    assert any("Correspondence: a@x" in l for l in lines)


def test_author_block_is_empty_without_authors():
    assert build_pdf.author_block({}) == []


def test_zenodo_record_exists_and_is_valid_json():
    assert ZENODO.exists(), ".zenodo.json is missing"
    data = json.loads(ZENODO.read_text(encoding="utf-8"))
    assert data["title"]
    assert data["upload_type"] == "software"
    assert data["creators"]


def test_zenodo_record_still_flags_its_own_placeholders():
    """It is prepared, not finished; the notes field must say so."""
    data = json.loads(ZENODO.read_text(encoding="utf-8"))
    raw = ZENODO.read_text(encoding="utf-8")
    if "TODO" in raw:
        assert "TODO" in data.get("notes", ""), \
            "placeholders remain but the notes do not warn about them"


def test_the_manuscript_banner_matches_the_metadata_state():
    """The paper must not say it is ready while fields are TODO, or vice versa."""
    md = (ROOT / "paper" / "manuscript.md").read_text(encoding="utf-8")
    banner = "Submission metadata still required" in md
    incomplete = bool(build_pdf.missing_metadata(build_pdf.load_metadata(META)))
    assert banner == incomplete, (
        "the manuscript's front-matter banner disagrees with metadata.yaml: "
        f"banner present={banner}, fields outstanding={incomplete}")
