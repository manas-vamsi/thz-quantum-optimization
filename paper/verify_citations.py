"""Resolve every DOI and arXiv identifier in the manuscript against its registry.

A reference list is a set of factual claims, and an unverified one is the
easiest place for an error to hide. This resolves each DOI through CrossRef and
each arXiv identifier through arXiv, prints the author and title the registry
actually returns, and exits non-zero if anything fails to resolve -- so the
claim "these references are real" is checkable rather than asserted.

    python paper/verify_citations.py

It caught a genuine misattribution during drafting: arXiv:2206.11040 was
credited to the wrong authors, which the registry record disproved.

Network access is required. arXiv's API is blocked from some networks and
returns HTTP 406; those entries are reported as unchecked rather than as
failures, since an unreachable host is not evidence of a bad citation.
"""

from __future__ import annotations

import io
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
MANUSCRIPT = HERE / "manuscript.md"

UA = {"User-Agent": "thz-qubo-citation-check/1.0 (mailto:noreply@example.com)",
      "Accept": "*/*"}


def fetch(url: str, tries: int = 3) -> bytes:
    last: Exception | None = None
    for k in range(tries):
        try:
            return urllib.request.urlopen(
                urllib.request.Request(url, headers=UA), timeout=30).read()
        except urllib.error.HTTPError as e:
            if e.code == 429 and k < tries - 1:      # rate limited: back off
                time.sleep(4 * (k + 1))
                last = e
                continue
            raise
        except Exception as e:                        # transient network trouble
            last = e
            time.sleep(2 * (k + 1))
    raise last if last else RuntimeError("unreachable")


def check_doi(doi: str) -> tuple:
    url = f"https://api.crossref.org/works/{urllib.parse.quote(doi, safe='')}"
    m = json.loads(fetch(url))["message"]
    authors = m.get("author") or []
    first = authors[0].get("family", "?") if authors else "?"
    year = (m.get("issued", {}).get("date-parts") or [[None]])[0][0]
    title = (m.get("title") or ["<no title>"])[0]
    return f"{first} {year}", title


def check_arxiv(ident: str) -> tuple:
    xml = fetch(f"http://export.arxiv.org/api/query?id_list={ident}&max_results=1")
    text = xml.decode("utf-8", "replace")
    title = re.search(r"<entry>.*?<title>(.*?)</title>", text, re.S)
    author = re.search(r"<entry>.*?<author>\s*<name>(.*?)</name>", text, re.S)
    if not title:
        raise ValueError("no entry returned")
    return (author.group(1) if author else "?"), " ".join(title.group(1).split())


def main() -> None:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    md = MANUSCRIPT.read_text(encoding="utf-8")

    # DOIs may contain parentheses; stop at whitespace or punctuation that cannot
    # belong to one, and drop a trailing sentence period.
    dois = sorted({d.rstrip(".") for d in re.findall(r"10\.\d{4,9}/[^\s,;<>]+", md)})
    arxivs = sorted(set(re.findall(r"arXiv:([\w.\-/]+\d)", md)))

    failed, unchecked = [], []

    print(f"{len(dois)} DOIs\n")
    for d in dois:
        try:
            who, title = check_doi(d)
            print(f"  ok        {d}\n            {who}: {title[:88]}")
        except Exception as e:
            print(f"  FAILED    {d}  ({type(e).__name__})")
            failed.append(d)

    print(f"\n{len(arxivs)} arXiv identifiers\n")
    for a in arxivs:
        try:
            who, title = check_arxiv(a)
            print(f"  ok        arXiv:{a}\n            {who}: {title[:88]}")
        except urllib.error.HTTPError as e:
            if e.code == 406:
                print(f"  unchecked arXiv:{a}  (arXiv API blocked from this network)")
                unchecked.append(a)
            else:
                print(f"  FAILED    arXiv:{a}  HTTP {e.code}")
                failed.append(a)
        except Exception as e:
            print(f"  FAILED    arXiv:{a}  ({type(e).__name__})")
            failed.append(a)

    total = len(dois) + len(arxivs)
    print(f"\n{total - len(failed) - len(unchecked)} of {total} resolved, "
          f"{len(unchecked)} unchecked, {len(failed)} failed")
    if failed:
        print("failed: " + ", ".join(failed))
    raise SystemExit(1 if failed else 0)


if __name__ == "__main__":
    main()
