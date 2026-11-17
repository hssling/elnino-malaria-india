"""Append the spatial/trend keys (papers/A/revision/spatial_trend_keys.json, written by
spatial_trends.py) to papers/A/results.json after `Rscript A/07_revision.R --step=results`.
Text append: every existing byte of results.json is kept; re-running replaces only the appended block."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RES = ROOT / "papers/A/results.json"
SP = ROOT / "papers/A/revision/spatial_trend_keys.json"
SCRIPT = "papers/A/revision/spatial_trends.py"


def main() -> int:
    text = RES.read_text(encoding="utf-8")
    current = json.loads(text)
    new = json.loads(SP.read_text(encoding="utf-8"))
    present = [k for k in new if k in current]
    if present:
        if len(present) != len(new):
            sys.exit("partial spatial keys present; rerun 07_revision.R --step=results first")
        # strip the previously appended block (it starts at the first spatial key)
        cut = text.index('  "' + next(iter(new)) + '"')
        text = text[:cut].rstrip().rstrip(",") + "\n}\n"
    body = text.rstrip()
    assert body.endswith("}")
    body = body[:-1].rstrip()
    parts = []
    for k, v in new.items():
        entry = {"value": v["value"], "display": v["display"], "script": SCRIPT,
                 "description": v["description"]}
        block = json.dumps({k: entry}, indent=2, ensure_ascii=False)[2:-2]
        parts.append(block)
    out = body + ",\n" + ",\n".join(parts) + "\n}\n"
    merged = json.loads(out)
    assert len(merged) == len(json.loads(text)) + len(new)
    RES.write_text(out, encoding="utf-8", newline="\n")
    print(f"appended {len(new)} keys; total {len(merged)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
