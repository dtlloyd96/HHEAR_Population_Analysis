# A CAS Registry Number ends in a check digit: sum the preceding digits weighted
# 1,2,3,... from the right, mod 10. This validates every entry offline and
# catches transcription errors deterministically -- no lookup required.
import csv, re, sys

def cas_ok(c):
    m = re.fullmatch(r'(\d{2,7})-(\d{2})-(\d)', c)
    if not m: return None                      # not CAS-shaped
    body = m.group(1) + m.group(2)
    total = sum(int(d) * (i + 1) for i, d in enumerate(reversed(body)))
    return total % 10 == int(m.group(3))

rows = list(csv.DictReader(open("Output/Hero_Numbers/SuppTable2_FINAL.csv")))
good = bad = missing = malformed = 0
bad_rows, mal_rows = [], []
for r in rows:
    raw = (r["CAS"] or "").strip()
    if raw in ("", "NA"):
        missing += 1; continue
    # PCB congener pairs are stored as two CAS joined by underscore
    parts = raw.split("_")
    results = [cas_ok(p) for p in parts]
    if any(x is None for x in results):
        malformed += 1; mal_rows.append((r["Code"], raw, r["Category"])); continue
    if all(results): good += 1
    else:
        bad += 1; bad_rows.append((r["Code"], r["Full_Name"][:40], raw, r["Category"]))

print(f"rows: {len(rows)}   checksum OK: {good}   FAILS: {bad}   malformed: {malformed}   missing: {missing}\n")
if bad_rows:
    print("=== CHECKSUM FAILURES (digit is wrong somewhere) ===")
    for c,n,v,cat in bad_rows: print(f"  {c:<16} {v:<26} {cat:<16} {n}")
if mal_rows:
    print("\n=== NOT CAS-SHAPED ===")
    for c,v,cat in mal_rows: print(f"  {c:<16} {v:<30} {cat}")
