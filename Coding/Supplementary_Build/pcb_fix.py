# Rebuild every PCB CAS from its congener number via PubChem ("PCB 99" resolves
# directly). The stored values are misaligned across the congener-pair rows --
# several carry a neighbouring pair's numbers -- so patching individually would
# just move the collision around.
import csv, json, re, subprocess, time, urllib.parse
def ok(c):
    m=re.fullmatch(r'(\d{2,7})-(\d{2})-(\d)',c)
    if not m: return False
    b=m.group(1)+m.group(2)
    return sum(int(d)*(i+1) for i,d in enumerate(reversed(b)))%10==int(m.group(3))
CACHE={}
def cas_for(cong):
    if cong in CACHE: return CACHE[cong]
    u=("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
       +urllib.parse.quote(f"PCB {cong}")+"/synonyms/JSON")
    r=subprocess.run(["curl","-s","--max-time","20",u],capture_output=True,text=True)
    out=None
    try:
        inf=json.loads(r.stdout).get("InformationList",{}).get("Information",[])
        if inf:
            cas=[s for s in inf[0].get("Synonym",[]) if re.fullmatch(r'\d{2,7}-\d{2}-\d',s) and ok(s)]
            if cas: out=cas[0]
    except Exception: pass
    CACHE[cong]=out; time.sleep(0.25); return out

p="Input/Harmonized_Datasets/Targeted_Analyte_Category_Mapping.csv"
rows=list(csv.DictReader(open(p)))
changed=unresolved=same=0
for r in rows:
    c=r["Code"]
    if not c.startswith("PCB"): continue
    congs=re.findall(r'\d+', c)          # PCB90_101_89 -> [90,101,89]
    if not congs: continue
    got=[cas_for(x) for x in congs]
    if any(g is None for g in got):
        if r["CAS"] not in ("","NA"):
            print(f"  BLANK {c:<16} {r['CAS']:<34} (congener(s) unresolved)")
            r["CAS"], r["CAS_Source"] = "", "unconfirmed"; unresolved+=1
        continue
    new="_".join(got)
    if r["CAS"]!=new:
        print(f"  FIX   {c:<16} {r['CAS']:<34} -> {new}")
        r["CAS"], r["CAS_Source"] = new, "pubchem_congener"; changed+=1
    else: same+=1
w=csv.DictWriter(open(p,"w",newline=""),fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
print(f"\nPCB rows: {changed} corrected, {same} already correct, {unresolved} blanked")
