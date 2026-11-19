import json, re, subprocess, urllib.parse

def ok(c):
    m = re.fullmatch(r'(\d{2,7})-(\d{2})-(\d)', c)
    if not m: return False
    b = m.group(1) + m.group(2)
    return sum(int(d)*(i+1) for i,d in enumerate(reversed(b))) % 10 == int(m.group(3))

def pubchem(name):
    url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
           + urllib.parse.quote(name) + "/synonyms/JSON")
    r = subprocess.run(["curl","-s","--max-time","25",url], capture_output=True, text=True)
    try: d = json.loads(r.stdout)
    except Exception: return None, "no JSON"
    info = d.get("InformationList", {}).get("Information", [])
    if not info: return None, d.get("Fault",{}).get("Code","no match")
    syn = info[0].get("Synonym", [])
    cas = [s for s in syn if re.fullmatch(r'\d{2,7}-\d{2}-\d', s) and ok(s)]
    return (info[0].get("CID"), cas[:5]), None

QUERIES = [
 ("TDCCA", "55401-03-06", ["trans-3-(2,2-dichlorovinyl)-2,2-dimethylcyclopropanecarboxylic acid","trans-DCCA"]),
 ("MCIOP", "895844-09-7", ["mono-carboxyisooctyl phthalate","monocarboxyisooctyl phthalate"]),
 ("MHXP",  "245539-57-9", ["monohexyl phthalate","mono-n-hexyl phthalate"]),
 ("PTU",   "96-45-7",     ["propylenethiourea","4-methyl-2-imidazolidinethione"]),
 ("ETU",   "96-45-7",     ["ethylenethiourea"]),
 ("MNBP",  "34-74-2",     ["monobutyl phthalate","mono-n-butyl phthalate"]),
]
for code, stored, names in QUERIES:
    print(f"\n=== {code}  stored {stored}  (checksum {'OK' if ok(stored) else 'FAIL'}) ===")
    for n in names:
        res, err = pubchem(n)
        if err: print(f"   {n!r}: {err}"); continue
        cid, cas = res
        print(f"   {n!r}\n     CID {cid}  CAS: {', '.join(cas) if cas else '(none listed)'}")
