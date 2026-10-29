import json, re, os, time, csv
from collections import defaultdict, Counter

d = json.load(open(r"C:\invscan\inventory.json", encoding="utf-8"))
F = d["files"]
OUT = os.path.join(os.getcwd(), "inventory_data"); os.makedirs(OUT, exist_ok=True)

# fix: "ALE-1_1.fastq (1).gz" style names are FASTQ, not archives
for f in F:
    if f["cat"] == "Archives" and re.search(r"\.f(ast)?q", f["name"].lower()):
        f["cat"], f["type"] = "Sequencing reads", "FASTQ"

PLATFORM = [("m84", "PacBio Revio (HiFi)"), ("LH", "Illumina NovaSeq X"), ("A0", "Illumina NovaSeq 6000"),
            ("VH", "Illumina NextSeq 2000"), ("ST-E", "Illumina HiSeq X"), ("D00", "Illumina HiSeq 2500"), ("SRR", "Public (SRA)")]
def platform(h):
    h = (h or "").lstrip("@")
    for p, n in PLATFORM:
        if h.startswith(p): return n
    return ""

# curated datasets: (drive, path regex relative to drive root, id, title, organism, assay)
DS = [
 ("I:", r"^Chlamy WGS\\", "chlamy-wgs", "Chlamy WGS: CC125 vs 4-9", "Chlamydomonas reinhardtii", "WGS resequencing"),
 ("I:", r"^chlamy_wgs_variant_calling\\", "chlamy-varcall", "Chlamy WGS variant calling", "Chlamydomonas reinhardtii", "Variant calling (SNV/SV/CNV)"),
 ("I:", r"^Chlamy mRNA seq\\", "chlamy-rnaseq", "Chlamy mRNA-seq: WT vs 4-9", "Chlamydomonas reinhardtii", "RNA-seq"),
 ("I:", r"^HS2 vs\. ALE_long read seq_260621\\", "hs2-ale-hifi", "HS2 vs ALE HiFi long-read genomes", "Chlorella sp. HS2", "Long-read genome + SV"),
 ("I:", r"^Chlorella\\1\. ALE vs\. HS2", "hs2-ale-2304", "ALE vs HS2 WGS (2023.04)", "Chlorella sp. HS2", "WGS resequencing"),
 ("I:", r"^Chlorella\\2\. HMF vs\. HS2", "hs2-hmf-2311", "HMF vs HS2 WGS (2023.11)", "Chlorella sp. HS2", "WGS resequencing"),
 ("I:", r"^Chlorella\\3\. HMF vs\. HS2 vs\. ALE", "hs2-hmf-ale-2410", "HMF vs HS2 vs ALE WGS (2024.10)", "Chlorella sp. HS2", "WGS resequencing"),
 ("I:", r"^Chlorella\\ALE Chlorella sp\. HS2 in hydrolysate", "hs2-ale-hydro", "ALE HS2 in hydrolysate WGS", "Chlorella sp. HS2", "WGS resequencing"),
 ("I:", r"^Chlorella\\조대현", "hs2-omics-legacy", "HS2 omics (2016–17 genome + RNA-seq)", "Chlorella sp. HS2", "Genome + RNA-seq"),
 ("I:", r"^HS19 vs\. HS2\\", "hs19-hs2", "HS19 vs HS2 WGS", "Chlorella sp. HS2", "WGS resequencing"),
 ("I:", r"^HS2_RNA-seq_Hetero\\", "hs2-hetero", "HS2 heterotrophy RNA-seq (2019)", "Chlorella sp. HS2", "RNA-seq"),
 ("I:", r"^26\.10\.16_", "pld-lipidomics", "PLD KO vs WT4 lipidomics (vendor)", "Schizochytrium sp.", "Lipidomics (LC-MS)"),
 ("I:", r"^발표자료 수업자료 위주 백업\\.*HS6 work_Folding", "hs6-docking-bak", "HS6 folding & docking (backup)", "Schizochytrium sp.", "Structures + variants"),
 ("I:", r"^발표자료 수업자료 위주 백업\\", "hs6-genome-bak", "HS6 genome work (backup)", "Schizochytrium sp.", "Genome assembly + variants"),
 ("C:", r"^Desktop\\HS6 vs\. 4_RNA-seq_260720\\", "hs6-rnaseq", "HS6 vs #4 RNA-seq + genome", "Schizochytrium sp.", "RNA-seq + genome"),
 ("C:", r"^Desktop\\Published manuscripts\\Paper with SB", "hs2-hsv-tx", "HS-V transcriptome (published)", "Chlorella sp. HS2", "RNA-seq"),
 ("C:", r"^Downloads\\(4-9|CC125)_", "dl-chlamy-wgs", "Chlamy WGS copy in Downloads", "Chlamydomonas reinhardtii", "WGS resequencing"),
 ("C:", r"^Downloads\\HS19_", "dl-hs19", "HS19 read copy in Downloads", "Chlorella sp. HS2", "WGS resequencing"),
 ("C:", r"^Downloads\\(DEG|All_DEGs|HS6_vs_WT4)", "dl-deg", "DEG result tables in Downloads", "Mixed", "RNA-seq results"),
 ("C:", r"^Downloads\\GSE146789", "dl-gse146789", "GSE146789 HS2 Trinity transcriptome", "Chlorella sp. HS2", "Public transcriptome"),
 ("C:", r"^Desktop\\20240619 ", "pigment-2406", "Pigment HPLC runs (2024.06)", "Microalgae", "Pigment HPLC"),
 ("C:", r"^Desktop\\Pigment analysis results_221005", "pigment-2210", "Pigment HPLC runs (2022.10)", "Microalgae", "Pigment HPLC"),
 ("C:", r"^Documents\\MEGA X", "mega-phylo", "MEGA X alignments & trees", "Mixed", "Phylogenetics"),
 ("C:", r"^Desktop\\Schizochytrium (MS|transformation)", "schizo-constructs", "Schizochytrium constructs & sequences", "Schizochytrium sp.", "Plasmid / sequences"),
]
SKIP_PROJ = ("GBox", "openLCA", "features", "plugins", "Documents\\R")

def rel(f):
    p = f["path"]
    if f["drive"] == "I:": return p[3:]
    return p[len(r"C:\Users\user") + 1:]

def assign(f):
    r = rel(f)
    for drv, rx, i, *_ in DS:
        if drv == f["drive"] and re.search(rx, r): return i
    return None

def is_omics_archive(f):
    zc = f.get("zip_contents") or {}
    return any(k not in ("Other", "Archives", "Tabular data") for k in zc)

def sample_of(name):
    n = name.lower()
    n = re.sub(r"\s*\(\d+\)", "", n)
    n = re.sub(r"\.(fastq|fq)(\.gz)?$", "", n)
    n = re.sub(r"_trim", "", n)
    n = re.sub(r"(_r?[12])$", "", n)
    return n

meta = {x[2]: dict(id=x[2], dataset=x[3], organism=x[4], assay=x[5], drive=x[0]) for x in DS}
agg = defaultdict(lambda: dict(files=0, bytes=0, raw_bytes=0, samples=set(), lens=Counter(), plats=Counter(), types=Counter(),
                               first=9e12, last=0, folders=Counter(), has_raw=False, has_proc=False, has_arch=False))
files_rows = []
for f in F:
    if any(s in f["project"] for s in SKIP_PROJ): continue
    omics = f["cat"] not in ("Archives", "Disk images", "Tabular data") or (f["cat"] == "Archives" and is_omics_archive(f))
    did = assign(f)
    if did is None and not omics: continue
    if did is None and f["cat"] == "Tabular data": continue
    if did and f["cat"] in ("Tabular data", "Disk images"): continue
    if did is None:
        did = "misc:" + f["drive"] + f["project"]
        meta.setdefault(did, dict(id=did, dataset=(f["project"].split("\\")[0] + " loose files") if f["project"].endswith("(root files)") else f["project"].split("\\")[-1],
                                  organism="—", assay="Other omics files", drive=f["drive"]))
    a = agg[did]
    a["files"] += 1; a["bytes"] += f["bytes"]; a["types"][f["type"]] += 1
    a["first"] = min(a["first"], f["mtime"]); a["last"] = max(a["last"], f["mtime"])
    a["folders"][os.path.dirname(f["path"])] += 1
    if f["type"] == "FASTQ":
        a["has_raw"] = True; a["raw_bytes"] += f["bytes"]
        if not f["name"].lower().startswith(("r1.", "r2.")) and f["bytes"] > 5e7:
            a["samples"].add(sample_of(f["name"]))
        if f.get("read_len"): a["lens"][f["read_len"] if f["read_len"] < 1000 else "HiFi"] += 1
        p = platform(f.get("header"))
        if p: a["plats"][p] += 1
    elif f["cat"] == "Archives": a["has_arch"] = True
    else: a["has_proc"] = True
    if f["cat"] == "Structures": continue  # summarised, not listed one by one
    detail = ""
    if f.get("read_len"): detail = f"{'HiFi ~' + format(f['read_len'], ',') + ' bp' if f['read_len'] > 1000 else str(f['read_len']) + ' bp reads'} · {platform(f.get('header'))}"
    elif f.get("n_samples") is not None: detail = f"{f['n_samples']} sample(s): {', '.join(f['samples'][:4])}" + (f" · {f['source']}" if f.get("source") else "")
    elif f.get("n_seqs"): detail = f"{f['n_seqs']:,} seqs · {f['total_bp']/1e6:.1f} Mb"
    elif f.get("zip_contents"):
        zc = {k: v[0] for k, v in f["zip_contents"].items() if k not in ("Other",)}
        detail = "zip: " + ", ".join(f"{v} {k.lower()}" for k, v in sorted(zc.items(), key=lambda x: -x[1]))
    files_rows.append(dict(path=f["path"], file=f["name"], type=f["type"], category=f["cat"], dataset=meta[did]["dataset"], dataset_id=did,
                           drive=f["drive"], folder=os.path.dirname(f["path"]), gb=round(f["bytes"] / 1e9, 3),
                           modified=time.strftime("%Y-%m-%d", time.localtime(f["mtime"])), detail=detail))

ds_rows = []
for did, a in agg.items():
    m = meta[did]
    if a["bytes"] < 1e6 and did.startswith("misc:"): continue
    state = " + ".join(x for x, ok in (("raw reads", a["has_raw"]), ("processed", a["has_proc"]), ("archived", a["has_arch"])) if ok)
    top = a["types"].most_common(3)
    loc = os.path.commonpath(list(a["folders"])) if a["folders"] else ""
    lens = a["lens"].most_common(1)
    ds_rows.append(dict(id=did, dataset=m["dataset"], organism=m["organism"], assay=m["assay"], drive=m["drive"],
                        platform=" / ".join(p for p, _ in a["plats"].most_common(2)),
                        read_length=(str(lens[0][0]) + (" bp" if lens[0][0] != "HiFi" else "")) if lens else "",
                        samples=len(a["samples"]) or None, files=a["files"], gb=round(a["bytes"] / 1e9, 2), raw_gb=round(a["raw_bytes"] / 1e9, 2),
                        state=state, main_types=", ".join(f"{t} ({n})" for t, n in top),
                        first=time.strftime("%Y-%m-%d", time.localtime(a["first"])), last=time.strftime("%Y-%m-%d", time.localtime(a["last"])),
                        location=loc, curated=not did.startswith("misc:")))
ds_rows.sort(key=lambda r: -r["gb"])

# duplicates: same name + size in different folders
groups = defaultdict(list)
for f in F:
    if f["bytes"] > 20e6 and f["cat"] not in ("Disk images",) and not any(s in f["project"] for s in SKIP_PROJ):
        groups[(re.sub(r"\s*\(\d+\)", "", f["name"].lower()), f["bytes"])].append(f)
dup_rows = []
for (n, b), fs in groups.items():
    folders = sorted({os.path.dirname(x["path"]) for x in fs})
    if len(fs) < 2: continue
    if fs[0]["cat"] in ("Archives", "Tabular data") and not is_omics_archive(fs[0]): continue
    dup_rows.append(dict(file=fs[0]["name"], type=fs[0]["type"], gb_each=round(b / 1e9, 3), copies=len(fs),
                         reclaimable_gb=round(b * (len(fs) - 1) / 1e9, 3), locations=" | ".join(x["path"] for x in fs)))
dup_rows.sort(key=lambda r: -r["reclaimable_gb"])

# file-type summary (omics only)
tsum = defaultdict(lambda: [0, 0])
for r in files_rows: tsum[(r["category"], r["type"])][0] += 1; tsum[(r["category"], r["type"])][1] += r["gb"]
struct = [f for f in F if f["cat"] == "Structures" and not any(s in f["project"] for s in SKIP_PROJ)]
tsum[("Structures", "PDB/mmCIF")] = [len(struct), sum(x["bytes"] for x in struct) / 1e9]
type_rows = [dict(category=c, type=t, files=v[0], gb=round(v[1], 2)) for (c, t), v in sorted(tsum.items(), key=lambda x: -x[1][1])]

# drive context
drv = []
for name, used, free in (("C:", 295.3, 635.2), ("I:", 946.5, 2779.3)):
    om = sum(r["gb"] for r in ds_rows if r["drive"] == name)
    drv.append(dict(drive=name, used_gb=used, free_gb=free, omics_gb=round(om, 1)))
wsl = [f for f in F if f["name"].lower() == "ext4.vhdx"]

def dump(name, rows):
    json.dump(rows, open(os.path.join(OUT, name), "w", encoding="utf-8"), ensure_ascii=False)
    print(name, len(rows), os.path.getsize(os.path.join(OUT, name)))
dump("datasets.json", ds_rows); dump("files.json", files_rows); dump("duplicates.json", dup_rows); dump("filetypes.json", type_rows)
dump("drives.json", [dict(**x, wsl_vhdx_gb=round(sum(w["bytes"] for w in wsl if w["drive"] == x["drive"]) / 1e9, 1)) for x in drv])
for r in ds_rows: print(r["drive"], r["gb"], r["samples"], r["read_length"], r["platform"], "|", r["dataset"], "|", r["state"])
print("dups", len(dup_rows), round(sum(r["reclaimable_gb"] for r in dup_rows), 1))
for r in dup_rows[:15]: print(r["copies"], r["reclaimable_gb"], r["file"])

# ---- lipidomics: expand the vendor package (PLD KO vs WT4) ----
import zipfile
LZ = r"I:\26.10.16_고마바이오_WTvsPLD.zip"
lip_rows = []
if os.path.exists(LZ):
    z = zipfile.ZipFile(LZ)
    base = next(n for n in z.namelist() if n.endswith("PLD_KO_WT4_DiffScreening.xls"))
    lines = z.read(base).decode("utf-8", "replace").splitlines()
    hdr = lines[0].split("\t"); ix = {h: i for i, h in enumerate(hdr)}
    cls = defaultdict(lambda: {"n": 0, "up": 0, "down": 0})
    for ln in lines[1:]:
        c = ln.split("\t")
        if len(c) < len(hdr): continue
        g = c[ix["Sub Class"]] if c[ix["Sub Class"]] not in ("", "NA") else "Unassigned"
        st = c[ix["State"]].lower()
        cls[g]["n"] += 1
        if st.startswith("up"): cls[g]["up"] += 1
        elif st.startswith("down"): cls[g]["down"] += 1
    lip_rows = [dict(lipid_class=k, lipids=v["n"], up=v["up"], down=v["down"]) for k, v in sorted(cls.items(), key=lambda x: -x[1]["n"])]
    samples = [l.split("\t") for l in z.read(next(n for n in z.namelist() if n.endswith("SampleList.xls"))).decode().splitlines()[1:] if l.strip()]
    for r in ds_rows:
        if r["id"] == "pld-lipidomics":
            r.update(samples=len(samples), platform="LC-MS untargeted, pos + neg (BGI Dr.Tom)",
                     read_length=f"{len(lines) - 1} lipids", state="raw data + processed + archived",
                     main_types="LC-MS raw (2), result tables, report PDFs")
    pick = re.compile(r"(rawdata/(pos|neg)/Rawdata_1\.tar\.gz|SampleList\.xls|compare_group\.xls|DiffScreening\.xls|identiNumber\.xls|PCAscores\.xls|_Volcano\.xls|C_en\.pdf)$")
    zt = time.strftime("%Y-%m-%d", time.localtime(os.path.getmtime(LZ)))
    for i in z.infolist():
        if not pick.search(i.filename): continue
        m = re.search(r"rawdata/(pos|neg)/", i.filename)
        typ = f"LC-MS raw ({m.group(1)})" if m else ("Vendor report (PDF)" if i.filename.endswith(".pdf") else "Lipid result table")
        det = {"SampleList.xls": f"{len(samples)} samples: " + ", ".join(f"{s[1]} ({s[2]})" for s in samples),
               "DiffScreening.xls": f"{len(lines) - 1} lipids with fold change, p/q values, class, chain length",
               "identiNumber.xls": f"{len(lines) - 1} lipids identified"}
        files_rows.append(dict(path=LZ + "::" + i.filename, file=os.path.basename(i.filename), type=typ, category="Proteomics & metabolomics",
                               dataset="PLD KO vs WT4 lipidomics (vendor)", dataset_id="pld-lipidomics", drive="I:", folder=LZ + " (inside zip)",
                               gb=round(i.file_size / 1e9, 3), modified=zt,
                               detail=next((v for k, v in det.items() if i.filename.endswith(k)), "MS run archive" if m else "")))
    dump("datasets.json", ds_rows); dump("files.json", files_rows); dump("lipidomics.json", lip_rows)
    print(lip_rows[:12])

# ---- lineage (taxonomy, NCBI-style ranks) and strains per dataset ----
LINEAGE = {
 "Chlorella sp. HS2":         ["Archaeplastida", "Chlorophyta", "Trebouxiophyceae", "Chlorellales", "Chlorellaceae", "Chlorella sp. HS2"],
 "Chlamydomonas reinhardtii": ["Archaeplastida", "Chlorophyta", "Chlorophyceae", "Chlamydomonadales", "Chlamydomonadaceae", "Chlamydomonas reinhardtii"],
 "Schizochytrium sp.":        ["SAR (Stramenopiles)", "Bigyra", "Labyrinthulomycetes", "Thraustochytrida", "Thraustochytriidae", "Schizochytrium sp."],
}
ORG_FIX = {r"misc:C:Desktop\Trang-HS2 RNAseq": "Chlorella sp. HS2"}
STRAINS = {
 "hs2-ale-hifi": "HS2, ALE", "hs2-omics-legacy": "HS2", "hs2-hsv-tx": "HS2 wild type, HS-LZ mutants", "hs2-ale-hydro": "HS2 (W), ALE (M)",
 "hs2-ale-2304": "HS2, ALE", "hs2-hmf-2311": "HS2, HMF", "hs2-hmf-ale-2410": "HS2, HMF, ALE", "hs19-hs2": "HS2, HS19", "dl-hs19": "HS19",
 "hs2-hetero": "HS2", "dl-gse146789": "HS2 (public)", r"misc:C:Desktop\Trang-HS2 RNAseq": "HS2",
 "chlamy-rnaseq": "WT, 4-9 (ref. CC-4532)", "chlamy-wgs": "CC-125, 4-9", "dl-chlamy-wgs": "CC-125, 4-9", "chlamy-varcall": "CC-125, 4-9",
 "hs6-rnaseq": "HS6, #4 (WT)", "hs6-genome-bak": "HS6, #4 (ref. Aurantiochytrium sp. 9LR)", "hs6-docking-bak": "HS6, #4",
 "pld-lipidomics": "PLD KO, WT4", "schizo-constructs": "Schizochytrium sp. (ref. Aurantiochytrium sp. KH105)",
}
RANKS = ["supergroup", "phylum", "class", "order", "family", "species"]
for r in ds_rows:
    r["organism"] = ORG_FIX.get(r["id"], r["organism"])
    lin = LINEAGE.get(r["organism"])
    if lin:
        r.update(dict(zip(RANKS, lin)))
    else:
        label = "Microalgae, species not recorded" if r["organism"] == "Microalgae" else "Several organisms or none"
        r.update(supergroup="Unplaced", phylum="", **{"class": ""}, order="", family="", species=label)
    r["strains"] = STRAINS.get(r["id"], "")
dump("datasets.json", ds_rows)
for r in ds_rows: print(r["species"], "|", r["strains"], "|", r["dataset"])
