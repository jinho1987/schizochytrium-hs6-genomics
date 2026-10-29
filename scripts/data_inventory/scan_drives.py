import os, sys, json, time, gzip, zipfile, re
from collections import defaultdict

ROOTS = ["I:\\",
         r"C:\Users\user\Documents", r"C:\Users\user\Desktop", r"C:\Users\user\Downloads",
         r"C:\Users\user\HS2 ALE 실험 데이터", r"C:\Users\user\OneDrive"]
SKIP = {"$recycle.bin", "system volume information", "appdata", "anaconda3", "miniconda3", "envs",
        "pip_cache", "hf_cache", ".git", "node_modules", "site-packages", "__pycache__", "recovery",
        "seagate", "start_here_mac.app", ".ipynb_checkpoints", "%ee36cc", "lib", ".conda", "pkgs"}

# (category, subtype) by suffix; longest match first
EXT = [
 ("Sequencing reads", "FASTQ", [".fastq.gz", ".fq.gz", ".fastq", ".fq", ".fastq.bz2", ".fq.bz2", ".fastq.zst"]),
 ("Sequencing reads", "Long-read (PacBio/ONT)", [".subreads.bam", ".hifi_reads.bam", ".ccs.bam", ".fast5", ".pod5", ".bax.h5", ".bas.h5"]),
 ("Alignments", "BAM/CRAM/SAM", [".bam", ".cram", ".sam"]),
 ("Alignments", "Index", [".bai", ".csi", ".crai"]),
 ("Variants", "VCF/BCF", [".vcf.gz", ".vcf", ".bcf", ".gvcf.gz", ".g.vcf.gz"]),
 ("Assemblies & sequences", "FASTA", [".fasta", ".fa", ".fna", ".faa", ".ffn", ".fas", ".fasta.gz", ".fa.gz", ".fna.gz", ".faa.gz", ".scaf", ".contigs"]),
 ("Assemblies & sequences", "GenBank/EMBL", [".gb", ".gbk", ".gbff", ".embl"]),
 ("Assemblies & sequences", "Assembly graph", [".gfa", ".gfa.gz"]),
 ("Annotation", "GFF/GTF/BED", [".gff", ".gff3", ".gtf", ".gff.gz", ".gff3.gz", ".gtf.gz", ".bed", ".bed.gz", ".narrowpeak", ".broadpeak"]),
 ("Annotation", "Functional (eggNOG/KEGG/GO)", [".emapper.annotations", ".annotations", ".ko", ".kegg"]),
 ("Expression & counts", "Count matrices", [".counts", ".featurecounts", ".genes.results", ".isoforms.results", ".quant.sf", ".sf", ".abundance.tsv", ".mtx", ".mtx.gz"]),
 ("Expression & counts", "Single-cell / HDF5", [".h5ad", ".loom", ".h5", ".hdf5", ".zarr", ".h5seurat"]),
 ("Coverage tracks", "bigWig/bedGraph", [".bw", ".bigwig", ".bedgraph", ".bdg", ".tdf"]),
 ("Proteomics & metabolomics", "MS raw/processed", [".mzml", ".mzxml", ".mgf", ".raw", ".wiff", ".d", ".mzid", ".pepxml", ".cdf", ".qgd", ".lcd"]),
 ("Structures", "PDB/mmCIF", [".pdb", ".cif", ".mmcif", ".pdbqt", ".sdf", ".mol2"]),
 ("Phylogeny & alignments", "MSA/trees", [".aln", ".clustal", ".sto", ".stockholm", ".phy", ".phylip", ".nwk", ".newick", ".tree", ".treefile", ".meg", ".mas", ".hmm", ".nex"]),
 ("Reference DBs", "BLAST/index", [".nhr", ".nin", ".nsq", ".phr", ".pin", ".psq", ".ndb", ".not", ".ntf", ".nto", ".njs", ".bt2", ".ht2", ".mmi", ".fai", ".dict", ".sa", ".amb", ".ann", ".bwt", ".pac"]),
 ("QC reports", "FastQC/MultiQC", ["_fastqc.zip", "_fastqc.html", "multiqc_report.html", ".qc.html"]),
 ("qPCR & instrument", "qPCR/plate", [".eds", ".pcrd", ".rdml", ".sds", ".lc96p", ".ixo"]),
 ("Tabular data", "CSV/TSV/Excel", [".csv", ".tsv", ".tab", ".xlsx", ".xls", ".parquet", ".feather", ".rds", ".rdata", ".txt.gz"]),
 ("Archives", "zip/tar/egg", [".zip", ".tar.gz", ".tgz", ".tar", ".7z", ".egg", ".rar", ".gz", ".alz"]),
 ("Disk images", "VM/ISO", [".vhdx", ".iso", ".vmdk"]),
]
LOOKUP = sorted([(s, c, t) for c, t, ss in EXT for s in ss], key=lambda x: -len(x[0]))
OMICS = {"Sequencing reads", "Alignments", "Variants", "Assemblies & sequences", "Annotation", "Expression & counts",
         "Coverage tracks", "Proteomics & metabolomics", "Phylogeny & alignments", "Structures", "Reference DBs", "QC reports", "qPCR & instrument"}

def classify(name):
    n = name.lower()
    if re.search(r"polymorphic_(snp|indel|sv)", n) and n.endswith((".xlsx", ".xls", ".csv", ".pdf")):
        return "Variants", "Variant tables (vendor)"
    if n.startswith("data_statistics") and n.endswith((".xlsx", ".xls")):
        return "QC reports", "Run statistics (vendor)"
    if re.search(r"analysis_report|wgs_analysis|_report_v\d", n) and n.endswith(".pdf"):
        return "QC reports", "Vendor report (PDF)"
    if re.search(r"(^|[_\s-])degs?[_\s.-]", n) and n.endswith((".zip", ".xlsx", ".xls", ".csv", ".tsv")):
        return "Expression & counts", "Count/DE tables"
    if n.endswith(".txt") and any(k in n for k in ("count", "featurecount", "deseq", "edger", "tpm", "fpkm")):
        return "Expression & counts", "Count/DE tables"
    if n.endswith((".csv", ".tsv", ".xlsx", ".txt")) and any(k in n for k in ("deg", "deseq", "diff_exp", "differential", "_de_", "volcano")):
        return "Expression & counts", "Count/DE tables"
    for s, c, t in LOOKUP:
        if n.endswith(s):
            return c, t
    return None, None

def fastq_peek(path):
    try:
        op = gzip.open if path.lower().endswith(".gz") else open
        with op(path, "rt", errors="replace") as f:
            lens, hdr = [], None
            for i, line in enumerate(f):
                if i % 4 == 0 and hdr is None: hdr = line.strip()[:120]
                if i % 4 == 1: lens.append(len(line.strip()))
                if i >= 4000: break
        if not lens: return {}
        return {"read_len": round(sum(lens) / len(lens)), "max_len": max(lens), "header": hdr}
    except Exception as e:
        return {"err": str(e)[:80]}

def vcf_peek(path):
    try:
        op = gzip.open if path.lower().endswith(".gz") else open
        with op(path, "rt", errors="replace") as f:
            contigs, src = 0, None
            for line in f:
                if line.startswith("##contig"): contigs += 1
                elif line.startswith("##source=") and not src: src = line.strip()[9:80]
                elif line.startswith("#CHROM"):
                    s = line.rstrip("\n").split("\t")[9:]
                    return {"samples": s[:30], "n_samples": len(s), "contigs": contigs, "source": src}
                elif not line.startswith("#"): break
    except Exception as e:
        return {"err": str(e)[:80]}
    return {}

def fasta_peek(path, size):
    if size > 3e9: return {}
    try:
        op = gzip.open if path.lower().endswith(".gz") else open
        n, total, first = 0, 0, None
        with op(path, "rt", errors="replace") as f:
            for line in f:
                if line.startswith(">"):
                    n += 1
                    if first is None: first = line.strip()[:100]
                else: total += len(line.strip())
        return {"n_seqs": n, "total_bp": total, "first": first}
    except Exception as e:
        return {"err": str(e)[:80]}

def zip_peek(path):
    try:
        with zipfile.ZipFile(path) as z:
            cats = defaultdict(lambda: [0, 0]); names = []
            for i in z.infolist():
                if i.is_dir(): continue
                c, t = classify(os.path.basename(i.filename))
                k = c or "Other"
                cats[k][0] += 1; cats[k][1] += i.file_size
                if c in OMICS and len(names) < 12: names.append(i.filename[-90:])
            return {"zip_contents": {k: v for k, v in cats.items()}, "zip_omics_examples": names}
    except Exception as e:
        return {"err": str(e)[:80]}

def project_of(path, root):
    rel = os.path.relpath(path, root)
    parts = rel.split(os.sep)
    drive = root.rstrip("\\")
    if len(parts) == 1: return drive, "(root files)", ""
    sub = parts[1] if len(parts) > 2 else ""
    return drive, parts[0], sub

files = []; dirstats = defaultdict(lambda: {"n": 0, "bytes": 0, "omics_n": 0, "omics_bytes": 0, "cats": defaultdict(lambda: [0, 0]), "newest": 0, "oldest": 9e12})
t0 = time.time(); total_seen = 0
for root in ROOTS:
    if not os.path.exists(root): continue
    for dp, dns, fns in os.walk(root):
        dns[:] = [d for d in dns if d.lower() not in SKIP and not d.startswith("$")]
        # treat .d folders (Agilent) and .zarr as files
        for d in list(dns):
            if d.lower().endswith((".d", ".zarr")):
                fns.append(d); dns.remove(d)
        for fn in fns:
            p = os.path.join(dp, fn); total_seen += 1
            try:
                st = os.stat(p)
                size = st.st_size if os.path.isfile(p) else sum(os.path.getsize(os.path.join(a, b)) for a, _, bs in os.walk(p) for b in bs)
            except Exception:
                continue
            c, t = classify(fn)
            drive, proj, sub = project_of(p, root)
            if root.startswith("C:"): drive, proj = "C:", os.path.basename(root.rstrip("\\")) + "\\" + proj
            k = (drive, proj)
            ds = dirstats[k]
            ds["n"] += 1; ds["bytes"] += size
            ds["newest"] = max(ds["newest"], st.st_mtime); ds["oldest"] = min(ds["oldest"], st.st_mtime)
            if c:
                ds["cats"][c][0] += 1; ds["cats"][c][1] += size
                if c in OMICS:
                    ds["omics_n"] += 1; ds["omics_bytes"] += size
            keep = c in OMICS or (c in ("Archives", "Disk images") and size > 50e6) or (c == "Archives" and re.search(r"analysis_result|assembly|omics|result", fn.lower()) is not None) or (c == "Tabular data" and size > 5e6)
            if keep:
                rec = {"path": p, "name": fn, "drive": drive, "project": proj, "sub": sub, "cat": c, "type": t, "bytes": size, "mtime": st.st_mtime}
                if c == "Reference DBs" or (c == "Alignments" and t == "Index"):
                    rec["minor"] = True
                files.append(rec)
print("walked", total_seen, "files in", round(time.time() - t0), "s; kept", len(files), file=sys.stderr)

# peeks
for r in files:
    p, n = r["path"], r["name"].lower()
    if r["type"] == "FASTQ": r.update(fastq_peek(p))
    elif r["type"] == "VCF/BCF" and not n.endswith(".bcf"): r.update(vcf_peek(p))
    elif r["type"] == "FASTA": r.update(fasta_peek(p, r["bytes"]))
    elif n.endswith(".zip") and r["bytes"] < 30e9: r.update(zip_peek(p))
print("peeked in", round(time.time() - t0), "s", file=sys.stderr)

out = {"generated": time.strftime("%Y-%m-%d %H:%M"), "roots": ROOTS, "files_seen": total_seen,
       "projects": [{"drive": k[0], "project": k[1], **{kk: vv for kk, vv in v.items() if kk != "cats"},
                     "cats": {c: v2 for c, v2 in v["cats"].items()}} for k, v in dirstats.items()],
       "files": files}
json.dump(out, open(sys.argv[1], "w", encoding="utf-8"), ensure_ascii=False)
