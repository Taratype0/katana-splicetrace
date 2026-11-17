# katana-motiftrace (`v2.3.0`)

[![CI](https://github.com/Taratype0/katana-motiftrace/actions/workflows/ci.yml/badge.svg)](https://github.com/Taratype0/katana-motiftrace/actions/workflows/ci.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

Standalone pure-Python package and CLI for **calibrated bulk RNA splicing RBP regulator discovery** and **mature-mRNA splicing motif consequence tracing** (`SE`, `A3SS`, `A5SS`, `MXE`, `RI`), bundling the decontaminated **Curated Splicing & Epitranscriptome Library `v3.0` (`169` RBPs / `621` experimental PWMs)** offline out-of-the-box with built-in `--species human|mouse` ortholog mapping and a zero-code **5-format extensible motif database engine**.

---

## Key Capabilities
- **Null-Calibrated RBP Regulator Discovery (`max-T` FWER)**: 1:5 (`log_depth, baseline_PSI, event_type`)-matched case-control & continuous rescue engine with `CpG` + `GC` + sequence-length FWL orthogonalization (`0%` false positives on real `GSE206660` negative controls).
- **Mature-mRNA Motif Consequence Tracing**: Reconstructs inclusion (`I`) vs exclusion (`E`) spliceforms across `SE`, `A3SS`, `A5SS`, `MXE`, and `RI` events to identify exact RBP binding sites skipped, gained, or shifted across exon-exon junctions.
- **Bundled Library `v3.0` (`169` RBPs / `621` PWMs)**: Covers canonical splicing regulators, nuclear $\text{m}^6\text{A}$ readers/writers (`YTHDC1` `[U/G]GGACU`, `YTHDC2`, `YTHDF1-3`, `RBM15/15B`, `METTL3`), and ENCODE RBNS/eCLIP disease factors (`RBM20`, `RBM10`, `RBM39`, `U2AF1`, `SF3B4`, `SUGP1`, `SRRM4`, `PUF60`).
- **Zero-Code Extensible Motif Database**: Directly import or permanently register `.meme`, `.motif` (HOMER/eCLIP), `.jaspar`/`.pfm`, `.tsv`/`.csv` (IUPAC/$k$-mer tables), or `.json.gz` packs.

---

## Installation

```bash
# Install directly from GitHub
pip install git+https://github.com/Taratype0/katana-motiftrace.git

# Or install pre-built wheel from GitHub Releases
pip install https://github.com/Taratype0/katana-motiftrace/releases/download/v2.3.0/katana_motiftrace-2.3.0-py3-none-any.whl
```

---

## Standalone CLI Usage (No Nextflow / R Required)

```bash
# 1. Direct rMATS directory or .tar.gz archive (human hg38 or mouse mm10/mm39)
katana-motiftrace \
  --rmats ./rmats_out \
  --fasta hg38.fa \
  --species human \
  --output-dir ./motiftrace_results

# 2. Mouse dataset with custom external motif files (.meme / .motif / .jaspar / .tsv)
katana-motiftrace \
  --rmats ./Setd2_rmats \
  --fasta mm10.fa \
  --species mouse \
  --extra-motifs ./custom_clip.motif ./my_lab_kmers.tsv \
  --output-dir ./Setd2_motiftrace_out

# 3. Inspect bundled + installed RBPs or permanently install a custom motif pack
katana-motiftrace --list-rbps --filter-rbp YTH
katana-motiftrace --install-pack ./encode_extra.meme --pack-name encode_extra
```

---

## Python API Usage

```python
import katana_motiftrace as kmt

lib = kmt.load_bundled_human_rbp_library(species="mouse")
events = kmt.load_rmats_measurable_events("Setd2_rmats/", "mm10.fa")
res = kmt.run_bulk_continuous_rescue_v22(events, "mm10.fa", lib)
print(res.rbp_ranking.head(10))
```

---

## Katana Standalone Package Suite
| Package | Scope | Repository |
| :--- | :--- | :--- |
| **`katana-motiftrace`** | Bulk RBP regulator discovery & mature-mRNA motif tracing | [Taratype0/katana-motiftrace](https://github.com/Taratype0/katana-motiftrace) |
| **`katana-cosplice`** | Bulk count-debiased ExcessBioVar + Signed TOM co-splicing networks | [Taratype0/katana-cosplice](https://github.com/Taratype0/katana-cosplice) |
| **`katana-sc-motiftrace`** | Smart-seq3 single-cell RBP discovery, motif tracing & PSI modeling | [Taratype0/katana-sc-motiftrace](https://github.com/Taratype0/katana-sc-motiftrace) |
| **`katana-sc-cosplice`** | Smart-seq3 single-cell Hybrid-Manifold co-splicing module discovery | [Taratype0/katana-sc-cosplice](https://github.com/Taratype0/katana-sc-cosplice) |
