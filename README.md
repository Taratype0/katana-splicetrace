# katana-splicetrace (`v2.4.3`) — KatanaAS Suite

[![CI](https://github.com/Taratype0/katana-splicetrace/actions/workflows/ci.yml/badge.svg)](https://github.com/Taratype0/katana-splicetrace/actions/workflows/ci.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

Part of the **KatanaAS** 4-package Bulk & Single-Cell (Smart-seq3) alternative splicing suite. Supports **two first-class invocation modes**:
1. **Standalone Python Package & CLI** (`pip install` -> Python API / `katana-splicetrace` / `splicetrace` / `python -m katana_splicetrace`)
2. **Nextflow `nf-core` DSL2 Module & Subworkflow** (`process KATANA_SPLICETRACE` with native `--emit-versions-yml`)

---

## Mode 1: Standalone Package Installation & Invocation

### Installation
```bash
# Direct one-line wheel install from GitHub Release
pip install https://github.com/Taratype0/katana-splicetrace/releases/download/v2.4.3/katana_splicetrace-2.4.3-py3-none-any.whl

# Or install from source with server Numba acceleration & Matplotlib visualization extras
pip install "git+https://github.com/Taratype0/katana-splicetrace.git#egg=katana-splicetrace[server,viz]"
```

### CLI Usage (`katana-splicetrace`, `splicetrace`, or `python -m katana_splicetrace`)
```bash
katana-splicetrace \
  --rmats ./rmats_out \
  --fasta GRCm39.primary_assembly.genome.fa \
  --gtf gencode.vM34.annotation.gtf \
  --splicenet-modules ./splicenet_out/module_assignments_v2.tsv \
  --species mouse --threads 16 \
  --output-dir ./splicetrace_out
```

### Python API Usage
```python
import katana_splicetrace as kst

kst.configure_server_threads(16)
lib = kst.load_bundled_human_rbp_library(species='mouse')
events = kst.load_rmats_measurable_events('./rmats_out', 'GRCm39.fa')
res = kst.run_bulk_continuous_rescue_v22(events, 'GRCm39.fa', lib, threads=16)
prot_df = kst.trace_protein_and_nmd_consequences(
    res.event_contrasts[res.event_contrasts['_is_fg']],
    'gencode.vM34.annotation.gtf', 'GRCm39.fa',
    indexed_fasta_cls=kst.IndexedFasta, species='mouse'
)
```

---

## Mode 2: Nextflow `nf-core` DSL2 Integration

Every KatanaAS CLI natively implements `--emit-versions-yml versions.yml --process-name "${task.process}"` and SLURM multi-core `--threads ${task.cpus}`.

```groovy
include { KATANA_SPLICETRACE } from './modules/local/katana_splicetrace/main'
```

---

## KatanaAS 4-Package Suite Architecture

| Package | Version | Modality | Responsibility |
| :--- | :---: | :---: | :--- |
| **[katana-splicetrace](https://github.com/Taratype0/katana-splicetrace)** | `v2.4.3` | Bulk RNA-seq | RBP discovery (`169` RBPs / `621` PWMs), mature-mRNA motifs, `SpliceNet` bridge & `GTF/UniProt` ORF/NMD + 5 plots |
| **[katana-splicenet](https://github.com/Taratype0/katana-splicenet)** | `v2.3.2` | Bulk RNA-seq | Pure-Python `ExcessBioVar` + Adaptive Scale-Free Signed TOM co-splicing network & GO enrichment |
| **[katana-sc-splicetrace](https://github.com/Taratype0/katana-sc-splicetrace)** | `v2.4.1` | Smart-seq3 SC | Single-cell RBP discovery, mature-mRNA motifs, `SC-SpliceNet` bridge & `GTF/UniProt` ORF/NMD + 5 plots |
| **[katana-sc-splicenet](https://github.com/Taratype0/katana-sc-splicenet)** | `v2.2.2` | Smart-seq3 SC | Hybrid-Manifold `ExcessBioVar` + Empirical-Bayes Fisher-z shrunk Signed TOM single-cell modules |
