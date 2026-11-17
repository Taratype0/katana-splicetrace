# katana-splicetrace (`v2.4.0`) — KatanaAS Suite

[![CI](https://github.com/Taratype0/katana-splicetrace/actions/workflows/ci.yml/badge.svg)](https://github.com/Taratype0/katana-splicetrace/actions/workflows/ci.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

KatanaAS high-speed bulk alternative splicing RBP regulator discovery, mature-mRNA motif tracing, and automated GTF ORF/NMD/domain consequence engine.

## Highlights in `v2.4.0`
- **Multi-Core Server Acceleration (`Numba @njit` + `mmap` FASTA + `float32` BLAS)**: Scales across all SLURM / HPC CPU cores (`--threads`), accelerating full 169-RBP / 621-PWM scans by **`15x–40x`**.
- **Automated GTF `ORF / Frameshift / NMD / Protein Domain` Consequence Engine (`--gtf`)**: Reconstructs full-length Inclusion vs Exclusion host mRNAs from standard GENCODE/Ensembl `.gtf[.gz]`, evaluates the Nagy-Maquat `50-nt` EJC **NMD decay / Poison Exon** rule, and detects disrupted conserved catalytic/binding domains, transmembrane helices, NLS & IDRs — zero manual configuration required.
- **Bundled Library `v3.0` (`169` RBPs / `621` PWMs) + 5-Format Extensible Plugin Engine**: Includes `YTHDC1` (`UGGACU`), `SRRM4`, `RBM20`, `RBM10`, `RBM39`, `SON`, `--species human|mouse`, and `.meme` / `.motif` / `.jaspar` / `.tsv` / `.json.gz` custom motif expansion.

## CLI Quickstart
```bash
# Install with server multi-core Numba acceleration
pip install "git+https://github.com/Taratype0/katana-splicetrace.git#egg=katana-splicetrace[server]"

# Run end-to-end RBP discovery + mature-mRNA motif + protein ORF/NMD/domain tracing
katana-splicetrace --rmats ./rmats_out --fasta GRCm39.fa --gtf gencode.vM34.annotation.gtf --species mouse --threads 16 --output-dir ./splicetrace_out
```

---

## KatanaAS Standalone Alternative Splicing Suite
| Package | Scope | Repository |
| :--- | :--- | :--- |
| **`katana-splicetrace`** | Bulk RBP regulator discovery, mature-mRNA motif & ORF/NMD/domain consequence tracing | [Taratype0/katana-splicetrace](https://github.com/Taratype0/katana-splicetrace) |
| **`katana-splicenet`** | Bulk count-debiased `ExcessBioVar` + Signed TOM co-splicing module networks | [Taratype0/katana-splicenet](https://github.com/Taratype0/katana-splicenet) |
| **`katana-sc-splicetrace`** | Smart-seq3 single-cell RBP discovery, motif & ORF/NMD tracing, and PSI modeling | [Taratype0/katana-sc-splicetrace](https://github.com/Taratype0/katana-sc-splicetrace) |
| **`katana-sc-splicenet`** | Smart-seq3 single-cell Hybrid-Manifold co-splicing module discovery | [Taratype0/katana-sc-splicenet](https://github.com/Taratype0/katana-sc-splicenet) |
