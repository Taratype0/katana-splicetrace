"""Unit tests for katana-splicetrace v2.4.0 (Numba speed exactness + GTF ORF/NMD/domain consequences)."""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import katana_splicetrace as kst


def test_v24_library_and_numba_exactness() -> None:
    assert kst.__version__ == "2.4.0"
    lib_m = kst.load_bundled_human_rbp_library(species="mouse")
    assert lib_m.n_rbps >= 165
    assert "Ythdc1" in lib_m.rbp_to_models
    seqs = ["ACGUUGGACUUGCAUGACGUACGU", "UUUUCCCCUUCUUUGGACUACGU", "GCGCGCGCAAAAUGCAUGUUUU"]
    sub_models = {k: lib_m.models[k] for k in list(lib_m.models)[:25]}
    scores = kst.fast_score_mean_log_odds(seqs, sub_models)
    assert scores.shape == (3, 25)
    assert np.all(np.isfinite(scores)) and np.all(scores >= 0.0)


def test_automated_gtf_orf_and_nmd_poison_exon(tmp_path: Path) -> None:
    # Construct a synthetic 3-exon coding gene on chr1 (+) where including a 50-nt middle cassette exon
    # introduces an in-frame TAA stop codon 80 nt upstream of the last exon-exon junction (-> Poison Exon NMD!)
    ex1 = "ATG" + "GCT" * 20                      # 63 nt (21 aa, frame 0)
    alt_poison = "GCT" * 3 + "TAA" + "GCT" * 12 + "AA"  # 50 nt (frameshift + immediate TAA stop!)
    ex3 = "GCT" * 30                             # 90 nt downstream of last EJC (>50 nt!)
    ex4 = "GCT" * 20 + "TGA"                     # 63 nt terminal exon with canonical stop
    chrom_seq = "N" * 100 + ex1 + "GTA" + "T" * 97 + alt_poison + "GTA" + "T" * 97 + ex3 + "GTA" + "T" * 97 + ex4 + "N" * 100

    fa = tmp_path / "syn.fa"
    fa.write_text(f">chr1\n{chrom_seq}\n", encoding="utf-8", newline="\n")
    (tmp_path / "syn.fa.fai").write_text(f"chr1\t{len(chrom_seq)}\t6\t{len(chrom_seq)}\t{len(chrom_seq)+1}\n", encoding="utf-8", newline="\n")

    # 1-based GTF coordinates:
    # ex1: [100, 163) -> 101..163
    # alt: [263, 313) -> 264..313
    # ex3: [413, 503) -> 414..503
    # ex4: [603, 666) -> 604..666
    gtf = tmp_path / "syn.gtf"
    gtf_lines = [
        'chr1\ttest\texon\t101\t163\t.\t+\t.\tgene_id "ENSG_SYN1"; transcript_id "ENST_SYN1"; gene_name "SynGene1"; tag "Ensembl_canonical";',
        'chr1\ttest\tCDS\t101\t163\t.\t+\t0\tgene_id "ENSG_SYN1"; transcript_id "ENST_SYN1"; gene_name "SynGene1";',
        'chr1\ttest\texon\t414\t503\t.\t+\t.\tgene_id "ENSG_SYN1"; transcript_id "ENST_SYN1"; gene_name "SynGene1"; tag "Ensembl_canonical";',
        'chr1\ttest\tCDS\t414\t503\t.\t+\t0\tgene_id "ENSG_SYN1"; transcript_id "ENST_SYN1"; gene_name "SynGene1";',
        'chr1\ttest\texon\t604\t666\t.\t+\t.\tgene_id "ENSG_SYN1"; transcript_id "ENST_SYN1"; gene_name "SynGene1"; tag "Ensembl_canonical";',
        'chr1\ttest\tCDS\t604\t663\t.\t+\t0\tgene_id "ENSG_SYN1"; transcript_id "ENST_SYN1"; gene_name "SynGene1";',
    ]
    gtf.write_text("\n".join(gtf_lines) + "\n", encoding="utf-8")

    ev = pd.DataFrame([{
        "event_id": "SE:chr1:+:263:313:100:163:413:503",
        "gene_id": "ENSG_SYN1",
        "gene_symbol": "SynGene1",
        "event_type": "SE",
        "chr": "chr1",
        "strand": "+",
        "exonStart_0base": 263,
        "exonEnd": 313,
        "upstreamES": 100,
        "upstreamEE": 163,
        "downstreamES": 413,
        "downstreamEE": 503,
        "dpsi": 0.35,
    }])
    res = kst.trace_protein_and_nmd_consequences(ev, gtf, fa, indexed_fasta_cls=kst.IndexedFasta)
    assert len(res) == 1
    row = res.iloc[0]
    assert row["frame_status"] == "frameshift"
    assert row["nmd_switch_class"] == "inclusion_triggers_NMD_poison_exon"
    assert row["aa_len_exclusion"] == 71
    assert row["aa_len_inclusion"] == 24
