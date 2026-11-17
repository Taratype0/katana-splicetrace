"""Unit tests for katana-splicetrace v2.4.1 (Numba exactness, GTF phase/start_codon, UniProt 194k catalog, & 5 plots)."""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import katana_splicetrace as kst
from katana_splicetrace.cli import main as cli_main


def test_v241_library_uniprot_and_cli() -> None:
    assert kst.__version__ == "2.4.1"
    lib_m = kst.load_bundled_human_rbp_library(species="mouse")
    assert lib_m.n_rbps >= 165
    assert "Ythdc1" in lib_m.rbp_to_models
    dom_m = kst.load_bundled_uniprot_domains(species="mouse")
    assert len(dom_m) >= 14000
    assert "macf1" in dom_m and dom_m["macf1"]["uniprot_len"] > 5000
    seqs = ["ACGUUGGACUUGCAUGACGUACGU", "UUUUCCCCUUCUUUGGACUACGU", "GCGCGCGCAAAAUGCAUGUUUU"]
    sub_models = {k: lib_m.models[k] for k in list(lib_m.models)[:25]}
    scores = kst.fast_score_mean_log_odds(seqs, sub_models)
    assert scores.shape == (3, 25)
    assert cli_main(["--list-rbps", "--filter-rbp", "YTHDC1"]) == 0


def test_gtf_phase_offset_and_nmd_poison_exon_and_viz(tmp_path: Path) -> None:
    # Construct gene where CDS starts with phase=2 (2 leading overhang bases "TT" before "ATG" frame!)
    ex1 = "TT" + "ATG" + "GCT" * 20                # 65 nt (phase 2 -> skips "TT", translates 21 aa in-frame!)
    alt_poison = "GCT" * 3 + "TAA" + "GCT" * 12 + "AA"  # 50 nt poison exon
    ex3 = "GCT" * 30                               # 90 nt (>50 nt EJC rule)
    ex4 = "GCT" * 20 + "TGA"                       # 63 nt terminal exon
    chrom_seq = "N" * 100 + ex1 + "GTA" + "T" * 97 + alt_poison + "GTA" + "T" * 97 + ex3 + "GTA" + "T" * 97 + ex4 + "N" * 100

    fa = tmp_path / "syn.fa"
    fa.write_text(f">chr1\n{chrom_seq}\n", encoding="utf-8", newline="\n")
    (tmp_path / "syn.fa.fai").write_text(f"chr1\t{len(chrom_seq)}\t6\t{len(chrom_seq)}\t{len(chrom_seq)+1}\n", encoding="utf-8", newline="\n")

    gtf = tmp_path / "syn.gtf"
    gtf_lines = [
        'chr1\ttest\texon\t101\t165\t.\t+\t.\tgene_id "ENSG_SYN1"; transcript_id "ENST_SYN1"; gene_name "Rbfox2"; transcript_type "protein_coding"; tag "Ensembl_canonical";',
        'chr1\ttest\tCDS\t101\t165\t.\t+\t2\tgene_id "ENSG_SYN1"; transcript_id "ENST_SYN1"; gene_name "Rbfox2"; transcript_type "protein_coding";',
        'chr1\ttest\texon\t416\t505\t.\t+\t.\tgene_id "ENSG_SYN1"; transcript_id "ENST_SYN1"; gene_name "Rbfox2"; transcript_type "protein_coding";',
        'chr1\ttest\tCDS\t416\t505\t.\t+\t0\tgene_id "ENSG_SYN1"; transcript_id "ENST_SYN1"; gene_name "Rbfox2"; transcript_type "protein_coding";',
        'chr1\ttest\texon\t606\t668\t.\t+\t.\tgene_id "ENSG_SYN1"; transcript_id "ENST_SYN1"; gene_name "Rbfox2"; transcript_type "protein_coding";',
        'chr1\ttest\tCDS\t606\t665\t.\t+\t0\tgene_id "ENSG_SYN1"; transcript_id "ENST_SYN1"; gene_name "Rbfox2"; transcript_type "protein_coding";',
    ]
    gtf.write_text("\n".join(gtf_lines) + "\n", encoding="utf-8")

    ev = pd.DataFrame([{
        "event_id": "SE:chr1:+:265:315:100:165:415:505",
        "gene_id": "ENSG_SYN1",
        "gene_symbol": "Rbfox2",
        "event_type": "SE",
        "chr": "chr1",
        "strand": "+",
        "exonStart_0base": 265,
        "exonEnd": 315,
        "upstreamES": 100,
        "upstreamEE": 165,
        "downstreamES": 415,
        "downstreamEE": 505,
        "dpsi": 0.35,
    }])
    res = kst.trace_protein_and_nmd_consequences(ev, gtf, fa, indexed_fasta_cls=kst.IndexedFasta, species="mouse")
    assert len(res) == 1
    row = res.iloc[0]
    assert row["frame_status"] == "frameshift"
    assert row["nmd_switch_class"] == "inclusion_triggers_NMD_poison_exon"
    assert row["aa_len_exclusion"] == 71
    assert row["aa_len_inclusion"] == 24

    png1, _ = kst.plot_consequence_landscape_v24({"c1": res}, "c1", tmp_path / "f1")
    png2, _ = kst.plot_dual_track_nmd_architecture_v24(res, {}, tmp_path / "f2", "c1")
    assert png1.is_file() and png2.is_file()
