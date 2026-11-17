"""Katana-SpliceTrace 2.4.0 — High-speed bulk AS RBP discovery, mature-mRNA motif & ORF/NMD/domain consequence tracing."""
from katana_splicetrace.continuous_rescue import MotifRescueResultV22, run_bulk_continuous_rescue_v22
from katana_splicetrace.extensible_db import (
    export_katana_motif_pack,
    iupac_to_pwm,
    load_external_motif_entries,
    merge_motif_entries_into_library,
)
from katana_splicetrace.fast_accel import HAS_NUMBA, configure_server_threads, fast_score_mean_log_odds
from katana_splicetrace.library import (
    CuratedMotifLibrary,
    load_bundled_human_rbp_library,
    load_curated_human_rbp_library,
)
from katana_splicetrace.matching import match_events_v2
from katana_splicetrace.protein_consequences import load_gtf_coding_transcripts, trace_protein_and_nmd_consequences
from katana_splicetrace.rmats_io import load_rmats_measurable_events
from katana_splicetrace.rna_forms import IndexedFasta, extract_event_rna_v2, trace_motif_instances_v2
from katana_splicetrace.scanner import run_motiftrace_v2_scan

__version__ = "2.4.0"

__all__ = [
    "__version__",
    "HAS_NUMBA",
    "configure_server_threads",
    "fast_score_mean_log_odds",
    "CuratedMotifLibrary",
    "MotifRescueResultV22",
    "IndexedFasta",
    "load_bundled_human_rbp_library",
    "load_curated_human_rbp_library",
    "load_external_motif_entries",
    "merge_motif_entries_into_library",
    "export_katana_motif_pack",
    "iupac_to_pwm",
    "load_rmats_measurable_events",
    "match_events_v2",
    "extract_event_rna_v2",
    "trace_motif_instances_v2",
    "load_gtf_coding_transcripts",
    "trace_protein_and_nmd_consequences",
    "run_motiftrace_v2_scan",
    "run_bulk_continuous_rescue_v22",
]
