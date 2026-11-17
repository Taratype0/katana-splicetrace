from __future__ import annotations
import katana_motiftrace as pkg

def test_v3_library_and_ythdc1() -> None:
    assert pkg.__version__ == "2.3.0"
    lib_h = pkg.load_bundled_human_rbp_library(species="human")
    assert lib_h.n_rbps >= 165
    assert "YTHDC1" in lib_h.rbp_to_models
    assert "RBM20" in lib_h.rbp_to_models
    assert "SRRM4" in lib_h.rbp_to_models
    lib_m = pkg.load_bundled_human_rbp_library(species="mouse")
    assert "Ythdc1" in lib_m.rbp_to_models
