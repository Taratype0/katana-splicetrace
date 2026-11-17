"""Automated GTF-driven ORF, Reading-Frame, NMD (50-nt EJC rule), and Protein Domain Consequence Tracer.

Zero custom user configuration required:
  - Given standard reference `.gtf` / `.gtf.gz` + indexed genome `.fa`, automatically maps each
    alternative splicing event (`SE`, `A3SS`, `A5SS`, `MXE`, `RI`) onto its canonical protein-coding
    host transcript, splices full-length Inclusion vs Exclusion mRNAs, translates both ORFs, evaluates
    the Nagy & Maquat (50-nt upstream of last exon-exon junction) NMD mRNA decay switch, and scans
    for disrupted conserved protein domains & intrinsic structural elements (TM helix, NLS, IDR, PTM).
"""
from __future__ import annotations

from dataclasses import dataclass
import gzip
from pathlib import Path
import re
import numpy as np
import pandas as pd

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

KD_HYDROPATHY = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9, "A": 1.8, "G": -0.4,
    "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3, "P": -1.6, "H": -3.2, "E": -3.5,
    "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
}

# Conserved PROSITE / Pfam-core catalytic & binding signatures (100% species-agnostic)
CONSERVED_DOMAIN_SIGNATURES: list[tuple[str, re.Pattern[str]]] = [
    ("P-loop_NTPase_WalkerA", re.compile(r"[AG].{4}GK[ST]")),
    ("Protein_Kinase_Catalytic_Core", re.compile(r"[LIVMFYC].[HY].D[LIVMFY]K.{2}N[LIVMFYCT]{3}")),
    ("DEAD_DEAH_Box_RNA_Helicase", re.compile(r"D[EA]A[DH]")),
    ("C2H2_Zinc_Finger", re.compile(r"C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H")),
    ("RNA_Recognition_Motif_RNP1", re.compile(r"[RK]G[FY][GA][FY][ILV].?[FY]")),
    ("EF_Hand_Calcium_Binding", re.compile(r"D.{1}[DNS].{1}[ILVFYW].{1}[DENSTG][DNQGHRK].{4}[LIVMC]")),
    ("WD40_Beta_Propeller_Repeat", re.compile(r"[LIVMSTAC][FY].{1}GH.{3}[LIVMFY].{2}[ST][GSA].{3}D")),
    ("SH3_Proline_Recognition", re.compile(r"W.{10,22}G[WFY].{2}[NSA]Y")),
    ("WW_Interaction_Domain", re.compile(r"W.{20,23}W.{2}P")),
    ("SNARE_Vesicle_Fusion_Motif", re.compile(r"[LVIM].{6}[LVIM].{6}[QE].{6}[LVIM].{6}[LVIM]")),
    ("RS_Arginine_Serine_Splicing_Domain", re.compile(r"(?:RS|SR){4,}")),
    ("RGG_RG_RNA_Phase_Separation_Box", re.compile(r"(?:RGG|RG.{0,2}RGG){2,}")),
]

_ATTR_GENE_ID = re.compile(r'gene_id\s+"?([^";\s]+)')
_ATTR_GENE_NAME = re.compile(r'gene_name\s+"?([^";\s]+)')
_ATTR_TX_ID = re.compile(r'transcript_id\s+"?([^";\s]+)')
_ATTR_APPRIS = re.compile(r'tag\s+"?(appris_principal[^";\s]*|Ensembl_canonical|MANE_Select)')


@dataclass
class HostTranscriptModel:
    transcript_id: str
    gene_id: str
    gene_name: str
    chrom: str
    strand: str
    exons: list[tuple[int, int]]       # sorted ascending 0-based half-open genomic intervals
    cds_blocks: list[tuple[int, int]]  # sorted ascending 0-based half-open genomic intervals
    is_canonical: bool


def load_gtf_coding_transcripts(
    gtf_path: Path | str,
    target_genes: set[str] | None = None,
) -> dict[str, list[HostTranscriptModel]]:
    """Fast single-pass streaming GTF/GTF.gz parser indexed by clean gene_id AND lowercase gene_symbol."""
    p = Path(gtf_path)
    opener = gzip.open if p.name.endswith(".gz") else open

    clean_targets: set[str] | None = None
    byte_tokens: list[bytes] | None = None
    if target_genes:
        clean_targets = {g.split(".")[0].strip() for g in target_genes if g} | {g.lower().strip() for g in target_genes if g}

    tx_meta: dict[str, tuple[str, str, str, str, bool]] = {}
    tx_exons: dict[str, list[tuple[int, int]]] = {}
    tx_cds: dict[str, list[tuple[int, int]]] = {}

    with opener(p, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line[0] == "#":
                continue
            parts = line.split("\t", 8)
            if len(parts) < 9:
                continue
            ftype = parts[2]
            if ftype not in ("exon", "CDS"):
                continue
            attrs = parts[8]
            mg = _ATTR_GENE_ID.search(attrs)
            if not mg:
                continue
            gid = mg.group(1).split(".")[0]
            mn = _ATTR_GENE_NAME.search(attrs)
            gname = mn.group(1) if mn else gid

            if clean_targets is not None and (gid not in clean_targets and gname.lower() not in clean_targets):
                continue

            mt = _ATTR_TX_ID.search(attrs)
            if not mt:
                continue
            tid = mt.group(1).split(".")[0]
            start0 = int(parts[3]) - 1
            end0 = int(parts[4])
            if end0 <= start0:
                continue
            is_can = bool(_ATTR_APPRIS.search(attrs))
            if tid not in tx_meta:
                tx_meta[tid] = (gid, gname, parts[0], parts[6], is_can)
            elif is_can and not tx_meta[tid][4]:
                tx_meta[tid] = (gid, gname, parts[0], parts[6], True)

            if ftype == "exon":
                tx_exons.setdefault(tid, []).append((start0, end0))
            else:
                tx_cds.setdefault(tid, []).append((start0, end0))

    by_gene: dict[str, list[HostTranscriptModel]] = {}
    for tid, (gid, gname, chrom, strand, is_can) in tx_meta.items():
        ex = sorted(set(tx_exons.get(tid, ())))
        cds = sorted(set(tx_cds.get(tid, ())))
        if not ex:
            continue
        model = HostTranscriptModel(
            transcript_id=tid,
            gene_id=gid,
            gene_name=gname,
            chrom=chrom,
            strand=strand,
            exons=ex,
            cds_blocks=cds,
            is_canonical=is_can,
        )
        by_gene.setdefault(gid, []).append(model)
        by_gene.setdefault(gname.lower(), []).append(model)

    # Sort each gene's transcripts: protein-coding first, canonical/APPRIS first, longest CDS, longest mRNA
    for k in by_gene:
        by_gene[k].sort(
            key=lambda m: (
                len(m.cds_blocks) > 0,
                m.is_canonical,
                sum(b - a for a, b in m.cds_blocks),
                sum(b - a for a, b in m.exons),
            ),
            reverse=True,
        )
    return by_gene


def _merge_sorted_blocks(blocks: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not blocks:
        return []
    out: list[tuple[int, int]] = []
    for a, b in sorted(blocks):
        if b <= a:
            continue
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def _build_full_isoform_exons(host_exons: list[tuple[int, int]], row: dict) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    """Splice exact Inclusion vs Exclusion local event coordinates into the full-length host transcript exon chain."""
    kind = str(row["event_type"])
    strand = str(row["strand"])
    if kind == "SE":
        u0, u1 = int(row["upstreamES"]), int(row["upstreamEE"])
        t0, t1 = int(row["exonStart_0base"]), int(row["exonEnd"])
        d0, d1 = int(row["downstreamES"]), int(row["downstreamEE"])
        win_lo, win_hi = u1, d0
        kept = [(a, b) for a, b in host_exons if b <= win_lo or a >= win_hi]
        # Ensure flanking exon edges match u1 and d0
        if not any(b == u1 for _, b in kept):
            kept.append((u0, u1))
        if not any(a == d0 for a, _ in kept):
            kept.append((d0, d1))
        inc = _merge_sorted_blocks(kept + [(t0, t1)])
        exc = _merge_sorted_blocks(kept)
        return inc, exc

    if kind == "MXE":
        u0, u1 = int(row["upstreamES"]), int(row["upstreamEE"])
        a0, a1 = int(row["1stExonStart_0base"]), int(row["1stExonEnd"])
        b0, b1 = int(row["2ndExonStart_0base"]), int(row["2ndExonEnd"])
        d0, d1 = int(row["downstreamES"]), int(row["downstreamEE"])
        kept = [(x, y) for x, y in host_exons if y <= u1 or x >= d0]
        if not any(y == u1 for _, y in kept):
            kept.append((u0, u1))
        if not any(x == d0 for x, _ in kept):
            kept.append((d0, d1))
        inc_ex, exc_ex = ((a0, a1), (b0, b1)) if strand == "+" else ((b0, b1), (a0, a1))
        return _merge_sorted_blocks(kept + [inc_ex]), _merge_sorted_blocks(kept + [exc_ex])

    if kind == "RI":
        u0, u1 = int(row["upstreamES"]), int(row["upstreamEE"])
        d0, d1 = int(row["downstreamES"]), int(row["downstreamEE"])
        kept = [(x, y) for x, y in host_exons if y <= u0 or x >= d1]
        inc = _merge_sorted_blocks(kept + [(u0, d1)])
        exc = _merge_sorted_blocks(kept + [(u0, u1), (d0, d1)])
        return inc, exc

    # A3SS / A5SS
    l0, l1 = int(row["longExonStart_0base"]), int(row["longExonEnd"])
    s0, s1 = int(row["shortES"]), int(row["shortEE"])
    f0, f1 = int(row["flankingES"]), int(row["flankingEE"])
    span_lo, span_hi = min(l0, f0), max(l1, f1)
    kept = [(x, y) for x, y in host_exons if y <= span_lo or x >= span_hi]
    inc = _merge_sorted_blocks(kept + [(l0, l1), (f0, f1)])
    exc = _merge_sorted_blocks(kept + [(s0, s1), (f0, f1)])
    return inc, exc


def _genomic_pos_to_tx_offset(exons_asc: list[tuple[int, int]], strand: str, gpos: int) -> int | None:
    ordered = exons_asc if strand == "+" else exons_asc[::-1]
    cursor = 0
    for a, b in ordered:
        if a <= gpos < b:
            return cursor + (gpos - a if strand == "+" else (b - 1 - gpos))
        cursor += b - a
    return None


def _translate_and_evaluate_nmd(
    dna_5to3: str,
    exon_lengths_5to3: list[int],
    start_tx_pos: int | None,
) -> dict:
    """Translate spliced mRNA from `start_tx_pos` and evaluate Nagy-Maquat 50-nt EJC NMD rule."""
    total_len = len(dna_5to3)
    # Exon-exon junction positions in 5'->3' mature mRNA coordinates
    ejc_positions = list(np.cumsum(exon_lengths_5to3)[:-1]) if len(exon_lengths_5to3) >= 2 else []
    last_ejc = int(ejc_positions[-1]) if ejc_positions else None

    if start_tx_pos is None or start_tx_pos < 0 or start_tx_pos + 3 > total_len:
        # Fallback: find longest canonical ATG ORF if host start codon itself was spliced out
        best_orf = ""
        best_stop_end = None
        for m in re.finditer(r"ATG", dna_5to3):
            st = m.start()
            aa_chars = []
            stop_end = None
            for p in range(st, total_len - 2, 3):
                aa = CODON_TABLE.get(dna_5to3[p : p + 3], "X")
                if aa == "*":
                    stop_end = p + 3
                    break
                aa_chars.append(aa)
            if len(aa_chars) > len(best_orf):
                best_orf = "".join(aa_chars)
                best_stop_end = stop_end
        nmd = bool(best_stop_end is not None and last_ejc is not None and (last_ejc - best_stop_end) >= 50)
        dist_ejc = int(last_ejc - best_stop_end) if (best_stop_end is not None and last_ejc is not None) else 0
        return {
            "protein": best_orf,
            "aa_len": len(best_orf),
            "has_stop": best_stop_end is not None,
            "nmd_sensitive": nmd,
            "ptc_to_last_ejc_nt": dist_ejc,
            "start_preserved": False,
        }

    aa_list: list[str] = []
    stop_end: int | None = None
    for p in range(start_tx_pos, total_len - 2, 3):
        aa = CODON_TABLE.get(dna_5to3[p : p + 3], "X")
        if aa == "*":
            stop_end = p + 3
            break
        aa_list.append(aa)
    prot = "".join(aa_list)
    dist_ejc = int(last_ejc - stop_end) if (stop_end is not None and last_ejc is not None) else 0
    nmd = bool(stop_end is not None and last_ejc is not None and dist_ejc >= 50)
    return {
        "protein": prot,
        "aa_len": len(prot),
        "has_stop": stop_end is not None,
        "nmd_sensitive": nmd,
        "ptc_to_last_ejc_nt": dist_ejc,
        "start_preserved": True,
    }


def _scan_protein_features(aa_seq: str) -> set[str]:
    if not aa_seq or len(aa_seq) < 8:
        return set()
    feats: set[str] = set()
    for name, pat in CONSERVED_DOMAIN_SIGNATURES:
        if pat.search(aa_seq):
            feats.add(name)
    # Transmembrane helix (19-aa window mean Kyte-Doolittle hydropathy >= 1.6)
    if len(aa_seq) >= 19:
        hyd = np.array([KD_HYDROPATHY.get(c, 0.0) for c in aa_seq], dtype=float)
        win = np.convolve(hyd, np.ones(19) / 19.0, mode="valid")
        if float(win.max()) >= 1.6:
            feats.add("Transmembrane_Helix")
    # Classical Monopartite / Bipartite Nuclear Localization Signal (NLS)
    if re.search(r"K[KR].[KR]|[KR]{2}.{10,12}[KR]{3,5}", aa_seq):
        feats.add("Nuclear_Localization_Signal_NLS")
    # Proline/Serine-rich Intrinsically Disordered / Phosphorylation Tract (IDR)
    if re.search(r"[SP]{6,}|[Q]{6,}|[ED]{6,}", aa_seq):
        feats.add("Intrinsically_Disordered_LowComplexity_Tract")
    if len(re.findall(r"[ST]P", aa_seq)) >= 4:
        feats.add("CDK_MAPK_Phospho_Hotspot_(SP/TP)")
    return feats


def trace_protein_and_nmd_consequences(
    events_df: pd.DataFrame,
    gtf_path: Path | str,
    fasta_path: Path | str,
    *,
    indexed_fasta_cls: type,
) -> pd.DataFrame:
    """Automatically annotate splicing events with ORF reading-frame, NMD decay switch, and protein domain impacts."""
    if events_df.empty:
        return pd.DataFrame()

    target_genes = set(events_df["gene_id"].dropna().astype(str)) | set(
        events_df.get("gene_symbol", pd.Series(dtype=str)).dropna().astype(str)
    )
    tx_index = load_gtf_coding_transcripts(gtf_path, target_genes=target_genes)
    reader = indexed_fasta_cls(Path(fasta_path))

    rows: list[dict] = []
    with reader.open_stream() as stream:
        for rec in events_df.to_dict("records"):
            eid = str(rec["event_id"])
            gid = str(rec.get("gene_id", "")).split(".")[0]
            gsym = str(rec.get("gene_symbol", gid))
            dpsi = float(rec.get("dpsi", 0.0))
            etype = str(rec["event_type"])
            strand = str(rec["strand"])

            candidates = tx_index.get(gid) or tx_index.get(gsym.lower()) or []
            # Pick best coding transcript overlapping the event span
            ev_coords = [
                int(rec[c])
                for c in ("exonStart_0base", "exonEnd", "longExonStart_0base", "longExonEnd", "1stExonStart_0base", "2ndExonEnd", "riExonStart_0base", "riExonEnd")
                if c in rec and pd.notna(rec[c])
            ]
            ev_lo, ev_hi = (min(ev_coords), max(ev_coords)) if ev_coords else (0, 0)
            host: HostTranscriptModel | None = None
            for cand in candidates:
                tx_lo, tx_hi = cand.exons[0][0], cand.exons[-1][1]
                if tx_lo <= ev_hi and tx_hi >= ev_lo:
                    host = cand
                    break
            if host is None and candidates:
                host = candidates[0]

            if host is None or not host.cds_blocks:
                rows.append({
                    "event_id": eid,
                    "gene_id": gid,
                    "gene_symbol": gsym,
                    "event_type": etype,
                    "dpsi": dpsi,
                    "host_transcript_id": host.transcript_id if host else "unmatched",
                    "coding_region_location": "noncoding_host",
                    "frame_status": "noncoding_transcript",
                    "nmd_switch_class": "noncoding",
                    "aa_len_inclusion": 0,
                    "aa_len_exclusion": 0,
                    "delta_aa": 0,
                    "disrupted_protein_features": "none",
                    "functional_consequence_summary": "Non-coding transcript isoform change",
                })
                continue

            try:
                inc_exons, exc_exons = _build_full_isoform_exons(host.exons, rec)
                ord_inc = inc_exons if strand == "+" else inc_exons[::-1]
                ord_exc = exc_exons if strand == "+" else exc_exons[::-1]
                dna_inc = reader.fetch_blocks(stream, host.chrom, strand, ord_inc)
                dna_exc = reader.fetch_blocks(stream, host.chrom, strand, ord_exc)
            except Exception:
                continue

            # Canonical start codon genomic coordinate (5'-most base of CDS)
            start_gpos = host.cds_blocks[0][0] if strand == "+" else (host.cds_blocks[-1][1] - 1)
            cds_lo, cds_hi = host.cds_blocks[0][0], host.cds_blocks[-1][1]
            overlaps_cds = bool(ev_lo < cds_hi and ev_hi > cds_lo)
            if not overlaps_cds:
                is_5utr = (ev_hi <= cds_lo) if strand == "+" else (ev_lo >= cds_hi)
                loc = "5prime_UTR" if is_5utr else "3prime_UTR"
            else:
                loc = "CDS_coding"

            st_inc = _genomic_pos_to_tx_offset(inc_exons, strand, start_gpos)
            st_exc = _genomic_pos_to_tx_offset(exc_exons, strand, start_gpos)

            t_inc = _translate_and_evaluate_nmd(dna_inc, [b - a for a, b in ord_inc], st_inc)
            t_exc = _translate_and_evaluate_nmd(dna_exc, [b - a for a, b in ord_exc], st_exc)

            delta_nt = len(dna_inc) - len(dna_exc)
            if t_inc["protein"] == t_exc["protein"]:
                f_status = "utr_cis_regulatory_only"
            elif not (t_inc["start_preserved"] and t_exc["start_preserved"]):
                f_status = "translation_start_altered"
            elif abs(delta_nt) % 3 == 0:
                f_status = "in_frame_indel"
            else:
                f_status = "frameshift"

            if t_inc["nmd_sensitive"] and not t_exc["nmd_sensitive"]:
                nmd_cls = "inclusion_triggers_NMD_poison_exon"
            elif t_exc["nmd_sensitive"] and not t_inc["nmd_sensitive"]:
                nmd_cls = "skipping_triggers_NMD_decay"
            elif not t_inc["nmd_sensitive"] and not t_exc["nmd_sensitive"]:
                nmd_cls = "productive_protein_isoform_switch" if t_inc["protein"] != t_exc["protein"] else "coding_silent_utr_switch"
            else:
                nmd_cls = "dual_NMD_background"

            feats_inc = _scan_protein_features(t_inc["protein"])
            feats_exc = _scan_protein_features(t_exc["protein"])
            diff_feats = sorted(feats_inc ^ feats_exc)
            # Also scan the differential peptide window directly for short local signals (NLS / phospho / IDR)
            if not diff_feats and t_inc["protein"] != t_exc["protein"]:
                p1, p2 = t_inc["protein"], t_exc["protein"]
                pre = 0
                while pre < min(len(p1), len(p2)) and p1[pre] == p2[pre]:
                    pre += 1
                suf = 0
                while suf < (min(len(p1), len(p2)) - pre) and p1[len(p1) - 1 - suf] == p2[len(p2) - 1 - suf]:
                    suf += 1
                seg = p1[max(0, pre - 6) : len(p1) - suf + 6] + p2[max(0, pre - 6) : len(p2) - suf + 6]
                diff_feats = sorted(_scan_protein_features(seg))

            feat_str = ";".join(diff_feats) if diff_feats else "none"
            delta_aa = int(t_inc["aa_len"] - t_exc["aa_len"])

            # Concise, high-impact mechanistic summary
            if nmd_cls == "inclusion_triggers_NMD_poison_exon":
                fav = "Promotes NMD decay" if dpsi > 0 else "Rescues productive protein from NMD"
                summary = f"Poison exon ({f_status}): inclusion introduces PTC ({t_inc['ptc_to_last_ejc_nt']} nt upstream of last EJC) -> {fav}"
            elif nmd_cls == "skipping_triggers_NMD_decay":
                fav = "Triggers NMD mRNA decay" if dpsi < 0 else "Maintains full-length ORF"
                summary = f"Essential coding exon ({f_status}): skipping introduces PTC ({t_exc['ptc_to_last_ejc_nt']} nt upstream of last EJC) -> {fav}"
            elif f_status == "utr_cis_regulatory_only":
                summary = f"ORF intact ({t_inc['aa_len']} aa); remodels {loc} ({delta_nt:+d} nt)"
            else:
                dom_note = f"; alters {feat_str}" if diff_feats else ""
                summary = f"Productive {f_status}: {t_inc['aa_len']} vs {t_exc['aa_len']} aa ({delta_aa:+d} aa){dom_note}"

            rows.append({
                "event_id": eid,
                "gene_id": gid,
                "gene_symbol": gsym,
                "event_type": etype,
                "dpsi": dpsi,
                "host_transcript_id": host.transcript_id,
                "coding_region_location": loc,
                "frame_status": f_status,
                "nmd_switch_class": nmd_cls,
                "aa_len_inclusion": int(t_inc["aa_len"]),
                "aa_len_exclusion": int(t_exc["aa_len"]),
                "delta_aa": delta_aa,
                "disrupted_protein_features": feat_str,
                "functional_consequence_summary": summary,
            })
    return pd.DataFrame(rows)
