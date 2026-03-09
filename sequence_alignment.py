from __future__ import annotations
import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List
import time

@dataclass
class ScoringScheme:
    match:    int = 2
    mismatch: int = -1
    gap_open: int = -2
    gap_extend: int = 0

    def score(self, a: str, b: str) -> int:
        return self.match if a == b else self.mismatch


@dataclass
class AlignmentResult:
    seq1_aligned:   str
    seq2_aligned:   str
    score:          float
    identity:       float
    algorithm:      str
    elapsed_ms:     float
    start1:         int = 0
    start2:         int = 0
    matrix:         Optional[np.ndarray] = field(default=None, repr=False)

    def __str__(self) -> str:
        mid = "".join("|" if a == b and a != "-" else " "
                      for a, b in zip(self.seq1_aligned, self.seq2_aligned))
        return (
            f"Algorithm : {self.algorithm}\n"
            f"Score     : {self.score}\n"
            f"Identity  : {self.identity:.1%}\n"
            f"Time      : {self.elapsed_ms:.2f} ms\n\n"
            f"  {self.seq1_aligned}\n"
            f"  {mid}\n"
            f"  {self.seq2_aligned}\n"
        )


@dataclass
class MSAResult:
    """Result of a Multiple Sequence Alignment via progressive alignment."""
    sequences:      List[str]          # original (unaligned) sequences
    aligned:        List[str]          # aligned sequences (same length, padded with '-')
    labels:         List[str]          # sequence labels (e.g. "Seq1", "Seq2", ...)
    guide_tree:     List[tuple]        # list of (i, j, distance) merge steps
    avg_identity:   float              # mean pairwise identity over the MSA columns
    elapsed_ms:     float
    algorithm:      str = "Progressive MSA (NW guide tree)"


@dataclass
class BlastHit:
    start1: int
    start2: int
    length: int
    score:  float
    seq1_region: str
    seq2_region: str


_DIAG = 1
_UP   = 2
_LEFT = 3
_STOP = 0


def needleman_wunsch(
    seq1: str,
    seq2: str,
    scheme: ScoringScheme = ScoringScheme(),
    store_matrix: bool = False,
) -> AlignmentResult:
    t0 = time.perf_counter()

    m, n = len(seq1), len(seq2)
    GAP = scheme.gap_open

    if m == 0 and n == 0:
        elapsed = (time.perf_counter() - t0) * 1000
        return AlignmentResult(
            seq1_aligned="", seq2_aligned="",
            score=0, identity=0.0,
            algorithm="Needleman-Wunsch (global)",
            elapsed_ms=elapsed, matrix=None,
        )
    if m == 0:
        elapsed = (time.perf_counter() - t0) * 1000
        return AlignmentResult(
            seq1_aligned="-" * n, seq2_aligned=seq2,
            score=GAP * n, identity=0.0,
            algorithm="Needleman-Wunsch (global)",
            elapsed_ms=elapsed, matrix=None,
        )
    if n == 0:
        elapsed = (time.perf_counter() - t0) * 1000
        return AlignmentResult(
            seq1_aligned=seq1, seq2_aligned="-" * m,
            score=GAP * m, identity=0.0,
            algorithm="Needleman-Wunsch (global)",
            elapsed_ms=elapsed, matrix=None,
        )

    s1 = np.frompyfunc(lambda c: c, 1, 1)(np.array(list(seq1)))
    s2 = np.frompyfunc(lambda c: c, 1, 1)(np.array(list(seq2)))
    sub = np.where(
        s1[:, None] == s2[None, :],
        scheme.match,
        scheme.mismatch,
    ).astype(np.int32)

    dp  = np.zeros((m + 1, n + 1), dtype=np.int32)
    ptr = np.zeros((m + 1, n + 1), dtype=np.uint8)

    dp[0, :] = GAP * np.arange(n + 1)
    dp[:, 0] = GAP * np.arange(m + 1)
    ptr[0, 1:] = _LEFT
    ptr[1:, 0] = _UP

    for d in range(m + n - 1):
        r_start = max(0, d - n + 1)
        r_end   = min(m - 1, d) + 1
        rows = np.arange(r_start, r_end, dtype=np.int32)
        cols = d - rows

        diag = dp[rows, cols] + sub[rows - 1, cols - 1]
        up   = dp[rows,     cols + 1] + GAP
        left = dp[rows + 1, cols    ] + GAP

        best = np.maximum(np.maximum(diag, up), left)
        dp[rows + 1, cols + 1] = best

        p = np.where(best == diag, _DIAG,
            np.where(best == up,   _UP,   _LEFT)).astype(np.uint8)
        ptr[rows + 1, cols + 1] = p

    a1, a2 = [], []
    i, j = m, n
    while i > 0 or j > 0:
        p = ptr[i, j]
        if p == _DIAG:
            a1.append(seq1[i - 1]); a2.append(seq2[j - 1]); i -= 1; j -= 1
        elif p == _UP:
            a1.append(seq1[i - 1]); a2.append("-"); i -= 1
        else:
            a1.append("-"); a2.append(seq2[j - 1]); j -= 1

    a1 = "".join(reversed(a1))
    a2 = "".join(reversed(a2))
    matches = sum(x == y and x != "-" for x, y in zip(a1, a2))
    identity = matches / max(len(a1), 1)

    elapsed = (time.perf_counter() - t0) * 1000
    return AlignmentResult(
        seq1_aligned=a1,
        seq2_aligned=a2,
        score=int(dp[m, n]),
        identity=identity,
        algorithm="Needleman-Wunsch (global)",
        elapsed_ms=elapsed,
        matrix=dp if store_matrix else None,
    )


def smith_waterman(
    seq1: str,
    seq2: str,
    scheme: ScoringScheme = ScoringScheme(),
    store_matrix: bool = False,
) -> AlignmentResult:
    t0 = time.perf_counter()

    m, n = len(seq1), len(seq2)
    GAP = scheme.gap_open

    if m == 0 or n == 0:
        elapsed = (time.perf_counter() - t0) * 1000
        return AlignmentResult(
            seq1_aligned="", seq2_aligned="",
            score=0, identity=0.0,
            algorithm="Smith-Waterman (local)",
            elapsed_ms=elapsed, matrix=None,
        )

    dp  = np.zeros((m + 1, n + 1), dtype=np.int32)
    ptr = np.zeros((m + 1, n + 1), dtype=np.uint8)

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            s = scheme.match if seq1[i - 1] == seq2[j - 1] else scheme.mismatch
            diag = dp[i - 1, j - 1] + s
            up   = dp[i - 1, j    ] + GAP
            left = dp[i,     j - 1] + GAP
            best = max(diag, up, left, 0)
            dp[i, j] = best
            if best == 0:
                ptr[i, j] = _STOP
            elif best == diag:
                ptr[i, j] = _DIAG
            elif best == up:
                ptr[i, j] = _UP
            else:
                ptr[i, j] = _LEFT

    max_score = int(dp.max())
    i, j = map(int, np.unravel_index(dp.argmax(), dp.shape))

    a1, a2 = [], []
    while i > 0 and j > 0 and ptr[i, j] != _STOP:
        p = ptr[i, j]
        if p == _DIAG:
            a1.append(seq1[i - 1]); a2.append(seq2[j - 1]); i -= 1; j -= 1
        elif p == _UP:
            a1.append(seq1[i - 1]); a2.append("-"); i -= 1
        else:
            a1.append("-"); a2.append(seq2[j - 1]); j -= 1

    a1 = "".join(reversed(a1))
    a2 = "".join(reversed(a2))
    matches = sum(x == y and x != "-" for x, y in zip(a1, a2))
    identity = matches / max(len(a1), 1)

    elapsed = (time.perf_counter() - t0) * 1000
    return AlignmentResult(
        seq1_aligned=a1,
        seq2_aligned=a2,
        score=max_score,
        identity=identity,
        algorithm="Smith-Waterman (local)",
        elapsed_ms=elapsed,
        start1=i,
        start2=j,
        matrix=dp if store_matrix else None,
    )


# ── Multiple Sequence Alignment (Progressive / ClustalW-style) ───────────────

def _pairwise_distance(seq1: str, seq2: str, scheme: ScoringScheme) -> float:
    """
    Distance = 1 - identity of a Needleman-Wunsch pairwise alignment.
    Used to build the guide tree.
    """
    res = needleman_wunsch(seq1, seq2, scheme)
    return 1.0 - res.identity


def _build_guide_tree(sequences: List[str], scheme: ScoringScheme) -> List[tuple]:
    """
    UPGMA-style greedy guide tree.
    Returns a list of merge steps: (index_a, index_b, distance).
    The cheapest (most similar) pair is merged first.
    """
    n = len(sequences)
    # distance matrix
    dist = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d = _pairwise_distance(sequences[i], sequences[j], scheme)
            dist[i][j] = dist[j][i] = d

    active = list(range(n))
    steps: List[tuple] = []

    while len(active) > 1:
        best_d = float("inf")
        best_pair = (0, 1)
        for idx_a in range(len(active)):
            for idx_b in range(idx_a + 1, len(active)):
                i, j = active[idx_a], active[idx_b]
                if dist[i][j] < best_d:
                    best_d = dist[i][j]
                    best_pair = (i, j)

        i, j = best_pair
        steps.append((i, j, best_d))
        active.remove(j)   # merge j into i

    return steps


def _align_profiles(profile_a: List[str], profile_b: List[str],
                    scheme: ScoringScheme) -> tuple[List[str], List[str]]:
    """
    Align two groups of already-aligned sequences (profiles) using the
    consensus of each profile as a representative.

    A gap in either profile is propagated to ALL sequences in that profile.
    This is the 'once a gap, always a gap' rule used by ClustalW.
    """
    def consensus(seqs: List[str]) -> str:
        if not seqs:
            return ""
        out = []
        for col in zip(*seqs):
            non_gap = [c for c in col if c != "-"]
            out.append(non_gap[0] if non_gap else "-")
        return "".join(out)

    rep_a = consensus(profile_a)
    rep_b = consensus(profile_b)

    res = needleman_wunsch(rep_a, rep_b, scheme)
    col_a, col_b = res.seq1_aligned, res.seq2_aligned

    def expand_profile(profile: List[str], aligned_rep: str,
                       original_rep: str) -> List[str]:
        """
        Re-insert gaps introduced by NW into every sequence of the profile.
        We walk along `aligned_rep`; whenever we see a '-' that was NOT
        in `original_rep`, we insert a '-' at that position in every seq.
        """
        orig_iter = iter(range(len(original_rep)))
        new_seqs = [""] * len(profile)
        orig_idx = 0

        for ch in aligned_rep:
            if ch == "-":
                # a new gap column — insert '-' everywhere in this profile
                for k in range(len(profile)):
                    new_seqs[k] += "-"
            else:
                # a real character — copy the original column from every seq
                for k in range(len(profile)):
                    new_seqs[k] += profile[k][orig_idx] if orig_idx < len(profile[k]) else "-"
                orig_idx += 1

        return new_seqs

    new_a = expand_profile(profile_a, col_a, rep_a)
    new_b = expand_profile(profile_b, col_b, rep_b)
    return new_a, new_b


def progressive_msa(
    sequences: List[str],
    labels:    Optional[List[str]] = None,
    scheme:    ScoringScheme = ScoringScheme(),
) -> MSAResult:
    """
    Progressive Multiple Sequence Alignment (ClustalW-style).

    Algorithm:
      1. Compute all pairwise NW distances  → O(N² · L²)
      2. Build a greedy UPGMA guide tree    → O(N²)
      3. Merge profiles in guide-tree order → O(N · L²)

    Parameters
    ----------
    sequences : list of raw (unaligned) DNA/protein strings
    labels    : optional sequence names; defaults to ["Seq1", "Seq2", ...]
    scheme    : ScoringScheme shared with pairwise steps

    Returns
    -------
    MSAResult with .aligned containing all sequences padded to the same length.
    """
    t0 = time.perf_counter()

    if not sequences:
        raise ValueError("At least one sequence is required.")
    if len(sequences) == 1:
        elapsed = (time.perf_counter() - t0) * 1000
        lbl = labels or ["Seq1"]
        return MSAResult(
            sequences=sequences,
            aligned=list(sequences),
            labels=lbl,
            guide_tree=[],
            avg_identity=1.0,
            elapsed_ms=elapsed,
        )

    labels = labels or [f"Seq{i + 1}" for i in range(len(sequences))]
    seqs = [s.upper().strip() for s in sequences]

    # Step 1 + 2: guide tree
    guide_tree = _build_guide_tree(seqs, scheme)

    # Step 3: progressive merging
    # Each cluster is a list of (label_index, aligned_string) pairs.
    # We track which original index belongs to which cluster.
    clusters: List[List[int]]          = [[i] for i in range(len(seqs))]
    profiles: List[List[str]]          = [[s] for s in seqs]

    for (i, j, _dist) in guide_tree:
        new_a, new_b = _align_profiles(profiles[i], profiles[j], scheme)
        # merge j into i
        profiles[i] = new_a + new_b
        clusters[i] = clusters[i] + clusters[j]
        profiles[j] = []
        clusters[j] = []

    # The final aligned profile is whichever cluster survived
    final_profile = next(p for p in profiles if p)
    final_cluster = next(c for c in clusters if c)

    # Re-order back to original sequence order
    order = {orig_idx: pos for pos, orig_idx in enumerate(final_cluster)}
    aligned = [""] * len(seqs)
    for pos, orig_idx in enumerate(final_cluster):
        aligned[orig_idx] = final_profile[pos]

    # Pad to equal length (should already be equal, but defensive)
    max_len = max(len(a) for a in aligned)
    aligned = [a.ljust(max_len, "-") for a in aligned]

    # Compute average pairwise identity over MSA columns
    n = len(aligned)
    total_identity = 0.0
    pairs = 0
    for a in range(n):
        for b in range(a + 1, n):
            matches = sum(x == y and x != "-"
                          for x, y in zip(aligned[a], aligned[b]))
            col_len = max(len(aligned[a]), 1)
            total_identity += matches / col_len
            pairs += 1
    avg_identity = total_identity / max(pairs, 1)

    elapsed = (time.perf_counter() - t0) * 1000
    return MSAResult(
        sequences=sequences,
        aligned=aligned,
        labels=labels,
        guide_tree=guide_tree,
        avg_identity=avg_identity,
        elapsed_ms=elapsed,
    )


# ── BlastLike ─────────────────────────────────────────────────────────────────

class BlastLike:
    def __init__(self, k: int = 6, drop: int = 16, band: int = 30, scheme: ScoringScheme = ScoringScheme()):
        self.k      = k
        self.drop   = drop
        self.band   = band
        self.scheme = scheme
        self._db: str = ""
        self._index: dict[str, list[int]] = {}

    def build_index(self, db_sequence: str) -> None:
        self._db    = db_sequence
        self._index = {}
        k = self.k
        for i in range(len(db_sequence) - k + 1):
            kmer = db_sequence[i:i + k]
            self._index.setdefault(kmer, []).append(i)

    def _find_seeds(self, query: str) -> dict[int, list[tuple[int, int]]]:
        k = self.k
        diags: dict[int, list[tuple[int, int]]] = {}
        for q in range(len(query) - k + 1):
            kmer = query[q:q + k]
            for d in self._index.get(kmer, []):
                diag = d - q
                diags.setdefault(diag, []).append((q, d))
        return diags

    def _merge_seeds(self, seeds: list[tuple[int, int]]) -> list[tuple[int, int, int, int]]:
        seeds_sorted = sorted(seeds, key=lambda s: s[0])
        merged = []
        if not seeds_sorted:
            return merged

        qs, ds = seeds_sorted[0]
        qe, de = qs + self.k, ds + self.k

        for q, d in seeds_sorted[1:]:
            if q <= qe + self.drop:
                qe = max(qe, q + self.k)
                de = max(de, d + self.k)
            else:
                merged.append((qs, ds, qe, de))
                qs, ds, qe, de = q, d, q + self.k, d + self.k

        merged.append((qs, ds, qe, de))
        return merged

    def _extend(self, query: str, q_s: int, q_e: int, db_s: int, db_e: int) -> BlastHit:
        b = self.band
        q_lo  = max(0, q_s  - b);  q_hi  = min(len(query),    q_e  + b)
        db_lo = max(0, db_s - b);  db_hi = min(len(self._db), db_e + b)

        q_sub  = query[q_lo:q_hi]
        db_sub = self._db[db_lo:db_hi]

        res = smith_waterman(q_sub, db_sub, self.scheme)
        return BlastHit(
            start1=q_lo  + res.start1,
            start2=db_lo + res.start2,
            length=len(res.seq1_aligned),
            score=res.score,
            seq1_region=res.seq1_aligned,
            seq2_region=res.seq2_aligned,
        )

    def search(self, query: str, top_n: int = 5) -> list[BlastHit]:
        diags = self._find_seeds(query)
        hits: list[BlastHit] = []

        for diag, seeds in diags.items():
            for qs, ds, qe, de in self._merge_seeds(seeds):
                hit = self._extend(query, qs, qe, ds, de)
                if hit.score > 0:
                    hits.append(hit)

        hits.sort(key=lambda h: h.score, reverse=True)
        kept: list[BlastHit] = []
        used_db: list[tuple[int, int]] = []
        for h in hits:
            overlap = any(
                not (h.start2 + h.length <= s or h.start2 >= e)
                for s, e in used_db
            )
            if not overlap:
                kept.append(h)
                used_db.append((h.start2, h.start2 + h.length))
            if len(kept) >= top_n:
                break

        return kept


def _banner(title: str) -> None:
    print("\n" + "═" * 60)
    print(f"  {title}")
    print("═" * 60)


def demo() -> None:
    scheme = ScoringScheme(match=2, mismatch=-1, gap_open=-2)

    _banner("1. Needleman-Wunsch — Global Alignment")
    seq1 = "ACGTTGCATGCA"
    seq2 = "ACGTCATGCA"
    res  = needleman_wunsch(seq1, seq2, scheme)
    print(res)

    _banner("2. Smith-Waterman — Local Alignment")
    seq1 = "GGTTGACTAGCGGATCGATCGT"
    seq2 = "TTAATTGACAGCGGAA"
    res  = smith_waterman(seq1, seq2, scheme)
    print(res)

    _banner("3. Progressive MSA (3 sequences)")
    msa = progressive_msa(
        ["ACGTTGCATGCA", "ACGTCATGCA", "ACGCATGCA"],
        labels=["Seq1", "Seq2", "Seq3"],
        scheme=scheme,
    )
    print(f"Avg identity: {msa.avg_identity:.1%}  |  Time: {msa.elapsed_ms:.2f} ms")
    for lbl, aln in zip(msa.labels, msa.aligned):
        print(f"  {lbl}: {aln}")

    _banner("4. BLAST-like Heuristic Search")
    import random
    random.seed(42)
    genome = "".join(random.choices("ACGT", k=2000))
    query       = "ACGTACGTACGTACGT"
    insert_pos  = 700
    genome      = genome[:insert_pos] + query + genome[insert_pos + len(query):]

    blast = BlastLike(k=6, band=30, scheme=scheme)
    t0    = time.perf_counter()
    blast.build_index(genome)
    hits  = blast.search(query, top_n=3)
    elapsed = (time.perf_counter() - t0) * 1000

    print(f"Query  : {query}")
    print(f"Genome : {len(genome)} bp  |  Time: {elapsed:.2f} ms\n")
    for rank, h in enumerate(hits, 1):
        print(f"Hit #{rank}  score={h.score}  db_pos={h.start2}–{h.start2+h.length}")
        print(f"  Q: {h.seq1_region}")
        print(f"  G: {h.seq2_region}\n")


if __name__ == "__main__":
    demo()
