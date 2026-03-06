from __future__ import annotations
import numpy as np
from dataclasses import dataclass, field
from typing import Optional
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

    s1 = np.array(list(seq1))
    s2 = np.array(list(seq2))
    sub = np.where(s1[:, None] == s2[None, :], scheme.match, scheme.mismatch).astype(np.int32)

    dp  = np.zeros((m + 1, n + 1), dtype=np.int32)
    ptr = np.zeros((m + 1, n + 1), dtype=np.uint8)

    for d in range(m + n - 1):
        r_start = max(0, d - n + 1)
        r_end   = min(m - 1, d) + 1
        rows = np.arange(r_start, r_end, dtype=np.int32)
        cols = d - rows

        diag = dp[rows,     cols    ] + sub[rows - 1, cols - 1]
        up   = dp[rows,     cols + 1] + GAP
        left = dp[rows + 1, cols    ] + GAP

        best = np.maximum(np.maximum(np.maximum(diag, up), left), 0)
        dp[rows + 1, cols + 1] = best

        p = np.where(best == 0, _STOP,
            np.where(best == diag, _DIAG,
            np.where(best == up,   _UP,   _LEFT))).astype(np.uint8)
        ptr[rows + 1, cols + 1] = p

    max_score = dp.max()
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
        score=int(max_score),
        identity=identity,
        algorithm="Smith-Waterman (local)",
        elapsed_ms=elapsed,
        start1=i,
        start2=j,
        matrix=dp if store_matrix else None,
    )


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

    _banner("3. BLAST-like Heuristic Search")
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
