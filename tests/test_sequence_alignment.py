import pytest
import numpy as np
from sequence_alignment import (
    ScoringScheme,
    AlignmentResult,
    BlastHit,
    BlastLike,
    needleman_wunsch,
    smith_waterman,
)

SCHEME = ScoringScheme(match=2, mismatch=-1, gap_open=-2)



class TestScoringScheme:
    def test_match_returns_match_score(self):
        assert SCHEME.score("A", "A") == 2

    def test_mismatch_returns_mismatch_score(self):
        assert SCHEME.score("A", "C") == -1

    def test_default_values(self):
        s = ScoringScheme()
        assert s.match == 2
        assert s.mismatch == -1
        assert s.gap_open == -2



class TestNeedlemanWunsch:
    def test_identical_sequences_max_score(self):
        res = needleman_wunsch("ACGT", "ACGT", SCHEME)
        assert res.score == 8          
        assert res.identity == pytest.approx(1.0)
        assert res.seq1_aligned == "ACGT"
        assert res.seq2_aligned == "ACGT"

    def test_completely_different_sequences(self):
        res = needleman_wunsch("AAAA", "CCCC", SCHEME)
        assert res.score == -4          
        assert res.identity == pytest.approx(0.0)

    def test_gap_introduced_for_unequal_lengths(self):
        res = needleman_wunsch("ACGT", "AGT", SCHEME)
        assert "-" in res.seq1_aligned or "-" in res.seq2_aligned

    def test_aligned_strings_same_length(self):
        res = needleman_wunsch("ACGTTGCATGCA", "ACGTCATGCA", SCHEME)
        assert len(res.seq1_aligned) == len(res.seq2_aligned)

    def test_empty_both_sequences(self):
        res = needleman_wunsch("", "", SCHEME)
        assert res.score == 0
        assert res.identity == 0.0

    def test_one_empty_sequence_all_gaps(self):
        res = needleman_wunsch("", "ACGT", SCHEME)
        assert res.seq1_aligned == "----"
        assert res.seq2_aligned == "ACGT"
        assert res.score == -8         

    def test_algorithm_label(self):
        res = needleman_wunsch("ACGT", "ACGT", SCHEME)
        assert "Needleman" in res.algorithm

    def test_elapsed_ms_non_negative(self):
        res = needleman_wunsch("ACGT", "ACGT", SCHEME)
        assert res.elapsed_ms >= 0

    def test_matrix_stored_when_requested(self):
        res = needleman_wunsch("ACGT", "ACGT", SCHEME, store_matrix=True)
        assert isinstance(res.matrix, np.ndarray)
        assert res.matrix.shape == (5, 5)   # (m+1, n+1)

    def test_matrix_none_by_default(self):
        res = needleman_wunsch("ACGT", "ACGT", SCHEME)
        assert res.matrix is None

    def test_single_match(self):
        res = needleman_wunsch("A", "A", SCHEME)
        assert res.score == 2
        assert res.identity == pytest.approx(1.0)

    def test_single_mismatch(self):
        res = needleman_wunsch("A", "T", SCHEME)
        assert res.score == -1



class TestSmithWaterman:
    def test_finds_local_alignment(self):
        res = smith_waterman("TTACGTAGTT", "GGACGTCC", SCHEME)
        assert res.score == 8
        assert res.seq1_aligned == "ACGT"
        assert res.seq2_aligned == "ACGT"

    def test_unrelated_sequences_score_zero(self):
        res = smith_waterman("AAAA", "CCCC", SCHEME)
        assert res.score == 0

    def test_identical_sequences(self):
        res = smith_waterman("ACGT", "ACGT", SCHEME)
        assert res.score == 8
        assert res.identity == pytest.approx(1.0)

    def test_score_never_negative(self):
        res = smith_waterman("GGGG", "CCCC", SCHEME)
        assert res.score >= 0

    def test_empty_sequences(self):
        res = smith_waterman("", "", SCHEME)
        assert res.score == 0

    def test_algorithm_label(self):
        res = smith_waterman("ACGT", "ACGT", SCHEME)
        assert "Smith" in res.algorithm

    def test_aligned_strings_same_length(self):
        res = smith_waterman("ACGTTGCATGCA", "ACGTCATGCA", SCHEME)
        assert len(res.seq1_aligned) == len(res.seq2_aligned)

    def test_matrix_stored_when_requested(self):
        res = smith_waterman("ACGT", "ACGT", SCHEME, store_matrix=True)
        assert isinstance(res.matrix, np.ndarray)

    def test_local_score_equals_global_for_identical(self):
        nw = needleman_wunsch("ACGT", "ACGT", SCHEME)
        sw = smith_waterman("ACGT", "ACGT", SCHEME)
        assert sw.score == nw.score



class TestBlastLike:
    def test_finds_exact_insert(self):
        query = "ACGTACGT"
        genome = "TTTTTT" + query + "TTTTTT"
        blast = BlastLike(k=6, scheme=SCHEME)
        blast.build_index(genome)
        hits = blast.search(query, top_n=1)
        assert len(hits) == 1
        assert hits[0].score > 0

    def test_no_hits_on_unrelated_genome(self):
        query = "ACGTACGT"
        genome = "C" * 100
        blast = BlastLike(k=6, scheme=SCHEME)
        blast.build_index(genome)
        hits = blast.search(query, top_n=5)
        assert len(hits) == 0

    def test_top_n_limits_results(self):
        import random
        random.seed(0)
        genome = "".join(random.choices("ACGT", k=500))
        blast = BlastLike(k=3, scheme=SCHEME)
        blast.build_index(genome)
        hits = blast.search("ACGT", top_n=2)
        assert len(hits) <= 2

    def test_hits_sorted_by_score_descending(self):
        import random
        random.seed(1)
        genome = "".join(random.choices("ACGT", k=1000))
        blast = BlastLike(k=4, scheme=SCHEME)
        blast.build_index(genome)
        hits = blast.search("ACGTACGT", top_n=5)
        scores = [h.score for h in hits]
        assert scores == sorted(scores, reverse=True)

    def test_blast_hit_fields(self):
        query = "ACGTACGT"
        genome = "NNNN" + query + "NNNN"
        blast = BlastLike(k=6, scheme=SCHEME)
        blast.build_index(genome)
        hits = blast.search(query, top_n=1)
        if hits:
            h = hits[0]
            assert h.start1 >= 0
            assert h.start2 >= 0
            assert h.length > 0
            assert isinstance(h.seq1_region, str)
            assert isinstance(h.seq2_region, str)