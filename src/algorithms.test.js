import { needlemanWunsch, smithWaterman } from "./algorithms";

describe("needlemanWunsch", () => {
  test("identical sequences yield maximum score and 100% identity", () => {
    const result = needlemanWunsch("ACGT", "ACGT");
    expect(result.score).toBe(8); // 4 matches × match score 2
    expect(result.identity).toBeCloseTo(1.0);
    expect(result.seq1Aligned).toBe("ACGT");
    expect(result.seq2Aligned).toBe("ACGT");
  });

  test("completely different sequences produce negative score", () => {
    const result = needlemanWunsch("AAAA", "CCCC");
    expect(result.score).toBe(-4); // 4 mismatches × -1
    expect(result.identity).toBeCloseTo(0.0);
  });

  test("introduces a gap for unequal length sequences", () => {
    const result = needlemanWunsch("ACGT", "AGT");
    const hasgap =
      result.seq1Aligned.includes("-") || result.seq2Aligned.includes("-");
    expect(hasgap).toBe(true);
  });

  test("empty sequences return score 0", () => {
    const result = needlemanWunsch("", "");
    expect(result.score).toBe(0);
    expect(result.identity).toBe(0);
  });

  test("aligned sequences have equal length", () => {
    const result = needlemanWunsch("ACGTTGCATGCA", "ACGTCATGCA");
    expect(result.seq1Aligned.length).toBe(result.seq2Aligned.length);
  });
});

describe("smithWaterman", () => {
  test("finds local alignment within longer sequences", () => {
    const result = smithWaterman("TTACGTAGTT", "GGACGTCC");
    expect(result.score).toBe(8);
    expect(result.seq1Aligned).toBe("ACGT");
    expect(result.seq2Aligned).toBe("ACGT");
  });

  test("unrelated sequences return score 0", () => {
    const result = smithWaterman("AAAA", "CCCC");
    expect(result.score).toBe(0);
  });

  test("identical sequences align fully", () => {
    const result = smithWaterman("ACGT", "ACGT");
    expect(result.score).toBe(8);
    expect(result.identity).toBeCloseTo(1.0);
  });

  test("score is never negative", () => {
    const result = smithWaterman("GGGG", "CCCC");
    expect(result.score).toBeGreaterThanOrEqual(0);
  });

  test("empty sequences return score 0", () => {
    const result = smithWaterman("", "");
    expect(result.score).toBe(0);
  });
});
