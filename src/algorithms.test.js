import { needlemanWunsch, smithWaterman } from "./algorithms";

// ── Needleman-Wunsch ──────────────────────────────────────────────────────────

describe("needlemanWunsch", () => {
  test("identical sequences score at maximum", () => {
    const r = needlemanWunsch("ACGT", "ACGT", 2, -1, -2);
    expect(r.score).toBe(8);
    expect(r.identity).toBeCloseTo(1.0);
    expect(r.gaps).toBe(0);
  });

  test("completely different sequences produce negative score", () => {
    const r = needlemanWunsch("AAAA", "TTTT", 2, -1, -2);
    expect(r.score).toBeLessThan(0);
  });

  test("alignment length equals sum of matches + mismatches + gaps", () => {
    const r = needlemanWunsch("ACGTTGCATGCA", "ACGTCATGCA", 2, -1, -2);
    expect(r.a1.length).toBe(r.a2.length);
  });

  test("empty sequence returns score 0", () => {
    const r = needlemanWunsch("", "ACGT", 2, -1, -2);
    expect(r.score).toBe(0);
    expect(r.a1).toBe("");
  });

  test("single character match", () => {
    const r = needlemanWunsch("A", "A", 2, -1, -2);
    expect(r.score).toBe(2);
    expect(r.identity).toBe(1);
  });

  test("single character mismatch", () => {
    const r = needlemanWunsch("A", "T", 2, -1, -2);
    expect(r.score).toBe(-1);
  });
});

// ── Smith-Waterman ────────────────────────────────────────────────────────────

describe("smithWaterman", () => {
  test("score is always non-negative", () => {
    const r = smithWaterman("AAAA", "TTTT", 2, -1, -2);
    expect(r.score).toBeGreaterThanOrEqual(0);
  });

  test("finds local alignment within longer sequences", () => {
    // GACAGCGG is shared between the two
    const r = smithWaterman("TTTTTGACAGCGGATTT", "TTAATTGACAGCGGAA", 2, -1, -2);
    expect(r.score).toBeGreaterThan(0);
    expect(r.a1).toContain("GACAGCGG");
  });

  test("identical sequences: local score equals global score", () => {
    const nw = needlemanWunsch("ACGT", "ACGT", 2, -1, -2);
    const sw = smithWaterman("ACGT", "ACGT", 2, -1, -2);
    expect(sw.score).toBe(nw.score);
  });

  test("alignment strings are equal length", () => {
    const r = smithWaterman("ACGTTGCATGCA", "ACGTCATGCA", 2, -1, -2);
    expect(r.a1.length).toBe(r.a2.length);
  });

  test("maxI and maxJ are within bounds", () => {
    const r = smithWaterman("ACGT", "ACGT", 2, -1, -2);
    expect(r.maxI).toBeLessThanOrEqual(4);
    expect(r.maxJ).toBeLessThanOrEqual(4);
  });
});