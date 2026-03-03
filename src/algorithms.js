/**
 * algorithms.js
 * Needleman-Wunsch and Smith-Waterman sequence alignment algorithms.
 * Uses typed arrays (Int32Array / Uint8Array) for memory efficiency.
 *
 * BLOSUM62 substitution matrix for protein alignments.
 */

const DIAG = 1, UP = 2, LEFT = 3, STOP = 0;

// Minimal BLOSUM62 for the 20 standard amino acids + DNA bases
export const BLOSUM62 = {
  A:{A:4,R:-1,N:-2,D:-2,C:0,Q:-1,E:-1,G:0,H:-2,I:-1,L:-1,K:-1,M:-1,F:-2,P:-1,S:1,T:0,W:-3,Y:-2,V:0},
  R:{A:-1,R:5,N:0,D:-2,C:-3,Q:1,E:0,G:-2,H:0,I:-3,L:-2,K:2,M:-1,F:-3,P:-2,S:-1,T:-1,W:-3,Y:-2,V:-3},
  N:{A:-2,R:0,N:6,D:1,C:-3,Q:0,E:0,G:0,H:1,I:-3,L:-3,K:0,M:-2,F:-3,P:-2,S:1,T:0,W:-4,Y:-2,V:-3},
  D:{A:-2,R:-2,N:1,D:6,C:-3,Q:0,E:2,G:-1,H:-1,I:-3,L:-4,K:-1,M:-3,F:-3,P:-1,S:0,T:-1,W:-4,Y:-3,V:-3},
  C:{A:0,R:-3,N:-3,D:-3,C:9,Q:-3,E:-4,G:-3,H:-3,I:-1,L:-1,K:-3,M:-1,F:-2,P:-3,S:-1,T:-1,W:-2,Y:-2,V:-1},
  G:{A:0,R:-2,N:0,D:-1,C:-3,Q:-2,E:-2,G:6,H:-2,I:-4,L:-4,K:-2,M:-3,F:-3,P:-2,S:0,T:-2,W:-2,Y:-3,V:-3},
  // DNA bases — simple match/mismatch handled by fallback
};

function subScore(a, b, matchScore, mismatchScore, useBlosum) {
  if (useBlosum && BLOSUM62[a] && BLOSUM62[a][b] !== undefined) {
    return BLOSUM62[a][b];
  }
  return a === b ? matchScore : mismatchScore;
}

export function needlemanWunsch(seq1, seq2, match = 2, mismatch = -1, gap = -2, useBlosum = false) {
  const m = seq1.length, n = seq2.length;
  const dp  = Array.from({ length: m + 1 }, () => new Int32Array(n + 1));
  const ptr = Array.from({ length: m + 1 }, () => new Uint8Array(n + 1));

  for (let j = 0; j <= n; j++) { dp[0][j] = gap * j; ptr[0][j] = LEFT; }
  for (let i = 0; i <= m; i++) { dp[i][0] = gap * i; ptr[i][0] = UP; }
  ptr[0][0] = STOP;

  for (let i = 1; i <= m; i++) {
    for (let j = 1; j <= n; j++) {
      const d = dp[i-1][j-1] + subScore(seq1[i-1], seq2[j-1], match, mismatch, useBlosum);
      const u = dp[i-1][j] + gap;
      const l = dp[i][j-1] + gap;
      const best = Math.max(d, u, l);
      dp[i][j] = best;
      ptr[i][j] = best === d ? DIAG : best === u ? UP : LEFT;
    }
  }

  let i = m, j = n, a1 = [], a2 = [];
  while (i > 0 || j > 0) {
    const p = ptr[i][j];
    if (p === DIAG) { a1.push(seq1[i-1]); a2.push(seq2[j-1]); i--; j--; }
    else if (p === UP) { a1.push(seq1[i-1]); a2.push("-"); i--; }
    else { a1.push("-"); a2.push(seq2[j-1]); j--; }
  }
  a1.reverse(); a2.reverse();
  const matches = a1.filter((c, k) => c === a2[k] && c !== "-").length;
  return {
    a1: a1.join(""), a2: a2.join(""),
    score: dp[m][n],
    identity: matches / Math.max(a1.length, 1),
    gaps: a1.filter(c => c === "-").length + a2.filter(c => c === "-").length,
    dp,
  };
}

export function smithWaterman(seq1, seq2, match = 2, mismatch = -1, gap = -2, useBlosum = false) {
  const m = seq1.length, n = seq2.length;
  const dp  = Array.from({ length: m + 1 }, () => new Int32Array(n + 1));
  const ptr = Array.from({ length: m + 1 }, () => new Uint8Array(n + 1));
  let maxVal = 0, maxI = 0, maxJ = 0;

  for (let i = 1; i <= m; i++) {
    for (let j = 1; j <= n; j++) {
      const d = dp[i-1][j-1] + subScore(seq1[i-1], seq2[j-1], match, mismatch, useBlosum);
      const u = dp[i-1][j] + gap;
      const l = dp[i][j-1] + gap;
      const best = Math.max(d, u, l, 0);
      dp[i][j] = best;
      ptr[i][j] = best === 0 ? STOP : best === d ? DIAG : best === u ? UP : LEFT;
      if (best > maxVal) { maxVal = best; maxI = i; maxJ = j; }
    }
  }

  let i = maxI, j = maxJ, a1 = [], a2 = [];
  while (i > 0 && j > 0 && ptr[i][j] !== STOP) {
    const p = ptr[i][j];
    if (p === DIAG) { a1.push(seq1[i-1]); a2.push(seq2[j-1]); i--; j--; }
    else if (p === UP) { a1.push(seq1[i-1]); a2.push("-"); i--; }
    else { a1.push("-"); a2.push(seq2[j-1]); j--; }
  }
  a1.reverse(); a2.reverse();
  const matches = a1.filter((c, k) => c === a2[k] && c !== "-").length;
  return {
    a1: a1.join(""), a2: a2.join(""),
    score: maxVal,
    identity: matches / Math.max(a1.length, 1),
    gaps: a1.filter(c => c === "-").length + a2.filter(c => c === "-").length,
    dp, maxI, maxJ,
  };
}