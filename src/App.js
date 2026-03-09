import { useState, useCallback } from "react";
import "./App.css";
import { needlemanWunsch, smithWaterman } from "./algorithms";

// ── MatrixHeatmap ─────────────────────────────────────────────────────────────
function MatrixHeatmap({ dp, seq1, seq2, maxI, maxJ, algorithm }) {
  if (!dp) return (
    <div className="matrix-placeholder">Run alignment to see the matrix~</div>
  );
  if (seq1.length > 30 || seq2.length > 30) return (
    <div style={{ color: "#c4a0b0", textAlign: "center", padding: "2rem", fontFamily: "monospace", fontSize: 13 }}>
      ✿ Matrix shown for sequences ≤ 30bp. Your sequences are too long to display.
    </div>
  );

  const m = seq1.length, n = seq2.length;
  let min = Infinity, max = -Infinity;
  for (let i = 0; i <= m; i++)
    for (let j = 0; j <= n; j++) {
      min = Math.min(min, dp[i][j]);
      max = Math.max(max, dp[i][j]);
    }

  const color = (v) => {
    const t = max === min ? 0.5 : (v - min) / (max - min);
    return `rgb(${Math.round(220 + t * 20)},${Math.round(160 - t * 60)},${Math.round(180 - t * 30)})`;
  };
  const cellSize = Math.min(32, Math.floor(460 / Math.max(m, n)));
  const fs = cellSize < 24 ? 9 : 11;

  return (
    <div style={{ overflowX: "auto" }}>
      <table style={{ borderCollapse: "collapse", fontFamily: "monospace", margin: "0 auto" }}>
        <thead>
          <tr>
            <td style={{ width: cellSize, height: cellSize }} />
            <td style={{ width: cellSize, height: cellSize }} />
            {seq2.split("").map((c, j) => (
              <td key={j} style={{ width: cellSize, height: cellSize, textAlign: "center", color: "#b06080", fontSize: fs, fontWeight: 700 }}>{c}</td>
            ))}
          </tr>
        </thead>
        <tbody>
          {Array.from({ length: m + 1 }, (_, i) => (
            <tr key={i}>
              <td style={{ textAlign: "center", color: "#b06080", fontSize: fs, fontWeight: 700, paddingRight: 4 }}>
                {i === 0 ? "" : seq1[i - 1]}
              </td>
              {Array.from({ length: n + 1 }, (_, j) => {
                const isMax = algorithm === "sw" && i === maxI && j === maxJ;
                return (
                  <td key={j} style={{
                    width: cellSize, height: cellSize,
                    background: color(dp[i][j]),
                    textAlign: "center", fontSize: fs, color: "#fff",
                    outline: isMax ? "2px solid #e8608a" : "none",
                    borderRadius: 3,
                  }}>
                    {dp[i][j]}
                  </td>
                );
              })}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}

// ── AlignedPair ────────────────────────────────────────────────────────────────
function AlignedPair({ a1, a2 }) {
  if (!a1) return null;
  const mid = a1.split("").map((c, i) =>
    c === a2[i] && c !== "-" ? "♡" : c !== "-" && a2[i] !== "-" ? "·" : " "
  );
  const chunk = 40;
  const segments = [];
  for (let s = 0; s < a1.length; s += chunk) {
    segments.push({
      s1:  a1.slice(s, s + chunk),
      s2:  a2.slice(s, s + chunk),
      mid: mid.slice(s, s + chunk).join(""),
      pos: s,
    });
  }
  return (
    <div style={{ fontFamily: "'Courier New', monospace", fontSize: 13, lineHeight: 1.9 }}>
      {segments.map(({ s1, s2, mid, pos }, idx) => (
        <div key={idx} style={{ marginBottom: 14 }}>
          <span style={{ color: "#c4a0b0", fontSize: 11, marginRight: 8 }}>{pos + 1}</span>
          <div style={{ display: "inline-block" }}>
            <div>{s1.split("").map((c, i) => (
              <span key={i} style={{ color: c === "-" ? "#e8608a" : s2[i] === c ? "#a855a0" : "#d4904a" }}>{c}</span>
            ))}</div>
            <div style={{ color: "#d4a0b8", letterSpacing: "1px" }}>{mid}</div>
            <div>{s2.split("").map((c, i) => (
              <span key={i} style={{ color: c === "-" ? "#e8608a" : s1[i] === c ? "#a855a0" : "#d4904a" }}>{c}</span>
            ))}</div>
          </div>
        </div>
      ))}
      <div style={{ marginTop: 10, display: "flex", gap: 20, fontSize: 12, color: "#c4a0b0" }}>
        <span><span style={{ color: "#a855a0" }}>■</span> Match</span>
        <span><span style={{ color: "#d4904a" }}>■</span> Mismatch</span>
        <span><span style={{ color: "#e8608a" }}>■</span> Gap</span>
      </div>
    </div>
  );
}

// ── MSAView ────────────────────────────────────────────────────────────────────
// Renders the multi-sequence MSA result returned from the API.
function MSAView({ msaResult }) {
  if (!msaResult) return null;
  const { aligned, labels, avg_identity, guide_tree, elapsed_ms } = msaResult;
  const chunk = 40;

  // Build consensus row: show the character if all non-gap cols agree, else '*'
  const colCount = aligned[0]?.length ?? 0;
  const consensus = Array.from({ length: colCount }, (_, col) => {
    const chars = aligned.map(row => row[col]).filter(c => c !== "-");
    if (chars.length === 0) return "-";
    return chars.every(c => c === chars[0]) ? chars[0] : "·";
  });

  const segments = [];
  for (let s = 0; s < colCount; s += chunk) {
    segments.push({
      rows: aligned.map(row => row.slice(s, s + chunk)),
      cons: consensus.slice(s, s + chunk),
      pos:  s,
    });
  }

  const labelWidth = Math.max(...labels.map(l => l.length), 4) + 2;

  return (
    <div>
      {/* Stats bar */}
      <div style={{ display: "flex", gap: 20, marginBottom: 16, fontSize: 12, color: "#ad6080" }}>
        <span>Avg identity: <strong style={{ color: "#880e4f" }}>{(avg_identity * 100).toFixed(1)}%</strong></span>
        <span>Sequences: <strong style={{ color: "#880e4f" }}>{aligned.length}</strong></span>
        <span>Length: <strong style={{ color: "#880e4f" }}>{colCount} cols</strong></span>
        <span>Time: <strong style={{ color: "#880e4f" }}>{elapsed_ms.toFixed(1)} ms</strong></span>
      </div>

      {/* Guide tree summary */}
      <div style={{ marginBottom: 12, fontSize: 11, color: "#c4a0b0", fontStyle: "italic" }}>
        Guide tree merges: {guide_tree.map(([i, j, d], idx) => (
          <span key={idx}>
            {idx > 0 && " → "}}>
            ({labels[i]} + {labels[j]}, dist={d.toFixed(2)})
          </span>
        ))}
      </div>

      {/* Alignment blocks */}
      <div style={{ fontFamily: "'Courier New', monospace", fontSize: 13, lineHeight: 1.8 }}>
        {segments.map(({ rows, cons, pos }, segIdx) => (
          <div key={segIdx} style={{ marginBottom: 18 }}>
            <div style={{ color: "#c4a0b0", fontSize: 10, marginBottom: 4 }}>pos {pos + 1}–{pos + rows[0].length}</div>
            {rows.map((row, ri) => (
              <div key={ri} style={{ display: "flex", alignItems: "center" }}>
                <span style={{
                  color: "#ad6080", fontSize: 11, width: `${labelWidth}ch`,
                  flexShrink: 0, fontWeight: 600,
                }}>
                  {labels[ri]}
                </span>
                <span>
                  {row.split("").map((c, ci) => {
                    const isGap = c === "-";
                    const isMatch = cons[ci] !== "·" && !isGap;
                    return (
                      <span key={ci} style={{
                        color: isGap ? "#e8608a" : isMatch ? "#a855a0" : "#d4904a",
                      }}>{c}</span>
                    );
                  })}
                </span>
              </div>
            ))}
            {/* Consensus row */}
            <div style={{ display: "flex", alignItems: "center", marginTop: 2 }}>
              <span style={{ color: "#c4a0b0", fontSize: 11, width: `${labelWidth}ch`, flexShrink: 0 }}>
                consensus
              </span>
              <span style={{ color: "#c4a0b0", letterSpacing: "0px" }}>
                {cons.map((c, ci) => (
                  <span key={ci} style={{ color: c === "·" ? "#d4904a" : c === "-" ? "#e8608a" : "#a855a0" }}>{c}</span>
                ))}
              </span>
            </div>
          </div>
        ))}
      </div>

      <div style={{ marginTop: 10, display: "flex", gap: 20, fontSize: 12, color: "#c4a0b0" }}>
        <span><span style={{ color: "#a855a0" }}>■</span> Conserved</span>
        <span><span style={{ color: "#d4904a" }}>■</span> Variable</span>
        <span><span style={{ color: "#e8608a" }}>■</span> Gap</span>
      </div>
    </div>
  );
}

// ── App ────────────────────────────────────────────────────────────────────────
const PRESETS = [
  { label: "✿ Similar",   s1: "ACGTTGCATGCA",           s2: "ACGTCATGCA"         },
  { label: "✿ Divergent", s1: "GCTAGCTAGCTA",           s2: "ATGATGATGATG"       },
  { label: "✿ Local hit", s1: "TTTTTGACAGCGGATTT",      s2: "TTAATTGACAGCGGAA"   },
  { label: "✿ Protein",   s1: "MKTLLLTLVVVTIVCLDLGAVK", s2: "MKTLLVTIVCLDLGAVKD" },
];

const MSA_PRESETS = [
  {
    label: "✿ DNA trio",
    seqs: ["ACGTTGCATGCA", "ACGTCATGCA", "ACGCATGCA"],
    labels: ["Seq1", "Seq2", "Seq3"],
  },
  {
    label: "✿ Protein trio",
    seqs: ["MKTLLLTLVVVTIVCLDLGAVK", "MKTLLVTIVCLDLGAVKD", "MKTLVTIVCLDL"],
    labels: ["Human", "Mouse", "Zebrafish"],
  },
];

const API_BASE = process.env.REACT_APP_API_URL || "";

export default function App() {
  // ── pairwise state ──
  const [seq1,     setSeq1]     = useState("ACGTTGCATGCA");
  const [seq2,     setSeq2]     = useState("ACGTCATGCA");
  const [algo,     setAlgo]     = useState("nw");
  const [match,    setMatch]    = useState(2);
  const [mismatch, setMismatch] = useState(-1);
  const [gap,      setGap]      = useState(-2);
  const [blosum,   setBlosum]   = useState(false);
  const [result,   setResult]   = useState(null);
  const [tab,      setTab]      = useState("alignment");

  // ── MSA state ──
  const [mode,      setMode]      = useState("pairwise");  // "pairwise" | "msa"
  const [msaSeqs,   setMsaSeqs]   = useState(["ACGTTGCATGCA", "ACGTCATGCA", "ACGCATGCA"]);
  const [msaLabels, setMsaLabels] = useState(["Seq1", "Seq2", "Seq3"]);
  const [msaResult, setMsaResult] = useState(null);
  const [msaError,  setMsaError]  = useState(null);
  const [msaLoading, setMsaLoading] = useState(false);

  // ── pairwise run (client-side, instant) ──
  const run = useCallback(() => {
    const s1 = seq1.toUpperCase().replace(/[^A-Z]/g, "");
    const s2 = seq2.toUpperCase().replace(/[^A-Z]/g, "");
    if (!s1 || !s2) return;
    const fn = algo === "nw" ? needlemanWunsch : smithWaterman;
    const t0 = performance.now();
    const res = fn(s1, s2, match, mismatch, gap, blosum);
    res.elapsed = (performance.now() - t0).toFixed(2);
    res.algo = algo;
    setResult(res);
    setTab("alignment");
  }, [seq1, seq2, algo, match, mismatch, gap, blosum]);

  // ── MSA run (server-side via Flask API) ──
  const runMSA = useCallback(async () => {
    const sequences = msaSeqs.map(s => s.toUpperCase().replace(/[^A-Z]/g, "")).filter(Boolean);
    if (sequences.length < 2) {
      setMsaError("Please enter at least 2 non-empty sequences.");
      return;
    }
    setMsaLoading(true);
    setMsaError(null);
    setMsaResult(null);
    try {
      const resp = await fetch(`${API_BASE}/api/align/msa`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          sequences,
          labels: msaLabels.slice(0, sequences.length),
          match,
          mismatch,
          gap_open: gap,
        }),
      });
      const data = await resp.json();
      if (!resp.ok) {
        setMsaError(data.error || "Server error.");
      } else {
        setMsaResult(data);
      }
    } catch (err) {
      setMsaError(`Network error: ${err.message}`);
    } finally {
      setMsaLoading(false);
    }
  }, [msaSeqs, msaLabels, match, mismatch, gap]);

  const updateMsaSeq = (idx, val) => {
    setMsaSeqs(prev => prev.map((s, i) => i === idx ? val.toUpperCase() : s));
  };
  const updateMsaLabel = (idx, val) => {
    setMsaLabels(prev => prev.map((l, i) => i === idx ? val : l));
  };
  const addSequence = () => {
    if (msaSeqs.length >= 10) return;
    setMsaSeqs(prev => [...prev, ""]);
    setMsaLabels(prev => [...prev, `Seq${prev.length + 1}`]);
  };
  const removeSequence = (idx) => {
    if (msaSeqs.length <= 2) return;
    setMsaSeqs(prev => prev.filter((_, i) => i !== idx));
    setMsaLabels(prev => prev.filter((_, i) => i !== idx));
  };

  const exportResult = () => {
    if (!result) return;
    const text = [
      `Algorithm : ${result.algo === "nw" ? "Needleman-Wunsch (Global)" : "Smith-Waterman (Local)"}`,
      `Score     : ${result.score}`,
      `Identity  : ${(result.identity * 100).toFixed(1)}%`,
      `Gaps      : ${result.gaps}`,
      `Length    : ${result.a1.length} cols`,
      `Time      : ${result.elapsed} ms`,
      `BLOSUM62  : ${blosum ? "yes" : "no"}`,
      ``,
      `Seq1: ${result.a1}`,
      `      ${result.a1.split("").map((c, i) => c === result.a2[i] && c !== "-" ? "|" : " ").join("")}`,
      `Seq2: ${result.a2}`,
    ].join("\n");
    const blob = new Blob([text], { type: "text/plain" });
    const url  = URL.createObjectURL(blob);
    const a    = document.createElement("a");
    a.href = url; a.download = "alignment_result.txt"; a.click();
    URL.revokeObjectURL(url);
  };

  const exportMSA = () => {
    if (!msaResult) return;
    const lines = msaResult.labels.map((lbl, i) =>
      `>${lbl}\n${msaResult.aligned[i]}`
    );
    lines.push(`>consensus\n${msaResult.aligned[0].split("").map((_, col) => {
      const chars = msaResult.aligned.map(r => r[col]).filter(c => c !== "-");
      return chars.length && chars.every(c => c === chars[0]) ? chars[0] : "*";
    }).join("")}`);
    const blob = new Blob([lines.join("\n")], { type: "text/plain" });
    const url  = URL.createObjectURL(blob);
    const a    = document.createElement("a");
    a.href = url; a.download = "msa_result.fasta"; a.click();
    URL.revokeObjectURL(url);
  };

  return (
    <div className="app">
      <div className="deco" style={{ top: "5%",  left: "-2%" }}>🌸</div>
      <div className="deco" style={{ bottom: "10%", right: "-2%" }}>🌷</div>

      {/* Header */}
      <div className="header">
        <div className="logo">🧬</div>
        <div>
          <div style={{ fontFamily: "Playfair Display, serif", fontSize: 24, fontWeight: 600, color: "#880e4f" }}>
            Sequence Alignment Studio
          </div>
          <div style={{ fontSize: 12, color: "#ad6080", marginTop: 3, fontStyle: "italic" }}>
            Needleman-Wunsch · Smith-Waterman · BLOSUM62 · Progressive MSA
          </div>
        </div>

        {/* Mode selector */}
        <div style={{ marginLeft: "auto", display: "flex", gap: 10, alignItems: "center" }}>
          {['pairwise', 'msa'].map(m => (
            <button key={m}
              className={`algo-btn ${mode === m ? "active" : ""}`}
              onClick={() => setMode(m)}>
              {m === "pairwise" ? "Pairwise" : "Multiple (MSA)"}
            </button>
          ))}
        </div>
      </div>

      <div className="main">

        {/* ── PAIRWISE MODE ── */}
        {mode === "pairwise" && (
          <>
            {/* Algorithm selector */}
            <div style={{ display: "flex", gap: 10, marginBottom: 12 }}>
              {['NW', 'SW'].map(a => (
                <button key={a}
                  className={`algo-btn ${algo === a.toLowerCase() ? "active" : ""}`}
                  onClick={() => setAlgo(a.toLowerCase())}>
                  {a === "NW" ? "Global (NW)" : "Local (SW)"}
                </button>
              ))}
            </div>

            {/* Presets */}
            <div style={{ display: "flex", gap: 8, marginBottom: 18, flexWrap: "wrap", alignItems: "center" }}>
              <span style={{ fontSize: 12, color: "#c4a0b0", fontStyle: "italic" }}>try a preset:</span>
              {PRESETS.map(p => (
                <button key={p.label} className="preset-btn"
                  onClick={() => { setSeq1(p.s1); setSeq2(p.s2); }}>
                  {p.label}
                </button>
              ))}
            </div>

            {/* Input card */}
            <div className="card">
              <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 16, marginBottom: 16 }}>
                {[['Sequence 1 ✦', seq1, setSeq1], ['Sequence 2 ✦', seq2, setSeq2]].map(([label, val, setter]) => (
                  <div key={label}>
                    <div className="section-label">{label}</div>
                    <textarea value={val} onChange={e => setter(e.target.value.toUpperCase())} rows={2} />
                  </div>
                ))}
              </div>
              <div style={{ display: "flex", gap: 20, alignItems: "center", flexWrap: "wrap" }}>
                {[['Match ♡', match, setMatch], ['Mismatch', mismatch, setMismatch], ['Gap', gap, setGap]].map(([label, val, setter]) => (
                  <label key={label} style={{ fontSize: 12, color: "#ad6080", display: "flex", alignItems: "center", gap: 8 }}>
                    {label}
                    <input type="number" value={val} onChange={e => setter(Number(e.target.value))} className="number-input" />
                  </label>
                ))}
                <label className="blosum-toggle">
                  <input type="checkbox" checked={blosum} onChange={e => setBlosum(e.target.checked)} />
                  BLOSUM62 (proteins)
                </label>
                <button onClick={run} className="run-btn">Align ✦</button>
              </div>
            </div>

            {/* Stats row */}
            {result && (
              <div style={{ display: "flex", gap: 10, marginBottom: 16, flexWrap: "wrap", alignItems: "center" }}>
                {[['Score',     result.score], ['Identity',  (result.identity * 100).toFixed(1) + "%"], ['Gaps',      result.gaps], ['Length',    result.a1.length + " cols"], ['Time',      result.elapsed + " ms"]].map(([k, v]) => (
                  <div key={k} className="stat-card">
                    <div className="stat-label">{k}</div>
                    <div className="stat-value">{v}</div>
                  </div>
                ))}
                <button onClick={exportResult} className="export-btn">⬇ Export .txt</button>
              </div>
            )}

            {/* Tabs */}
            <div style={{ display: "flex", borderBottom: "1px solid #f8bbd0" }}>
              {['alignment', 'matrix'].map(t => (
                <button key={t} className={`tab-btn ${tab === t ? "active" : ""}`} onClick={() => setTab(t)}>
                  {t === "alignment" ? "✦ Alignment View" : "✦ DP Matrix"}
                </button>
              ))}
            </div>

            <div className="content-box">
              {result && tab === "alignment" && <AlignedPair a1={result.a1} a2={result.a2} />}
              {result && tab === "matrix"    && (
                <MatrixHeatmap dp={result.dp} seq1={seq1.toUpperCase()} seq2={seq2.toUpperCase()}
                  maxI={result.maxI} maxJ={result.maxJ} algorithm={result.algo} />
              )}
              {!result && (
                <div style={{ color: "#d4a0b8", textAlign: "center", paddingTop: 60, fontStyle: "italic" }}>
                  ✿ press Align to begin ✿
                </div>
              )}
            </div>

            {/* Info box */}
            <div className="info-box">
              {algo === "nw" ? (
                <>
                  <strong style={{ color: "#880e4f", fontFamily: "Playfair Display, serif" }}>Needleman-Wunsch (Global)</strong>
                  {" "}— fills an (m+1)×(n+1) DP table. Every cell scores max(diagonal+sub, up+gap, left+gap).
                  Traceback from bottom-right to top-left. Complexity: <code style={{ color: "#e91e96" }}>O(mn)</code> time &amp; space.
                </>
              ) : (
                <>
                  <strong style={{ color: "#880e4f", fontFamily: "Playfair Display, serif" }}>Smith-Waterman (Local)</strong>
                  {" "}— same recurrence but scores are floored at 0, allowing alignment to restart anywhere.
                  Traceback begins at the highest-scoring cell and stops at 0. Complexity: <code style={{ color: "#e91e96" }}>O(mn)</code> time &amp; space.
                </>
              )}
            </div>
          </>
        )}

        {/* ── MSA MODE ── */}
        {mode === "msa" && (
          <>
            {/* MSA Presets */}
            <div style={{ display: "flex", gap: 8, marginBottom: 18, flexWrap: "wrap", alignItems: "center" }}>
              <span style={{ fontSize: 12, color: "#c4a0b0", fontStyle: "italic" }}>try a preset:</span>
              {MSA_PRESETS.map(p => (
                <button key={p.label} className="preset-btn"
                  onClick={() => { setMsaSeqs(p.seqs); setMsaLabels(p.labels); setMsaResult(null); setMsaError(null); }}>
                  {p.label}
                </button>
              ))}
            </div>

            {/* MSA Input card */}
            <div className="card">
              <div style={{ marginBottom: 12 }}>
                {msaSeqs.map((seq, idx) => (
                  <div key={idx} style={{ display: "grid", gridTemplateColumns: "120px 1fr auto", gap: 10, marginBottom: 10, alignItems: "center" }}>
                    <input
                      value={msaLabels[idx] || `Seq${idx + 1}`}
                      onChange={e => updateMsaLabel(idx, e.target.value)}
                      style={{ fontSize: 12, padding: "4px 8px", border: "1px solid #f8bbd0", borderRadius: 6, color: "#880e4f", fontWeight: 600 }}
                      placeholder={`Label ${idx + 1}`}
                    />
                    <textarea
                      value={seq}
                      onChange={e => updateMsaSeq(idx, e.target.value)}
                      rows={1}
                      placeholder={`Sequence ${idx + 1} (DNA or protein)`}
                      style={{ resize: "vertical" }}
                    />
                    <button
                      onClick={() => removeSequence(idx)}
                      disabled={msaSeqs.length <= 2}
                      style={{ background: "none", border: "1px solid #f8bbd0", borderRadius: 6, color: "#e8608a", cursor: "pointer", padding: "4px 10px", fontSize: 16 }}
                      title="Remove sequence"
                    >✕</button>
                  </div>
                ))}
              </div>

              <div style={{ display: "flex", gap: 12, alignItems: "center", flexWrap: "wrap" }}>
                <button onClick={addSequence} disabled={msaSeqs.length >= 10} className="preset-btn">
                  + Add sequence
                </button>
                {[['Match ♡', match, setMatch], ['Mismatch', mismatch, setMismatch], ['Gap', gap, setGap]].map(([label, val, setter]) => (
                  <label key={label} style={{ fontSize: 12, color: "#ad6080", display: "flex", alignItems: "center", gap: 8 }}>
                    {label}
                    <input type="number" value={val} onChange={e => setter(Number(e.target.value))} className="number-input" />
                  </label>
                ))}
                <button onClick={runMSA} className="run-btn" disabled={msaLoading}>
                  {msaLoading ? "Aligning…" : "Align ✦"}
                </button>
              </div>
            </div>

            {/* MSA error */}
            {msaError && (
              <div style={{ color: "#e8608a", background: "#fff0f3", border: "1px solid #f8bbd0", borderRadius: 8, padding: "10px 16px", marginBottom: 12, fontSize: 13 }}>
                ✕ {msaError}
              </div>
            )}

            {/* MSA result */}
            {msaResult && (
              <>
                <div style={{ display: "flex", justifyContent: "flex-end", marginBottom: 8 }}>
                  <button onClick={exportMSA} className="export-btn">⬇ Export .fasta</button>
                </div>
                <div className="content-box">
                  <MSAView msaResult={msaResult} />
                </div>
              </>
            )}

            {!msaResult && !msaError && !msaLoading && (
              <div className="content-box">
                <div style={{ color: "#d4a0b8", textAlign: "center", paddingTop: 60, fontStyle: "italic" }}>
                  ✿ enter sequences and press Align ✿
                </div>
              </div>
            )}

            {/* MSA Info box */}
            <div className="info-box">
              <strong style={{ color: "#880e4f", fontFamily: "Playfair Display, serif" }}>Progressive MSA (ClustalW-style)</strong>
              {" "}— builds a pairwise distance matrix (NW identity), constructs a greedy UPGMA guide tree,
              then merges sequence profiles in tree order using the <em>"once a gap, always a gap"</em> rule.
              Complexity: <code style={{ color: "#e91e96" }}>O(N² · L²)</code> where N = # sequences, L = length.
            </div>
          </>
        )}

      </div>
    </div>
