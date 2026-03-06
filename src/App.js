import { useState, useCallback } from "react";
import "./App.css";
import { needlemanWunsch, smithWaterman } from "./algorithms";

// ── MatrixHeatmap ───────────────────────────────────────────────────────────
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

// ── AlignedPair ───────────────────────────────────────────────────────────
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

// ── App ───────────────────────────────────────────────────────────────────
const PRESETS = [
  { label: "✿ Similar",   s1: "ACGTTGCATGCA",           s2: "ACGTCATGCA"         },
  { label: "✿ Divergent", s1: "GCTAGCTAGCTA",           s2: "ATGATGATGATG"       },
  { label: "✿ Local hit", s1: "TTTTTGACAGCGGATTT",      s2: "TTAATTGACAGCGGAA"   },
  { label: "✿ Protein",   s1: "MKTLLLTLVVVTIVCLDLGAVK", s2: "MKTLLVTIVCLDLGAVKD" },
];

export default function App() {
  const [seq1,     setSeq1]     = useState("ACGTTGCATGCA");
  const [seq2,     setSeq2]     = useState("ACGTCATGCA");
  const [algo,     setAlgo]     = useState("nw");
  const [match,    setMatch]    = useState(2);
  const [mismatch, setMismatch] = useState(-1);
  const [gap,      setGap]      = useState(-2);
  const [blosum,   setBlosum]   = useState(false);
  const [result,   setResult]   = useState(null);
  const [tab,      setTab]      = useState("alignment");

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
            Needleman-Wunsch · Smith-Waterman · BLOSUM62
          </div>
        </div>
        <div style={{ marginLeft: "auto", display: "flex", gap: 10 }}>
          {["NW", "SW"].map(a => (
            <button key={a}
              className={`algo-btn ${algo === a.toLowerCase() ? "active" : ""}`}
              onClick={() => setAlgo(a.toLowerCase())}>
              {a === "NW" ? "Global (NW)" : "Local (SW)"}
            </button>
          ))}
        </div>
      </div>

      <div className="main">
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
            {[["Sequence 1 ✦", seq1, setSeq1], ["Sequence 2 ✦", seq2, setSeq2]].map(([label, val, setter]) => (
              <div key={label}>
                <div className="section-label">{label}</div>
                <textarea value={val} onChange={e => setter(e.target.value.toUpperCase())} rows={2} />
              </div>
            ))}
          </div>
          <div style={{ display: "flex", gap: 20, alignItems: "center", flexWrap: "wrap" }}>
            {[["Match ♡", match, setMatch], ["Mismatch", mismatch, setMismatch], ["Gap", gap, setGap]].map(([label, val, setter]) => (
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
            {[
              ["Score",     result.score],
              ["Identity",  (result.identity * 100).toFixed(1) + "%"],
              ["Gaps",      result.gaps],
              ["Length",    result.a1.length + " cols"],
              ["Time",      result.elapsed + " ms"],
            ].map(([k, v]) => (
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
          {["alignment", "matrix"].map(t => (
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
      </div>
    </div>
  );
}
