import { useState, useEffect } from "react";

const DIAG = 1, UP = 2, LEFT = 3, STOP = 0;

function buildSub(s1, s2, match, mismatch) {
  return s1.split("").map(a => s2.split("").map(b => (a === b ? match : mismatch)));
}

function needlemanWunsch(seq1, seq2, match = 2, mismatch = -1, gap = -2) {
  const m = seq1.length, n = seq2.length;
  const dp  = Array.from({ length: m + 1 }, () => new Int32Array(n + 1));
  const ptr = Array.from({ length: m + 1 }, () => new Uint8Array(n + 1));
  const sub = buildSub(seq1, seq2, match, mismatch);
  for (let j = 0; j <= n; j++) { dp[0][j] = gap * j; ptr[0][j] = LEFT; }
  for (let i = 0; i <= m; i++) { dp[i][0] = gap * i; ptr[i][0] = UP; }
  ptr[0][0] = STOP;
  for (let i = 1; i <= m; i++) {
    for (let j = 1; j <= n; j++) {
      const d = dp[i-1][j-1] + sub[i-1][j-1];
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
  return { a1: a1.join(""), a2: a2.join(""), score: dp[m][n], identity: matches / Math.max(a1.length, 1), dp };
}

function smithWaterman(seq1, seq2, match = 2, mismatch = -1, gap = -2) {
  const m = seq1.length, n = seq2.length;
  const dp  = Array.from({ length: m + 1 }, () => new Int32Array(n + 1));
  const ptr = Array.from({ length: m + 1 }, () => new Uint8Array(n + 1));
  const sub = buildSub(seq1, seq2, match, mismatch);
  let maxVal = 0, maxI = 0, maxJ = 0;
  for (let i = 1; i <= m; i++) {
    for (let j = 1; j <= n; j++) {
      const d = dp[i-1][j-1] + sub[i-1][j-1];
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
  return { a1: a1.join(""), a2: a2.join(""), score: maxVal, identity: matches / Math.max(a1.length, 1), dp, maxI, maxJ };
}

function MatrixHeatmap({ dp, seq1, seq2, maxI, maxJ, algorithm }) {
  if (!dp || seq1.length > 18 || seq2.length > 18) {
    return (
      <div style={{ color: "#c4a0b0", textAlign: "center", padding: "2rem", fontFamily: "monospace", fontSize: 13 }}>
        {seq1.length > 18 || seq2.length > 18
          ? "✿ Matrix hidden for sequences longer than 18bp. Use shorter sequences to visualize ✿"
          : "Run alignment to see the matrix~"}
      </div>
    );
  }
  const m = seq1.length, n = seq2.length;
  let min = Infinity, max = -Infinity;
  for (let i = 0; i <= m; i++) for (let j = 0; j <= n; j++) {
    min = Math.min(min, dp[i][j]);
    max = Math.max(max, dp[i][j]);
  }
  const color = (v) => {
    const t = max === min ? 0.5 : (v - min) / (max - min);
    const r = Math.round(220 + t * 20);
    const g = Math.round(160 - t * 60);
    const b = Math.round(180 - t * 30);
    return `rgb(${r},${g},${b})`;
  };
  const cellSize = Math.min(32, Math.floor(320 / Math.max(m, n)));
  const fontSize = cellSize < 24 ? 9 : 11;
  return (
    <div style={{ overflowX: "auto" }}>
      <table style={{ borderCollapse: "collapse", fontFamily: "monospace", margin: "0 auto" }}>
        <thead>
          <tr>
            <td style={{ width: cellSize, height: cellSize }} />
            <td style={{ width: cellSize, height: cellSize }} />
            {seq2.split("").map((c, j) => (
              <td key={j} style={{ width: cellSize, height: cellSize, textAlign: "center", color: "#b06080", fontSize, fontWeight: 700 }}>{c}</td>
            ))}
          </tr>
        </thead>
        <tbody>
          {Array.from({ length: m + 1 }, (_, i) => (
            <tr key={i}>
              <td style={{ textAlign: "center", color: "#b06080", fontSize, fontWeight: 700, paddingRight: 4 }}>
                {i === 0 ? "" : seq1[i - 1]}
              </td>
              {Array.from({ length: n + 1 }, (_, j) => {
                const isMax = algorithm === "sw" && i === maxI && j === maxJ;
                const val = dp[i][j];
                return (
                  <td key={j} style={{
                    width: cellSize, height: cellSize,
                    background: color(val),
                    textAlign: "center",
                    fontSize,
                    color: "#fff",
                    outline: isMax ? "2px solid #e8608a" : "none",
                    borderRadius: 3,
                  }}>
                    {val}
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

function AlignedPair({ a1, a2 }) {
  if (!a1) return null;
  const mid = a1.split("").map((c, i) => c === a2[i] && c !== "-" ? "♡" : c !== "-" && a2[i] !== "-" ? "·" : " ");
  const chunk = 40;
  const segments = [];
  for (let s = 0; s < a1.length; s += chunk) {
    segments.push({ s1: a1.slice(s, s + chunk), s2: a2.slice(s, s + chunk), mid: mid.slice(s, s + chunk).join("") });
  }
  return (
    <div style={{ fontFamily: "'Courier New', monospace", fontSize: 13, lineHeight: 1.9 }}>
      {segments.map(({ s1, s2, mid }, idx) => (
        <div key={idx} style={{ marginBottom: 14 }}>
          <div>{s1.split("").map((c, i) => (
            <span key={i} style={{ color: c === "-" ? "#e8608a" : s2[i] === c ? "#a855a0" : "#d4904a" }}>{c}</span>
          ))}</div>
          <div style={{ color: "#d4a0b8", letterSpacing: "1px" }}>{mid}</div>
          <div>{s2.split("").map((c, i) => (
            <span key={i} style={{ color: c === "-" ? "#e8608a" : s1[i] === c ? "#a855a0" : "#d4904a" }}>{c}</span>
          ))}</div>
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

const styles = `
  @import url('https://fonts.googleapis.com/css2?family=Playfair+Display:ital,wght@0,400;0,600;1,400&family=DM+Sans:wght@300;400;500&display=swap');
  * { box-sizing: border-box; margin: 0; padding: 0; }
  .app {
    min-height: 100vh;
    background: radial-gradient(ellipse at 0% 0%, #fce4ec 0%, transparent 50%),
                radial-gradient(ellipse at 100% 100%, #f3e5f5 0%, transparent 50%),
                radial-gradient(ellipse at 50% 50%, #fdf0f5 0%, #fce4ec 100%);
    font-family: 'DM Sans', sans-serif;
    color: #5a2d4a;
    padding-bottom: 60px;
  }
  .header {
    background: linear-gradient(135deg, #f8bbd0 0%, #e1bee7 100%);
    border-bottom: 1px solid #f48fb1;
    padding: 28px 36px;
    display: flex; align-items: center; gap: 18px;
    box-shadow: 0 4px 20px rgba(233,30,99,0.1);
  }
  .logo {
    width: 48px; height: 48px; border-radius: 50%;
    background: linear-gradient(135deg, #e91e96, #ab47bc);
    display: flex; align-items: center; justify-content: center;
    font-size: 22px; box-shadow: 0 4px 15px rgba(233,30,99,0.3);
  }
  .algo-btn {
    padding: 8px 20px; border-radius: 20px; border: 2px solid #f48fb1;
    cursor: pointer; font-size: 13px; font-weight: 500;
    transition: all 0.25s; font-family: 'DM Sans', sans-serif;
  }
  .algo-btn.active {
    background: linear-gradient(135deg, #e91e96, #ab47bc);
    border-color: transparent; color: #fff;
    box-shadow: 0 4px 15px rgba(233,30,99,0.35);
  }
  .algo-btn:not(.active) { background: rgba(255,255,255,0.6); color: #ad6080; }
  .main { max-width: 900px; margin: 0 auto; padding: 28px 24px 0; }
  .card {
    background: rgba(255,255,255,0.7); border: 1px solid #f8bbd0;
    border-radius: 16px; padding: 20px;
    box-shadow: 0 4px 20px rgba(233,30,99,0.08);
  }
  .section-label {
    font-size: 10px; color: #ad6080; text-transform: uppercase;
    letter-spacing: 2px; margin-bottom: 6px; font-weight: 500;
  }
  textarea {
    width: 100%; padding: 12px 14px;
    background: rgba(255,255,255,0.8); border: 1.5px solid #f8bbd0;
    border-radius: 10px; color: #5a2d4a; font-size: 13px;
    font-family: 'Courier New', monospace; resize: vertical; outline: none;
    transition: border-color 0.2s;
  }
  textarea:focus { border-color: #e91e96; box-shadow: 0 0 0 3px rgba(233,30,99,0.1); }
  .preset-btn {
    padding: 5px 14px; border-radius: 20px; border: 1.5px solid #f48fb1;
    background: rgba(255,255,255,0.6); color: #ad6080; font-size: 12px;
    cursor: pointer; transition: all 0.2s; font-family: 'DM Sans', sans-serif;
  }
  .preset-btn:hover { background: #fce4ec; border-color: #e91e96; color: #880e4f; }
  .number-input {
    width: 56px; padding: 6px 10px; border-radius: 8px;
    border: 1.5px solid #f8bbd0; background: rgba(255,255,255,0.8);
    color: #5a2d4a; font-size: 13px; font-family: 'Courier New', monospace;
    outline: none; transition: border-color 0.2s;
  }
  .number-input:focus { border-color: #e91e96; }
  .run-btn {
    margin-left: auto; padding: 11px 32px; border-radius: 25px; border: none;
    background: linear-gradient(135deg, #e91e96, #ab47bc);
    color: #fff; font-weight: 600; font-size: 14px; cursor: pointer;
    box-shadow: 0 6px 20px rgba(233,30,99,0.4);
    transition: transform 0.15s, box-shadow 0.15s; font-family: 'DM Sans', sans-serif;
  }
  .run-btn:hover { transform: translateY(-2px); box-shadow: 0 8px 25px rgba(233,30,99,0.5); }
  .stat-card {
    flex: 1; background: rgba(255,255,255,0.7); border: 1px solid #f8bbd0;
    border-radius: 12px; padding: 12px 16px;
  }
  .stat-label { font-size: 10px; color: #ad6080; text-transform: uppercase; letter-spacing: 1.5px; }
  .stat-value { font-size: 17px; font-weight: 600; color: #e91e96; margin-top: 3px; font-family: 'Playfair Display', serif; }
  .tab-btn {
    padding: 11px 22px; border: none; cursor: pointer; background: transparent;
    font-size: 13px; font-weight: 500; transition: all 0.2s; font-family: 'DM Sans', sans-serif;
  }
  .tab-btn.active { color: #e91e96; border-bottom: 2px solid #e91e96; }
  .tab-btn:not(.active) { color: #c4a0b0; border-bottom: 2px solid transparent; }
  .content-box {
    background: rgba(255,255,255,0.6); border: 1px solid #f8bbd0;
    border-top: none; border-radius: 0 0 14px 14px; padding: 24px; min-height: 200px;
  }
  .info-box {
    margin-top: 20px; background: rgba(255,255,255,0.5); border: 1px solid #f8bbd0;
    border-radius: 14px; padding: 18px 22px; font-size: 12.5px; color: #ad6080; line-height: 1.8;
  }
  .deco {
    position: fixed; pointer-events: none; opacity: 0.07;
    font-size: 120px; user-select: none;
  }
`;

export default function App() {
  const [seq1, setSeq1] = useState("ACGTTGCATGCA");
  const [seq2, setSeq2] = useState("ACGTCATGCA");
  const [algo, setAlgo] = useState("nw");
  const [match, setMatch] = useState(2);
  const [mismatch, setMismatch] = useState(-1);
  const [gap, setGap] = useState(-2);
  const [result, setResult] = useState(null);
  const [tab, setTab] = useState("alignment");

  const presets = [
    { label: "✿ Similar",   s1: "ACGTTGCATGCA",        s2: "ACGTCATGCA"       },
    { label: "✿ Divergent", s1: "GCTAGCTAGCTA",        s2: "ATGATGATGATG"     },
    { label: "✿ Local hit", s1: "TTTTTGACAGCGGATTT",  s2: "TTAATTGACAGCGGAA" },
    { label: "✿ Protein",   s1: "MKTLLLTLVVVTIVCLDLGAVK", s2: "MKTLLVTIVCLDLGAVKD" },
  ];

  const run = () => {
    const s1 = seq1.toUpperCase().replace(/[^A-Z]/g, "");
    const s2 = seq2.toUpperCase().replace(/[^A-Z]/g, "");
    if (!s1 || !s2) return;
    const fn  = algo === "nw" ? needlemanWunsch : smithWaterman;
    const t0  = performance.now();
    const res = fn(s1, s2, match, mismatch, gap);
    res.elapsed = (performance.now() - t0).toFixed(2);
    res.algo    = algo;
    setResult(res);
    setTab("alignment");
  };

useEffect(() => { run(); }, []);

  return (
    <>
      <style>{styles}</style>
      <div className="app">
        <div className="deco" style={{ top: "5%", left: "-2%" }}>🌸</div>
        <div className="deco" style={{ bottom: "10%", right: "-2%" }}>🌷</div>
        <div className="deco" style={{ top: "40%", right: "3%", fontSize: 80 }}>✿</div>

        <div className="header">
          <div className="logo">🧬</div>
          <div>
            <div style={{ fontFamily: "Playfair Display, serif", fontSize: 24, fontWeight: 600, color: "#880e4f" }}>
              Sequence Alignment Studio
            </div>
            <div style={{ fontSize: 12, color: "#ad6080", marginTop: 3, fontStyle: "italic" }}>
              Needleman-Wunsch · Smith-Waterman · BLAST-like heuristic
            </div>
          </div>
          <div style={{ marginLeft: "auto", display: "flex", gap: 10 }}>
            {["NW", "SW"].map(a => (
              <button key={a} className={`algo-btn ${algo === a.toLowerCase() ? "active" : ""}`}
                onClick={() => setAlgo(a.toLowerCase())}>
                {a === "NW" ? "Global (NW)" : "Local (SW)"}
              </button>
            ))}
          </div>
        </div>

        <div className="main">
          <div style={{ display: "flex", gap: 8, marginBottom: 18, flexWrap: "wrap", alignItems: "center" }}>
            <span style={{ fontSize: 12, color: "#c4a0b0", fontStyle: "italic" }}>try a preset:</span>
            {presets.map(p => (
              <button key={p.label} className="preset-btn"
                onClick={() => { setSeq1(p.s1); setSeq2(p.s2); }}>
                {p.label}
              </button>
            ))}
          </div>

          <div className="card" style={{ marginBottom: 14 }}>
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
              <button onClick={run} className="run-btn">Align ✦</button>
            </div>
          </div>

          {result && (
            <div style={{ display: "flex", gap: 10, marginBottom: 16, flexWrap: "wrap" }}>
              {[
                ["Score",     result.score],
                ["Identity",  (result.identity * 100).toFixed(1) + "%"],
                ["Length",    result.a1.length + " cols"],
                ["Time",      result.elapsed + " ms"],
                ["Algorithm", result.algo === "nw" ? "Needleman-Wunsch" : "Smith-Waterman"],
              ].map(([k, v]) => (
                <div key={k} className="stat-card">
                  <div className="stat-label">{k}</div>
                  <div className="stat-value">{v}</div>
                </div>
              ))}
            </div>
          )}

          <div style={{ display: "flex", borderBottom: "1px solid #f8bbd0" }}>
            {["alignment", "matrix"].map(t => (
              <button key={t} className={`tab-btn ${tab === t ? "active" : ""}`} onClick={() => setTab(t)}>
                {t === "alignment" ? "✦ Alignment View" : "✦ DP Matrix"}
              </button>
            ))}
          </div>

          <div className="content-box">
            {result && tab === "alignment" && <AlignedPair a1={result.a1} a2={result.a2} />}
            {result && tab === "matrix" && (
              <MatrixHeatmap dp={result.dp} seq1={seq1.toUpperCase()} seq2={seq2.toUpperCase()}
                maxI={result.maxI} maxJ={result.maxJ} algorithm={result.algo} />
            )}
            {!result && (
              <div style={{ color: "#d4a0b8", textAlign: "center", paddingTop: 60, fontStyle: "italic" }}>
                ✿ press Align to begin ✿
              </div>
            )}
          </div>

          <div className="info-box">
            {algo === "nw" ? (
              <>
                <strong style={{ color: "#880e4f", fontFamily: "Playfair Display, serif" }}>Needleman-Wunsch (Global)</strong>
                {" "}— fills an (m+1)×(n+1) DP table. Every cell scores max(diagonal+sub, up+gap, left+gap).
                Traceback from bottom-right to top-left. Complexity: <code style={{ color: "#e91e96" }}>O(mn)</code> time & space.
              </>
            ) : (
              <>
                <strong style={{ color: "#880e4f", fontFamily: "Playfair Display, serif" }}>Smith-Waterman (Local)</strong>
                {" "}— same recurrence but scores are floored at 0, allowing alignment to restart anywhere.
                Traceback begins at the highest-scoring cell and stops at 0. Complexity: <code style={{ color: "#e91e96" }}>O(mn)</code> time & space.
              </>
            )}
          </div>
        </div>
      </div>
    </>
  );
}