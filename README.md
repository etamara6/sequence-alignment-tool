
# 🧬 Sequence Alignment Studio

![CI](https://github.com/etamara6/sequence-alignment-tool/actions/workflows/ci.yml/badge.svg)
[![Python](https://img.shields.io/badge/python-3.8%2B-3776ab?logo=python&logoColor=white)](https://www.python.org/)
[![React](https://img.shields.io/badge/react-18-61dafb?logo=react&logoColor=black)](https://reactjs.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)

> An interactive, full-stack genomic sequence alignment tool — visualise Needleman-Wunsch and Smith-Waterman in real time, powered by a Python/Flask backend and a React frontend.

---

[![Live Demo](https://img.shields.io/badge/live%20demo-UI-46E3B7?logo=render)](https://sequence-alignment-ui.onrender.com)
[![API](https://img.shields.io/badge/api-live-46E3B7?logo=render)](https://sequence-alignment-api.onrender.com/api/health)

---

## 🧠 What This Project Demonstrates

| Skill | How it shows up |
|---|---|
| Algorithm implementation | NW & SW from scratch in both Python (NumPy) and JavaScript (TypedArrays) |
| Full-stack integration | React ↔ Flask REST API, CORS, JSON contract |
| Software engineering | Dataclasses, type hints, edge-case handling, clean separation of concerns |
| Testing discipline | pytest + pytest-cov on Python; Jest on JS; CI on every push |
| Bioinformatics knowledge | BLOSUM62, gap penalties, traceback, identity calculation |

---

## 🔬 Algorithm Deep Dive

### Needleman-Wunsch — Global Alignment

Globally aligns two sequences end-to-end. Optimal when sequences are of similar length and believed to be homologous throughout.

**Recurrence:**
```
dp[i][j] = max(
    dp[i-1][j-1] + sub(seq1[i], seq2[j]),   # match/mismatch
    dp[i-1][j]   + gap,                       # deletion
    dp[i][j-1]   + gap                        # insertion
)
```
Initialisation: `dp[i][0] = i × gap`, `dp[0][j] = j × gap`  
Traceback: from `dp[m][n]` → `dp[0][0]`

### Smith-Waterman — Local Alignment

Finds the highest-scoring local sub-alignment. Best for spotting conserved domains between divergent sequences.

**Key difference from NW:**
```
dp[i][j] = max(0, diag, up, left)   # floor at zero
```
Traceback: from the **maximum cell** → stops at first `0`

### Scoring Parameters

| Parameter | Default | Effect |
|---|---|---|
| Match | `+2` | Reward for identical residues |
| Mismatch | `-1` | Penalty for substitution |
| Gap open | `-2` | Penalty per gap character |

---

## 📥 Example

**Input**
```
Sequence 1:  ACGTTGCATGCA
Sequence 2:  ACGTCATGCA
Algorithm:   Needleman-Wunsch
```

**Output**
```
Algorithm : Needleman-Wunsch (Global)
Score     : 16
Identity  : 83.3%
Time      : 0.41 ms

  ACGTTGCATGCA
  |||| |||||||
  ACGT-CATGCA-
```

---

## 🏗️ Architecture

```
┌──────���──────────────────┐        HTTP/JSON        ┌──────────────────────────┐
│   React Frontend        │ ──── POST /api/align ──▶ │   Flask Backend          │
│   src/App.js            │ ◀─── AlignmentResult ─── │   api.py                 │
│   src/algorithms.js     │                          │   sequence_alignment.py  │
│   (JS fallback impl.)   │                          │   (NumPy-powered)        │
└─────────────────────────┘                          └──────────────────────────┘
```

- The **React app** handles all UI, input validation, and DP matrix heatmap rendering.
- The **Flask API** exposes `POST /api/align` — the Python implementation runs on the server.
- The **JS implementation** in `src/algorithms.js` serves as a client-side fallback and is independently tested with Jest.

---

## 🚀 Running Locally

### Prerequisites
- Node.js ≥ 16 & npm ≥ 8
- Python 3.8+

### Step 1 — Clone
```bash
git clone https://github.com/etamara6/sequence-alignment-tool.git
cd sequence-alignment-tool
```

### Step 2 — Start the Python API
```bash
pip install -r requirements.txt
python api.py
# → http://localhost:5000
```

### Step 3 — Start the React app
```bash
npm install
npm start
# → http://localhost:3000  (proxies /api/* to Flask automatically)
```

### Step 4 — (Optional) Configure API URL

By default the React app proxies `/api/*` to `http://localhost:5000` via the `proxy` field in `package.json`.  
To point to a different backend (e.g. a deployed API), set the environment variable before starting:

```bash
REACT_APP_API_URL=https://your-api-url.onrender.com npm start
```

### Step 5 — (Optional) Python CLI demo
```bash
python sequence_alignment.py
```

---

## 🧪 Testing

### Python
```bash
pytest --cov=. --cov-report=term-missing
```

### JavaScript
```bash
npm test
```

CI runs both suites automatically on every push via GitHub Actions.

---

## 📁 Project Structure

```
sequence-alignment-tool/
├── .github/workflows/ci.yml        # CI — JS + Python, with coverage
├── public/                         # Static HTML shell
├── src/
│   ├── App.js                      # Main React component + UI logic
│   ├── algorithms.js               # JS NW + SW (TypedArrays)
│   ├── algorithms.test.js          # Jest tests for JS algorithms
│   ├── validation.js               # Sequence input validation
│   └── App.test.js
├── tests/
│   └── test_sequence_alignment.py  # pytest suite — NW, SW, edge cases
├── api.py                          # Flask REST API
├── sequence_alignment.py           # Python NW + SW + BLAST-like engine
├── requirements.txt
└── package.json
```

---

## 📄 License

MIT © [etamara6](https://github.com/etamara6)
