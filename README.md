
# 🧬 Sequence Alignment Tool

An interactive genomic sequence alignment tool implementing three core bioinformatics algorithms from scratch, with a React visualizer frontend.

---

## ✦ Algorithms Implemented

- **Needleman-Wunsch** — global pairwise alignment (exact, O(mn))
- **Smith-Waterman** — local pairwise alignment (exact, O(mn))
- **BLAST-like heuristic** — seed-and-extend approximate search (fast, ~O(n))

All algorithms are implemented from scratch with NumPy anti-diagonal vectorisation for performance optimization.

---

## ✦ Project Structure

```
sequence-alignment-tool/
│
├── sequence_alignment.py       # Python backend — all three algorithms
│
├── src/
│   └── App.js                  # React visualizer (interactive frontend)
│
├── public/
│   └── index.html              # HTML entry point
│
├── .gitignore                  # Git ignore rules
├── README.md                   # This file
├── package.json                # React dependencies
└── package-lock.json           # Locked dependency versions
└── requirements.txt
└── LICENSE   
```

---

## ✦ How to Run

### Python (algorithms only)
Make sure you have Python 3.9+ and NumPy installed:
```
pip install numpy
python sequence_alignment.py
```

### React Visualizer (interactive UI)
Make sure you have Node.js installed:
```
npm install
npm start
```
Then open your browser at `http://localhost:3001`

---

## ✦ Features

- Interactive alignment of custom DNA/protein sequences
- Toggle between Global (NW) and Local (SW) alignment
- Adjustable scoring parameters (match, mismatch, gap)
- Color-coded alignment view with match ♡ indicators
- DP matrix heatmap visualization
- Preset example sequences to try
- Performance timing for each alignment

---

## ✦ Built With

- Python 3.9+
- NumPy
- React (Create React App)
