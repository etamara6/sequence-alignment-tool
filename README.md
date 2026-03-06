# Sequence Alignment Studio

A high-performance genomic sequence alignment tool built to demonstrate modern web development and algorithmic skills. This application provides a rich user interface for running both global (Needleman-Wunsch) and local (Smith-Waterman) alignments. It features dual implementations in Python and JavaScript, showcasing versatility across different technology stacks.


## 🌟 Features

*   **Dual Algorithm Support:** Implements both Needleman-Wunsch for global alignments and Smith-Waterman for local alignments.
*   **Interactive UI:** A polished and responsive interface built with React for setting sequences and scoring parameters.
*   **Dynamic Programming Matrix Visualization:** Displays the complete DP matrix as a heatmap, providing insight into the algorithm's process.
*   **BLOSUM62 Scoring:** Includes support for the BLOSUM62 substitution matrix for protein sequence alignments.
*   **Performance-Optimized:** The JavaScript implementation uses TypedArrays for memory efficiency, while the Python version uses NumPy for vectorized calculations.
*   **Text Export:** Results, including the alignment, score, and identity, can be exported to a `.txt` file.

## 🛠️ Tech Stack

*   **Frontend:** React.js, HTML5, CSS3
*   **Backend/Algorithm Logic:**
    *   JavaScript (ES6+)
    *   Python 3 with NumPy
*   **Testing:** Jest, React Testing Library
*   **CI/CD:** GitHub Actions

## 🚀 Getting Started

### Prerequisites

*   Node.js (v16 or later)
*   npm (v8 or later)

### Running the React Application

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/etamara6/sequence-alignment-tool.git
    cd sequence-alignment-tool
    ```

2.  **Install dependencies:**
    ```bash
    npm install
    ```

3.  **Start the development server:**
    The application will be available at `http://localhost:3000`.
    ```bash
    npm start
    ```

### Running the Python Implementation

The Python script `sequence_alignment.py` contains the backend logic and a command-line demo.

1.  **Prerequisites:**
    *   Python 3.8+
    *   NumPy

2.  **Install NumPy:**
    ```bash
    pip install numpy
    ```

3.  **Run the demo:**
    ```bash
    python sequence_alignment.py
    ```

## 🧪 Testing

The project is configured with Jest and React Testing Library. To run the test suite:

```bash
npm test
```

All tests run automatically on every push to the `main` branch via the GitHub Actions CI workflow.
