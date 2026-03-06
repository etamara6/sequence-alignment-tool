
from __future__ import annotations

from flask import Flask, jsonify, request
from flask_cors import CORS

from sequence_alignment import (
    AlignmentResult,
    ScoringScheme,
    needleman_wunsch,
    smith_waterman,
)

app = Flask(__name__)
CORS(app)

_ALGORITHMS = {
    "needleman_wunsch": needleman_wunsch,
    "smith_waterman": smith_waterman,
}

_MAX_LEN = 2_000


@app.get("/api/health")
def health():
    return jsonify({"status": "ok"}), 200


@app.post("/api/align")
def align():
    data = request.get_json(silent=True)
    if not data:
        return jsonify({"error": "Request body must be JSON."}), 400

    seq1 = data.get("seq1", "").strip().upper()
    seq2 = data.get("seq2", "").strip().upper()
    algorithm = data.get("algorithm", "needleman_wunsch")

    if algorithm not in _ALGORITHMS:
        return jsonify({"error": f"Unknown algorithm '{algorithm}'."}), 400

    if len(seq1) > _MAX_LEN or len(seq2) > _MAX_LEN:
        return jsonify({"error": f"Sequences must be ≤ {_MAX_LEN} characters."}), 400

    try:
        scheme = ScoringScheme(
            match=int(data.get("match", 2)),
            mismatch=int(data.get("mismatch", -1)),
            gap_open=int(data.get("gap_open", -2)),
        )
    except (ValueError, TypeError) as exc:
        return jsonify({"error": f"Invalid scoring parameter: {exc}"}), 400

    result: AlignmentResult = _ALGORITHMS[algorithm](seq1, seq2, scheme)

    return jsonify({
        "seq1_aligned": result.seq1_aligned,
        "seq2_aligned": result.seq2_aligned,
        "score":        result.score,
        "identity":     round(result.identity, 4),
        "algorithm":    result.algorithm,
        "elapsed_ms":   round(result.elapsed_ms, 3),
    }), 200


if __name__ == "__main__":
    app.run(debug=True, port=5000)
    return jsonify({"status": "ok"}), 200


@app.post("/api/align")
def align():
    data = request.get_json(silent=True)
    if not data:
        return jsonify({"error": "Request body must be JSON."}), 400

    seq1 = data.get("seq1", "").strip().upper()
    seq2 = data.get("seq2", "").strip().upper()
    algorithm = data.get("algorithm", "needleman_wunsch")

    if algorithm not in _ALGORITHMS:
        return jsonify({"error": f"Unknown algorithm '{algorithm}'."}), 400

    if len(seq1) > _MAX_LEN or len(seq2) > _MAX_LEN:
        return jsonify({"error": f"Sequences must be ≤ {_MAX_LEN} characters."}), 400

    try:
        scheme = ScoringScheme(
            match=int(data.get("match", 2)),
            mismatch=int(data.get("mismatch", -1)),
            gap_open=int(data.get("gap_open", -2)),
        )
    except (ValueError, TypeError) as exc:
        return jsonify({"error": f"Invalid scoring parameter: {exc}"}), 400

    result: AlignmentResult = _ALGORITHMS[algorithm](seq1, seq2, scheme)

    return jsonify({
        "seq1_aligned": result.seq1_aligned,
        "seq2_aligned": result.seq2_aligned,
        "score":        result.score,
        "identity":     round(result.identity, 4),
        "algorithm":    result.algorithm,
        "elapsed_ms":   round(result.elapsed_ms, 3),
    }), 200


if __name__ == "__main__":

    app.run(debug=True, port=5000)
