from __future__ import annotations

import os

from flask import Flask, jsonify, request
from flask_cors import CORS

from sequence_alignment import (
    AlignmentResult,
    MSAResult,
    ScoringScheme,
    needleman_wunsch,
    smith_waterman,
    progressive_msa,
)

app = Flask(__name__)
CORS(app)

_ALGORITHMS = {
    "needleman_wunsch": needleman_wunsch,
    "smith_waterman": smith_waterman,
}

_MAX_LEN = 2_000
_MAX_SEQS = 10   # safety cap for MSA


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
        return jsonify({"error": f"Sequences must be <= {_MAX_LEN} characters."}), 400

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


@app.post("/api/align/msa")
def align_msa():
    """
    Progressive Multiple Sequence Alignment endpoint.

    Expected JSON body:
    {
        "sequences": ["ACGT...", "ACGT...", "ACGT..."],   // 2–10 sequences
        "labels":    ["Seq1", "Seq2", "Seq3"],             // optional
        "match":     2,                                    // optional
        "mismatch":  -1,                                   // optional
        "gap_open":  -2                                    // optional
    }

    Response:
    {
        "aligned":       ["ACGT-...", "A-GTAC...", ...],
        "labels":        ["Seq1", "Seq2", ...],
        "guide_tree":    [[0, 1, 0.25], [0, 2, 0.40]],
        "avg_identity":  0.82,
        "algorithm":     "Progressive MSA (NW guide tree)",
        "elapsed_ms":    12.3
    }
    """
    data = request.get_json(silent=True)
    if not data:
        return jsonify({"error": "Request body must be JSON."}), 400

    raw_seqs = data.get("sequences", [])
    if not isinstance(raw_seqs, list) or len(raw_seqs) < 2:
        return jsonify({"error": "Provide at least 2 sequences in 'sequences'."}), 400
    if len(raw_seqs) > _MAX_SEQS:
        return jsonify({"error": f"At most {_MAX_SEQS} sequences are supported."}), 400

    sequences = [str(s).strip().upper() for s in raw_seqs]
    for s in sequences:
        if len(s) > _MAX_LEN:
            return jsonify({"error": f"Each sequence must be <= {_MAX_LEN} characters."}), 400
        if not s:
            return jsonify({"error": "Empty sequences are not allowed."}), 400

    labels = data.get("labels", None)
    if labels is not None:
        if not isinstance(labels, list) or len(labels) != len(sequences):
            return jsonify({"error": "'labels' must be a list matching the number of sequences."}), 400
        labels = [str(l) for l in labels]

    try:
        scheme = ScoringScheme(
            match=int(data.get("match", 2)),
            mismatch=int(data.get("mismatch", -1)),
            gap_open=int(data.get("gap_open", -2)),
        )
    except (ValueError, TypeError) as exc:
        return jsonify({"error": f"Invalid scoring parameter: {exc}"}), 400

    result: MSAResult = progressive_msa(sequences, labels=labels, scheme=scheme)

    return jsonify({
        "aligned":       result.aligned,
        "labels":        result.labels,
        "guide_tree":    [[i, j, round(d, 4)] for i, j, d in result.guide_tree],
        "avg_identity":  round(result.avg_identity, 4),
        "algorithm":     result.algorithm,
        "elapsed_ms":    round(result.elapsed_ms, 3),
    }), 200


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=int(os.environ.get("PORT", 5000)))



