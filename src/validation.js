/**
 * Validates a nucleotide or protein sequence.
 * Returns { valid: boolean, error: string|null, cleaned: string }
 */
export function validateSequence(seq) {
  // Empty input
  if (!seq || seq.trim().length === 0) {
    return { valid: false, error: "Sequence cannot be empty.", cleaned: "" };
  }

  const cleaned = seq.trim().toUpperCase().replace(/\s+/g, "");

  // Too long — O(m×n) will freeze the browser tab
  if (cleaned.length > 1000) {
    return {
      valid: false,
      error: `Sequence is too long (${cleaned.length} bp). Maximum is 1000 characters.`,
      cleaned,
    };
  }

  // Invalid characters — only DNA bases or standard amino acids allowed
  const validDNA     = /^[ACGTN\-]+$/;
  const validProtein = /^[ACDEFGHIKLMNPQRSTVWY\-]+$/;

  if (!validDNA.test(cleaned) && !validProtein.test(cleaned)) {
    return {
      valid: false,
      error: "Invalid characters detected. Use DNA bases (A, C, G, T) or standard amino acid letters.",
      cleaned,
    };
  }

  return { valid: true, error: null, cleaned };
}