// Utility functions accesible to all scripts

/**
 * Returns true if the string is an accession number
 */
function isAccn(string) {
  const accn_pat = /^EPI_ISL_[0-9]+$/i;  // case-insensitive
  return accn_pat.test(string);
}

/**
 * Retruns true if the string is a lineage 
 */
function isLineage(string) {
  const lin_pat = /[A-Z]\.[0-9]+/i;
  return lin_pat.test(string);
}

/**
 * Returns a string in an ISO8601 format
 *
 * @param {Date} date: The date to be formated
 */
function formatDate(date) {
  return d3.utcFormat("%Y-%m-%d")(date);
}
