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
  // const lin_pat = /[A-Z]\.[0-9]+/i;
  // return lin_pat.test(string);
  return lineage_to_cid[string.toUpperCase()] != undefined;
}

/**
 * Returns a string in an ISO8601 format
 *
 * @param {Date} date: The date to be formated
 */
function formatDate(date) {
  return d3.utcFormat("%Y-%m-%d")(date);
}

/**
 * Returns true if the date is in the correct format (YYYY-MM-DD)
 * @param {String} date 
 */
function isDate(date) {
  const date_pat = /^\d{4}\-(0?[1-9]|1[012])\-(0?[1-9]|[12][0-9]|3[01])$/;
  return date_pat.test(date);
}

/**
 * Returns a date in UTC 
 *
 * @param {String} date: The date to be converted
 */
function utcDate(date) {
   const dateObj = new Date(date);
   return new Date(dateObj.getTime() + dateObj.getTimezoneOffset() * 60000);
}
