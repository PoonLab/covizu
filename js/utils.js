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

/**
 * Returns true if the date is in the correct format (YYYY-MM-DD)
 * @param {String} date 
 */
function isDate(date) {
  const date_pat = /^\d{4}\-(0?[1-9]|1[012])\-(0?[1-9]|[12][0-9]|3[01])$/;
  return date_pat.test(date);
}


/**
 *
 * Retrives sParam if it exists in url
 * @param {string} sParam
*/
function getUrlParameter(sParam) {
    var sPageURL = window.location.search.substring(1),
        sURLVariables = sPageURL.split('&'),
        sParameterName,
        i;

    for (i = 0; i < sURLVariables.length; i++) {
        sParameterName = sURLVariables[i].split('=');

        if (sParameterName[0] === sParam) {
            return typeof sParameterName[1] === undefined ? true : decodeURIComponent(sParameterName[1]);
        }
    }
    return false;
};


/**
/* Returns encoded URI string from JSON object
/*
*/
function jsonToURI(json){ return encodeURIComponent(JSON.stringify(json)); }

