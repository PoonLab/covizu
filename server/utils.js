/**
 * Returns unique elements in given array.
 * @param {Array} arr
 * @returns {string[]}
 */
const unique = (arr) => {
    var key, history = {};
    for (var i = 0; i < arr.length; i++) {
        key = arr[i];
        if (history[key] === undefined) {
            history[key] = 1;
        }
    }
    return (Object.keys(history));
}


/**
 * Returns most common element of array.  If there is a tie, then
 * the function returns the right-most value.
 * @param {Array} arr:  array of elements to sort
 * @return most common element
 */
const mode = (arr) => {
    if (arr.length === 0) {
        return undefined;
    }
    var counts = {},
        key, max_key = arr[0], max_count = 1;
    for (var i = 0; i < arr.length; i++) {
        key = arr[i];
        if (counts[key] == null) {
            counts[key] = 1
            continue
        }
        counts[key]++;
        if (counts[key] > max_count) {
            max_count = counts[key];
            max_key = key;
        }
    }
    return (max_key);
}




/**
 * Merge the counts for two or more tables.
 * @param {Array} tables: Objects with integer values to be combined by key.
 * @returns {Object}
 */
function merge_tables(tables) {
    var total = {};
    for (tab of tables) {
        if (tab === null) {
            continue;
        }
        for (key of Object.keys(tab)) {
            if (total[key] === undefined) {
                total[key] = 0;
            }
            total[key] += tab[key];
        }
    }
    return (total);
}

/**
 * Tabulate values in array.
 * @param {Array} arr:  Array of values to tabulate
 * @returns {{}} Associative list of unique value: count pairs
 */
function tabulate(arr) {
    var val, counts = {};
    for (var i = 0; i < arr.length; i++) {
        val = arr[i];
        if (val === null) {
            continue;
        }
        if (counts[val] === undefined) {
            counts[val] = 0;
        }
        counts[val]++;
    }
    return (counts);
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

module.exports = {
    unique,
    mode,
    tabulate,
    merge_tables,
    utcDate
};