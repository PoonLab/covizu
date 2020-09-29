/**
 *  Configure SVG to display beadplots
 */
var marginB = {top: 50, right: 10, bottom: 50, left: 10},
    widthB = 550 - marginB.left - marginB.right,
    heightB = 1000 - marginB.top - marginB.bottom;

// set up plotting scales
var xValueB = function(d) { return d.x },
  xValue1B = function(d) { return d.x1; },
  xValue2B =  function(d) { return d.x2; },
  xScaleB = d3.scaleLinear().range([0, widthB]),
  xMap1B = function(d) { return xScaleB(d.x1); },
  xMap2B = function(d) { return xScaleB(d.x2); },
  xMapB = function(d) { return xScaleB(xValueB(d)); };

var yValueB = function(d) { return d.y; },
  yValue1B = function(d) { return d.y1; },
  yScaleB = d3.scaleLinear().range([0, heightB]),  // inversion
  yMap1B = function(d) { return yScaleB(d.y1); },
  yMap2B = function(d) { return yScaleB(d.y2); },
  yMapB = function(d) { return yScaleB(yValueB(d)); };


var visB = d3.select("div#svg-cluster")
  .append("svg")
  .attr("width", widthB + marginB.left + marginB.right)
  .attr("height", heightB + marginB.top + marginB.bottom)
  .append("g");

var visBaxis = d3.select("div#svg-clusteraxis")
  .append("svg")
  .attr("width", widthB + marginB.left + marginB.right)
  .attr("height", 25)
  .append("g");


// regular expression to remove redundant sequence name components
const pat = /^hCoV-19\/(.+\/.+)\/20[0-9]{2}$/gi;


/**
 * Returns unique elements in given array.
 * @param {Array} arr
 * @returns {string[]}
 */
function unique(arr) {
  var key, history = {};
  for (var i = 0; i< arr.length; i++) {
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
function mode(arr) {
  if (arr.length === 0) {
    return undefined;
  }
  var counts = {},
      key, max_key=arr[0], max_count = 1;
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
  return(max_key);
}


/**
 * Entabulate values in array.
 * @param {Array} arr:  Array of values to entabulate
 * @returns {{}} Associative list of value: count pairs
 */
function table(arr) {
  var val, counts = {};
  for (var i=0; i<arr.length; i++) {
    val = arr[i];
    if (counts[val] === undefined) {
      counts[val] = 0;
    }
    counts[val]++;
  }
  return(counts);
}


/**
 * Parse node and edge data from clusters JSON to a format that is
 * easier to map to SVG.
 * @param {Object} clusters:
 */
function parse_clusters(clusters) {
  var cluster, variant, coldates, samples, regions, country, labels,
      variants,  // horizontal line segment data + labels
      edgelist,  // vertical line segment data
      points,  // the "beads"
      beaddata = [];  // return value

  for (const cidx in clusters) {
    cluster = clusters[cidx];
    if (cluster.nodes.length === 1) {
      console.log('skip '+ cluster.nodes);
      beaddata.push({'variants': [], 'edgelist': [], 'points': []})
      continue;
    }

    // deconvolute edge list to get node list in preorder
    var nodelist = unique(cluster.edges.map(x => x.slice(0,2)).flat());
    //cluster.edges.map(x => x.slice(0,2)).flat().filter(onlyUnique);

    // extract the date range for each variant in cluster
    var y = 1;
    variants = [];
    points = [];
    for (const accn of nodelist) {
      variant = cluster.nodes[accn];
      coldates = variant.map(x => x.coldate);
      coldates.sort();

      country = variant.map(x => x.country);

      // parse samples within variant
      var isodates = unique(coldates),
          isodate;

      variants.push({
        'accession': accn,
        'label': variant[0].label1.replace(pat, "$1"),
        'x1': new Date(coldates[0]),  // min date
        'x2': new Date(coldates[coldates.length-1]),  // max date
        'count': coldates.length,
        'country': country,
        'region': country.map(x => countries[x]),
        'y1': y,  // horizontal line segment
        'y2': y,
        'numBeads': isodates.length,
        'parent': null,
        'dist': 0
      });

      for (var i=0; i<isodates.length; i++) {
        isodate = isodates[i];
        samples = variant.filter(x => x.coldate === isodate);
        country = samples.map(x => x.country);
        regions = country.map(x => countries[x]);

        // warn developers if no region for country
        if (regions.includes(undefined)) {
          console.log("Developer msg, need to update countries.json:");
          for (const j in regions.filter(x => x===undefined)) {
            console.log(`"${samples[j].country}"`);
          }
        }

        points.push({
          cidx,
          'x': new Date(isodate),
          'y': y,
          'count': samples.length,
          'accessions': samples.map(x => x.accession),
          'labels': samples.map(x => x.label1.replace(pat, "$1")),
          'region1': mode(regions),
          'region': regions,
          'country': country,
          'parent': null,
          'dist': 0
        })
      }
      y++;
    }

    // map earliest collection date of child node to vertical edges
    var edge, parent, child, dist;
    edgelist = [];
    for (var e = 0; e < cluster.edges.length; e++) {
      edge = cluster.edges[e];
      parent = variants.filter(x => x.accession === edge[0])[0];
      child = variants.filter(x => x.accession === edge[1])[0];
      dist = parseFloat(edge[2]);
      edgelist.push({
        'y1': parent.y1,
        'y2': child.y1,
        'x1': child.x1,  // vertical line segment
        'x2': child.x1,
        'parent': parent.label,
        'child': child.label,
        'dist': dist
      });

      // Assign parent and genomic distance of each variant
      let childvariants = variants.filter(x => x.y1 === child.y1);
      for (let v = 0; v < childvariants.length; v++) {
        childvariants[v].parent = parent.label;
        childvariants[v].dist = dist;
      }

      // Assign the parent and genomic distance of each point
      let childpoints = points.filter(x => x.y === child.y1);
      for (let c = 0; c < childpoints.length; c++) {
        childpoints[c].parent = parent.label;
        childpoints[c].dist = dist;
      }

      // update variant time range
      if (parent.x1 > child.x1) {
        parent.x1 = child.x1;
      }
      if (parent.x2 < child.x1) {
        parent.x2 = child.x1;
      }
    }

    // calculate consensus region for cluster
    beaddata.push({
      'variants': variants,
      'edgelist': edgelist,
      'points': points
    });

    // collect all region Arrays for all samples, all variants
    regions = points.map(x => x.region).flat();
    cluster['region'] = mode(regions);
    cluster['allregions'] = table(regions);

    // concatenate all sample labels within cluster for searching
    labels = points.map(x => x.labels).flat();
    cluster['searchtext'] = labels.join();
    cluster['label1'] = labels[0];

    // collect all countries
    cluster['country'] = variants.map(x => x.country).flat();
  }

  return(beaddata);
}

function create_selection(selected_obj) {
  d3.select("div#svg-cluster").selectAll("line").attr("stroke-opacity", 0.3);
  d3.select("div#svg-cluster")
    .selectAll("circle:not(.SelectedBead):not(.selectionH)")
    .attr("class", "not_SelectedBead");
  selected_obj.attr("class","SelectedBead");
  d3.selectAll("circle.SelectedBead").each(function(r) {
    visB.append("circle").lower()
    .attr('class', "selectionH")
    .attr("cx", xScaleB(xValueB(r)))
    .attr("cy", yScaleB(yValueB(r)))
    .attr("r", 4*Math.sqrt(r.count)+4)
    .attr("fill", "none")
    .attr("stroke", "grey")
    .attr("fill-opacity", 1)
    .attr("stroke-width", 5);
  });
}

function clear_selection() {
  $('#search-input').val('');
  d3.select("#svg-cluster").selectAll("line").attr("stroke-opacity", 1);
  d3.selectAll("circle:not(.selectionH)").attr("class", "default");
  d3.select('#svg-timetree').selectAll("rect:not(.clicked):not(.clickedH)").attr("class", "default");
  d3.selectAll("circle.selectionH").remove();
}
/**
 * Draw the beadplot for a specific cluster (identified through its
 * integer index <cid>) in the SVG.
 * @param {Number} cid:  integer index of cluster to draw as beadplot
 */
function beadplot(cid) {

  var variants = beaddata[cid].variants,
      edgelist = beaddata[cid].edgelist,
      points = beaddata[cid].points;

  // set plotting domain
  var mindate = d3.min(variants, xValue1B),
      maxdate = d3.max(variants, xValue2B),
      spandate = maxdate-mindate,  // in milliseconds
      min_y = d3.min(variants, yValue1B),
      max_y = d3.max(variants, yValue1B);

  // Create a div for the tooltip
  let bTooltip = d3.select("#tooltipContainer")
      .style("opacity", 0);

  // update vertical range for consistent spacing between variants
  heightB = max_y * 10;
  $("#svg-cluster > svg").attr("height", heightB + marginB.top + marginB.bottom);
  yScaleB = d3.scaleLinear().range([40, heightB]);
  xScaleB.domain([mindate - 0.5*spandate, maxdate]);
  yScaleB.domain([min_y, max_y]);

  // clear SVG
  visB.selectAll('*').remove();
  visBaxis.selectAll('*').remove();

  // draw horizontal line segments that represent variants in cluster
  visB.selectAll("lines")
      .data(variants)
      .enter().append("line")
      .attr(  "class", "lines")
      .attr("x1", xMap1B)
      .attr("x2", xMap2B)
      .attr("y1", yMap1B)
      .attr("y2", yMap2B)
      .attr("stroke-width", 3)
      .attr("stroke", "#777")
      .on("mouseover", function(d) {
        d3.select(this)
            .attr("stroke-width", 5);

        bTooltip.transition()       // Show tooltip
            .duration(50)
            .style("opacity", 0.9);

        let tooltipText = "";
        if (d.parent || d.dist) {
          tooltipText += `<b>Parent:</b> ${d.parent}<br><b>Genomic distance:</b> ${d.dist}<br><br>`;
        }

        tooltipText += region_to_string(table(d.region));
        tooltipText += `<br><b>Unique collection dates:</b> ${d.numBeads}<br>`;
        let formatDate = d3.timeFormat("%Y-%m-%d");
        tooltipText += `<br><b>Collection dates:</b><br>${formatDate(new Date(d.x1))} / ${formatDate(new Date(d.x2))}<br>`;

        // Tooltip appears 10 pixels left of the cursor
        bTooltip.html(tooltipText)
            .style("left", (d3.event.pageX + 10) + "px")
            .style("top", (d3.event.pageY + "px"));

      })
      .on("mouseout", function() {
        d3.select(this)
            .attr("stroke-width", 3);
        bTooltip.transition()     // Hide tooltip
            .duration(50)
            .style("opacity", 0);
      })
      .on("click", function(d) {
        gentable(d);
        draw_region_distribution(table(d.region));
        let var_samples = points.filter(x => x.y === d.y1);
        gen_details_table(var_samples);
      });

  // label variants with earliest sample name
  visB.selectAll("text")
      .data(variants)
      .enter().append("text")
      .style("font-size", "10px")
      .attr("text-anchor", "end")
      .attr("alignment-baseline", "middle")
      .attr("x", function(d) { return(xScaleB(d.x1)-5); })
      .attr("y", function(d) { return(yScaleB(d.y1)); })
      .text(function(d) { return(d.label); });

  // draw vertical line segments that represent edges in minimum spanning tree
  visB.selectAll("lines")
      .data(edgelist)
      .enter().append("line")
      .attr("class", "lines")
      .attr("x1", xMap1B)
      .attr("x2", xMap2B)
      .attr("y1", yMap1B)
      .attr("y2", yMap2B)
      .attr("stroke-width", 1)
      .attr("stroke", function(d) {
        if (d.dist < 1.5) {
          return("#55b7");
        } else if (d.dist < 2.5) {
          return("#77d7");
        } else {
          return("#99f7");
        }
      })
      .on("mouseover", function(d) {
        d3.select(this).attr("stroke-width", 3);

        bTooltip.transition()       // Show tooltip
            .duration(50)
            .style("opacity", 0.9);

        let tooltipText = `<b>Parent:</b> ${d.parent}<br><b>Child:</b> ${d.child}<br>`;
        tooltipText += `<b>Genomic distance:</b> ${d.dist}<br><br>`;

        let formatDate = d3.timeFormat("%Y-%m-%d");
        tooltipText += `<b>Collection date:</b> ${formatDate(new Date(d.x2))}`;

        // Tooltip appears 10 pixels left of the cursor
        bTooltip.html(tooltipText)
            .style("left", (d3.event.pageX + 10) + "px")
            .style("top", (d3.event.pageY + "px"));

      })
      .on("mouseout", function() {
        d3.select(this).attr("stroke-width", 1);

        bTooltip.transition()     // Hide tooltip
            .duration(50)
            .style("opacity", 0);
      })
      .on("click", function(d) {
        $("#text-node").text(`Parent: ${d.parent}\nChild: ${d.child}\nGenomic distance: ${d.dist}`);
      });

  // draw "beads" to represent samples per collection date
  visB.selectAll("circle")
      .data(points)
      .enter().append("circle")
      .attr("r", function(d) { return (4*Math.sqrt(d.count)); })
      .attr("cx", xMapB)
      .attr("cy", yMapB)
      .attr("class", "default")
      .attr("fill", function(d) {
        return(country_pal[d.region1]);
      })
      .attr("stroke", function(d) {
        return(country_pal[d.region1]);
      })
      .on("mouseover", function(d) {
        d3.select(this).attr("stroke-width", 2)
            .attr("r", 4*Math.sqrt(d.count)+3);

        bTooltip.transition()       // Show tooltip
            .duration(50)
            .style("opacity", 0.9);

        let tooltipText = "";
        if (d.parent || d.dist) {
          tooltipText += `<b>Parent:</b> ${d.parent}<br><b>Genomic distance:</b> ${d.dist}<br><br>`;
        }

        tooltipText += region_to_string(table(d.region));
        let formatDate = d3.timeFormat("%Y-%m-%d");
        tooltipText += `<br><b>Collection date:</b> ${formatDate(new Date(d.x))}<br>`;

        // Tooltip appears 10 pixels left of the cursor
        bTooltip.html(tooltipText)
            .style("left", (d3.event.pageX + 10) + "px")
            .style("top", (d3.event.pageY + "px"));
      })
      .on("mouseout", function(d) {
          d3.select(this).attr("stroke-width", 1)
              .attr("r", 4*Math.sqrt(d.count));
          bTooltip.transition()     // Hide tooltip
              .duration(50)
              .style("opacity", 0);
      })
      .on("click", function(d) {

        clear_selection();

        var cur_obj = d3.select(this);
        create_selection(cur_obj);

        d3.select("#svg-timetree").selectAll("line").attr("stroke-opacity", 1);

        gentable(d);
        draw_region_distribution(table(d.region));
        gen_details_table(d);
      });

  // draw x-axis
  visBaxis.append("g")
      .attr("transform", "translate(0,20)")
      .call(
        d3.axisTop(xScaleB)
          .ticks(5)
          .tickFormat(d3.timeFormat("%Y-%m-%d"))
      );

}

/**
 * Writes information about the region distribution to a string
 * @param {{}} my_regions: associative list of region and case count pairs
 * @returns {string} regStr: a string representation of the region distribution
 */
function region_to_string(my_regions) {
  // Display region distribution in tooltip
  let regStr = `<b>Number of cases</b><br>`,
      total = 0;

  for (let [r_key, r_value] of Object.entries(my_regions)) {
    regStr += `${r_key}: ${r_value}<br>`;
    total += r_value;
  }

  // Display total number of cases if variants are from multiple countries
  if(Object.keys(my_regions).length > 1) {
    regStr += `Total: ${total}<br>`
  }

  return regStr;
}

/**
 * Creates a table that displays sequence details (Sequence name, GISAID accession number, collection date)
 * @param obj: JS object or an array of JS Objects
 */
function gen_details_table(obj) {
  var details = [];
  let formatDate = d3.timeFormat("%Y-%m-%d");

  // Check for a list of samples
  if (Array.isArray(obj)) {

    // Iterate over each sample in the list
    for (let j = 0; j < obj.length; j++) {

      // "zip" the sequence details of each sample
      for (let i = 0; i < obj[j].accessions.length; i++) {
        let sample_details = [obj[j].accessions[i], obj[j].labels[i], formatDate(new Date(obj[j].x))];
        details.push(sample_details);
      }
    }
  }

  else {
    // "zip" the sequence details of each sample
    for (let i = 0; i < obj.accessions.length; i++) {
      let sample_details = [obj.accessions[i], obj.labels[i], formatDate(new Date(obj.x))];
      details.push(sample_details);
    }
  }


  thead.html(""); // clear table headers
  var sort_ascending = true;

  var headers = thead.append('tr')
      .selectAll('th')
      .data(seq_theaders)
      .enter()
      .append('th')
      .on('click', function (x, i) {
        // Reset sorting arrows
        thead.selectAll('a').classed('hide-before', false);
        thead.selectAll('a').classed('hide-after', false);

        // Sort columns
        if (sort_ascending) {
          t_rows.sort(function(a, b) { return d3.ascending(b[i], a[i]); });
          sort_ascending = false;
          d3.select(this).select('a').classed('hide-after', true);
          d3.select(this).select('a').classed('hide-before', false);
        } else {
          t_rows.sort(function(a, b) { return d3.descending(b[i], a[i]); });
          sort_ascending = true;
          d3.select(this).select('a').classed('hide-before', true);
          d3.select(this).select('a').classed('hide-after', false);
        }
      })
      .append('a')
      .text(function (x) { return x; })
      .classed('sort-by', true);

  // Create a row for each sample
  seq_tbody.html("");
  var t_rows = seq_tbody.selectAll('tr')
      .data(details)
      .enter()
      .append('tr')
      .on("mouseover", function () {
        d3.select(this).style("background-color", "grey");  // Highlight on mouseover
      })
      .on("mouseout", function () {            // Remove highlighting on mouseout
        d3.select(this).style("background-color", null);
      });

  // Create a cell for every row in the column
  t_rows.selectAll('td')
      .data(function (r) { return r; })
      .enter()
      .append('td')
      .text(function (x) { return x; })
      .style("font", "0.875em/1.2 Lato, sans-serif");
}


/**
 * Function to generate table from my_countries object on bead click
 * @param {Object} obj:  JS Object with country attribute
 */
function gentable(obj) {
  var my_countries = Object.entries(table(obj.country)),
      row, region;

  // annotate with region (continent)
  for (const i in my_countries) {
    row = my_countries[i];
    region = countries[row[0]];
    my_countries[i] = [region].concat(row);
  }

  country_table.selectAll('thead').remove()

  // https://stackoverflow.com/questions/32871044/how-to-update-d3-table
  //create a row <tr> for each item in my_countries array, set mouse over over color, then create one column <td> for each item in array
  var rows = country_tbody.selectAll("tr")
    .data(my_countries);
  rows.exit().remove();
  rows.enter()
    .append("tr")
    .on("mouseover", function(){
      d3.select(this).style("background-color", "grey");}) //mouseover highlight
    .on("mouseout", function(){
      d3.select(this).style("background-color", null);})  //mouseout unhighlight
    .selectAll("td")
    .data(function(d) { return d; })
    .enter()
    .append("td")
    .text(function(d) { return d; });

  var rows = country_tbody.selectAll("tr") //reselect all rows

  var cells = rows.selectAll("td")
    .data(function(d) { return d; })
    cells.text(function(d) { return d; })

  var sort_ascending = true;

  var headers = country_table.append("thead").append("tr")
      .selectAll("th")
      .data(theaders)
      .enter()
      .append("th")
      .on('click', function (x, i) {
        // Reset sorting arrows
        country_table.select('thead').selectAll('a').classed('hide-before', false);
        country_table.select('thead').selectAll('a').classed('hide-after', false);

        // Sort columns
        if (sort_ascending) {
          rows.sort(function(a, b) { return d3.ascending(b[i], a[i]); });
          sort_ascending = false;
          d3.select(this).select('a').classed('hide-after', true);
          d3.select(this).select('a').classed('hide-before', false);
        } else {
          rows.sort(function(a, b) { return d3.descending(b[i], a[i]); });
          sort_ascending = true;
          d3.select(this).select('a').classed('hide-before', true);
          d3.select(this).select('a').classed('hide-after', false);
        }
      })
      .append('a')
      .text(function (x) { return x; })
      .classed('sort-by', true);
}


/**
 * Draws a bar chart of the distribution of cases across regions
 * @param {{}} my_regions: associative list of region and case count pairs
 */
function draw_region_distribution(my_regions) {
  const regions = unique(Object.values(countries)).sort();
  var counts = [], count;

  regions.forEach(function(r) {
    count = my_regions[r];
    if (count === undefined) count = 0;
    counts.push({'region': r, 'count': count})
  });

  // Set the margins
  const margin = {top: 15, right: 10, bottom: 75, left: 55},
      width = 250 - margin.right - margin.left,
      height = 200 - margin.top - margin.bottom;

  // Create the barchart
  const svg = d3.select("#barplot")
      .html("")
      .append("svg")
      .attr("width", 250)
      .attr("height", 200);

  const chart = svg.append("g")
      .attr("transform", `translate(${margin.left}, ${margin.top})`);

  // Set the scale of the x-axis
  const xScale = d3.scaleBand()
      .range([0, width])
      .domain(counts.map((r) => r.region))
      .padding(0.1);

  chart.append("g")
      .attr("transform", `translate(0, ${height})`)
      .call(d3.axisBottom(xScale))
      .selectAll("text")
      .style("text-anchor", "end")
      .attr("transform", "rotate(-45)");

  // Set the scale of the y-axis
  const max_count = counts.map(x=>x.count).reduce(
      function(a,b) { return Math.max(a,b) }
      );

  const yScale = d3.scaleLinear()
      .range([height, 0])
      .domain([0, (max_count < 5) ? 5 : max_count])
      .nice();

  chart.append("g")
       .call(d3.axisLeft(yScale).ticks(3));

  // Create the chart
  const regionBars = chart.selectAll()
      .data(counts)
      .enter()
      .append("g");

  // Draw bars on the chart
  regionBars
      .append("rect")
      .attr("x", (r) => xScale(r.region))
      .attr("y", (r) => yScale(r.count))
      .attr("height", (r) => height - yScale(r.count))
      .attr("width", xScale.bandwidth())
      .attr("fill", function(d) { return(country_pal[d.region]); });

    // Write the case count above each bar
    regionBars.append("text")
        .style("font", "0.7em/1.2 Lato, sans-serif")
        .attr("x", (r) => xScale(r.region) + xScale.bandwidth() / 2)
        .attr("y", (r) => yScale(r.count) - 5)
        .attr("text-anchor", "middle")
        .text((r) => `${r.count}`);

    // Add axis labels
    svg.append("text")
        .style("font", "0.8em/1.2 Lato, sans-serif")
        .attr("transform", "rotate(-90)")
        .attr("x", -(height / 2) - margin.top)
        .attr("y", 0)
        .attr("dy", "1em")
        .attr("text-anchor", "middle")
        .text("Number of Cases");

    svg.append("text")
        .style("font", "0.8em/1.2 Lato, sans-serif")
        .attr("text-anchor", "middle")
        .attr("x", (width / 2) + margin.left)
        .attr("y", height + 85)
        .text("Region");
}
