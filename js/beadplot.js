/**
 *  Configure SVG to display beadplots
 */
var marginB = {top: 50, right: 10, bottom: 50, left: 10},
    widthB = 600 - marginB.left - marginB.right,
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
  var cluster, variant, coldates, samples, regions, country,
      variants,  // horizontal line segment data + labels
      edgelist,  // vertical line segment data
      points,  // the "beads"
      beaddata = [];  // return value

  for (var cidx = 0; cidx < clusters.length; cidx++) {
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

      variants.push({
        'accession': accn,
        'label': variant[0].label1.replace(pat, "$1"),
        'x1': new Date(coldates[0]),  // min date
        'x2': new Date(coldates[coldates.length-1]),  // max date
        'count': coldates.length,
        'country': table(variant.map(x => x.country)),
        'y1': y,  // horizontal line segment
        'y2': y
      });

      var isodates = unique(coldates),
          isodate;

      for (var i=0; i<isodates.length; i++) {
        isodate = isodates[i];
        samples = variant.filter(x => x.coldate === isodate);
        country = samples.map(x => x.country);
        regions = country.map(x => countries[x]);
        
        // warn developers if no region for country
        if (regions.includes(undefined)) {
          console.log("Developer msg, need to update countries.json:")
          for (const j in regions.filter(x => x===undefined)) {
            console.log(samples[j].country);
          }
        }
        
        points.push({
          'x': new Date(isodate),
          'y': y,
          'count': samples.length,
          'accessions': samples.map(x => x.accession),
          'labels': samples.map(x => x.label1.replace(pat, "$1")),
          'region1': mode(regions),
          'region': regions,
          'country': country
        })
      }
      y++;
    }

    // map earliest collection date of child node to vertical edges
    var edge, parent, child;
    edgelist = [];
    for (var e = 0; e < cluster.edges.length; e++) {
      edge = cluster.edges[e];
      parent = variants.filter(x => x.accession === edge[0])[0];
      child = variants.filter(x => x.accession === edge[1])[0];
      edgelist.push({
        'y1': parent.y1,
        'y2': child.y1,
        'x1': child.x1,  // vertical line segment
        'x2': child.x1,
        'parent': parent.label,
        'child': child.label,
        'dist': parseFloat(edge[2])
      });

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
    cluster['region'] = mode(points.map(x => x.region).flat());
    cluster['allregions'] = points.map(x => x.region).flat();
  }

  return(beaddata);
}


/**
 * Draw the beadplot for a specific cluster (identified through its
 * integer index <cid>) in the SVG.
 * @param {Number} cid:  integer index of cluster to draw as beadplot
 */
function beadplot(cid) {
  console.log(cid);

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
  var bTooltip = d3.select("body")
      .append("div")
      .attr("class", "tooltip")
      .style("opacity", 0);

  // update vertical range for consistent spacing between variants
  heightB = max_y * 10;
  $("#svg-cluster > svg").attr("height", heightB + marginB.top + marginB.bottom);
  yScaleB = d3.scaleLinear().range([40, heightB]);
  xScaleB.domain([mindate - 0.5*spandate, maxdate]);
  yScaleB.domain([min_y, max_y]);

  // clear SVG
  visB.selectAll('*').remove();

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
      .on("mouseover", function() {
        d3.select(this).attr("stroke-width", 5);
      })
      .on("mouseout", function() {
        d3.select(this).attr("stroke-width", 3);
      })
      .on("click", function(d) {
        var mystr = "";
        for (let [key, value] of Object.entries(d.country)) {
          mystr += `${key}: ${value}\n`;
        }
        $("#text-node").text(mystr);
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
      .on("mouseover", function() {
        d3.select(this).attr("stroke-width", 3);
      })
      .on("mouseout", function() {
        d3.select(this).attr("stroke-width", 1);
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
      .attr("fill", function(d) {
        return(country_pal[d.region1]);
      })
      .attr("stroke", "black")
      .on("mouseover", function(d) {
        d3.select(this).attr("stroke-width", 2)
            .attr("r", 4*Math.sqrt(d.count)+3);
        bTooltip.transition()       // Show tooltip
            .duration(200)
            .style("opacity", 0.9);

        // Display region distribution in tooltip
        let my_regions = table(d.region),
            tooltipText = `<b>Number of cases</b><br>`;

        for (let [r_key, r_value] of Object.entries(my_regions)) {
          tooltipText += `${r_key}: ${r_value}<br>`
        }
        // Display total number of cases if variants are from multiple countries
        if(Object.keys(my_regions).length > 1) {
          tooltipText += `Total: ${d.count}<br>`
        }
        // Display the sample date
        let formatDate = d3.timeFormat("%Y-%m-%d");
        tooltipText += `<br><b>Sample Date:</b> ${formatDate(new Date(d.x))}<br>`
        bTooltip.html(tooltipText)
            .style("left", (d3.event.pageX + 10) + "px")    // Tooltip appears 10 pixels left of the cursor
            .style("top", (d3.event.pageY + "px"));
      })
      .on("mouseout", function(d) {
        if (!selected.includes(this)) {
          d3.select(this).attr("stroke-width", 1)
              .attr("r", 4*Math.sqrt(d.count));
          bTooltip.transition()     // Hide tooltip
              .duration(500)
              .style("opacity", 0);
        }
      })
      .on("click", function(d) {
        // TODO: display first 3, collapsed text
        console.log(d.labels);

        // TODO: incorporate the following into tool-tip
        var my_countries = table(d.country)
        var mystr = gentable(my_countries);
	console.log(mystr)
        $("#text-node").html(mystr);
	console.log(mystr)

        draw_region_distribution(table(d.region));
      });

  // draw x-axis
  visB.append("g")
      .attr("transform", "translate(0,20)")
      .call(
        d3.axisTop(xScaleB)
          .ticks(5)
          .tickFormat(d3.timeFormat("%Y-%m-%d"))
      );

}
/**
 * Function to generate table from my_countries object on bead click
 * @param {json object} my_countries: json object containing key (countries) value (cases count) pairs
 */
function gentable(my_countries){
	tablehtml = '<table><tr><th id="Countryheader">Country</th><th id="ccheader">Case Count</th></tr>';
	for (let [key, value] of Object.entries(my_countries)) {
		tablehtml += '<tr><td>' + `${key}` + '</td><td>' + `${value}` + '</td></tr>';
		console.log(key,value)
	}
	return tablehtml+= '</table>';
}


/**
 * Draws a bar chart of the distribution of cases across regions
 * @param {{}} my_regions: associative list of region and case count pairs
 */
function draw_region_distribution(my_regions) {

    // Include all regions
    const regions = ["Africa", "Asia", "China", "Europe", "North America", "Oceania", "South America"];
    regions.forEach(function(reg) {
        if (!my_regions.hasOwnProperty(reg)) {my_regions[reg] = 0;}
    });

    // Store {region, count} pairs in a list, sorted in decreasing order
    const regKeys = Object.keys(my_regions),
        regionOrder = regKeys.sort(function(a, b) {return my_regions[b] - my_regions[a]}),
        sortedRegions = [];
    regionOrder.forEach(function(reg, idx) {
      sortedRegions.push({"region":regionOrder[idx], "count": my_regions[reg]})
    });

    // Set the margins
    const margin = 60,
        width = 600 - 2 * margin,
        height = 400 - 2 * margin;

    // Create the barchart
    const svg = d3.select("#text-node")
        .append("svg")
        .attr("width", 600)
        .attr("height", 400);
    const chart = svg.append("g")
        .attr("transform", `translate(${margin}, ${margin})`);

    // Set the scale of the x-axis
    const xScale = d3.scaleBand()
        .range([0, width])
        .domain(sortedRegions.map((r) => r.region))
        .padding(0.2);
    chart.append("g")
        .attr("transform", `translate(0, ${height})`)
        .call(d3.axisBottom(xScale));

    // Set the scale of the y-axis
    const yScale = d3.scaleLinear()
        .range([height, 0])
        .domain([0, sortedRegions[0].count])
        .nice()
    chart.append("g")
        .call(d3.axisLeft(yScale).ticks((sortedRegions[0].count < 10) ? sortedRegions[0].count : 10));

    const regionBars = chart.selectAll()
        .data(sortedRegions)
        .enter()
        .append("g")

    // Draw bars on the bar chart
    regionBars
        .append("rect")
        .attr("x", (r) => xScale(r.region))
        .attr("y", (r) => yScale(r.count))
        .attr("height", (r) => height - yScale(r.count))
        .attr("width", xScale.bandwidth())
        .attr("fill", "darkgray");

    // Write the case count above each bar
    regionBars.append("text")
        .style("font", "0.8em/1.2 Lato, sans-serif")
        .attr("x", (r) => xScale(r.region) + xScale.bandwidth() / 2)
        .attr("y", (r) => yScale(r.count) - 10)
        .attr("text-anchor", "middle")
        .text((r) => `${r.count}`);

    // Add axis labels
    svg.append("text")
        .attr("x", width / 2 + margin)
        .attr("y", 40)
        .attr("text-anchor", "middle")
        .text("Cases by Region");

    svg.append("text")
        .attr("x", -(height / 2) - margin)
        .attr("y", margin/ 3)
        .attr("transform", "rotate(-90)")
        .attr("text-anchor", "middle")
        .text("Number of Cases");

    svg.append("text")
        .attr("x", (width / 2) + margin)
        .attr("y", height + margin * 1.8)
        .attr("text-anchor", "middle")
        .text("Region")
}
