/**
 *  Configure SVG to display beadplots
 */
var marginB = {top: 50, right: 10, bottom: 50, left: 10},
    widthB = document.getElementById("svg-cluster").clientWidth - marginB.left - marginB.right,
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

var region_map = {};  // map country to region

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
 * Tabulate values in array.
 * @param {Array} arr:  Array of values to tabulate
 * @returns {{}} Associative list of unique value: count pairs
 */
function tabulate(arr) {
  var val, counts = {};
  for (var i=0; i<arr.length; i++) {
    val = arr[i];
    if (val === null) {
      continue;
    }
    if (counts[val] === undefined) {
      counts[val] = 0;
    }
    counts[val]++;
  }
  return(counts);
}


/**
 * Parse nodes that belong to the same variant.
 * A variant is a collection of genomes that are indistinguishable with respect to
 * (1) differences from the reference or (2) placement in the phylogenetic tree.
 * Visually, the variant is represented by a horizontal line segment spanning the
 * sample collection dates.
 * The samples that comprise a variant are gathered by collection date into "points"
 * along the horizontal line.  If there are no samples, then the variant is
 * "unsampled" and spans the entire width of the beadplot.
 *
 * @param {Object} variant:  Array, nested list of sample metadata
 * @param {number} y:  vertical position
 * @param {number} cidx:  cluster index for labeling points
 * @param {string} accn:  accession of baseline sample of variant
 * @param {Date} mindate:  used only for unsampled variants
 * @param {Date} maxdate:  used only for unsampled variants
 * @returns {{variants: [], points: []}}
 */
function parse_variant(variant, y, cidx, accn, mindate, maxdate) {
  var vdata, pdata = [];

  if (variant.length === 0) {
    // handle unsampled internal node
    if (mindate === null || maxdate === null) {
      // this would happen if a cluster comprised only one unsampled node - should be impossible!
      console.log("Error in parse_variants(): cannot draw unsampled variant without min and max dates");
    }
    vdata = {
      'accession': accn,
      'label': accn,
      'x1': new Date(mindate),  // cluster min date
      'x2': new Date(maxdate),  // cluster max date
      'y1': y,
      'y2': y,
      'count': 0,
      'country': null,
      'region': null,
      'numBeads': 0,
      'parent': null,
      'dist': 0,
      'unsampled': true
    };
  }
  else {
    // parse samples within variant, i.e., "beads"
    // each sample is an Array: [name, accession, location, coldate, gender, age, status]
    var label = variant[0][0].replace(pat, "$1"),
        coldates = variant.map(x => x[3]),
        isodate, samples, regions;

    coldates.sort();
    //Retrieving countries from variants?
    var geo = variant.map(x => x[2].split(' / ')),
        region = geo.map(x => x[0]),
        country = geo.map(x => x[1]),
        locale = geo.map(x => x[2]),
        isodates = unique(coldates);

    // remove underscores in country names
    //country = country.map(x => x.replace(/_/g," "));

    // update country to region map
    for (let i=0; i < country.length; i++) {
      let this_country = country[i],
          this_region = region[i];
      if (region_map[this_country] === undefined) {
        region_map[this_country] = this_region;
      }
    }

    vdata = {
      'accession': accn,
      'label': label,
      'x1': new Date(coldates[0]),  // min date
      'x2': new Date(coldates[coldates.length-1]),  // max date
      'y1': y,
      'y2': y,
      'count': coldates.length,
      'country': tabulate(country),
      'region': tabulate(region),
      'locale': locale,
      'numBeads': isodates.length,
      'parent': null,
      'dist': 0,
      'unsampled': false
    };

    for (let i=0; i<isodates.length; i++) {
      isodate = isodates[i];

      // FIXME: this is not very efficient
      samples = variant.filter(x => x[3] === isodate);
      geo = samples.map(x => x[2].split(' / '));
      region = geo.map(x => x[0]);
      country = geo.map(x => x[1]);
      locale = geo.map(x => x[2]);

      pdata.push({
        cidx,
        'variant': accn,
        'x': new Date(isodate),
        'y': y,
        'count': samples.length,
        'accessions': samples.map(x => x[1]),
        'labels': samples.map(x => x[0].replace(pat, "$1")),
        'region1': mode(region),  // most common region
        'region': tabulate(region),
        'country': tabulate(country),
        'parent': null,
        'dist': 0
      })
    }
  }

  return {'variant': vdata, 'points': pdata};
}


/**
 *
 * @param cluster
 * @param variants
 * @param points
 * @returns {[]}
 */
function parse_edgelist(cluster, variants, points) {
  // map earliest collection date of child node to vertical edges
  let edge, parent, child, dist, support,
      edgelist = [];

  // generate maps of variants and points keyed by accession
  let lookup_variant = {};
  variants.forEach(function(row) {
    lookup_variant[row.accession] = row;
  });

  let lookup_points = {},
      index_points = {};  // index by y-coordinate

  points.forEach(function(pt) {
    if (lookup_points[pt.variant] === undefined) {
      lookup_points[pt.variant] = [];
    }
    lookup_points[pt.variant].push(pt);
    if (index_points[pt.y] === undefined) {
      index_points[pt.y] = [];
    }
    index_points[pt.y].push(pt);
  })

  for (var e = 0; e < cluster.edges.length; e++) {
    edge = cluster.edges[e];
    //parent = variants.filter(x => x.accession === edge[0])[0];
    parent = lookup_variant[edge[0]];
    //child = variants.filter(x => x.accession === edge[1])[0];
    child = lookup_variant[edge[1]];

    if (parent === undefined || child === undefined) {
      // TODO: handle edge to unsampled node
      continue;
    }

    dist = parseFloat(edge[2]);
    if (edge[3] === null) {
      support = undefined;
    } else {
      support = parseFloat(edge[3]);
    }

    edgelist.push({
      'y1': parent.y1,
      'y2': child.y1,
      'x1': child.x1,  // vertical line segment
      'x2': child.x1,
      'parent': parent.label,
      'child': child.label,
      'dist': dist,
      'support': support
    });

    child.parent = parent.label;
    child.dist = dist;

    // Assign the parent and genomic distance of each point
    if (index_points[child.y1] !== undefined) {
      for (let pt of index_points[child.y1]) {
        pt.parent = parent.label;
        pt.dist = dist;
      }
    }
    /*
    for (let c = 0; c < points.length; c++) {
      if (points[c].y === child.y1) {
        points[c].parent = parent.label;
        points[c].dist = dist;
      }
    }
    */

    // update variant time range
    if (parent.x1 > child.x1) {
      parent.x1 = child.x1;
    }
    if (parent.x2 < child.x1) {
      parent.x2 = child.x1;
    }
  }
  return edgelist;
}


/**
 *
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
  return(total);
}

/**
 * Parse node and edge data from clusters JSON to a format that is
 * easier to map to SVG.
 * @param {Object} clusters:
 */
function parse_clusters(clusters) {
  var cluster, variant, coldates, regions, labels,
      accn, mindate, maxdate,
      result, vdata, pdata,
      variants,  // horizontal line segment data + labels
      edgelist,  // vertical line segment data
      points,  // the "beads"
      beaddata = [];  // return value

  for (const cidx in clusters) {
    cluster = clusters[cidx];

    variants = [];
    points = [];
    edgelist = [];

    // deal with edge case of cluster with only one variant, no edges
    if (Object.keys(cluster["nodes"]).length === 1) {
      variant = Object.values(cluster.nodes)[0];
      accn = Object.keys(cluster.nodes)[0];
      result = parse_variant(variant, 0, cidx, accn, null, null);
      vdata = result['variant'];
      variants.push(vdata);

      pdata = result['points'];
      points = points.concat(pdata);
    }
    else {
      // de-convolute edge list to get node list in preorder
      var nodelist = unique(cluster.edges.map(x => x.slice(0, 2)).flat());

      // date range of cluster
      coldates = nodelist.map(a => cluster.nodes[a].map(x => x[3])).flat();
      coldates.sort()
      mindate = coldates[0];
      maxdate = coldates[coldates.length - 1];

      // extract the date range for each variant in cluster
      var y = 1;
      for (const accn of nodelist) {
        // extract collection dates for all samples of this variant
        variant = cluster.nodes[accn];
        result = parse_variant(variant, y, cidx, accn, mindate, maxdate);
        variants.push(result['variant']);
        points = points.concat(result['points']);
        y++;
      }

      edgelist = parse_edgelist(cluster, variants, points);
    }

    beaddata.push({
      'variants': variants,
      'edgelist': edgelist,
      'points': points
    });

    // calculate consensus region for cluster
    // collect all region Arrays for all samples, all variants
    cluster['region'] = merge_tables(points.map(x => x.region));
    let max_freq = 0, max_region = '';
    for (let row of Object.entries(cluster['region'])) {
      if (row[1] > max_freq) {
        max_freq = row[1];
        max_region = row[0];
      }
    }
    cluster['region1'] = max_region;  // most common region

    // concatenate all sample labels within cluster for searching
    labels = points.map(x => x.labels).flat();

    // decompose labels and only keep unique substrings
    let uniq = new Set(labels.map(x => x.split('/')).flat());
    cluster['searchtext'] = Array.from(uniq).join();
    cluster['label1'] = labels[0];

    // collect all countries
    cluster['country'] = merge_tables(variants.map(x => x.country));
  }
  return beaddata;
}


function create_selection(selected_obj) {
  d3.select("div#svg-cluster").selectAll("line").attr("stroke-opacity", 0.3);
  d3.select("div#svg-cluster")
    .selectAll("circle:not(.SelectedBead):not(.selectionH)")
    .attr("class", "not_SelectedBead");

  selected_obj.attr("class", "SelectedBead");
}


function draw_halo(bead) {
  // remove other halos
  d3.selectAll(".selectionH").remove();

  visB.append("circle").lower()
      .attr('class', "selectionH")
      .attr("cx", xScaleB(xValueB(bead)))
      .attr("cy", yScaleB(yValueB(bead)))
      .attr("r", 4 * Math.sqrt(bead.count) + 4)
      .attr("fill", "none")
      .attr("stroke", "grey")
      .attr("fill-opacity", 1)
      .attr("stroke-width", 5);
}


/**
 * Reset beadplot display to default settings.
 */
function clear_selection() {
  // clear search field
  $('#search-input').val('');

  d3.select("#svg-cluster").selectAll("line")
      .attr("stroke-opacity", 1);
  d3.selectAll("circle:not(.selectionH)")
      .attr("class", "default");
  d3.select('#svg-timetree').selectAll("rect:not(.clicked):not(.clickedH)")
      .attr("class", "default");
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

  function redraw() {    
    currentWidth = document.getElementById("svg-cluster").clientWidth - marginB.left -  marginB.right;

    // set plotting domain
    var mindate = d3.min(variants, xValue1B),
        maxdate = d3.max(variants, xValue2B),
        spandate = maxdate-mindate,  // in milliseconds
        min_y = d3.min(variants, yValue1B),
        max_y = d3.max(variants, yValue1B);

    // update vertical range for consistent spacing between variants
    heightB = max_y * 10 + 40;
    $("#svg-cluster > svg").attr("height", heightB + marginB.top + marginB.bottom);
    yScaleB = d3.scaleLinear().range([20, heightB]);
    yScaleB.domain([min_y, max_y]);

    xScaleB = d3.scaleLinear().range([0, currentWidth]);
    // allocate space for text labels on left
    xScaleB.domain([mindate - 170*(spandate/currentWidth), maxdate]);

    // update horizontal range

    // clear SVG
    visB.selectAll('*').remove();
    visBaxis.selectAll('*').remove();

    $("#svg-cluster > svg").attr("width",currentWidth + marginB.left + marginB.right);
    $("#svg-clusteraxis > svg").attr("width", currentWidth + 20);


    // draw vertical line segments that represent edges in NJ tree
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
            return("#aad");
          } else if (d.dist < 2.5) {
            return("#bbd");
          } else {
            return("#ccf");
          }
        })
        .on("mouseover", function(d) {
          var edge = d3.select(this);
          edge.attr("stroke-width", 3);

          // Show tooltip
          cTooltip.transition()
              .duration(50)
              .style("opacity", 0.9);

          let tooltipText = `<b>Parent:</b> ${d.parent}<br/><b>Child:</b> ${d.child}<br/>`;
          tooltipText += `<b>Genomic distance:</b> ${Math.round(d.dist*100)/100}<br/>`;
          tooltipText += `<b>Support:</b> ${(d.support === undefined) ? 'n/a' : d.support}<br/>`
          tooltipText += `<b>Collection date:</b> ${formatDate(d.x2)}`;

          cTooltip.html(tooltipText)
              .style("left", function() {
                if (d3.event.pageX > window.innerWidth/2) {
                  return (
                    d3.event.pageX -
                    cTooltip.node().getBoundingClientRect().width +
                    "px"
                  );
                } else {
                  return d3.event.pageX + "px";
                }
              })
              .style("top", function() {
                if (d3.event.pageY > window.innerHeight/2) {
                  return (d3.event.pageY - cTooltip.node().getBoundingClientRect().height - 15) + "px";
                } else {
                  return d3.event.pageY + 15 + "px";
                }
              });
        })
        .on("mouseout", function() {
          d3.select(this).attr("stroke-width", 1);

          cTooltip.transition()     // Hide tooltip
              .duration(50)
              .style("opacity", 0);
        });

    // draw horizontal line segments that represent variants in cluster
    visB.selectAll("lines")
        .data(variants)
        .enter().append("line")
        .attr(  "class", "lines")
        .attr("x1", xMap1B)
        .attr("x2", xMap2B)
        .attr("y1", yMap1B)
        .attr("y2", yMap2B)
        .attr("stroke-width", function(d) {
          if (d.unsampled) {
            return 1;
          } else {
            return 3;
          }
        })
        .attr("stroke", function(d) {
          if (d.unsampled) {
            return "#ccc";
          } else {
            return "#777";
          }
        })
        .on("mouseover", function(d) {
          if (d.label !== null) {
            d3.select(this)
              .attr("stroke-width", 5);

            let tooltipText = "";
            if (d.parent || d.dist) {
              tooltipText += `<b>Parent:</b> ${d.parent}<br/><b>Genomic distance:</b> ${Math.round(d.dist*100)/100}<br/>`;
            }
            if (!d.unsampled) {
              tooltipText += region_to_string(tabulate(d.region));
              tooltipText += `<b>Unique collection dates:</b> ${d.numBeads}<br/>`;
              tooltipText += `<b>Collection dates:</b><br>${formatDate(d.x1)} / ${formatDate(d.x2)}`;
            }

            if (tooltipText.length > 0) {
              // Show tooltip
              cTooltip.transition()
                  .duration(50)
                  .style("opacity", 0.9);

              // Tooltip appears 10 pixels left of the cursor
              cTooltip.html(tooltipText)
                  .style("left", function() {
                    if (d3.event.pageX > window.innerWidth/2) {
                      return (
                        d3.event.pageX -
                        cTooltip.node().getBoundingClientRect().width +
                        "px"
                      );
                    } else {
                      return d3.event.pageX + "px";
                    }
                  })
                  .style("top", function() {
                    if (d3.event.pageY > window.innerHeight/2) {
                      return (d3.event.pageY - cTooltip.node().getBoundingClientRect().height - 15) + "px";
                    } else {
                      return d3.event.pageY + 10 + "px";
                    }
                  });
            }

          }
        })
        .on("mouseout", function() {
          d3.select(this)
              .attr("stroke-width", function(d) {
                // reset line width
                if (d.unsampled) {
                  return 1;
                } else {
                  return 3;
                }
              });
          cTooltip.transition()     // Hide tooltip
              .duration(50)
              .style("opacity", 0);
        })
        .on("click", function(d) {
          if (d.label !== null && d.country !== null) {
            gentable(d);
            draw_region_distribution(tabulate(d.region));
            let var_samples = points.filter(x => x.y === d.y1);
            gen_details_table(var_samples);
          }
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


    // draw "beads" to represent samples per collection date
    visB.selectAll("circle")
        .data(points)
        .enter().append("circle")
        .attr("r", function(d) { return (4*Math.sqrt(d.count)); })
        .attr("cx", xMapB)
        .attr("cy", yMapB)
        .attr("class", "default")
        .attr("id", function(d) { return d.accessions[0]; })
        .attr("search_hit", function(d) {
          if (search_results.get().beads.length == 0){
            return null;
          } else {
      return search_results.get().beads[d.accessions[0]] === undefined ? false : true;
          }
        })
        .attr("fill", function(d) {
          return(country_pal[d.region1]);
        })
        .attr("stroke", function(d) {
          return(country_pal[d.region1]);
        })
        .on("mouseover", function(d) {
          d3.select(this).attr("stroke-width", 2)
              .attr("r", 4*Math.sqrt(d.count)+3);

          cTooltip.transition()
              .duration(50)
              .style("opacity", 0.9);

          let tooltipText = "";
          if (d.parent || d.dist) {
            tooltipText += `<b>Parent:</b> ${d.parent}<br/><b>Genomic distance:</b> ${Math.round(d.dist*100)/100}<br/>`;
          }
          tooltipText += region_to_string(tabulate(d.region));
          tooltipText += `<b>Collection date:</b> ${formatDate(d.x)}`;

          // Tooltip appears 10 pixels left of the cursor
          cTooltip.html(tooltipText)
              .style("left", function() {
                if (d3.event.pageX > window.innerWidth/2) {
                  return (
                    d3.event.pageX -
                    cTooltip.node().getBoundingClientRect().width +
                    "px"
                  );
                } else {
                  return d3.event.pageX + "px";
                }
              })
              .style("top", function() {
                if (d3.event.pageY > window.innerHeight/2) {
                  return (d3.event.pageY - cTooltip.node().getBoundingClientRect().height - 15) + "px";
                } else {
                  return d3.event.pageY + 15 + "px";
                }
              });
        })
        .on("mouseout", function(d) {
            d3.select(this).attr("stroke-width", 1)
                .attr("r", 4*Math.sqrt(d.count));
            cTooltip//.transition()     // Hide tooltip
                //.duration(50)
                .style("opacity", 0);
        })
        .on("click", function(d) {
          clear_selection();
          draw_halo(d);
          gentable(d);
          draw_region_distribution(tabulate(d.region));
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
  redraw();
  window.addEventListener("resize", redraw);

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

  // Check for a list of samples
  if (Array.isArray(obj)) {

    // Iterate over each sample in the list
    for (let j = 0; j < obj.length; j++) {

      // "zip" the sequence details of each sample
      for (let i = 0; i < obj[j].accessions.length; i++) {
        let sample_details = [obj[j].accessions[i], obj[j].labels[i], formatDate(obj[j].x)];
        details.push(sample_details);
      }
    }
  }

  else {
    // "zip" the sequence details of each sample
    for (let i = 0; i < obj.accessions.length; i++) {
      let sample_details = [obj.accessions[i], obj.labels[i], formatDate(obj.x)];
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
        d3.select(this).style("background-color", "#e2e2e2");  // Highlight on mouseover
      })
      .on("mouseout", function () {            // Remove highlighting on mouseout
        d3.select(this).style("background-color", null);
      });

  // Create a cell for every row in the column
  t_rows.selectAll('td')
      .data(function (r) { return r; })
      .enter()
      .append('td')
      .append('span')
      .on("mouseover", function(x) {
        let gttx = d3.event.pageX,
            gtty = d3.event.pageY;

        if (x.startsWith("EPI_ISL")) {
          gisaid.getAcknowledgementData(x, function(gd) {
            cTooltip.transition()
                .duration(50)
                .style("opacity", 0.9);

            let tooltipText = "";
            tooltipText = `<b>Originating lab:</b> ${gd.covv_orig_lab}<br/>`;
            tooltipText += `<b>Submitting lab:</b> ${gd.covv_subm_lab}<br/>`;
            tooltipText += `<b>Authors:</b> ${gd.covv_authors}`;

            // Tooltip appears 10 pixels left of the cursor
            cTooltip.html(tooltipText)
                .style("left", (gttx - cTooltip.node().getBoundingClientRect().width - 20) + "px")
                .style("top", function() {
                  if (gtty > window.innerHeight/2) {
                    return (gtty - cTooltip.node().getBoundingClientRect().height) + "px";
                  } else {
                    return gtty + "px";
                  }
                });
          });
        }
      })
      .on("mouseout", function() {
        cTooltip.transition()
            .duration(50)
            .style("opacity", 0);
      })
      .text(function (x) { return x; })
      .style("font", "0.875em/1.2 Lato, sans-serif");
}


/**
 * Function to generate table from my_countries object on bead click
 * @param {Object} obj:  JS Object with country attribute
 */
function gentable(obj) {
  var my_countries = Object.entries(obj.country),  // returns Array of key/value pairs
      row, region;

  // annotate with region (continent)
  for (const i in my_countries) {
    row = my_countries[i];
    region = region_map[row[0]];
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
      d3.select(this).style("background-color", "#e2e2e2");}) //mouseover highlight
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
  const regions = ['Africa', 'Asia', 'Europe', 'North America', 'Oceania',
                   'South America'];
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


/**
 * Recursive function to serialize tree by postorder traversal of branches.
 * @param {str} parent:  unique node identifier
 * @param {Array} edgelist:  list of edges from beaddata Object
 * @returns {string}
 */
function serialize_branch(parent, edgelist) {
  var children = edgelist.filter(x => x.parent === parent),
      branch = edgelist.filter(x => x.child === parent),
      coldate;

  if (children.length === 0) {
    // terminal node, emit string
    // TODO: add count information (number of samples per variant)
    coldate = branch[0].x1.toISOString().split('T')[0];
    return branch[0].child+'|'+coldate+':'+branch[0].dist;  // node name : branch length
  }
  else {
    // follow branches toward tips
    let subtrees = children.map(x => serialize_branch(x.child, edgelist));
    let str = "("+subtrees.join(',')+")";
    if (branch.length > 0) {
      str += branch[0].support+':'+branch[0].dist;  // node support
    }
    return str;
  }
}

function serialize_beadplot(cidx) {
  var edgelist = beaddata[cidx].edgelist,
      root = edgelist[0].parent;
  return serialize_branch(root, edgelist)+';';
}
