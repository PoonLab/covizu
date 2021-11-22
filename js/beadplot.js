/**
 *  Configure SVG to display beadplots
 */
var marginB = {top: 50, right: 10, bottom: 50, left: 10},
    widthB = document.getElementById("svg-cluster").clientWidth - marginB.left - marginB.right,
    heightB = 1000 - marginB.top - marginB.bottom,
    pixelsPerDay = 10;

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
 * @param {Object} variant:  associative list member of cluster.nodes
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
      'x1': utcDate(mindate),  // cluster min date
      'x2': utcDate(maxdate),  // cluster max date
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
    var label = variant[0][2].replace(pat, "$1"),
        coldates = variant.map(x => x[0]),
        isodate, samples, regions;

    coldates.sort();
    //Retrieving countries from variants?
    var country = variant.map(x => x[2].split('/')[1]),
        isodates = unique(coldates);

    // remove underscores in country names
    country = country.map(x => x.replace(/_/g," "));

    vdata = {
      'accession': accn,
      'label': label,
      'x1': utcDate(coldates[0]),  // min date
      'x2': utcDate(coldates[coldates.length-1]),  // max date
      'y1': y,
      'y2': y,
      'count': coldates.length,
      'country': tabulate(country),
      'region': country.map(x => countries[x]),
      'numBeads': isodates.length,
      'parent': null,
      'dist': 0,
      'unsampled': false
    };

    for (var i=0; i<isodates.length; i++) {
      isodate = isodates[i];
      samples = variant.filter(x => x[0] === isodate);
      country = samples.map(x => x[2].split('/')[1]);
      country = country.map(x => x.replace(/_/g," "));
      regions = country.map(x => countries[x]);

      // warn developers if no region for country
      if (regions.includes(undefined)) {
        //console.log("Developer msg, need to update countries.json:");
        for (const j in regions.filter(x => x===undefined)) {
          let this_country = samples[j].country;
          // unsampled lineages have undefined country fields
          if (this_country !== undefined) {
            console.log(`Need to add "${samples[j].country}" to countries.json`);
          }
        }
      }

      pdata.push({
        cidx,
        'variant': accn,
        'x': utcDate(isodate),
        'y': y,
        'count': samples.length,
        'accessions': samples.map(x => x[1]),
        'labels': samples.map(x => x[2].replace(pat, "$1")),
        'region1': mode(regions),
        'region': regions,
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
      coldates = nodelist.map(a => cluster.nodes[a].map(x => x[0])).flat();
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
    regions = points.map(x => x.region).flat();
    cluster['region'] = mode(regions);
    cluster['allregions'] = tabulate(regions);

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

function parse_mutation_annotations(mut_annotations) {
  var mutations = [];

  // Sorts tips according to cluster_idx
  sorted_tips = [...tips]
  sorted_tips.sort(function(a,b) {
    return a.cluster_idx - b.cluster_idx
  });

  for (const cidx in tips) {
    let phenotype = [],
        mutations_list = [],
        frequency_list = [],
        tip = sorted_tips[cidx];

    if (!('mutations' in tip)) {
      tip['mutations'] = []
    }

    for (const [mutation, freq] of Object.entries(tip.mutations)) {
      mutations_list.push(mutation);
      frequency_list.push(parseFloat(freq).toFixed(2))
      if (mutation in mut_annotations) {
        phenotype.push(mut_annotations[mutation])
      }
      else {
        phenotype.push([])
      }
    }

    mutations.push({
      'mutation': mutations_list,
      'frequency' : frequency_list,
      'phenotype': phenotype
    })
  }

  return mutations
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
      .attr("stroke-width", 5)
      .attr("bead", bead.accessions[0]);
}


function draw_halo_front(bead) {
  d3.selectAll(".selectionH").remove();

  visB.append("circle")
    .data([{"y": bead.y, "count": bead.count}])
    .attr('class', "selectionH")
    .attr("cx", xScaleB(xValueB(bead)))
    .attr("cy", yScaleB(yValueB(bead)))
    .attr("r", 4 * Math.sqrt(bead.count) + 4)
    .attr("fill", "none")
    .attr("stroke", "grey")
    .attr("fill-opacity", 1)
    .attr("stroke-width", 5)
    .attr("bead", bead.accessions[0]);
}

/**
 * Reset beadplot display to default settings.
 */
function clear_selection() {
  // clear search field
  // $('#search-input').val('');

  // Update stats so that the "Previous" and "Next" buttons do not become enabled in the next search query
  const stats = search_results.update({
    current_point: -1,
    total_points: 0,
   });
   update_search_stats(stats); 

  // $('#error_message').text(``);

  // Changes opacity of currently clicked cluster
  d3.selectAll("rect.clicked").attr('class', "clicked");
  $("#loading").hide();
  $("#loading_text").text(``);

  d3.select("#svg-cluster").selectAll("line")
      .attr("stroke-opacity", 1);
  d3.selectAll("circle:not(.selectionH):not(.selectionLC)")
      .attr("class", "default");
  d3.select('#svg-timetree').selectAll("rect:not(.clicked):not(.clickedH)")
      .attr("class", "default");
  d3.selectAll("circle.selectionH").remove();

  // Clear vertical edge selection
  d3.selectAll(".selectionLH").attr("stroke-width", function(d) {
    if (d.unsampled) {
      return 1;
    } else {
      return 3;
    }
  });
  d3.selectAll(".selectionLC").attr("r", function(d) {
    if (this.classList.contains("selectionH")) {
      // draw_halo(d)
      return 4*Math.sqrt(d.count)+4;
    }
    return 4 * Math.sqrt(d.count);   
  });
  d3.selectAll(".selectionLC").attr("stroke-width", function(d) {
    if (this.classList.contains("selectionH")) {
      return 5;
    }
    return 1;   
  });
  d3.selectAll(".selectionL").attr("stroke-width", 1);

  d3.selectAll(".selectionLH").attr("class","lines");
  d3.selectAll(".selectionL").attr("class","lines");
  d3.selectAll(".selectionLC").attr("class","default");
}


/**
 * Draw the beadplot for a specific cluster (identified through its
 * integer index <cid>) in the SVG.
 * @param {Number} cid:  integer index of cluster to draw as beadplot
 */
function beadplot(cid) {
  // Update global cindex for SVG and NWK filenames
  cindex = cid;
  
  var variants = beaddata[cid].variants,
      edgelist = beaddata[cid].edgelist,
      points = beaddata[cid].points;

  // rescale slider
  let max_dist = Math.max(...edgelist.map(x => x.dist));
  if (isFinite(max_dist) == false) {
    max_dist = 2.0;
  }
  let slider = $("#vedge-slider");
  slider.slider("value", 2.0)
        .slider("option", "max", max_dist )
  move_arrow();
  $( "#custom-handle" ).text("2.0");

  // Beadplots aren't expanded by default
  if ($('#expand-option').attr('checked')) {
    $('.switch').trigger('click');
    $('#expand-option').removeAttr('checked')
  }

  function redraw() {
    // set plotting domain
    var mindate = d3.min(variants, xValue1B),
        maxdate = d3.max(variants, xValue2B),
        spandate = maxdate-mindate,  // in milliseconds
        min_y = d3.min(variants, yValue1B),
        max_y = d3.max(variants, yValue1B),
        numDays = d3.timeDay.count(mindate, maxdate),
        clientWidth = document.getElementById("svg-cluster").clientWidth;

    // Sets margin top to align vertical scrollbar with the beadplot
    $('#beadplot-vscroll').css('margin-top', document.getElementById("beadplot-title").clientHeight + document.getElementById("svg-clusteraxis").clientHeight +
                                            ($('#expand-option').attr('checked') ? $('#inner-hscroll').height() : 0));

    // Don't give the user the option to scroll horizontally if the beadplot cannot expand 
    if (numDays * pixelsPerDay > clientWidth) {
      $('.expand').show();
      $('.switch').show();
      if ($('#expand-option').attr('checked')) { $('#beadplot-hscroll').show(); }
    }
    else {
      $('.expand').hide();
      $('.switch').hide();
      $('#beadplot-hscroll').hide();
    }

    let currentWidth = null;
    if ($('#expand-option').attr('checked')) {
       currentWidth = (numDays * pixelsPerDay > clientWidth ? numDays * pixelsPerDay : clientWidth)
                      - marginB.left -  marginB.right;
    }
    else {
      currentWidth = clientWidth - marginB.left -  marginB.right;
    }

    // update vertical range for consistent spacing between variants
    heightB = max_y * 10 + 40;
    $("#svg-cluster > svg").attr("height", heightB + marginB.top + marginB.bottom);
    yScaleB = d3.scaleLinear().range([20, heightB]);
    yScaleB.domain([min_y, max_y]);

    xScaleB = d3.scaleLinear().range([0, currentWidth]);
    var xAxis = d3.scaleTime().range([0, currentWidth]);

    // allocate space for text labels on left
    xScaleB.domain([mindate - 170*(spandate/currentWidth), maxdate]);
    xAxis.domain([mindate - 170*(spandate/currentWidth), maxdate]);

    // update horizontal range

    // clear SVG
    visB.selectAll('*').remove();
    visBaxis.selectAll('*').remove();

    $("#svg-cluster > svg").attr("width",currentWidth + marginB.left + marginB.right);
    $("#inner-hscroll").css("width",currentWidth + marginB.left + marginB.right);
    
    // Add 15px to account for the scrollbar for the beadplot
    $("#svg-clusteraxis > svg").attr("width", currentWidth + marginB.left + marginB.right + 15);


    // draw vertical line segments that represent edges in NJ tree
    visB.selectAll("lines")
        .data(edgelist.filter(x => x.dist <= slider.slider("value")))
        .enter().append("line")
        .attr("class", "lines")
        .attr("x1", xMap1B)
        .attr("x2", xMap2B)
        .attr("y1", yMap1B)
        .attr("y2", yMap2B)
        .attr("stroke-width", 1)
        .attr("stroke", function(d) {
          return("#bbd");
        })
        .on("mouseover", function(d) {
          // visual feedback
          var edge = d3.select(this);
          edge.attr("stroke-width", 3);

          let parent_variant = d3.select(".lines#"+d.parent.replace('/', '-').replace(' ', '_')),
              child_variant = d3.select(".lines#"+d.child.replace('/', '-').replace(' ', '_'));

          if (!parent_variant.empty()) {
            if (parent_variant.datum().count > 0) {
              d3.selectAll("circle").filter(ci => ci.y === d.y1)
                  .attr("stroke-width", function() {
                    if (this.classList.contains("selectionH")) return 5;
                    return 1.5;
                  })
                  .attr("r", function(d) {
                    if (this.classList.contains("selectionH")) return 4*Math.sqrt(d.count)+7;
                    return 4 * Math.sqrt(d.count) + 3;
                  });
            }
            parent_variant.attr("stroke-width", 5);
          }

          if (!child_variant.empty()) {
            if (child_variant.datum().count > 0) {
              d3.selectAll("circle").filter(ci => ci.y === d.y2)
                  .attr("stroke-width", function() {
                    if (this.classList.contains("selectionH")) return 5;
                    return 1.5;
                  })
                  .attr("r", function(d) {
                    if (this.classList.contains("selectionH")) return 4*Math.sqrt(d.count)+7;
                    return 4*Math.sqrt(d.count)+3;
                  });
            }
            child_variant.attr("stroke-width", 5);
          }

          // Show tooltip
          cTooltip.transition()
              .duration(50)
              .style("opacity", 0.75);

          let tooltipText = `<b>${i18n_text.vedge_parent}:</b> ${d.parent}<br/><b>${i18n_text.vedge_child}:</b> ${d.child}<br/>`;
          tooltipText += `<b>${i18n_text.vedge_distance}:</b> ${Math.round(d.dist*100)/100}<br/>`;
          tooltipText += `<b>${i18n_text.vedge_support}:</b> ${(d.support === undefined) ? 'n/a' : d.support}<br/>`
          tooltipText += `<b>${i18n_text.vedge_coldate}:</b> ${formatDate(d.x2)}`;

          cTooltip.html(tooltipText)
              .style("left", function() {
                if (d3.event.pageX > window.innerWidth/2) {
                  return (
                    d3.event.pageX - 15 -
                    cTooltip.node().getBoundingClientRect().width +
                    "px"
                  );
                } else {
                  return d3.event.pageX + 15 + "px";
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
          if(this.classList.contains("selectionL")) {
            d3.select(this).attr("stroke-width", 3);
          }
          else {
            d3.select(this).attr("stroke-width", 1);
          }

          let parent_variant = d3.select(".lines#"+d.parent.replace('/', '-').replace(' ', '_')),
              child_variant = d3.select(".lines#"+d.child.replace('/', '-').replace(' ', '_'));

          if (!parent_variant.empty()) {
            if (parent_variant.datum().count > 0) {
              d3.selectAll("circle").filter(ci => ci.y === d.y1)
                  .attr("stroke-width", function() {
                    if (this.classList.contains("selectionH")) return 5;
                    if (this.classList.contains("selectionLC")) return 1.5;
                    return 1;
                  })
                  .attr("r", function(d) {
                    if (this.classList.contains("selectionH")) return 4*Math.sqrt(d.count)+4;
                    if (this.classList.contains("selectionLC")) return 4 * Math.sqrt(d.count) + 3;
                    return 4 * Math.sqrt(d.count);
                  });
            }

            if(this.classList.contains("selectionL") || parent_variant.node().classList.contains("selectionLH")) {
              parent_variant.attr("stroke-width", 5);
            }
            else {
              parent_variant.attr("stroke-width", 3);
            }
          }

          if (!child_variant.empty()) {
            if (child_variant.datum().count > 0) {
              d3.selectAll("circle").filter(ci => ci.y === d.y2)
                  .attr("stroke-width", function() {
                    if (this.classList.contains("selectionH")) return 5;
                    if (this.classList.contains("selectionLC")) return 1.5;
                    return 1;
                  })
                  .attr("r", function(d) {
                    if (this.classList.contains("selectionH")) return 4*Math.sqrt(d.count)+4;
                    if (this.classList.contains("selectionLC")) return 4*Math.sqrt(d.count)+3;
                    return 4 * Math.sqrt(d.count);
                  });
            }
            if(this.classList.contains("selectionL") || child_variant.node().classList.contains("selectionLH")) {
              child_variant.attr("stroke-width", 5);
            }
            else {
              child_variant.attr("stroke-width", 3);
            }
          }

          cTooltip.transition()     // Hide tooltip
              .duration(50)
              .style("opacity", 0);
        })
        .on("click", function(d) {
          // Clear previous selections 
          clear_selection();

          // Add attributes to keep track of line selection
          var edge = d3.select(this);
          edge.attr("class", "lines selectionL");

          let parent_variant = d3.select(".lines#"+d.parent.replace('/', '-').replace(' ', '_')),
              child_variant = d3.select(".lines#"+d.child.replace('/', '-').replace(' ', '_'));

          parent_variant.attr("class", "lines selectionLH");
          child_variant.attr("class", "lines selectionLH");

          if (!parent_variant.empty()) {
            if (parent_variant.datum().count > 0) {
              d3.selectAll("circle").filter(ci => ci.y === d.y1) 
                  .attr("class", function() {
                    if (this.classList.contains("selectionH")) return "selectionH selectionLC";
                    return "selectionLC";
                  })
            }
          }

          if (!child_variant.empty()) {
            if (child_variant.datum().count > 0) {
              d3.selectAll("circle").filter(ci => ci.y === d.y2)
                  .attr("class", function() {
                    if (this.classList.contains("selectionH")) return "selectionH selectionLC";
                    return "selectionLC";
                  })
            }
          }
        });

    // draw horizontal line segments that represent variants in cluster
    visB.selectAll("lines")
        .data(variants)
        .enter().append("line")
        .attr("class", "lines")
        .attr("id", function(d) {
          return d.label.replace('/', '-').replace(' ', '_');
        })
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
              tooltipText += `<b>${i18n_text.vedge_parent}:</b> ${d.parent}<br/><b>${i18n_text.vedge_distance}:</b> ${Math.round(d.dist*100)/100}<br/>`;
            }
            if (!d.unsampled) {
              tooltipText += region_to_string(tabulate(d.region));
              tooltipText += `<b>${i18n_text.hedge_unique_dates}:</b> ${d.numBeads}<br/>`;
              tooltipText += `<b>${i18n_text.hedge_coldates}:</b><br>${formatDate(d.x1)} / ${formatDate(d.x2)}`;
            }

            if (tooltipText.length > 0) {
              // Show tooltip
              cTooltip.transition()
                  .duration(50)
                  .style("opacity", 0.75);

              // Tooltip appears 10 pixels left of the cursor
              cTooltip.html(tooltipText)
                  .style("left", function() {
                    if (d3.event.pageX > window.innerWidth/2) {
                      return (
                        d3.event.pageX - 10 -
                        cTooltip.node().getBoundingClientRect().width +
                        "px"
                      );
                    } else {
                      return d3.event.pageX + 10 + "px";
                    }
                  })
                  .style("top", function() {
                    if (d3.event.pageY > window.innerHeight/2) {
                      return (d3.event.pageY - cTooltip.node().getBoundingClientRect().height - 10) + "px";
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
                if (this.classList.contains("selectionLH")) {
                  return 5;
                } else if (d.unsampled) {
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
        .attr("idx", function(d, i) { return i; })
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
              .style("opacity", 0.75);

          let tooltipText = "";
          if (d.parent || d.dist) {
            tooltipText += `<b>Parent:</b> ${d.parent}<br/><b>${i18n_text.vedge_distance}:</b> ${Math.round(d.dist*100)/100}<br/>`;
          }
          tooltipText += region_to_string(tabulate(d.region));
          tooltipText += `<b>${i18n_text.vedge_coldate}:</b> ${formatDate(d.x)}`;

          // Tooltip appears 10 pixels left of the cursor
          cTooltip.html(tooltipText)
              .style("left", function() {
                if (d3.event.pageX > window.innerWidth/2) {
                  return (
                    d3.event.pageX - 15 -
                    cTooltip.node().getBoundingClientRect().width +
                    "px"
                  );
                } else {
                  return d3.event.pageX + 15 + "px";
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
          if (this.classList.contains("selectionLC")) {
            d3.select(this).attr("stroke-width", 1.5)
                .attr("r", 4*Math.sqrt(d.count) + 3);
          } 
          else {
            d3.select(this).attr("stroke-width", 1)
                .attr("r", 4*Math.sqrt(d.count));
          }
            cTooltip//.transition()     // Hide tooltip
                //.duration(50)
                .style("opacity", 0);
        })
        .on("click", function(d) {
          clear_selection();
          d3.select(this).raise();
          draw_halo_front(d);

          gentable(d);
          draw_region_distribution(tabulate(d.region));
          gen_details_table(d);
          $('#search-input').val('');
          $('#end-date').val('');
          $('#start-date').val('');
          $('#error_message').text(``);

          // Disable buttons
          $('#search-button').attr("disabled", true);
          $('#clear_button').attr("disabled", true);
          $('#next_button').attr("disabled", true);
          $('#previous_button').attr("disabled", true);
          $('#search_stats').addClass("disabled_stats");
        });


    var tickCount = 0.005*currentWidth;
    if (tickCount <= 4) tickCount = 4;

    // Ensures that there aren't extra ticks in the time axis of the beadplot when the first and last coldates are within a few days
    var dateDiff = d3.timeDay.count(
                    d3.min(variants, function(d) {return d.x1}), 
                    d3.max(variants, function(d) {return d.x2}));
    if (dateDiff !=0 && dateDiff < tickCount) tickCount = dateDiff;

    // draw x-axis
    visBaxis.append("g")
        .attr("class", "treeaxis")
        .attr("transform", "translate(0,20)")
        .call(
          d3.axisTop(xAxis)
            .ticks(tickCount)
            .tickFormat(d3.timeFormat("%Y-%m-%d"))
        );

    // Sets the height of the vertical scrollbar element to have the same height as the beadplot
    $('#inner-vscroll').css('height', $('.beadplot-content > svg').height()); 
  }
  redraw();
  window.addEventListener("resize", redraw);

  // replace previous event listeners
  let listeners = $._data(slider[0], "events");
  if (listeners.slidechange !== undefined) {
    listeners.slidechange.pop();
  }
  if (listeners.slide !== undefined) {
    listeners.slide.pop();
  }

  slider.on(edgelist.length < 200 ? "slide" : "slidechange", function(event, ui) {
    if (ui.valueOf().value !== 2.0) {
      // avoid drawing beadplot twice on load

      // Update slider value before redrawing beadplot
      if (event.handleObj.type === "slide")
        slider.slider("value", ui.valueOf().value);

      redraw();
    }
  });
}


/**
 * Writes information about the region distribution to a string
 * @param {{}} my_regions: associative list of region and case count pairs
 * @returns {string} regStr: a string representation of the region distribution
 */
function region_to_string(my_regions) {
  // Display region distribution in tooltip
  let regStr = `<b>${i18n_text.tip_cases}:</b><br>`,
      total = 0;

  for (let [r_key, r_value] of Object.entries(my_regions)) {
    regStr += `&nbsp;&nbsp;${i18n_text.region_legend[r_key]}: ${r_value}<br>`;
    total += r_value;
  }

  // Display total number of cases if variants are from multiple countries
  if(Object.keys(my_regions).length > 1) {
    regStr += `${i18n_text.total}: ${total}<br>`
  }

  return regStr;
}

/**
 * Creates a table that displays mutations details (Mutation label, frequency, phenotypic effects)
 * @param obj: JS object or an array of JS Objects
 */

function gen_mut_table(obj) {
  var mutations = [];

  // Check for a list of samples
  if (Array.isArray(obj)) {
    // Iterate over each sample in the list
    for (let j = 0; j < obj.length; j++) {

      // "zip" the sequence details of each sample
      for (let i = 0; i < obj[j].mutations.length; i++) {
        let sample_details = [
          obj[j].mutations[i],
          obj[j].frequency[i],
          []
          // obj[j].phenotype[i]  
        ];
        mutations.push(sample_details);
      }
    } 
  }
  else {
    // "zip" the sequence details of each sample
    for (let i = 0; i < obj.mutations.length; i++) {
      let sample_details = [
        obj.mutations[i],
        obj.frequency[i],
        []
        // obj.phenotype[i]  
      ];
      mutations.push(sample_details);
    }
  }

  mut_thead.html(""); // clear table headers
  var sort_ascending = true;

  var headers = mut_thead.append('tr')
      .selectAll('th')
      .data(mut_theaders)
      .enter()
      .append('th')
      .on('click', function (x, i) {
        // Reset sorting arrows
        mut_thead.selectAll('a').classed('hide-before', false);
        mut_thead.selectAll('a').classed('hide-after', false);

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
      .classed('sort-by', true)

  // Create a row for each sample
  mut_tbody.html("");
  var t_rows = mut_tbody.selectAll('tr')
      .data(mutations)
      .enter()
      .append('tr')
      .on("mouseover", function (x) {
        // Highlight row on mouseover
        d3.select(this).style("background-color", "#e2e2e2");
      })
      .on("mouseout", function (x) {
        // Remove highlighting on mouseout
        d3.select(this).style("background-color", null);
      });

  // Create a cell for every row in the column
  t_rows.selectAll('td')
      .data(function (r) { return r.slice(0,3); })
      .enter()
      .append('td')
      .append('span')
      .on("mouseout", function() {
        cTooltip.transition()
            .duration(50)
            .style("opacity", 0);
      })
      .text(function (x) { return x; })
      .style("font", "0.875em/1.2 Lato, sans-serif");

  var t_cells = document.querySelectorAll("#tabs-3 table tbody tr")

  for (let j = 0; j < obj.phenotype.length; j++) {
    var phenotype = t_cells[j].children[2];
      for (let i = 0; i < obj.phenotype[j].length; i++) {
        if (obj.phenotype[j][i] == 'transmissibility') {
          var img = document.createElement("img");
          img.src = "../img/red_circle.png";
          img.classList.add('phenotype_icon')
          phenotype.appendChild(img);
        }
        if (obj.phenotype[j][i] == 'immunosuppression_variant_emergence') {
          var img = document.createElement("img");
          img.src = "../img/blue_triangle.png";
          img.classList.add('phenotype_icon')
          phenotype.appendChild(img);
        }
      }
  }    
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
        let sample_details = [
          obj[j].accessions[i],
          obj[j].labels[i],
          formatDate(obj[j].x),
          obj[j].accessions[0]  // variant labeled by accession
        ];
        details.push(sample_details);
      }
    }
  }

  else {
    // "zip" the sequence details of each sample
    for (let i = 0; i < obj.accessions.length; i++) {
      let sample_details = [
        obj.accessions[i],
        obj.labels[i],
        formatDate(obj.x),
        obj.accessions[0]
      ];
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
      .on("mouseover", function (x) {
        //console.log(x);
        let circle = d3.select("circle#"+x[3]).attr("stroke-width", 2);
        circle.attr("r", 4*Math.sqrt(circle.datum().count)+3);

        // Highlight row on mouseover
        d3.select(this).style("background-color", "#e2e2e2");
      })
      .on("mouseout", function (x) {
        let circle = d3.select("circle#"+x[3]).attr("stroke-width", 1);
        circle.attr("r", 4*Math.sqrt(circle.datum().count));
        // Remove highlighting on mouseout
        d3.select(this).style("background-color", null);
      });

  // Create a cell for every row in the column
  t_rows.selectAll('td')
      .data(function (r) { return r.slice(0,3); })
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
            tooltipText = `<b>${i18n_text.sample_orig_lab}:</b> ${gd.covv_orig_lab}<br/>`;
            tooltipText += `<b>${i18n_text.sample_subm_lab}:</b> ${gd.covv_subm_lab}<br/>`;
            tooltipText += `<b>${i18n_text.sample_authors}:</b> ${gd.covv_authors}`;

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
  var my_countries = Object.entries(obj.country),
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
      .text(i18n_text.number_cases);

  svg.append("text")
      .style("font", "0.8em/1.2 Lato, sans-serif")
      .attr("text-anchor", "middle")
      .attr("x", (width / 2) + margin.left)
      .attr("y", height + 85)
      .text(i18n_text.region);
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


function select_next_bead(next_node) {
  var select_bead = d3.selectAll('circle[id="'+next_node.id+'"]');

  // Moves the selected bead (circle) to the bottom of the list
  select_bead.raise();
  
  var working_bead = select_bead.nodes()[0];
  working_bead.scrollIntoView({block: "center"});
  update_table_individual_bead_front(d3.select(working_bead).datum());

  if (next_node.classList.contains("SelectedBead") && search_results.get().total_points > 0) {
    var search_index = search_results.get().beads[next_node.id];
    if (search_index !== undefined) {
      const stats = search_results.update({
        current_point: search_index
      });

      update_search_stats(stats);
    }
  }
}

/*********************** Moving Arrow Key Button ***********************/
function move_arrow() {
  var style = parseFloat($( "#custom-handle" ).css("left").split("px")[0]);

  if (style > 121) {
    $("#right-arrow").css("margin-left", style-121);
  }
  else
    $("#right-arrow").css("margin-left", 0);
}