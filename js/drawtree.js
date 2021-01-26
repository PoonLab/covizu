var margin = {top: 10, right: 20, bottom: 10, left: 0},
  width = 250 - margin.left - margin.right,
  height = 1200 - margin.top - margin.bottom;

// to store colour palettes
var sample_pal, coldate_pal, diverge_pal;

// set up plotting scales
var xValue = function(d) { return d.x; },
  xScale = d3.scaleLinear().range([0, width]),
  xMap1 = function(d) { return xScale(d.x1); },  // lines
  xMap2 = function(d) { return xScale(d.x2); },
  xWide = function(d) { return xScale(d.x2 - d.x1)};

var yValue = function(d) { return d.y; },
  yScale = d3.scaleLinear().range([height, 40]),  // inversion
  yMap1 = function(d) { return yScale(d.y1); },
  yMap2 = function(d) { return yScale(d.y2); },
  yMap = function(d) { return yScale(yValue(d)+0.4); };


var vis = d3.select("div#svg-timetree")
  .append("svg")
  .attr("width", width + margin.left + margin.right);
  //.attr("height", height + margin.top + margin.bottom);
  //.append("g");

var axis = d3.select("div#svg-timetreeaxis")
  .append("svg")
  .attr("width", width + margin.left + margin.right)
  .attr("height", 25)
  .append("g");

let cTooltip = d3.select("#tooltipContainer")
    .style("opacity", 0);

/**
 * Draw an open rectangle around the filled rectangle representing
 * a cluster/lineage in the time-scaled tree.
 * @param rect:  d3.Selection object
 */
function draw_cluster_box(rect) {
  var d = rect.datum();
  // draw a box around the cluster rectangle
  vis.append("rect")
    .attr('class', "clickedH")
    .attr("x", xMap1(d)-2)
    .attr("y", yMap(d)-2)
    .attr("width", xWide(d)+4)
    .attr("height", 14)
    .attr("fill", "white")
    .attr("stroke", "black")
    .attr("fill-opacity", 1)
    .attr("stroke-width", 2);

  // move the box to the background by promoting other objects
  rect.raise();
  d3.select("#svg-timetree")
      .selectAll("line")
      .raise();
  d3.select("#svg-timetree")
      .selectAll("text")
      .raise();
}

/**
 * Rectangular layout of tree, update nodes in place with x,y coordinates
 * @param {object} root
 */
function rectLayout(root) {
  // assign vertical positions to tips by postorder traversal
  var counter = 0;
  for (const node of traverse(root, 'postorder')) {
    if (node.children.length === 0) {
      // assign position to tip
      node.y = counter;
      counter++;
    } else {
      // ancestral node position is average of child nodes
      node.y = 0;
      for (var i = 0; i < node.children.length; i++) {
        var child = node.children[i];
        node.y += child.y;
      }
      node.y /= node.children.length;
    }
  }

  // assign horizontal positions by preorder traversal
  for (const node of traverse(root, 'preorder')) {
    if (node.parent === null) {
      // assign root to x=0
      node.x = 0.;
    } else {
      node.x = node.parent.x + node.branchLength;
    }
  }
}


/**
 * Draw time-scaled tree in SVG
 * @param {Object} timetree:  time-scaled phylogenetic tree imported as JSON
 * @returns {Array}  data frame
 */
function drawtree(timetree) {
  // generate tree layout (x, y coordinates
  rectLayout(timetree);

  var df = fortify(timetree),
    edgeset = edges(df, rectangular=true);

  // rescale SVG for size of tree
  var ntips = df.map(x => x.children.length === 0).reduce((x,y) => x+y);
  height = ntips*11 + margin.top + margin.bottom;
  vis.attr("height", height);
  yScale = d3.scaleLinear().range([height, 10]);  // add room for time axis

  // adjust d3 scales to data frame
  xScale.domain([
    d3.min(df, xValue)-0.05, d3.max(df, xValue)+0.2
  ]);
  yScale.domain([
    d3.min(df, yValue)-1, d3.max(df, yValue)+1
  ]);

  // draw lines
  vis.selectAll("lines")
    .data(edgeset)
    .enter().append("line")
    .attr("class", "lines")
    .attr("x1", xMap1)
    .attr("y1", yMap1)
    .attr("x2", xMap2)
    .attr("y2", yMap2)
    .attr("stroke-width", 1.5)
    .attr("stroke", "#777");

  return(df);
}


/**
 * Map cluster information to tips of the tree.
 * @param {Array} df: data frame extracted from time-scaled tree
 * @param {Array} clusters: data from clusters JSON
 * @returns {Array} subset of data frame annotated with cluster data
 */
function map_clusters_to_tips(df, clusters) {
  // extract accession numbers from phylogeny data frame
  var tips = df.filter(x => x.children.length===0),
      tip_labels = tips.map(x => x.thisLabel),  // accessions
      tip_stats;

  for (const cidx in clusters) {
    var cluster = clusters[cidx];
    if (cluster["nodes"].length === 1) {
      continue
    }

    // find variant in cluster that matches a tip label
    var labels = Object.keys(cluster["nodes"]),
        root = tip_labels.filter(value => value === cluster['lineage'])[0];
    if (root === undefined) {
      console.log("Failed to match cluster of index ", cidx, " to a tip in the tree");
      continue;
    }

    var root_idx = tip_labels.indexOf(root),  // row index in data frame
        root_xcoord = tips[root_idx].x;  // left side of cluster starts at end of tip

    // find most recent sample collection date
    var coldates = Array(),
        label, variant;

    for (var i=0; i<labels.length; i++) {
      label = labels[i];
      variant = cluster['nodes'][label];
      coldates = coldates.concat(variant.map(x => x.coldate));
    }
    coldates.sort();  // in place, ascending order

    var first_date = new Date(coldates[0]),
        last_date = new Date(coldates[coldates.length-1]);

    // augment data frame with cluster data
    tips[root_idx].cluster_idx = cidx;
    tips[root_idx].region = cluster.region;
    tips[root_idx].allregions = cluster.allregions;
    tips[root_idx].country = cluster.country;
    tips[root_idx].searchtext = cluster.searchtext;
    tips[root_idx].label1 = cluster["lineage"];
    tips[root_idx].count = coldates.length;
    tips[root_idx].varcount = labels.length;
    tips[root_idx].first_date = first_date;
    tips[root_idx].last_date = last_date;
    tips[root_idx].pdist = cluster.pdist;
    tips[root_idx].rdist = cluster.rdist;

    tips[root_idx].coldate = last_date;
    tips[root_idx].x1 = root_xcoord - ((last_date - first_date) / 3.154e10);
    tips[root_idx].x2 = root_xcoord;

    // map dbstats for lineage to tip
    tip_stats = dbstats["lineages"][cluster["lineage"]];
    tips[root_idx].max_ndiffs = tip_stats.max_ndiffs;
    tips[root_idx].mean_ndiffs = tip_stats.mean_ndiffs;
    tips[root_idx].nsamples = tip_stats.nsamples;
    tips[root_idx].residual = tip_stats.residual;
  }
  return tips;
}


/**
 * Convert x-coordinate of tree scale to Date scale by referring to a
 * tip in the tree.  Code from:
 * https://stackoverflow.com/questions/563406/add-days-to-javascript-date
 *
 * @param x:  float, distance from root in tree (years)
 * @param tip:  Object, representing a reference tip in the tree
 * @returns {string}  new date in ISO format (yyyy-mm-dd)
 */
function xaxis_to_date(x, tip) {
  var coldate = new Date(tip.coldate);  // collection date of reference tip
  coldate.setDate(coldate.getDate() + 365.25*(x - tip.x));
  return (coldate.toISOString().split('T')[0]);
}


/**
 * Add subtree objects to time-scaled tree.
 * @param {Array} tips, clusters that have been mapped to tips of tree
 */
function draw_clusters(tips) {
  axis.append("g")
    .attr("class", "treeaxis")
    .attr("transform", "translate(0,20)")
    .call(d3.axisTop(xScale)
      .ticks(3)
      .tickFormat(function(d) {
        return xaxis_to_date(d, tips[0])
      }));

  function mouseover(d) {
    d3.select("[cidx=cidx-" + d.cluster_idx + "]")
      .attr("txt_hover", "yes");

    cTooltip.transition()       // Show tooltip
            .duration(50)
            .style("opacity", 0.9);

    let ctooltipText = `<b>Mean pairwise distance:</b> ${d.pdist}<br/><b>Mean root distance:</b> ${d.rdist}<br>`;
    ctooltipText += region_to_string(d.allregions);
    ctooltipText += `<b>Number of variants:</b> ${d.varcount}<br>`;
    ctooltipText += `<b>Collection dates:</b><br>${formatDate(d.first_date)} / ${formatDate(d.last_date)}`;

    // Tooltip appears 10 pixels left of the cursor
    // Position tooltip based on the y-position of the cluster
    cTooltip.html(ctooltipText)
        .style("left", (d3.event.pageX + 15) + "px")
        .style("top", function(){
          if (d3.event.pageY > window.innerHeight/2) {
            return d3.event.pageY - cTooltip.node().getBoundingClientRect().height - 15 + "px";
          } else {
            return d3.event.pageY + 15 + "px";
          }
        });
  }

  vis.selectAll("rect")
    .data(tips)
    .enter()
    .lower()
    .append("rect")
    //.attr("selected", false)
    .attr("x", xMap1)
    .attr("y", yMap)
    .attr("width", xWide)
    .attr("height", 10)
    .attr("class", "default")
    .attr("cidx", function(d) { return "cidx-" + d.cluster_idx; })
    .attr("id", function(d, i) { return "id-" + i; })
    .on('mouseover', mouseover)
    .on("mouseout", function() {
      d3.select(this)
        .attr("txt_hover", null);

      cTooltip.transition()     // Hide tooltip
          .duration(50)
          .style("opacity", 0);
    })
    .on("click", function(d) {
      d3.selectAll("rect.clickedH").remove();

      beadplot(d.cluster_idx);

      if (isAccn($('#search-input').val())) {
        $('#search-input').val('');
      }

      search();

      // Reset search stats
      if ($('#search-input').val() !== "") {
	      var stats = search_stats.update({
                current_point: search_stats.get().start_idx[this.id],
                bead_indexer: 0,
        });
        update_search_stats(stats);
      }

      // reset all rectangles to high transparency
      if ($('#search-input').val() === "") {
        vis.selectAll("rect.clicked").attr('class', "default");
      }

      d3.select(this).attr("class", "clicked");

      draw_cluster_box(d3.select(this));

      $("#barplot").text(null);

      gentable(d);
      draw_region_distribution(d.allregions);

      // FIXME: this is the same div used for making barplot SVG
      $("#text-node").html(`Number of cases: ${d.count}<br/>Number of variants: ${d.varcount}<br/>`);

    });

  // generate colour palettes
  sample_pal = d3.scaleSequential(d3.interpolatePuBu)
      .domain(d3.extent(tips, function(d) { return d.nsamples; }));
  coldate_pal = d3.scaleSequential(d3.interpolateCividis)
      .domain(d3.extent(tips, function(d) { return d.last_date; }));
  diverge_pal = d3.scaleSequential(d3.interpolatePlasma)
      .domain(d3.extent(tips, function(d) { return d.residual; }));
  generate_legends();
  changeTreeColour();

  d3.select("#svg-timetree")
  .selectAll("line")
  .raise();

  vis.selectAll("text")
      .data(tips)
      .enter().append("text")
      .style("font-size", "10px")
      .attr("text-anchor", "start")
      .attr("alignment-baseline", "middle")
      .attr("cursor", "default")
      .attr("x", function(d) {
        return(xScale(d.x));
      })
      .attr("y", function(d) {
        return(yScale(d.y-0.15));
      })
      .text(function(d) { return(d.label1); })
      .on("mouseover", function(d) {
	mouseover(d);
      })
      .on("mouseout", function(d) {
        d3.select("[cidx=cidx-" + d.cluster_idx + "]").dispatch('mouseout');
      })
      .on("click", function(d) {
        d3.select("[cidx=cidx-" + d.cluster_idx + "]").dispatch('click');
      });
}


/**
 * Colour rect elements of tree to represent lineage attributes
 */
function changeTreeColour() {
  // hide legends, not knowing which one is showing
  $("#div-region-legend").hide();
  $("div#svg-sample-legend").hide();
  $("div#svg-coldate-legend").hide();
  $("div#svg-diverge-legend").hide();

  vis.selectAll("rect")
      .transition()
      .duration(300)
      .style("fill", function(d) {
        if (d !== undefined) {
          let opt = $("#select-tree-colours").val();
          if (opt === "Region") {
            $("#div-region-legend").show();
            return(country_pal[d.region]);
          }
          else if (opt === "No. samples") {
            $("div#svg-sample-legend").show();
            return(sample_pal(d.nsamples));  // placeholder values
          }
          else if (opt === "Collection date") {
            $("div#svg-coldate-legend").show();
            return(coldate_pal(d.last_date));
          }
          else {  // Divergence
            $("div#svg-diverge-legend").show();
            return(diverge_pal(d.residual));
          }
        }
      })
}

// bind to element
$("#select-tree-colours").change(function() {
  changeTreeColour();
});


// from https://observablehq.com/@d3/color-legend
function legend({
  color,
  title,
  tickSize = 6,
  width = 320,
  height = 44 + tickSize,
  marginTop = 18,
  marginRight = 0,
  marginBottom = 16 + tickSize,
  marginLeft = 0,
  ticks = width / 64,
  tickFormat,
  tickValues
} = {}) {
  const svg = d3.create("svg")
      .attr("width", width)
      .attr("height", height)
      .attr("viewBox", [0, 0, width, height])
      .style("overflow", "visible")
      .style("display", "block");

  let tickAdjust = g => g.selectAll(".tick line").attr("y1", marginTop + marginBottom - height);
  let x;

  // Continuous
  if (color.interpolate) {
    const n = Math.min(color.domain().length, color.range().length);

    x = color.copy().rangeRound(d3.quantize(d3.interpolate(marginLeft, width - marginRight), n));

    svg.append("image")
        .attr("x", marginLeft)
        .attr("y", marginTop)
        .attr("width", width - marginLeft - marginRight)
        .attr("height", height - marginTop - marginBottom)
        .attr("preserveAspectRatio", "none")
        .attr("xlink:href", ramp(color.copy().domain(d3.quantize(d3.interpolate(0, 1), n))).toDataURL());
  }

  // Sequential
  else if (color.interpolator) {
    x = Object.assign(color.copy()
        .interpolator(d3.interpolateRound(marginLeft, width - marginRight)),
        {range() { return [marginLeft, width - marginRight]; }});

    svg.append("image")
        .attr("x", marginLeft)
        .attr("y", marginTop)
        .attr("width", width - marginLeft - marginRight)
        .attr("height", height - marginTop - marginBottom)
        .attr("preserveAspectRatio", "none")
        .attr("xlink:href", ramp(color.interpolator()).toDataURL());

    // scaleSequentialQuantile doesn’t implement ticks or tickFormat.
    if (!x.ticks) {
      if (tickValues === undefined) {
        const n = Math.round(ticks + 1);
        tickValues = d3.range(n).map(i => d3.quantile(color.domain(), i / (n - 1)));
      }
      if (typeof tickFormat !== "function") {
        tickFormat = d3.format(tickFormat === undefined ? ",f" : tickFormat);
      }
    }
  }

  // Threshold
  else if (color.invertExtent) {
    const thresholds
        = color.thresholds ? color.thresholds() // scaleQuantize
        : color.quantiles ? color.quantiles() // scaleQuantile
        : color.domain(); // scaleThreshold

    const thresholdFormat
        = tickFormat === undefined ? d => d
        : typeof tickFormat === "string" ? d3.format(tickFormat)
        : tickFormat;

    x = d3.scaleLinear()
        .domain([-1, color.range().length - 1])
        .rangeRound([marginLeft, width - marginRight]);

    svg.append("g")
      .selectAll("rect")
      .data(color.range())
      .join("rect")
        .attr("x", (d, i) => x(i - 1))
        .attr("y", marginTop)
        .attr("width", (d, i) => x(i) - x(i - 1))
        .attr("height", height - marginTop - marginBottom)
        .attr("fill", d => d);

    tickValues = d3.range(thresholds.length);
    tickFormat = i => thresholdFormat(thresholds[i], i);
  }

  // Ordinal
  else {
    x = d3.scaleBand()
        .domain(color.domain())
        .rangeRound([marginLeft, width - marginRight]);

    svg.append("g")
      .selectAll("rect")
      .data(color.domain())
      .join("rect")
        .attr("x", x)
        .attr("y", marginTop)
        .attr("width", Math.max(0, x.bandwidth() - 1))
        .attr("height", height - marginTop - marginBottom)
        .attr("fill", color);

    tickAdjust = () => {};
  }

  svg.append("g")
      .attr("transform", `translate(0,${height - marginBottom})`)
      .call(d3.axisBottom(x)
        .ticks(ticks, typeof tickFormat === "string" ? tickFormat : undefined)
        .tickFormat(typeof tickFormat === "function" ? tickFormat : undefined)
        .tickSize(tickSize)
        .tickValues(tickValues))
      .call(tickAdjust)
      .call(g => g.select(".domain").remove())
      .call(g => g.append("text")
        .attr("x", marginLeft)
        .attr("y", marginTop + marginBottom - height - 6)
        .attr("fill", "currentColor")
        .attr("text-anchor", "start")
        .attr("font-weight", "bold")
        .attr("class", "title")
        .text(title));

  return svg.node();
}


function ramp(color, n = 256) {
  // https://stackoverflow.com/questions/60443356/legend-not-appearing-when-using-document-createelementcanvas
  const canvas = document.createElement('canvas');
  const context = canvas.getContext("2d");
  d3.select(canvas).attr("width", n)
    .attr("height", 1);
  for (let i = 0; i < n; ++i) {
    context.fillStyle = color(i / (n - 1));
    context.fillRect(i, 0, 1, 1);
  }
  return canvas;
}


function generate_legends() {
  // region legend with swatches
  let s = `<div style="display: flex; align-items: center; margin-left: 0px; padding-top: 6px; min-height: 33px; font: 10px sans-serif;">`;
  s += `<div style="width: 100%; columns: 60px;">`;
  for (const [key, value] of Object.entries(country_pal)) {
    s += `<div class="legend-item">`;
    s += `<div class="legend-swatch" style="background:${value};"></div>`;
    s += `<div class="legend-label" title="${key}">${key}</div>`;
    s += `</div>`;
  }
  s += `</div></div>`;
  $("#div-region-legend").html(s).hide();

  // sample size legend
  $("div#svg-sample-legend").html(legend({
    color: sample_pal,
    title: "Sample size",
    width: 240
  })).hide();

  // collection date legend
  var [millsec0, millsec1] = coldate_pal.domain(),
      days = (millsec1 - millsec0) / 8.64e7;
  $("div#svg-coldate-legend").html(legend({
    color: d3.scaleSequential([-days, 0], d3.interpolateCividis),
    title: "Last sample date (days before present)",
    width: 240
  })).hide();

  // divergence legend
  $("div#svg-diverge-legend").html(legend({
    color: diverge_pal,
    title: "Divergence (normalized residual from clock rate)",
    width: 240
  })).hide();
}
