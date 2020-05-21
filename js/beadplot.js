var marginB = {top: 50, right: 50, bottom: 50, left: 50},
    widthB = 900 - marginB.left - marginB.right,
    heightB = 900 - marginB.top - marginB.bottom;

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


// https://stackoverflow.com/questions/1960473/get-all-unique-values-in-a-javascript-array-remove-duplicates
function onlyUnique(value, index, self) {
  return self.indexOf(value) === index;
}

/**
 * Parse node and edge data from clusters JSON to a format that is
 * easier to map to SVG.
 * @param {Object} clusters:
 */
function parse_clusters(clusters) {
  var cluster, variant, coldates,
      variants,  // horizontal line segment data + labels
      edgelist,  // vertical line segment data
      points,  // the "beads"
      beaddata = [];  // return value

  for (var cidx = 0; cidx < clusters.length; cidx++) {
    cluster = clusters[cidx];
    if (cluster.nodes.length == 1) {
      console.log('skip '+ cluster.nodes);
      continue;
    }

    // deconvolute edge list to get node list in preorder
    var nodelist = cluster.edges.flat().filter(onlyUnique);

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
        'label': variant[0].label1,
        'x1': new Date(coldates[0]),  // min date
        'x2': new Date(coldates[coldates.length-1]),  // max date
        'count': coldates.length,
        'y1': y,  // horizontal line segment
        'y2': y
      });

      var isodates = coldates.filter(onlyUnique),
          isodate, count;

      for (var i=0; i<isodates.length; i++) {
        isodate = isodates[i];
        count = coldates.filter(x => x == isodate).length;
        // TODO: store labels, country data for each point
        points.push({
          'x': new Date(isodate),
          'y': y,
          'count': count
        })
      }
      y++;
    }

    // map earliest collection date of child node to vertical edges
    var edge, parent, child;
    edgelist = [];
    for (var e = 0; e < cluster.edges.length; e++) {
      edge = cluster.edges[e];
      parent = variants.filter(x => x.accession == edge[0])[0];
      child = variants.filter(x => x.accession == edge[1])[0];
      edgelist.push({
        'y1': parent.y1,
        'y2': child.y1,
        'x1': child.x1,  // vertical line segment
        'x2': child.x1
      });

      // update variant time range
      if (parent.x1 > child.x1) {
        parent.x1 = child.x1;
      }
      if (parent.x2 < child.x1) {
        parent.x2 = child.x1;
      }
    }

    beaddata.push({'variants': variants, 'edgelist': edgelist, 'points': points})
  }

  return(beaddata);
}


/**
 * Draw the beadplot for a specific cluster (identified through its
 * integer index <cid>) in the SVG.
 * @param cid
 */
function beadplot(cid) {
  console.log(cid);

  var variants = beaddata[cid].variants,
      edgelist = beaddata[cid].edgelist,
      points = beaddata[cid].points;

  // set plotting domain
  xScaleB.domain([
    d3.min(variants, xValue1B), d3.max(variants, xValue2B)
  ]);
  yScaleB.domain([
    d3.min(variants, yValue1B)-1, d3.max(variants, yValue1B)+1
  ]);

  visB.selectAll('*').remove();

  visB.selectAll("lines")
    .data(variants)
    .enter().append("line")
    .attr("class", "lines")
    .attr("x1", xMap1B)
    .attr("x2", xMap2B)
    .attr("y1", yMap1B)
    .attr("y2", yMap2B)
    .attr("stroke-width", 3)
    .attr("stroke", "#777");

  visB.selectAll("text")
    .data(variants)
    .enter().append("text")
    .style("font-size", "10px")
    .attr("text-anchor", "end")
    .attr("alignment-baseline", "middle")
    .attr("x", function(d) { return(xScaleB(d.x1)); })
    .attr("y", function(d) { return(yScaleB(d.y1)); })
    .text(function(d) { return(d.label); });

  visB.selectAll("lines")
    .data(edgelist)
    .enter().append("line")
    .attr("class", "lines")
    .attr("x1", xMap1B)
    .attr("x2", xMap2B)
    .attr("y1", yMap1B)
    .attr("y2", yMap2B)
    .attr("stroke-width", 1)
    .attr("stroke", "#777");

  visB.selectAll("circle")
    .data(points)
    .enter().append("circle")
    .attr("r", function(d) { return (4*Math.sqrt(d.count)); })
    .attr("cx", xMapB)
    .attr("cy", yMapB)
    .attr("fill", "white")
    .attr("stroke", "black");
}