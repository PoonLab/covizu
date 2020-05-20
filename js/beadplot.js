var marginB = {top: 10, right: 10, bottom: 10, left: 10},
    widthB = 900 - marginB.left - marginB.right,
    heightB = 500 - marginB.top - marginB.bottom;

// set up plotting scales
var xValueB = function(d) { return d.x; },
  xScaleB = d3.scaleLinear().range([0, widthB]),
  xMap1B = function(d) { return xScaleB(d.x1); },
  xMap2B = function(d) { return xScaleB(d.x2); };

var yValueB = function(d) { return d.y; },
  yScaleB = d3.scaleLinear().range([heightB, 0]),  // inversion
  yMap1B = function(d) { return yScaleB(d.y1); },
  yMap2B = function(d) { return yScaleB(d.y2); },
  yMapB = function(d) { return yScaleB(yValueB(d)); };


var visB = d3.select("div#svg-cluster")
  .append("svg")
  .attr("width", widthB + marginB.left + marginB.right)
  .attr("height", heightB + marginB.top + marginB.bottom)
  .append("g");


/**
 * Parse node and edge data from clusters JSON to a format that is
 * easier to map to SVG.
 * @param {Object} clusters:
 */
function parse_clusters(clusters) {
  var cluster, variant, coldates, variants, edgelist, beaddata = [];

  for (var cidx=0; cidx < clusters.length; cidx++) {
    variants = [];
    cluster = clusters[cidx];
    if (cluster.nodes.length == 1) {
      console.log('skip '+ cluster.nodes);
      continue;
    }

    // extract the date range for each variant in cluster
    var count = 1;
    for (const accn in cluster.nodes) {
      variant = cluster.nodes[accn];

      coldates = variant.map(x => new Date(x.coldate));
      coldates.sort();
      variants.push({
        'accession': accn,
        'label': variant[0].label1,
        'mindate': coldates[0],
        'maxdate': coldates[coldates.length-1],
        'y': count++  // assign then increment
      });
    }

    // map earliest collection dates to vertical edges
    var edge, parent, child;
    edgelist = [];
    for (var e = 0; e < cluster.edges.length; e++) {
      edge = cluster.edges[e];
      parent = variants.filter(x => x.accession == edge[0])[0];
      child = variants.filter(x => x.accession == edge[1])[0];
      edgelist.push({
        'parent': parent.accession,
        'child': child.accession,
        'mindate': parent.mindate < child.mindate ? parent.mindate : child.mindate
      });
    }

    beaddata.push({'variants': variants, 'edgelist': edgelist})
  }

  return(beaddata);
}
