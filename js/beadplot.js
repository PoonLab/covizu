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

    // extract the date range for each variant in cluster
    var y = 1;
    variants = [];
    points = [];
    for (const accn in cluster.nodes) {
      variant = cluster.nodes[accn];
      coldates = variant.map(x => x.coldate);
      coldates.sort();

      variants.push({
        'accession': accn,
        'label': variant[0].label1,
        'mindate': new Date(coldates[0]),  // x1
        'maxdate': new Date(coldates[coldates.length-1]),  // x2
        'count': coldates.length,
        'y': y
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
        'y1': parent.y,
        'y2': child.y,
        'mindate': child.mindate  // x-coordinate
      });
    }

    beaddata.push({'variants': variants, 'edgelist': edgelist, 'points': points})
  }

  return(beaddata);
}
