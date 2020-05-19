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
 *
 * @param {Number} cid:  integer index of cluster to draw
 */
function beadplot(cid) {
  var cluster = clusters[cid];
  if (cluster === undefined) {
    alert("Unknown cluster ID " + cid);
  }

  var variant, coldates, y = 1,
      node_to_y = Object();
  for (const accn in cluster.nodes) {
    variant = cluster.nodes[accn];
    coldates = variant.map(x => xfunc(x.coldate));
    coldates.sort();

  }
}
