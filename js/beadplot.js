/**
 *  Configure SVG to display beadplots
 */
var marginB = {top: 50, right: 10, bottom: 50, left: 10},
    widthB = document.getElementById("svg-cluster").clientWidth - marginB.left - marginB.right,
    heightB = 1000 - marginB.top - marginB.bottom,
    pixelsPerDay = 10,
    ccid = -1;

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


function parse_mutation_annotations(mut_annotations) {
  var mutations = [];

  // Sorts tips according to cluster_idx
  sorted_tips = [...tips, ...recombinant_tips]
  sorted_tips.sort(function(a,b) {
    return a.cluster_idx - b.cluster_idx
  });

  for (const cidx in sorted_tips) {
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

  d3.selectAll(".selectionLH").attr("class","lines hLine");
  d3.selectAll(".selectionL").attr("class","lines vLine");
  d3.selectAll(".selectionLC").attr("class","default");
}

/**
 * Updates beadplot when the expand switch is selected or when a resize event occurs
 */
function expand() {
  var mindate = d3.min(variants, xValue1B),
      maxdate = d3.max(variants, xValue2B),
      spandate = maxdate-mindate,  // in milliseconds
      numDays = d3.timeDay.count(mindate, maxdate),
      clientWidth = document.getElementById("svg-cluster").clientWidth;

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
    $('#svg-clusteraxis').css('padding-bottom', 0)
  }
  else {
    currentWidth = clientWidth - marginB.left -  marginB.right;
    $('#svg-clusteraxis').css('padding-bottom', $('#inner-hscroll').height())

  }

  xScaleB = d3.scaleLinear().range([0, currentWidth]);
  var xAxis = d3.scaleTime().range([0, currentWidth]);

  // allocate space for text labels on left
  xScaleB.domain([mindate - 170*(spandate/currentWidth), maxdate]);
  xAxis.domain([mindate - 170*(spandate/currentWidth), maxdate]);

  // update horizontal range
  $("#svg-cluster > svg").attr("width",currentWidth + marginB.left + marginB.right);
  $("#inner-hscroll").css("width",currentWidth + marginB.left + marginB.right);
  
  // Add 15px to account for the scrollbar for the beadplot
  $("#svg-clusteraxis > svg").attr("width", currentWidth + marginB.left + marginB.right + 15);

  // draw vertical line segments that represent edges in NJ tree
  visB.selectAll(".vLine")
      .transition().duration(700)
      .ease(d3.easeExp)
      .attr("x1", xMap1B)
      .attr("x2", xMap2B)
      .attr("y1", yMap1B)
      .attr("y2", yMap2B);

  // draw horizontal line segments that represent variants in cluster
  visB.selectAll(".hLine")
      .transition().duration(700)
      .ease(d3.easeExp)
      .attr("x1", xMap1B)
      .attr("x2", xMap2B)
      .attr("y1", yMap1B)
      .attr("y2", yMap2B);

  // label variants with earliest sample name
  visB.selectAll("text")
      .transition().duration(700)
      .ease(d3.easeExp)
      .attr("x", function(d) { return(xScaleB(d.x1)-5); })
      .attr("y", function(d) { return(yScaleB(d.y1)); })

  // draw "beads" to represent samples per collection date
  visB.selectAll("circle:not(.selectionH)")
      .transition().duration(700)
      .ease(d3.easeExp)
      .attr("cx", xMapB)
      .attr("cy", yMapB);

  var selectionH = d3.selectAll("circle.selectionH").nodes(),
      selectionLC = d3.selectAll("circle.selectionLC").nodes();

  // Add animation to selection on bead if it exists
  if (selectionH.length > 0) {
    visB.selectAll("circle.selectionH")
        .data(points.filter(x => x.accessions[0] === selectionH[0].attributes.bead.nodeValue))
        .transition().duration(700)
        .ease(d3.easeExp)
        .attr("cx", function(d) {
          return xScaleB(d.x)
        });
    setTimeout(function() {
      $('#svg-cluster').animate({
        scrollTop: selectionH[0].cy.baseVal.value - document.getElementById("svg-cluster").clientHeight/2,
        scrollLeft: selectionH[0].cx.baseVal.value - document.getElementById("svg-cluster").clientWidth/2,
      }, 200);
    }, 700);
  }

  if (selectionLC.length > 0) {
    setTimeout(function() {
      $('#svg-cluster').animate({
        scrollTop: selectionLC[0].cy.baseVal.value - document.getElementById("svg-cluster").clientHeight/2,
        scrollLeft: selectionLC[0].cx.baseVal.value - document.getElementById("svg-cluster").clientWidth/2,
      }, 200);
    }, 700);
  }

  var tickCount = 0.005*currentWidth;
  if (tickCount <= 4) tickCount = 4;

  // Ensures that there aren't extra ticks in the time axis of the beadplot when the first and last coldates are within a few days
  var dateDiff = d3.timeDay.count(
                  d3.min(variants, function(d) {return d.x1}), 
                  d3.max(variants, function(d) {return d.x2}));
  
  if (dateDiff !=0 && dateDiff < tickCount) tickCount = dateDiff;

  // draw x-axis
  visBaxis.selectAll(".treeaxis")
      .transition().duration(700)
      .ease(d3.easeExp)
      .call(
        d3.axisTop(xAxis)
          .ticks(tickCount)
          .tickFormat(d3.timeFormat("%Y-%m-%d"))
      );
}


/**
 * Draw the beadplot for a specific cluster (identified through its
 * integer index <cid>) in the SVG.
 * @param {Number} cid:  integer index of cluster to draw as beadplot
 */
 async function beadplot(cid) {
  // Update global cindex for SVG and NWK filenames
  if (cindex !== ccid) {
    cindex = cid;
    ccid = cindex
    edgelist = await getdata(`/api/edgelist/${cindex}`);
    edgelist.forEach(x => {
      x.x1 = utcDate(x.x1),
      x.x2 = utcDate(x.x2)
    });
    points = await getdata(`/api/points/${cindex}`);
    points.forEach(d => {
      d.x = utcDate(d.x)
    });
    variants = await getdata(`/api/variants/${cindex}`);
    variants.forEach(x => {
      x.x1 = utcDate(x.x1),
      x.x2 = utcDate(x.x2)
    });
    await fetch(`/api/lineage/${cindex}`)
    .then(response => response.text())
    .then(lin => lineage=lin)
  } 

  // rescale slider
  let max_dist = Math.max(...edgelist.map(x => x.dist));
  if (isFinite(max_dist) == false) {
    max_dist = 2.0;
  }
  let slider = $("#vedge-slider");
  slider.slider("option", "max", max_dist )
        .slider("value", max_dist < 2.0 ? max_dist : 2.0)
  move_arrow();
  $( "#custom-handle" ).text(max_dist < 2.0 ? max_dist.toString() : "2.0");

  // Beadplots aren't expanded by default
  if ($('#expand-option').attr('checked')) {
    $('.switch').trigger('click');
    $('#expand-option').removeAttr('checked')
  }

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
                                          ($('#expand-option').attr('checked') ? $('#inner-hscroll').height() : $('#inner-hscroll').height() * 2));

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
    $('#svg-clusteraxis').css('padding-bottom', 0)
  }
  else {
    currentWidth = clientWidth - marginB.left -  marginB.right;
    $('#svg-clusteraxis').css('padding-bottom', $('#inner-hscroll').height())

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
  var vLines = visB.selectAll("lines")
      .data(edgelist.filter(x => x.dist <= slider.slider("value")))
      .enter().append("line")

  draw_vertical_edges(vLines);

  // draw horizontal line segments that represent variants in cluster
  visB.selectAll("lines")
      .data(variants)
      .enter().append("line")
      .attr("class", "lines hLine")
      .attr("id", function(d) {
        return d.label.replace(/'/g, '-').replace(/\./g,'-').replace('/', '-').replace(' ', '_');
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
          draw_region_distribution(d.region);
          let var_samples = points.filter(x => x.y === d.y1);
          gen_details_table(var_samples);
        }

        // Clear previous selection
        clear_selection();

        // Visual feedback
        var edge = d3.select(this);
        edge.attr("class", "lines hLine selectionLH");

        var filtered_edgelist = edgelist.filter(x => x.dist <= slider.slider("value"))

        // Bold parent and child variants 
        for (i=0; i < filtered_edgelist.length; i++) {
          let edge_label = edge.data()[0].label
          if(filtered_edgelist[i].child == edge_label) {
            parent_variant = filtered_edgelist[i].parent.replace(/\./g,'-').replace('/', '-').replace(' ', '_')
            d3.select(`#${parent_variant}`).attr("stroke-width", 5).attr("class", "lines hLine selectionLH");
          }
          else if(filtered_edgelist[i].parent == edge_label) {
            child_variant = filtered_edgelist[i].child.replace(/\./g,'-').replace('/', '-').replace(' ', '_')
            d3.select(`#${child_variant}`).attr("stroke-width", 5).attr("class", "lines hLine selectionLH");
          }
        }

        // Bold edges connecting clicked variant to parent/children
        d3.selectAll(`.vLine[y1="${yScaleB(edge.data()[0].y1)}"]`).attr("stroke-width", 3).attr("class", "lines hLine selectionL")
        d3.selectAll(`.vLine[y2="${yScaleB(edge.data()[0].y1)}"]`).attr("stroke-width", 3).attr("class", "lines hLine selectionL")

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
        draw_region_distribution(d.region);
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

      slider_update();
    }
  });
}


/**
 * Filters vertical lines in the beadplot when the slider is moved
 */
async function slider_update() {
  var slider = $("#vedge-slider");

  visB.selectAll(".vLine").remove();

  var vLines = visB.selectAll("lines")
                .data(edgelist.filter(x => x.dist <= slider.slider("value")))
                .enter().insert("line", ".hLine");

  draw_vertical_edges(vLines)
}


function draw_vertical_edges(vLines) {
  vLines
    .attr("class", "lines vLine")
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

      let parent_variant = d3.select(".lines#"+d.parent.replace(/'/g, '-').replace(/\./g,'-').replace('/', '-').replace(' ', '_')),
          child_variant = d3.select(".lines#"+d.child.replace(/'/g, '-').replace(/\./g,'-').replace('/', '-').replace(' ', '_'));

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

      let parent_variant = d3.select(".lines#"+d.parent.replace(/'/g, '-').replace(/\./g,'-').replace('/', '-').replace(' ', '_')),
          child_variant = d3.select(".lines#"+d.child.replace(/'/g, '-').replace(/\./g,'-').replace('/', '-').replace(' ', '_'));

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
      edge.attr("class", "lines vLine selectionL");

      let parent_variant = d3.select(".lines#"+d.parent.replace(/'/g, '-').replace(/\./g,'-').replace('/', '-').replace(' ', '_')),
          child_variant = d3.select(".lines#"+d.child.replace(/'/g, '-').replace(/\./g,'-').replace('/', '-').replace(' ', '_'));

      parent_variant.attr("class", "lines hLine selectionLH");
      child_variant.attr("class", "lines hLine selectionLH");

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
      for (let i = 0; i < obj[j].mutation.length; i++) {
        let sample_details = [
          obj[j].mutation[i],
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
    for (let i = 0; i < obj.mutation.length; i++) {
      let sample_details = [
        obj.mutation[i],
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
  var phenotype_icons = {
    'vaccine_neutralization_efficacy': 'img/red_circle.png', 
    'anthropozoonotic_events': 'img/bat.png', 
    'gene_expression_increase': 'img/orange_star.png', 
    'ACE2_receptor_binding_affinity': 'img/purple_square.jpeg',
    'monoclonal_antibody_serial_passage_escape': 'img/antibody.png', 
    'convalescent_plasma_escape': 'img/green_pentagon.png', 
    'antibody_epitope_effects': 'img/blue_triangle.png'
  }

  for (let j = 0; j < obj.phenotype.length; j++) {
    var phenotype = t_cells[j].children[2];
      for (let i = 0; i < obj.phenotype[j].length; i++) {
        for (const [label, link] of Object.entries(phenotype_icons)) {
          if(obj.phenotype[j][i] == label) {
            var img = document.createElement("img");
            img.src = link;
            img.classList.add('phenotype_icon')
            phenotype.appendChild(img);
          }
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
  const regions = unique(Object.values(region_map)).sort()
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
  var root = edgelist[0].parent;
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
