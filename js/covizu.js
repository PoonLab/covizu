$( function() {
  $( "#tabs" ).tabs();
} );

$(document).tooltip({show: null});
$("#loading_text").text(``);
$("#loading").hide();

/*********************** Session ID check ***********************/
var sid = getUrlParameter('sid') // Load vars from url

payload = jsonToURI({"cmd":"state/session/validate",
 "api": {"version":1},
 "client_id":"cid-e9418c5b4b6e",
 "sid": sid})

$.post('https://gpsapi.epicov.org/epi3/gps_api?req='+ payload, function(data, status){
  //Not logged in
  if (data.rc != "ok"){
    var r = confirm('Unable to verify session credentials. Please access app through platform.gisaid.org. Press "OK" to redirect to GISAID homepage.')
    if (r == true){
      window.location.href = 'https://platform.gisaid.org'
      }
    throw new Error('Forbidden')
  }
});

/*********************** DIALOGS ***********************/


$( "#splash" ).dialog({
  autoOpen: true,
  width: 600,
  buttons: [{
    id: "splash-button",
    text: "Okay!",
    disabled: true,
    click: function() {
      $( this ).dialog( "close" );
    }
  }]
});

$("#splash-button").parents(".ui-dialog-buttonpane")
        .append("<div style=\"color: #4F2683\" id='splash-extra'>Loading JSON data from server (~10s)..." +
                "<img style=\"vertical-align: middle\" width=\"36px\" src=\"img/Loading_icon_GISAID.gif\"/>" +
                "</div>");

$( "#help-timetree" ).dialog({
  autoOpen: false,
  width: 600,
  buttons: {
    "Got it!": function() {
      $( this ).dialog( "close" );
    }
  }
});

$( "#help-beadplot" ).dialog({
  autoOpen: false,
  width: 600,
  buttons: {
    "Got it!": function() {
      $( this ).dialog( "close" );
    }
  }
});

$( "#help-search" ).dialog({
  autoOpen: false,
  width: 600,
  buttons: {
    "Got it!": function() {
      $( this ).dialog( "close" );
    }
  }
});

$( "#dialog" ).dialog({
  autoOpen: false,
  width: 600
});

$( "#ack-open" ).click(function() {
  $( "#dialog" ).dialog( "open" );
});


/*********************** LOAD JSON DATA ***********************/

// Forces the requested files to not be cached by the browser
$.ajaxSetup({ 
  cache: false 
});

// load database statistics
var dbstats, req;
req = $.getJSON("data/dbstats.json", function(data) {
  dbstats = data;
});
req.done(function() {
  $("#div-last-update").text(`Last update: ${dbstats.lastupdate}`);
  $("#div-number-genomes").text(`High-quality genomes: ${dbstats.noseqs}`);
});

// colors from https://personal.sron.nl/~pault/#sec:qualitative
var country_pal = {
  "Africa": "#EEDD88",
  "Asia": "#EE8866",
  "Europe": "#44BB99",
  "North America": "#99DDFF",
  "Oceania": "#FFAABB",
  "South America": "#AAAA00"
};

// load time-scaled phylogeny from server
var nwk, df, df_xbb, countries, region_map;
$.ajax({
  url: "data/timetree.nwk",
  success: function(data) {
    nwk = data;
  }
});

var clusters, tips, recombinant_tips, xbb_tips,
    accn_to_cid, cindex, lineage_to_cid, lineage;
var edgelist = [], points = [], variants = []
var map_cidx_to_id = [], id_to_cidx = [], display_id = {};

req = $.when(
  $.getJSON("/epicov_api/tips", function(data) {
    tips = data;
    tips.forEach(x => {
      x.first_date = new Date(x.first_date)
      x.last_date = new Date(x.last_date)
      x.coldate = new Date(x.coldate)
      x.mcoldate = new Date(x.mcoldate)
    });
  }),
  $.getJSON("/epicov_api/recombtips", function(data) {
    recombinant_tips = data;
    recombinant_tips.forEach(x => {
      x.first_date = new Date(x.first_date)
      x.last_date = new Date(x.last_date)
      x.coldate = new Date(x.coldate)
      x.mcoldate = new Date(x.mcoldate)
    });
  }),
  $.getJSON("/epicov_api/df", function(data) {
    df = data;
    df.forEach(x => {
      x.first_date = x.first_date ? new Date(x.first_date) : undefined
      x.last_date = x.last_date ? new Date(x.last_date) : undefined
      x.coldate = x.coldate ? new Date(x.coldate) : undefined
      x.mcoldate = x.coldate ? new Date(x.mcoldate) : undefined
    });
  }),
  $.getJSON("/epicov_api/xbb", function(data) {
    df_xbb = data;
    df_xbb.forEach(x => {
      x.first_date = x.first_date ? new Date(x.first_date) : undefined
      x.last_date = x.last_date ? new Date(x.last_date) : undefined
      x.coldate = x.coldate ? new Date(x.coldate) : undefined
      x.mcoldate = x.coldate ? new Date(x.mcoldate) : undefined
    });
  }),
  $.getJSON("/epicov_api/regionmap", function(data) {
    region_map = data;
  })
);


req.done( async function() {
  $("#splash-button").button("enable");
  $("#splash-extra").html("");  // remove loading animation

  // Maps id to a cidx
  const reverse_recombinant_tips = [...recombinant_tips].reverse()
  var i = 0, first, last;
  first = i;
  for (const index in reverse_recombinant_tips) {
    id_to_cidx[i++] = 'cidx-' + reverse_recombinant_tips[index].cluster_idx;
  }
  last = i - 1;
  display_id["other_recombinants"] = {"first": first, "last": last};

  first = i;
  xbb_tips = df_xbb.filter(x=>x.isTip);
  for (const index in xbb_tips) {
    id_to_cidx[i++] = 'cidx-' + xbb_tips[index].cluster_idx;
  }
  last = i - 1;
  display_id["xbb"] = {"first": first, "last": last};

  first = i;
  for (const index in tips) {
    id_to_cidx[i++] = 'cidx-' + tips[index].cluster_idx;
  }
  last = i - 1;
  display_id["non_recombinants"] = {"first": first, "last": last};

  // Maps cidx to an idx
  const reverseMapping = o => Object.keys(o).reduce((r, k) => Object.assign(r, { [o[k]]: (r[o[k]] || parseInt(k)) }), {})
  map_cidx_to_id = reverseMapping(id_to_cidx)

  switch($("#display-tree").val()) {
    case "XBB Lineages":
      drawtree(df_xbb);
      draw_clusters(tips);
      break;
    case "Other Recombinants":
    default:
      drawtree(df);
      draw_clusters(tips);
  }

  var rect = d3.selectAll("#svg-timetree > svg > rect"),
      node = rect.nodes()[rect.size()-1];

  // initial display
  // d3.select(node).dispatch("click");
  cindex = node.__data__.cluster_idx;
  d3.select(node).attr("class", "clicked");
  await beadplot(node.__data__.cluster_idx);
  $("#barplot").text(null);
  gentable(node.__data__);
  draw_region_distribution(node.__data__.region);
  gen_details_table(points);  // update details table with all samples
  draw_cluster_box(d3.select(node));

  /*
  rect = d3.selectAll("#svg-cluster > svg > g > circle");
  node = rect.nodes()[0];
  d3.select(node).dispatch("click");//.dispatch("mouseover");
   */

  // Maps lineage to a cidx
  await fetch(`/epicov_api/lineagetocid`)
  .then(response => response.json())
  .then(data => lineage_to_cid = data)

  $('#search-input').autocomplete({
    source: function(req, res) {
      $.ajax({
        url: `/epicov_api/getHits/${req.term}`,
        dataType: "json",
        type: "GET",
        data: {
          term: req.term
        },
        success: function(data) {
          res(data)
        },
        error: function(xhr) {
          console.log(xhr.statusText)
        }
      })
    },
    minLength: 1,
    delay: 0,
    select: function( event, ui ) {
        const accn = ui.item.value;
        //search(accn);
    }
  });


  /***********************  SEARCH INTERFACE ***********************/

  // Enable and Disable "Search", "Clear" "Next" and "Previous" buttons when needed
  function disable_buttons() {
    $('#search-button').attr("disabled", true);
    $('#clear_button').attr("disabled", true);
    $('#next_button').attr("disabled", true);
    $('#previous_button').attr("disabled", true);
    $('#search_stats').addClass("disabled_stats");
  }

  function enable_buttons() {
    if (search_results.get().total_points > 0) {
      $('#next_button').removeAttr("disabled");
      $('#previous_button').removeAttr("disabled");
      $('#search_stats').removeClass("disabled_stats");
    } else {
      $('#next_button').attr("disabled", true);
      $('#previous_button').attr("disabled", true);
      $('#search_stats').addClass("disabled_stats");
    }
  }

  disable_buttons();

  // Enables "search" and "clear" buttons if the input fields are not empty
  $('#search-input').on('change keyup search', function() {
    if ($('#search-input').val() != "" || $('#start-date').val() != "" || $('#end-date').val() != "") {
      $('#search-button').removeAttr("disabled");
      $('#clear_button').removeAttr("disabled");
    }
    else {
      clear_selection();
      $('#search_stats').text(`0 of 0 points`);
      disable_buttons();
    }
  });

  $('#search-input, #start-date, #end-date').on('keydown', async function(e) {
    $('#error_message').text(``);
    // Only resets search results if the backspace key is pressed
    if (search_results.get().total_points > 0 && (e.keyCode == 8)) {
      clear_selection();
      $('#search_stats').text(`0 of 0 points`);
      disable_buttons();
    }
    if (e.keyCode == 13 && ($('#search-input').val() != "" || $('#start-date').val() != "" || $('#end-date').val() != "")) {
      // type <enter> to run search
      // run_search();
      $('#error_message').text(``);
      $("#loading").show();
      $("#loading_text").text(`Loading. Please Wait...`);
      await wrap_search();
      enable_buttons();
      $("#loading").hide();
      $("#loading_text").text(``);
    }
  });

  // Removes the error message when the user clicks on the date picker
  $('#start-date').on('mousedown', function() {
    $('#error_message').text(``);
  });

  $('#end-date').on('mousedown', function() {
    $('#error_message').text(``);
  });

  $('#start-date').on('change keyup', function() {
    if ($('#start-date').val() != "" || $('#search-input').val() != "" || $('#end-date').val() != "") {
      $('#search-button').removeAttr("disabled");
      $('#clear_button').removeAttr("disabled");
    }
    else {
      clear_selection();
      $('#search_stats').text(`0 of 0 points`);
      disable_buttons();
    }
  });

  $('#end-date').on('change keyup', function() {
    if ($('#end-date').val() != "" || $('#search-input').val() != "" || $('#start-date').val() != "") {
      $('#search-button').removeAttr("disabled");
      $('#clear_button').removeAttr("disabled");
    }
    else {
      clear_selection();
      $('#search_stats').text(`0 of 0 points`);
      disable_buttons();
    }
  })

  const dateFormat = 'yy-mm-dd'; // ISO_8601
  $('#start-date').datepicker({
    dateFormat,
    onSelect: function(date_text){
      const start = utcDate(date_text);
      if ($('#start-date').val() != "") {
        $('#search-button').removeAttr("disabled");
        $('#clear_button').removeAttr("disabled");
      }
      else {
        clear_selection();
        $('#search_stats').text(`0 of 0 points`);
        disable_buttons();
      }
    }
  });

  $('#end-date').datepicker({
    dateFormat,
    onSelect: function(date_text){
      const end = utcDate(date_text);
      if ($('#end-date').val() != "") {
        $('#search-button').removeAttr("disabled");
        $('#clear_button').removeAttr("disabled");
      }
      else {
        clear_selection();
        $('#search_stats').text(`0 of 0 points`);
        disable_buttons();
      }
    }
  });

  $('#search-button').click(async function() {
    $('#error_message').text(``);
    $("#loading").show();
    $("#loading_text").text(`Loading. Please Wait...`);
    await wrap_search();
    enable_buttons();
    $("#loading").hide();
    $("#loading_text").text(``);
  });

  // Clear search
  $('#clear_button').click(function(){
    clear_selection();
    // search();
    $('#search-input').val('');
    $('#end-date').val('');
    $('#start-date').val('');
    $('#error_message').text(``);
    disable_buttons();
  });

  $('#next_button').click(async function() {
    var curr_bead = search_results.get().current_point;
    var bead_hits = search_results.get().beads;
    var bead_id_to_accession = Object.keys(bead_hits);
    var hit_ids = search_results.get().hit_ids;

    // Edge case: User clicks next from a cluster above the first cluster
    if (curr_bead == 0 && (parseInt(d3.selectAll("rect.clicked").nodes()[0].id.substring(3)) > hit_ids[hit_ids.length - 1])) {
      $("#loading").show();
      $("#loading_text").text(`Loading. Please Wait...`);
      await select_next_prev_bead(bead_id_to_accession, curr_bead);
      select_working_bead(bead_id_to_accession, curr_bead);
      $("#loading").hide();
      $("#loading_text").text(``);
    }
    else if (curr_bead + 1 < search_results.get().total_points) {
      var curr_cid, next_cid;
      await fetch(`/epicov_api/cid/${bead_id_to_accession[curr_bead]}`)
      .then(response => response.text())
      .then(data => curr_cid = data);

      await fetch(`/epicov_api/cid/${bead_id_to_accession[curr_bead + 1]}`)
      .then(response => response.text())
      .then(data => next_cid = data);

      if (curr_cid !== next_cid) {
        $("#loading").show();
        $("#loading_text").text("Loading...");
        await select_next_prev_bead(bead_id_to_accession, curr_bead+1);
        $("#loading").hide();
        $("#loading_text").text(``);
        select_working_bead(bead_id_to_accession, curr_bead + 1);
      }
      else
        select_working_bead(bead_id_to_accession, curr_bead + 1);

      const stats = search_results.update({
        current_point: curr_bead + 1
      });
      
      update_search_stats(stats);
    }

  });


  $('#previous_button').click(async function() {
    var curr_bead = search_results.get().current_point;
    var bead_hits = search_results.get().beads;
    var bead_id_to_accession = Object.keys(bead_hits);
    var hit_ids = search_results.get().hit_ids;

    var current_selection = d3.selectAll("rect.clicked").nodes()[0];
    if (current_selection.className.baseVal !== "SelectedCluster clicked") {
      if(parseInt(current_selection.id.substring(3)) < hit_ids[hit_ids.length - 1]) {
        $("#loading").show();
        $("#loading_text").text("Loading...");
        await select_next_prev_bead(bead_id_to_accession, curr_bead);
        select_working_bead(bead_id_to_accession, curr_bead);
        $("#loading").hide();
        $("#loading_text").text(``);
      }
    }
    else if (curr_bead - 1 >= 0) {
      var curr_cid, prev_cid;
      await fetch(`/epicov_api/cid/${bead_id_to_accession[curr_bead]}`)
      .then(response => response.text())
      .then(data => curr_cid = data);

      await fetch(`/epicov_api/cid/${bead_id_to_accession[curr_bead - 1]}`)
      .then(response => response.text())
      .then(data => prev_cid = data);

      // If the previous bead is not in the same cluster, selection of cluster needs to be modified
      if (curr_cid !== prev_cid) {
        $("#loading").show();
        $("#loading_text").text("Loading...");
        await select_next_prev_bead(bead_id_to_accession, curr_bead-1);
        $("#loading").hide();
        $("#loading_text").text(``);
        select_working_bead(bead_id_to_accession, curr_bead-1);
      }
      else
        select_working_bead(bead_id_to_accession, curr_bead - 1);

      const stats = search_results.update({
        current_point: curr_bead - 1
      });
      
      update_search_stats(stats);
    }
  });


  $(document).on('keydown', function(e) {
    // Ignore event if its inside an input field
    if (e.target.matches('input')) {
      return;
    }

    // User presses the left arrow key (37) or right arrow key (39)
    if (e.keyCode == 37 || e.keyCode == 39) {
      var selected_bead = d3.selectAll(".selectionH").nodes();

      if (selected_bead.length == 0) {
        var points_ui = d3.selectAll("#svg-cluster > svg > g > circle").nodes()[0];
        var working_bead = d3.selectAll('circle[id="'+points_ui.__data__.accessions[0]+'"]').nodes()[0];
        working_bead.scrollIntoView({block: "center"});
        update_table_individual_bead(d3.select(working_bead).datum());
      }
      else {
        var selected_accession = selected_bead[0].attributes.bead.nodeValue;
        var bead_node = d3.selectAll('circle[id="'+selected_accession+'"]');
        var bead_id = parseInt(bead_node.nodes()[0].attributes.idx.nodeValue);
        var total_nodes = d3.selectAll("#svg-cluster > svg > g > circle").nodes().length;

        if (e.keyCode == 37) {
          if (bead_id - 1 >= 0) {
            var prev_node = d3.selectAll('circle[idx="'+(bead_id-1).toString()+'"]').nodes()[0];
            select_next_bead(prev_node);
          }
        }
        else if (e.keyCode == 39) {
          if (bead_id + 1 < total_nodes - 1) {
            var next_node = d3.selectAll('circle[idx="'+(bead_id+1).toString()+'"]').nodes()[0];
            select_next_bead(next_node);
          }
        }
      }
    }
  });
});

/*********************** UPDATE TABLES ***********************/
// populate countries table
var country_table = d3.select("#country-table").append('table');
var theaders = ["Region", "Country", "Count"];

// to be populated in beadplot.js
var country_tbody = country_table.append("tbody");

// Populate sequence details table
var seq_table = d3.select("#seq-table").append('table')
    .attr("class", "details");
var thead = seq_table.append('thead');
var seq_theaders = ["Accession", "Name", "Date", "Locale", "Age", "Gender", "Status"];
var seq_tbody = seq_table.append('tbody');

// implement acknowledgements dialog
$( "#dialog" ).dialog({ autoOpen: false });


/*********************** FILE SAVING ***********************/

// implement save buttons
var blob;
function save_timetree() {
  blob = new Blob([nwk], {type: "text/plain;charset=utf-8"});
  saveAs(blob, "timetree.nwk");
}

function save_beadplot() {
  blob = new Blob([serialize_beadplot(cindex)],
      {type: "text/plain;charset=utf-8"});
  saveAs(blob, lineage + ".nwk");
}

/*  AP: pending testing on dev branch
$('#save_svg').click(function(){
  var config = {
    filename: 'customFileName',
  }
  d3_save_svg.save(d3.select('#svg-cluster>svg').node(), config);

});
*/
function export_svg() {
  var config = {filename: lineage};

  // Creates a duplicate of the beadplot
  var svg_beadplot = d3.select('#svg-cluster>svg').clone(true);
  var svg_axis = d3.select('#svg-clusteraxis>svg').clone(true);
  svg_beadplot.select("g").attr("transform", "translate(0, 30)");
  svg_axis.node().appendChild(svg_beadplot.selectAll("g").node())
  svg_axis.attr("height", svg_beadplot.attr("height"))
  
  d3_save_svg.save(svg_axis.node(), config);

  svg_axis.remove();
  svg_beadplot.remove()
}
