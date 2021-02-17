$( function() {
  $( "#tabs" ).tabs();
} );

$(document).tooltip({show: null});


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
                "<img style=\"vertical-align: middle\" width=\"33px\" src=\"img/Loading_icon_cropped.gif\"/>" +
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

// load database statistics
var dbstats, req;
req = $.getJSON("data/dbstats.json", function(data) {
  dbstats = data;
});
req.done(function() {
  $("#div-last-update").text(`Last update: ${dbstats.lastupdate}`);
  $("#div-number-genomes").text(`Number of genomes: ${dbstats.noseqs}`);
});

var country_pal = {
  "Africa": "#EEDD88",
  "Asia": "#BBCC33",
  "China": "#EE8866",
  "Europe": "#44BB99",
  "North America": "#99DDFF",
  "Oceania": "#FFAABB",
  "South America": "#77AADD"
};

// load time-scaled phylogeny from server
var nwk, df, countries;
$.ajax({
  url: "data/timetree.nwk",
  success: function(data) {
    nwk = data;
    df = readTree(data);
  }
});
$.getJSON("data/countries.json", function(data) {
  countries = data;
});


var clusters, beaddata, tips,
    accn_to_cid, cindex, lineage_to_cid;
var map_cidx_to_id = [], id_to_cidx = [];

// load cluster data from server
req = $.getJSON("data/clusters.json", function(data) {
  clusters = data;
});

req.done(function() {
  $("#splash-button").button("enable");
  $("#splash-extra").html("");  // remove loading animation

  beaddata = parse_clusters(clusters);
  tips = map_clusters_to_tips(df, clusters);
  //spinner.stop();
  draw_clusters(tips);

  var rect = d3.selectAll("#svg-timetree > svg > rect"),
      node = rect.nodes()[rect.size()-1];

  // initial display
  d3.select(node).dispatch("click");

  /*
  rect = d3.selectAll("#svg-cluster > svg > g > circle");
  node = rect.nodes()[0];
  d3.select(node).dispatch("click");//.dispatch("mouseover");
   */

  accn_to_cid = index_accessions(clusters);

  $('#search-input').autocomplete({
    source: get_autocomplete_source_fn(accn_to_cid),
    select: function( event, ui ) {
        const accn = ui.item.value;
        //search(accn);
    }
  });

  // Maps lineage to a cidx
  lineage_to_cid = index_lineage(clusters);

  // Maps cidx to an id
  var key;
  var rect = d3.selectAll('#svg-timetree > svg > rect:not(.clickedH)').nodes();
  for (var i = 0; i < rect.length; i++) {
	  key = d3.select(rect[i]).attr("cidx");
	  map_cidx_to_id[key] = parseInt(d3.select(rect[i]).attr("id").substring(3));
  }

  // Maps id to a cidx
  const reverseMapping = o => Object.keys(o).reduce((r, k) => Object.assign(r, { [o[k]]: (r[o[k]] || []).concat(k) }), {})
  id_to_cidx = reverseMapping(map_cidx_to_id);


  /***********************  SEARCH INTERFACE ***********************/
  // function run_search() {

  //     var query = $('#search-input').val();

  //     if (query !== "") {
  //       // revert selections
  //       d3.selectAll("rect.clicked").attr('class', "default");
  //       d3.selectAll("rect.clickedH").remove();
  //     }

  //     // Create new search stats
  //     const points = find_beads_points(beaddata)
  //             .filter(point => point.labels.some(label => label.includes(query)));

  //     // Map cluster index to id
  //     var map_to_id = [], key;
  //     var rect = d3.selectAll('#svg-timetree > svg > rect')
  //             .nodes()
  //             .sort((x, y) => d3.ascending(
  //                     parseInt(x.id.substring(3)),
  //                     parseInt(y.id.substring(3))
  //             ));
  //     for (var i = rect.length - 1; i >= 0; i--) {
  //       key = rect[i].id;
  //       map_to_id[key] = d3.select(rect[i]).attr("cidx");
  //     }

  //     // Count the number of hits in each cluster
  //     var count_hits_per_cluster = {};
  //     for (var i = 0; i < points.length; i++) {
  //       key = 'cidx-' + points[i].cidx;
  //       if (count_hits_per_cluster[key] == null) {
  //         count_hits_per_cluster[key] = 1
  //         continue
  //       }
  //       count_hits_per_cluster[key]++;
  //     }

  //     // First index of each cluster
  //     var start_idx = [], start_index = 0;
  //     for (const [key, value] of Object.entries(map_to_id)) {
  //       start_idx[key] = start_index;
  //       if (count_hits_per_cluster[value] != null) {
  //         start_index = start_index +  count_hits_per_cluster[value];
  //       }
  //     }

  //     const stats = search_stats.update({
  //       query,
  //       current_point: 0,
  //       total_points: points.length,
  //       points: points,
  //       bead_indexer: 0,
  //       start_idx: start_idx,
  //     });
  //     update_search_stats(stats);
  //     search();
  //     enable_buttons();
  // }

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
      disable_buttons();
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
  })

  $('#search-input').on('keydown', function(e) {
    $('#error_message').text(``);
    if (e.keyCode == 13 && ($('#search-input').val() != "" || $('#start-date').val() != "" || $('#end-date').val() != "")) {
      // type <enter> to run search
      // run_search();
      wrap_search();
      enable_buttons();
    }
  });

  // Removes the error message when the user clicks on the date picker
  $('#start-date').on('mousedown', function() {
    $('#error_message').text(``);
  });

  $('#end-date').on('mousedown', function() {
    $('#error_message').text(``);
  });

  $('#start-date').on('change keyup search', function() {
    if ($('#start-date').val() != "") {
      $('#search-button').removeAttr("disabled");
      $('#clear_button').removeAttr("disabled");
    }
    else {
      clear_selection();
      $('#search_stats').text(`0 of 0 points`);
      disable_buttons();
    }
  });

  $('#end-date').on('change keyup search', function() {
    if ($('#end-date').val() != "") {
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
      const start = new Date(date_text);
      if ($('#start-date').val() != "") {
        $('#search-button').removeAttr("disabled");
        $('#clear_button').removeAttr("disabled");
      }
      else {
        clear_selection();
        $('#search_stats').text(`0 of 0 points`);
        disable_buttons();
      }
      // const stats = search_stats.update({
      //   start,
      // });
      // update_search_stats(stats);
      // if (start && stats.end) {
      //   search_by_dates(beaddata, start, stats.end);
      //   run_search();
      // }
    }
  });

  $('#end-date').datepicker({
    dateFormat,
    onSelect: function(date_text){
      const end = new Date(date_text);
      if ($('#end-date').val() != "") {
        $('#search-button').removeAttr("disabled");
        $('#clear_button').removeAttr("disabled");
      }
      else {
        clear_selection();
        $('#search_stats').text(`0 of 0 points`);
        disable_buttons();
      }
      // search_by_dates(search_stats.get().start, end);
      // const stats = search_stats.update({
      //   end,
      // });
      // if (stats.start && end) {
      //   search_by_dates(beaddata, stats.start, end);
      //   run_search();
      // }
    }
  });

  $('#search-button').click(function() {
    wrap_search();
    enable_buttons();
  });

  // Clear search
  $('#clear_button').click(function(){
    clear_selection();
    // search();
    $('#end-date').val('');
    $('#start-date').val('');
    disable_buttons();
  });

  $('#next_button').click(function() {
    var curr_bead = search_results.get().current_point;
    var bead_hits = search_results.get().beads;
    var bead_id_to_accession = Object.keys(bead_hits);

    // Edge case: User clicks next from a cluster above the first cluster
    if (curr_bead == 0 && d3.selectAll("rect.clicked").nodes()[0].className.baseVal !== "SelectedCluster clicked") {
      select_next_prev_bead(bead_id_to_accession, curr_bead);
      var working_bead = d3.selectAll('circle[id="'+bead_id_to_accession[curr_bead]+'"]').nodes()[0];
      working_bead.scrollIntoView({block: "center"});
      update_table_individual_bead(d3.select(working_bead).datum());
    }
    else if (curr_bead + 1 < search_results.get().total_points) {
      if (accn_to_cid[bead_id_to_accession[curr_bead]] != accn_to_cid[bead_id_to_accession[curr_bead + 1]]) {
        select_next_prev_bead(bead_id_to_accession, curr_bead+1);
      }
      var working_bead = d3.selectAll('circle[id="'+bead_id_to_accession[curr_bead + 1]+'"]').nodes()[0];
      working_bead.scrollIntoView({block: "center"});
      update_table_individual_bead(d3.select(working_bead).datum());

      const stats = search_results.update({
        current_point: curr_bead + 1
      });
      
      update_search_stats(stats);
    }
  });


  $('#previous_button').click(function(){
    var curr_bead = search_results.get().current_point;
    var bead_hits = search_results.get().beads;
    var bead_id_to_accession = Object.keys(bead_hits);
    var hit_ids = search_results.get().hit_ids;

    var current_selection = d3.selectAll("rect.clicked").nodes()[0];
    if (current_selection.className.baseVal !== "SelectedCluster clicked") {
      if(parseInt(current_selection.id.substring(3)) < hit_ids[hit_ids.length - 1]) {
        select_next_prev_bead(bead_id_to_accession, curr_bead);
        var working_bead = d3.selectAll('circle[id="'+bead_id_to_accession[curr_bead]+'"]').nodes()[0];
        working_bead.scrollIntoView({block: "center"});
        update_table_individual_bead(d3.select(working_bead).datum());
      }
    }
    else if (curr_bead - 1 >= 0) {
      // If the previous bead is not in the same cluster, selection of cluster needs to be modified
      if (accn_to_cid[bead_id_to_accession[curr_bead]] != accn_to_cid[bead_id_to_accession[curr_bead - 1]]) {
        select_next_prev_bead(bead_id_to_accession, curr_bead-1);
      }
      var working_bead = d3.selectAll('circle[id="'+bead_id_to_accession[curr_bead - 1]+'"]').nodes()[0];
      working_bead.scrollIntoView({block: "center"});
      update_table_individual_bead(d3.select(working_bead).datum());

      const stats = search_results.update({
        current_point: curr_bead - 1
      });
      
      update_search_stats(stats);
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
var seq_table = d3.select("#seq-table").append('table');
var thead = seq_table.append('thead');
var seq_theaders = ["Accession", "Name", "Date"];
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
  saveAs(blob, clusters[cindex].lineage + ".nwk");
}