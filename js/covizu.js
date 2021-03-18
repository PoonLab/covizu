$( function() {
  $( "#tabs" ).tabs();

  var handle = $( "#custom-handle" );

  $("#vedge-slider").slider({
    create: function( event, ui ) {
      handle.text( $( this ).slider( "value" ) );
    },
    slide: function( event, ui ) {
      handle.text( ui.value );
    },
    min: 0.5,
    step: 0.01
  });
} );

$(document).tooltip({show: null});
$("#loading_text").text(``);
$("#loading").hide();


/*********************** DIALOGS ***********************/

$( "#splash" ).dialog({
  autoOpen: true,
  width: 600,
  buttons: [{
    id: "splash-button",
    text: i18n_text.okay,
    disabled: true,
    click: function() {
      $( this ).dialog( "close" );
    }
  }]
});

$("#splash-button").parents(".ui-dialog-buttonpane")
        .append("<div style=\"color: #4F2683\" id='splash-extra'>" + i18n_text.loading_json +
                "<img style=\"vertical-align: middle\" width=\"33px\" src=\"img/Loading_icon_cropped.gif\"/>" +
                "</div>");

$( "#help-timetree" ).dialog({
  autoOpen: false,
  width: 600,
  buttons: [{
    text: i18n_text.got_it,
    click: function() {
      $( this ).dialog( "close" );
    }
  }]
});

$( "#help-beadplot" ).dialog({
  autoOpen: false,
  width: 600,
  buttons: [{
    text: i18n_text.got_it,
    click: function() {
      $( this ).dialog( "close" );
    }
  }]
});

$( "#help-search" ).dialog({
  autoOpen: false,
  width: 600,
  buttons: [{
    text: i18n_text.got_it,
    click: function() {
      $( this ).dialog( "close" );
    }
  }]
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
  dbstats.nlineages = Object.keys(dbstats.lineages).length;
});
req.done(function() {
  $("#div-last-update").text(`${i18n_text.last_update}: ${dbstats.lastupdate}`);
  $("#div-number-genomes").text(`${i18n_text.number_genomes}: ${dbstats.noseqs}`);
  $("#div-number-lineages").text(`${i18n_text.number_lineages}: ${dbstats.nlineages}`);
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
  // d3.select(node).dispatch("click");
  cindex = node.__data__.cluster_idx;
  d3.select(node).attr("class", "clicked");
  beadplot(node.__data__.cluster_idx);
  $("#barplot").text(null);
  gentable(node.__data__);
  draw_region_distribution(node.__data__.allregions);
  gen_details_table(beaddata[node.__data__.cluster_idx].points);  // update details table with all samples
  draw_cluster_box(d3.select(node));

  /*
  rect = d3.selectAll("#svg-cluster > svg > g > circle");
  node = rect.nodes()[0];
  d3.select(node).dispatch("click");//.dispatch("mouseover");
   */

  accn_to_cid = index_accessions(clusters);

  // Maps lineage to a cidx
  lineage_to_cid = index_lineage(clusters);

  $('#search-input').autocomplete({
    source: get_autocomplete_source_fn(accn_to_cid, lineage_to_cid),
    select: function( event, ui ) {
        const accn = ui.item.value;
        //search(accn);
    }
  });


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
    if ($('#search-input').val() !== "" || $('#start-date').val() !== "" ||
        $('#end-date').val() !== "") {
      $('#search-button').removeAttr("disabled");
      $('#clear_button').removeAttr("disabled");
    }
    else {
      clear_selection();
      $('#search_stats').text(i18n_text.zero_points);
      disable_buttons();
    }
  });

  $('#search-input, #start-date, #end-date').on('keydown', function(e) {
    $('#error_message').text(``);
    if (search_results.get().total_points > 0) {
      clear_selection();
      $('#search_stats').text(i18n_text.zero_points);
      disable_buttons();
    }
    if (e.keyCode === 13 && ($('#search-input').val() !== "" || $('#start-date').val() !== "" ||
        $('#end-date').val() !== "")) {
      // type <enter> to run search
      // run_search();
      $('#error_message').text(``);
      $("#loading").show();
      $("#loading_text").text(i18n_text.loading);
      setTimeout(function() {
        wrap_search();
        enable_buttons();
        $("#loading").hide();
        $("#loading_text").text(``);
      }, 20);
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
      $('#search_stats').text(i18n_text.zero_points);
      disable_buttons();
    }
  });

  $('#end-date').on('change keyup', function() {
    if ($('#end-date').val() !== "" || $('#search-input').val() !== "" ||
        $('#start-date').val() !== "") {
      $('#search-button').removeAttr("disabled");
      $('#clear_button').removeAttr("disabled");
    }
    else {
      clear_selection();
      $('#search_stats').text(i18n_text.zero_points);
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
        $('#search_stats').text(i18n_text.zero_points);
        disable_buttons();
      }
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
        $('#search_stats').text(i18n_text.zero_points);
        disable_buttons();
      }
    }
  });

  $('#search-button').click(function() {
    $('#error_message').text(``);
    $("#loading").show();
    $("#loading_text").text(i18n_text.loading);
    setTimeout(function() {
      wrap_search();
      enable_buttons();
      $("#loading").hide();
      $("#loading_text").text(``);
    }, 20);
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

  $('#next_button').click(function() {
    var curr_bead = search_results.get().current_point;
    var bead_hits = search_results.get().beads;
    var bead_id_to_accession = Object.keys(bead_hits);
    var hit_ids = search_results.get().hit_ids;

    // console.log(first_bead.id);
    // Edge case: User clicks next from a cluster above the first cluster
    if (curr_bead == 0 && (parseInt(d3.selectAll("rect.clicked").nodes()[0].id.substring(3)) > hit_ids[hit_ids.length - 1])) {
      $("#loading").show();
      $("#loading_text").text(i18n_text.loading);
      setTimeout(function() {
        select_next_prev_bead(bead_id_to_accession, curr_bead);
        select_working_bead(bead_id_to_accession, curr_bead);
        $("#loading").hide();
        $("#loading_text").text(``);
      }, 20);
    }
    else if (curr_bead + 1 < search_results.get().total_points) {
      if (accn_to_cid[bead_id_to_accession[curr_bead]] !==
          accn_to_cid[bead_id_to_accession[curr_bead + 1]]) {
        $("#loading").show();
        $("#loading_text").text(i18n_text.loading);
        setTimeout(function() {
          select_next_prev_bead(bead_id_to_accession, curr_bead + 1);
          select_working_bead(bead_id_to_accession, curr_bead + 1);
          $("#loading").hide();
          $("#loading_text").text(``);
        }, 20);
      }
      else
        select_working_bead(bead_id_to_accession, curr_bead + 1);

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
        $("#loading").show();
        $("#loading_text").text(i18n_text.loading);
        setTimeout(function() {
          select_next_prev_bead(bead_id_to_accession, curr_bead);
          select_working_bead(bead_id_to_accession, curr_bead);
          $("#loading").hide();
          $("#loading_text").text(``);
        }, 20);
      }
    }
    else if (curr_bead - 1 >= 0) {
      // If the previous bead is not in the same cluster, selection of cluster needs to be modified
      if (accn_to_cid[bead_id_to_accession[curr_bead]] !==
          accn_to_cid[bead_id_to_accession[curr_bead - 1]]) {
        $("#loading").show();
        $("#loading_text").text(i18n_text.loading);
        setTimeout(function() {
          select_next_prev_bead(bead_id_to_accession, curr_bead-1);
          select_working_bead(bead_id_to_accession, curr_bead-1);
          $("#loading").hide();
          $("#loading_text").text(``);
        }, 20);
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
var theaders = i18n_text.country_theaders;

// to be populated in beadplot.js
var country_tbody = country_table.append("tbody");

// Populate sequence details table
var seq_table = d3.select("#seq-table").append('table');
var thead = seq_table.append('thead');
var seq_theaders = i18n_text.sample_theaders;
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

function export_svg() {
  var config = {filename: clusters[cindex].lineage};
  var svg_beadplot = d3.select('#svg-cluster>svg').node();
  d3_save_svg.save(svg_beadplot, config);
  svg_beadplot.removeAttribute("style");
  svg_beadplot.removeAttribute("version");
  svg_beadplot.removeAttribute("xmlns");
  svg_beadplot.removeAttribute("xmlns:xlink");
}
