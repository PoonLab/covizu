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
    accn_to_cid, cindex;

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

  /***********************  SEARCH INTERFACE ***********************/
  function run_search() {

      var query = $('#search-input').val();

      if (query !== "") {
        // revert selections
        d3.selectAll("rect.clicked").attr('class', "default");
        d3.selectAll("rect.clickedH").remove();
      }

      // Create new search stats
      const points = find_beads_points(beaddata)
              .filter(point => point.labels.some(label => label.includes(query)));

      // Map cluster index to id
      var map_to_id = [], key;
      var rect = d3.selectAll('#svg-timetree > svg > rect')
              .nodes()
              .sort((x, y) => d3.ascending(
                      parseInt(x.id.substring(3)),
                      parseInt(y.id.substring(3))
              ));
      for (var i = rect.length - 1; i >= 0; i--) {
        key = rect[i].id;
        map_to_id[key] = d3.select(rect[i]).attr("cidx");
      }

      // Count the number of hits in each cluster
      var count_hits_per_cluster = {};
      for (var i = 0; i < points.length; i++) {
        key = 'cidx-' + points[i].cidx;
        if (count_hits_per_cluster[key] == null) {
          count_hits_per_cluster[key] = 1
          continue
        }
        count_hits_per_cluster[key]++;
      }

      // First index of each cluster
      var start_idx = [], start_index = 0;
      for (const [key, value] of Object.entries(map_to_id)) {
        start_idx[key] = start_index;
        if (count_hits_per_cluster[value] != null) {
          start_index = start_index +  count_hits_per_cluster[value];
        }
      }

      const stats = search_stats.update({
        query,
        current_point: 0,
        total_points: points.length,
        points: points,
        bead_indexer: 0,
        start_idx: start_idx,
      });
      update_search_stats(stats);
      search();
      enable_buttons();
  }

  // Enable and Disable "Search", "Clear" "Next" and "Previous" buttons when needed
  function disable_buttons() {
    $('#search-button').attr("disabled", true);
    $('#clear_button').attr("disabled", true);
    $('#next_button').attr("disabled", true);
    $('#previous_button').attr("disabled", true);
    $('#search_stats').addClass("disabled_stats");
  }

  function enable_buttons() {
    if (search_stats.get().total_points > 0) {
      $('#next_button').removeAttr("disabled");
      $('#previous_button').removeAttr("disabled");
      $('#search_stats').removeClass("disabled_stats");
    } else {
      disable_buttons();
    }
  }

  disable_buttons();

  // Enables "search" and "clear" buttons if the input fields are not empty
  $('#search-input').on('change keyup', function() {
    if ($('#search-input').val() != "") {
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
    if (e.keyCode == 13 && $('#search-input').val() != "") {
      // type <enter> to run search
      run_search();
      wrap_search();
    }
  });

  const dateFormat = 'yy-mm-dd'; // ISO_8601
  $('#start-date').datepicker({
    dateFormat,
    onSelect: function(date_text){
      const start = new Date(date_text);
      const stats = search_stats.update({
        start,
      });
      update_search_stats(stats);
      if (start && stats.end) {
        search_by_dates(beaddata, start, stats.end);
        run_search();
      }
    }
  });

  $('#end-date').datepicker({
    dateFormat,
    onSelect: function(date_text){
      const end = new Date(date_text);
      search_by_dates(search_stats.get().start, end);
      const stats = search_stats.update({
        end,
      });
      if (stats.start && end) {
        search_by_dates(beaddata, stats.start, end);
        run_search();
      }
    }
  });

  $('#search-button').click(function() {
    run_search();
  });

  // Clear search
  $('#clear_button').click(function(){
    clear_selection();
    search();
    $('#end-date').val('');
    $('#start-date').val('');
    disable_buttons();
  });

  // Iterate results
  $('#next_button').click(function() {
    // retrieve current rect element
    var current_cluster = d3.selectAll(".clicked").node();
    var current_index = current_cluster.id;

    if (search_stats.get().current_point+1 < search_stats.get().total_points) {
      // increment bead index and current point display
      var stats = search_stats.update({
        current_point: search_stats.get().start_idx[current_index] +
                search_stats.get().bead_indexer + 1,
        bead_indexer: search_stats.get().bead_indexer + 1,
      });

      update_search_stats(stats);

      // Select bead hits in the current cluster
      var selected_points = d3.selectAll(".SelectedBead");

      // Move to the next cluster if next bead is in next cluster
      if (search_stats.get().bead_indexer > selected_points.nodes().length) {

        // Select cluster hits
        var selected_clusters = d3.selectAll(".SelectedCluster, .clicked").nodes();
        selected_clusters.sort((x, y) => d3.ascending(parseInt(x.id.substring(3)),
                parseInt(y.id.substring(3))));

        // Find current cluster
        var matched = selected_clusters.findIndex(function(d, i) {
          if (d.id === current_index)
            return i;
        })

        // Move to next cluster and reset bead indexer
        var next_cluster = selected_clusters[matched - 1];
        d3.select(next_cluster).dispatch('click');
        stats = search_stats.update({
                current_point: search_stats.get().start_idx[next_cluster.id] + search_stats.get().bead_indexer + 1,
                bead_indexer: search_stats.get().bead_indexer + 1,
        });
        update_search_stats(stats);

        // Select bead hits in next cluster
        selected_points = d3.selectAll(".SelectedBead");
      }

      // Scroll to the next bead
      var working_bead = selected_points.nodes()[search_stats.get().bead_indexer-1];
      working_bead.scrollIntoView({block: "center"});

      var selected_bead = d3.select(working_bead).datum();

      draw_halo(selected_bead);
      gentable(selected_bead);
      draw_region_distribution(tabulate(selected_bead.region));
      gen_details_table(selected_bead);
    }


  });  // end click next button

  $('#previous_button').click(function(){

    // Find current cluster
    var current_cluster = d3.selectAll(".clicked").node();
    var current_index = current_cluster.id;

    if (search_stats.get().current_point > 0) {
      // Find index of current bead
      var stats = search_stats.update({
        current_point: search_stats.get().start_idx[current_index] +
                search_stats.get().bead_indexer - 1,
        bead_indexer: search_stats.get().bead_indexer - 1,
      });
      update_search_stats(stats);

      // Select bead hits in the current cluster
      var selected_points = d3.selectAll(".SelectedBead");

      // Move to the previous cluster if previous bead is in previous cluster
      if (search_stats.get().bead_indexer < 1) {

        // Select cluster hits
        var selected_clusters = d3.selectAll(".SelectedCluster, .clicked").nodes();
        selected_clusters.sort((x, y) => d3.ascending(parseInt(x.id.substring(3)),
                parseInt(y.id.substring(3)))
        );

        // Find current cluster
        var matched = selected_clusters.findIndex(function(d, i) {
          if (d.id === current_index)
            return i;
        })

        // Move to next cluster and reset bead indexer
        var previous_cluster = selected_clusters[matched + 1];
        d3.select(previous_cluster).dispatch('click');

        // Select bead hits in next cluster
        var selected_points = d3.selectAll(".SelectedBead");

        stats = search_stats.update({
          current_point: search_stats.get().start_idx[previous_cluster.id] + selected_points.nodes().length,
          bead_indexer: selected_points.nodes().length,
        });
        update_search_stats(stats);

      }

      // Scroll to the next bead
      var working_bead = selected_points.nodes()[search_stats.get().bead_indexer-1];
      working_bead.scrollIntoView({block: "center"});

      var selected_bead = d3.select(working_bead).datum();
      draw_halo(selected_bead);
      gentable(selected_bead);
      draw_region_distribution(tabulate(selected_bead.region));
      gen_details_table(selected_bead);
    }

  });  // end click previous button

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
  saveAs(blob, clusters[cindex].lineage + ".nwk");
}