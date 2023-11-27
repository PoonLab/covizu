$( function() {
  $( "#tabs" ).tabs();

  var handle = $( "#custom-handle" );

  $("#vedge-slider").slider({
    create: function( event, ui ) {
      handle.text( $( this ).slider( "value" ) );
    },
    slide: function( event, ui ) {
      move_arrow();
      handle.text( ui.value );
    },
    stop: function (event, ui) {
      move_arrow();
    },
    min: 0.5,
    step: 0.01
  });

  // Prevents the default action when keydown event is detected
  handle.unbind('keydown')

} );

$(document).tooltip({show: null});
$("#loading_text").text(``);
$("#loading").hide();
$('#beadplot-hscroll').hide();
$("#loading-tree").hide();


/*********************** DIALOGS ***********************/

$( "#splash" ).dialog({
  autoOpen: true,
  width: 600,
  buttons: [{
    id: "tour-button",
    text: i18n_text.tour,
    disabled: true,
    click: function() {
      $( this ).dialog( "close" )

      // Walkthrough
      const driver = new Driver({
        allowClose: false,
        padding: 2.5,
        doneBtnText: i18n_text.done_bttn,
        closeBtnText: i18n_text.close_bttn,
        nextBtnText: i18n_text.next_bttn,
        prevBtnText: i18n_text.previous_bttn,
      });

      driver.defineSteps([
        {
          element: '#tree-axis',
          popover: {
            className: 'first-step-popover-class',
            title: i18n_text.tour_tree_title,
            description: i18n_text.tour_tree_desc,
            position: 'right'
          }
        },
        {
          element: '.clicked',
          popover: {
            className: 'first-step-popover-class',
            title: i18n_text.tour_lin_title,
            description: i18n_text.tour_lin_desc,
          }
        },
        {
          element: '.legend-container',
          popover: {
            title: i18n_text.tour_legend_title,
            description: i18n_text.tour_legend_desc,
          }
        },
        {
          element: '#beadplot-container',
          popover: {
            title: i18n_text.tour_beadplot_title,
            description: i18n_text.tour_beadplot_desc,
            position: 'right'
          }
        },
        {
          element: document.querySelector('#svg-cluster > svg > g > text'),
          popover: {
            title: i18n_text.tour_variant_title,
            description: i18n_text.tour_variant_desc,
            position: 'bottom-center'
          }
        },
        {
          element: document.querySelector('[idx="0"]'),
          popover: {
            title: i18n_text.tour_bead_title,
            description: i18n_text.tour_bead_desc,
            position: 'left-center'
          }
        },
        {
          element: '.search-bar-container',
          popover: {
            title: i18n_text.tour_search_title,
            description: i18n_text.tour_search_desc,
          }
        },
        {
          element: '#tabs',
          popover: {
            title: i18n_text.tour_tables_title,
            description: i18n_text.tour_tables_desc,
            position: 'left-center'
          }
        },
        {
          element: '#intro',
          popover: {
            title: i18n_text.tour_end_title,
            description: i18n_text.tour_end_desc,
            position: 'bottom'
          }
        },
      ]);
      driver.start()
    }
  },
  {
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
  console.log("dbstats = ", dbstats)
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

var phenotypes = {
  'Vaccine neutralization efficacy': 'img/red_circle.png', 
  'Anthropozoonotic events': 'img/bat.png', 
  'Gene expression increase': 'img/orange_star.png', 
  'ACE2 receptor binding affinity': 'img/purple_square.jpeg',
  'Monoclonal antibody serial passage escape': 'img/antibody.png', 
  'Convalescent plasma escape': 'img/green_pentagon.png', 
  'Antibody epitope effects': 'img/blue_triangle.png'
}

// load time-scaled phylogeny from server
var nwk, df, df_xbb, countries, mut_annotations, region_map;

// $.getJSON("data/mut_annotations.json", function(data) {
//   mut_annotations = data;
//   console.log("MUTATION ANNOTATION ",data)
// });

var clusters, beaddata, tips, recombinant_tips,
    accn_to_cid, cindex, lineage_to_cid, lineage;
var edgelist = [], points = [], variants = []
var map_cidx_to_id = [], id_to_cidx = [];

req = $.when(
  
  $.getJSON("data/mut_annotations.json", function(data) {
    mut_annotations = data;
    console.log("MUTATION ANNOTATION ",data)
  }),

  $.getJSON("/api/tips", function(data) {
    tips = data;
    tips.forEach(x => {
      x.first_date = new Date(x.first_date)
      x.last_date = new Date(x.last_date)
      x.coldate = new Date(x.coldate)
      x.mcoldate = new Date(x.mcoldate)
    });
  }),
  $.getJSON("/api/recombtips", function(data) {
    recombinant_tips = data;
    recombinant_tips.forEach(x => {
      x.first_date = new Date(x.first_date)
      x.last_date = new Date(x.last_date)
      x.coldate = new Date(x.coldate)
      x.mcoldate = new Date(x.mcoldate)
    });
  }),
  $.getJSON("/api/df", function(data) {
    df = data;
    df.forEach(x => {
      x.first_date = x.first_date ? new Date(x.first_date) : undefined
      x.last_date = x.last_date ? new Date(x.last_date) : undefined
      x.coldate = x.coldate ? new Date(x.coldate) : undefined
      x.mcoldate = x.coldate ? new Date(x.mcoldate) : undefined
    });
  }),
  $.getJSON("/api/xbb", function(data) {
    df_xbb = data;
    df_xbb.forEach(x => {
      x.first_date = x.first_date ? new Date(x.first_date) : undefined
      x.last_date = x.last_date ? new Date(x.last_date) : undefined
      x.coldate = x.coldate ? new Date(x.coldate) : undefined
      x.mcoldate = x.coldate ? new Date(x.mcoldate) : undefined
    });
  }),
  $.getJSON("/api/regionmap", function(data) {
    region_map = data;
  })
);

req.done(async function() {

  var urlParams = new URLSearchParams(window.location.search);
  var search = urlParams.get('search') || '';

  $("#splash-button").button("enable");
  $("#tour-button").button("enable");
  $("#splash-extra").html("");  // remove loading animation
  
  mutations = parse_mutation_annotations(mut_annotations);
  //spinner.stop();
  var curr_date = new Date();
  curr_date.setFullYear(curr_date.getFullYear() - 1);

  // Maps id to a cidx
  const reverse_recombinant_tips = [...recombinant_tips].reverse()
  var all_tips = [...tips, ...reverse_recombinant_tips]
  console.log("tips = ",tips)
  console.log("recombinant_tips = ",recombinant_tips)
  console.log("all_tips = ",all_tips)
  for (i in all_tips) {
    id_to_cidx[i] = 'cidx-' + all_tips[i].cluster_idx
  }

  // Maps cidx to an id
  const reverseMapping = o => Object.keys(o).reduce((r, k) => Object.assign(r, { [o[k]]: (r[o[k]] || parseInt(k)) }), {})
  map_cidx_to_id = reverseMapping(id_to_cidx)

  switch($("#display-tree").val()) {
    case "XBB Lineages":
      await redraw_tree(df_xbb, formatDate(curr_date), redraw=false);
      break;
    case "Recombinants":
    default:
      await redraw_tree(df, formatDate(curr_date), redraw=false);
  }

  //spinner.stop();
  var rect = d3.selectAll("#svg-timetree > svg > rect"),
      node = rect.nodes()[rect.size()-1];

  // initial display
  // d3.select(node).dispatch("click");
  cindex = node.__data__.cluster_idx;
  console.log("NODE = ",node);
  d3.select(node).attr("class", "clicked");
  window.addEventListener("resize", expand, true);

  /*
  rect = d3.selectAll("#svg-cluster > svg > g > circle");
  node = rect.nodes()[0];
  d3.select(node).dispatch("click");//.dispatch("mouseover");
   */

  // Maps lineage to a cidx
  await fetch(`/api/lineagetocid`)
  .then(response => response.json())
  .then(data => lineage_to_cid = data)
  .then(()=>{console.log("lineage_to_cid",lineage_to_cid)})

  $('#search-input').autocomplete({
    source: function(req, res) {
      $.ajax({
        url: `/api/getHits/${req.term}`,
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

  // Changes the value of the slider and triggers a change event
  function move_slider(direction) {
    var slider = $("#vedge-slider");

    if (direction === "LEFT")
      slider.slider("value", slider.slider("value") - slider.slider("option", "step"))
    else 
      slider.slider("value", slider.slider("value") + slider.slider("option", "step"))
    
    $("#custom-handle").text( slider.slider( "value" ) );
    move_arrow();
    slider_update();
  }

  disable_buttons();

  if (search !== '') {
    $('#search-input').val(search)
    $('#error_message').text(``);
    $('#search-button').removeAttr("disabled");
    $('#clear_button').removeAttr("disabled");
    reset_tree();
    wrap_search();
    enable_buttons();

    // Close the introductory text popup
    $( "#splash" ).dialog("close");

    if ($('#error_message').text() !== '') {
      $('.modal-text').text(`No matches for ${search}. Please try again.`);
      $('.modal').fadeIn(300);
    }
  }
  else {
    d3.select('#cidx-' + cindex).attr("class", "clicked")
    await beadplot(node.__data__.cluster_idx);
    $("#barplot").text(null);
    gentable(node.__data__);
    draw_region_distribution(node.__data__.allregions);
    gen_details_table(points);  // update details table with all samples
    console.log("Generating with ", mutations,cindex)
    gen_mut_table(mutations[cindex]);
    draw_cluster_box(d3.select(node));
  }

  // Calculate the height for the tree container and beadplot container
  set_height()

  // Enables "search" and "clear" buttons if the input fields are not empty
  $('#search-input').on('change keyup search', function() {
    if ($('#search-input').val() !== "" || $('#start-date').val() !== "" ||
        $('#end-date').val() !== "") {
      $('#search-button').removeAttr("disabled");
      $('#clear_button').removeAttr("disabled");
    }
    else {
      clear_selection();
      disable_buttons();
    }
  });

  $('#search-input, #start-date, #end-date').on('keydown', async function(e) {
    $('#error_message').text(``);
    // Only resets search results if the backspace key is pressed
    if (search_results.get().total_points > 0 && (e.keyCode == 8)) {
      clear_selection();
      disable_buttons();
    }
    if (e.keyCode === 13 && ($('#search-input').val() !== "" || $('#start-date').val() !== "" ||
        $('#end-date').val() !== "")) {
      // type <enter> to run search
      // run_search();
      $('#error_message').text(``);
      $("#loading").show();
      $("#loading_text").text(i18n_text.loading);
      await reset_tree(partial_redraw=true);
      await wrap_search();
      enable_buttons();

      $("#loading").hide();
      $("#loading_text").text(``);
      $("#tree-slider").slider({ disabled: true});
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
        disable_buttons();
      }
    }
  });

  $('#search-button').click(async function() {
    $('#error_message').text(``);
    $("#loading").show();
    $("#loading_text").text(i18n_text.loading);
    await reset_tree(partial_redraw=true);
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
    $("#tree-slider").slider({ disabled: false });
  });

  $('#next_button').click(async function() {
    var curr_bead = search_results.get().current_point;
    var bead_hits = search_results.get().beads;
    var bead_id_to_accession = Object.keys(bead_hits);
    var hit_ids = search_results.get().hit_ids;

    // console.log(first_bead.id);
    // Edge case: User clicks next from a cluster above the first cluster
    if (curr_bead == 0 && (parseInt(d3.selectAll("rect.clicked").nodes()[0].id.substring(3)) > hit_ids[hit_ids.length - 1])) {
      $("#loading").show();
      $("#loading_text").text(i18n_text.loading);
      await select_next_prev_bead(bead_id_to_accession, curr_bead);
      select_working_bead(bead_id_to_accession, curr_bead);
      $("#loading").hide();
      $("#loading_text").text(``);
    }
    else if (curr_bead + 1 < search_results.get().total_points) {
      var curr_cid, next_cid;
      await fetch(`/api/cid/${bead_id_to_accession[curr_bead]}`)
      .then(response => response.text())
      .then(data => curr_cid = data);

      await fetch(`/api/cid/${bead_id_to_accession[curr_bead + 1]}`)
      .then(response => response.text())
      .then(data => next_cid = data);

      if (curr_cid !== next_cid) {
        $("#loading").show();
        $("#loading_text").text(i18n_text.loading);
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

  $('#expand-option').on('change', function() {
    if (!$('#expand-option').attr('checked')) {
      $('#expand-option').attr('checked', 'checked');
      $('#beadplot-hscroll').show();
    }
    else {
      $('#expand-option').removeAttr('checked');
      $('#beadplot-hscroll').hide();
    }
    expand();
  });

  $(window).on('resize', set_height);

  // Sets the scrolling speed when scrolling through the beadplot
  const element = document.querySelector("#svg-cluster");

  element.addEventListener('wheel', (event) => {
    event.preventDefault();

    element.scrollBy({
      top: Math.abs(event.deltaY) == 0 ? 0 : event.deltaY * 2,
      left: Math.abs(event.deltaX) == 0 ? 0 : event.deltaX * 2
    });
  });

  // Sets the scrolling speed when scrolling through the time-scaled tree
  const timetree = document.querySelector("#svg-timetree");

  timetree.addEventListener('wheel', (event) => {
    event.preventDefault();

    timetree.scrollBy({
      top: Math.abs(event.deltaY) == 0 ? 0 : event.deltaY * 2,
    });
  });

  // Sets the display date and cutoff line to move along with the slider
  $(function() {

    var handle = $( "#tree-slider-handle" );
    var cutoff_date = $( "#cutoff-date" );
    var cutoff_line = $("#cutoff-line");
    var tree_cutoff = $("#tree-cutoff");

    const tree_multiplier = 100000000; 
    var min = Math.floor((d3.min(df, xValue)-0.05) * tree_multiplier);
    var max = Math.ceil(date_to_xaxis(d3.max(df, function(d) {return d.last_date})) * tree_multiplier);
    var start_value = date_to_xaxis(curr_date) * tree_multiplier;

    $("#tree-slider").slider({
      create: function( event, ui ) {
        cutoff_date.text(xaxis_to_date($( this ).slider( "value" )/tree_multiplier, tips[0], d3.min(tips, function(d) {return d.first_date}), d3.max(tips, function(d) {return d.last_date})));
        var cutoff_pos = handle.position().left;
        tree_cutoff.css('left', cutoff_pos);
      },
      slide: async function( event, ui ) {
        move_arrow();
        
        cutoff_date.text(xaxis_to_date(ui.value/tree_multiplier, tips[0], d3.min(tips, function(d) {return d.first_date}), d3.max(tips, function(d) {return d.last_date})));
        await handle.change()
          
        cutoff_line.css('visibility', 'visible');
        var cutoff_pos = handle.position().left;
        cutoff_line.css('left', cutoff_pos + 29);
        tree_cutoff.css('left', cutoff_pos);
      },
      change: async function (event, ui) {
        if (event.originalEvent) {
          cutoff_line.css('visibility', 'hidden');
          var cutoff_pos = handle.position().left;
          tree_cutoff.css('left', cutoff_pos);
  
          $("#loading").show();
          $("#loading_text").text(i18n_text.loading);
          switch($("#display-tree").val()) {
            case "XBB Lineages":
              await redraw_tree(df_xbb, cutoff_date.text());
              break;
            default:
              await redraw_tree(df, cutoff_date.text());
          }
          $("#loading").hide();
          $("#loading_text").text(``);
        }
      },
      min: min,
      max: max,
      value: (start_value > min && start_value < max) ? start_value : min
    });

    // Prevents the default action when keydown event is detected
    handle.unbind('keydown')

  } );

  // Sets the beadplot and time axis to move when the horizontal scrollbar is moved
  $('#beadplot-hscroll').scroll(function() {
    $("#svg-cluster").scrollLeft($(this).scrollLeft());
    $('#svg-clusteraxis').scrollLeft($(this).scrollLeft());
  });

  // Sets the beadplot and time axis to move when the vertical scrollbar is moved
  $('#beadplot-vscroll').scroll(function() {
    $("#svg-cluster").scrollTop($(this).scrollTop());
    $('#svg-clusteraxis').scrollTop($(this).scrollTop());
  });

  // Sets the time axis, vertical and horizontal scrollbar to move when scrolling through the beadplot
  $('#svg-cluster').scroll(function(e) {
    $("#beadplot-hscroll").scrollLeft($(this).scrollLeft());
    $('#svg-clusteraxis').scrollLeft($(this).scrollLeft());
    $("#beadplot-vscroll").scrollTop($(this).scrollTop());
  });

  // Vertical scrollbar for the time-scaled tree
  $('#tree-vscroll').scroll(function() {
    $("#svg-timetree").scrollTop($(this).scrollTop());
  });

  $('#svg-timetree').scroll(function() {
    $("#tree-vscroll").scrollTop($(this).scrollTop());
  });

  // Closing the modal
  $('.closeBtn').on('click', function () {
    $('.modal').fadeOut(300);
  });

  // Close the modal when clicking outside the modal
  $('.modal').on('click', function () {
    $('.modal').fadeOut(300);
  }).children().click(function () {
    return false;
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
        $("#loading_text").text(i18n_text.loading);
        await select_next_prev_bead(bead_id_to_accession, curr_bead);
        select_working_bead(bead_id_to_accession, curr_bead);
        $("#loading").hide();
        $("#loading_text").text(``);
      }
    }
    else if (curr_bead - 1 >= 0) {
      var curr_cid, prev_cid;
      await fetch(`/api/cid/${bead_id_to_accession[curr_bead]}`)
      .then(response => response.text())
      .then(data => curr_cid = data);

      await fetch(`/api/cid/${bead_id_to_accession[curr_bead - 1]}`)
      .then(response => response.text())
      .then(data => prev_cid = data);

      // If the previous bead is not in the same cluster, selection of cluster needs to be modified
      if (curr_cid !== prev_cid) {
        $("#loading").show();
        $("#loading_text").text(i18n_text.loading);
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

  // Moves the slider to the left/right when the arrow buttons are clicked
  $('#left-arrow').click(function() {
    $("#custom-handle").addClass("ui-slider-active")
    move_slider("LEFT");
  });

  $('#right-arrow').click(function() {
    $("#custom-handle").addClass("ui-slider-active")
    move_slider("RIGHT");
  });

  // Adds/Removes the selected/active state of the slider depending on where the user clicks
  $(document).on('mousedown', function(e) {
    if (e.target.matches("div.larrow")) 
      $('.larrow').addClass("selected")
    else if (e.target.matches("div.rarrow")) 
      $('.rarrow').addClass("selected")
    else {
      $('.rarrow').removeClass("selected")
      $('.larrow').removeClass("selected")
      $("#custom-handle").removeClass("ui-slider-active")
    }
  });

  // Listens to the keyup event to change the background of the arrow buttons to the default color
  $(document).on('keyup', function(e) {
    if(e.target.matches("div.ui-slider-handle") || $('.larrow').hasClass('selected') || $('.rarrow').hasClass('selected')) {
      if (e.keyCode == 37) 
        $('#left-arrow').removeClass("clicked");
      else if (e.keyCode == 39) 
        $('#right-arrow').removeClass("clicked");
      
      return;
    }
  });

  $(document).on('keydown', function(e) {
    // Ignore event if its inside an input field
    if (e.target.matches('input')) {
      return;
    }

    // If the slider or the arrow buttons are selected, then the arrow button's background color is changed and the slider is moved
    if(e.target.matches("div.ui-slider-handle") || $('.larrow').hasClass('selected') || $('.rarrow').hasClass('selected')) {
      if (e.keyCode == 37) {
        $('#left-arrow').addClass("clicked");
        move_slider("LEFT");
      }
      else if (e.keyCode == 39) {
        $('#right-arrow').addClass("clicked");
        move_slider("RIGHT");
      }
      return;
    }

    // User presses the left arrow key (37) or right arrow key (39)
    if (e.keyCode == 37 || e.keyCode == 39) {
      e.preventDefault() //issue #352
      var selected_bead = d3.selectAll(".selectionH").nodes();

      if (selected_bead.length == 0) {
        var points_ui = d3.selectAll("#svg-cluster > svg > g > circle").nodes()[0];
        var selected_bead = d3.selectAll('circle[id="'+points_ui.__data__.accessions[0]+'"]');
        selected_bead.raise();
        var working_bead = selected_bead.nodes()[0];
        working_bead.scrollIntoView({block: "center"});
        update_table_individual_bead_front(d3.select(working_bead).datum());
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

// Populate mutation details table

// Prepare the legend
var mut_table_legend = []

var num_labels = 0, phenotype_labels = [];
for (const [label, link] of Object.entries(phenotypes)) {
  if (num_labels != 0 && num_labels % 2 == 0) {
    mut_table_legend.push(phenotype_labels);
    phenotype_labels = []
  }
  phenotype_labels.push({
    label: i18n_text.phenotypes[label],
    src: link
  });
  num_labels++;
}

if (phenotype_labels.length !== 0) {
  mut_table_legend.push(phenotype_labels)
}

var mut_table = d3.select('#mut-table').append('table');
var mut_thead = mut_table.append('thead')
var mut_theaders = i18n_text.mutation_threaders;
var mut_tbody = mut_table.append("tbody")

var mut_legend_rows = d3.select("#muttable-legend")
                        .append("tbody")
                        .selectAll("tr")
                        .data(mut_table_legend)
                        .enter()
                        .append("tr")

var mut_legend_cells = mut_legend_rows
                        .selectAll("td")
                        .data(function (r) { return r.slice(0,2); })
                        .enter()
                        .append("td")

mut_legend_cells
  .append("img")
  .attr("src", function(d) { return d.src })
  .attr("class", "phenotype_icon")

mut_legend_cells
  .append("text")
  .text(function(d) { return d.label; })

// implement acknowledgements dialog
$( "#dialog" ).dialog({ autoOpen: false });


/*********************** FILE SAVING ***********************/
// implement save buttons
var blob;
function save_timetree() {
  var filename = $("#display-tree").val() === "XBB Lineages" ? "xbbtree.nwk" : "timetree.nwk"

  $.ajax({
    url: `data/${filename}`,
    success: function(data) {
      nwk = data;
      blob = new Blob([nwk], {type: "text/plain;charset=utf-8"});
      saveAs(blob, filename);
    }
  });
}

function save_beadplot() {
  // download beadplot as a Newick tree string
  blob = new Blob([serialize_beadplot(cindex)],
      {type: "text/plain;charset=utf-8"});
  saveAs(blob, lineage + ".nwk");
}

function export_svg() {
  // download beadplot as an SVG file
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

function export_csv() {
  var all_tips = [...tips, ...recombinant_tips, ...df_xbb];

  // write lineage-level information to CSV file for download
  var csvFile = 'lineage,mean.diffs,clock.residual,num.cases,num.variants,min.coldate,max.coldate,mean.coldate';
  var lineage_info = []
  for (tip of all_tips) {
    if (tip.isTip === undefined || tip.isTip)
      lineage_info.push([`${tip.thisLabel === undefined ? tip.label1 : tip.thisLabel},${Math.round(100*tip.mean_ndiffs)/100.},${Math.round(100*tip.residual)/100.},${tip.nsamples},${tip.varcount},${formatDate(tip.first_date)},${formatDate(tip.last_date)},${formatDate(tip.mcoldate)}`]);
  }
  csvFile = csvFile + "\n" + lineage_info.join("\n");
  blob = new Blob([csvFile], {type: "text/csv"});
  saveAs(blob, "lineage_stats.csv");
}

function export_muttable() {
  // download mutation frequencies (>=50%) for every lineage as CSV file
  var csvFile = 'lineage,mutation,frequency\n';
  var mutation_info = [], muts;
  for (lineage in dbstats['lineages']) {
    muts = dbstats['lineages'][lineage]['mutations']
    for (const [mut, count] of Object.entries(muts)) {
      mutation_info.push([`${lineage},${mut},${count}`]);
    }
  }
  csvFile = csvFile + mutation_info.join("\n") + "\n";
  blob = new Blob([csvFile], {type: "text/csv"});
  saveAs(blob, "mutations.csv");
}

function set_height() {
  // Calculate the height for the tree container and beadplot container
  $('#tree-vscroll').css('height',$(window).height() - $('#tree-vscroll').offset().top - 50)
  $('.tree-content').css('height',$(window).height() - $('.tree-content').offset().top - 50)
  $('#cutoff-line').css('height',$(window).height() - $('#cutoff-line').offset().top - 50)
  $('.beadplot-content').css('height',$(window).height() - $('.beadplot-content').offset().top - 50)
  $('#beadplot-vscroll').css('height',$(window).height() - $('#beadplot-vscroll').offset().top - 50)
}