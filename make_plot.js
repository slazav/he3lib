// default options (can be changed in the document if needed)

// for the whole plot
var plot_defs = {
    width:  640, // global width/height
    height: 480,
    lmarg:   60, // margins
    rmarg:   10,
    tmarg:   10,
    bmarg:   30,
    xlabel:  "X", // x/y labels
    ylabel:  "Y",
    xticks:  10,  // number of ticks
    yticks:  10,
    xmin: "auto", // plot range
    xmax: "auto",
    ymin: "auto",
    ymax: "auto",
    xfmt: "g",
    yfmt: "g",
    xlog: 0, // log scale base or 0 for linear scale
    ylog: 0,
    xmarg: 0, // data x margings
    ymarg: 0, // data y margings
    data: [] // plot data
}

// for a data plot
var plot_line_defs = {
    xcol:  0, // x/y columns in the data
    ycol:  1,
    xsc:   1,
    ysc:   1,
    color: "steelblue",
    width: 3,
    url:   "", // url
    class: "line"
}

// for points
var plot_point_defs = {
    xpt:   0,
    ypt:   0,
    color: "black",
    width: 6,
    class: "point"
}

// for text
var plot_text_defs = {
    color:  "blue",
    weight: "bold",
    align:  "start", // start/middle/end
    class:  "text"
}


/***************************************************************/
// plot_setdef -- set default parameters
// Usage:
//   opts={ a:1, b:2, c:3 };
//   defs={ a:'a1', d:'d1' };
//   plot_setdef(opts, defs);
// Result:
//   opts: {a:1, b:2, c:3, d:d1}

function plot_setdef(opts, defs){
  for (var par in defs){
    if (typeof opts[par] === 'undefined') opts[par] = defs[par];
  }
}

/***************************************************************/
// Update data range and data coordinates on the plot

function plot_update(plot){
  // run only once, after reading all the data
  if (plot.total != plot.ready) return;

  // calculate the ranges
  var xmin = plot.opts.xmin;
  var xmax = plot.opts.xmax;
  var ymin = plot.opts.ymin;
  var ymax = plot.opts.ymax;
  if (!isFinite(xmin)){
    var l = d3.min(plot.lines,  function(d) { return d3.min(d.xy, function(g) {return g.x;}); });
    var p = d3.min(plot.points, function(d) { return d.x; });
    if (l === undefined) l=p;
    if (p === undefined) p=l;
    xmin = Math.min(l, p) - plot.opts.xmarg;
  }
  if (!isFinite(xmax)){
    var l = d3.max(plot.lines,  function(d) { return d3.max(d.xy, function(g) {return g.x;}); });
    var p = d3.max(plot.points, function(d) { return d.x; });
    if (l === undefined) l=p;
    if (p === undefined) p=l;
    xmax = Math.max(l, p) + plot.opts.xmarg;
  }
  if (!isFinite(ymin)){
    var l = d3.min(plot.lines,  function(d) { return d3.min(d.xy, function(g) {return g.y;}); });
    var p = d3.min(plot.points, function(d) { return d.y; });
    if (l === undefined) l=p;
    if (p === undefined) p=l;
    ymin = Math.min(l, p) - plot.opts.ymarg;
  }
  if (!isFinite(ymax)){
    var l = d3.max(plot.lines,  function(d) { return d3.max(d.xy, function(g) {return g.y;}); });
    var p = d3.max(plot.points, function(d) { return d.y; });
    if (l === undefined) l=p;
    if (p === undefined) p=l;
    ymax = Math.max(l, p) + plot.opts.ymarg;
  }

  // update axes
  plot.x.domain([xmin,xmax]);
  plot.y.domain([ymin,ymax]);
  plot.svg_xax.call(plot.d3_xax);
  plot.svg_yax.call(plot.d3_yax);

  // Update coordinates of all objects after changing scales
  plot.line_layer.selectAll("path").data(plot.lines)
    .attr("d", function(d) { return plot.valueline(d.xy); });
  plot.point_layer.selectAll("circle").data(plot.points)
    .attr("cx", function(d) { return plot.x(d.x);})
    .attr("cy", function(d) { return plot.y(d.y);});
  plot.text_layer.selectAll("text").data(plot.texts)
    .attr("x", function(d) { return plot.x(d.x);})
    .attr("y", function(d) { return plot.y(d.y);});
}


/***************************************************************/
// plot_add_line -- add data to the plot

function plot_add_line(plot, text, opts){

  plot_setdef(opts, plot_line_defs); // set default options

  // extract data from a text
  data=[];
  var lines = text.split("\n");
  for (i=0; i<lines.length; i++){
    if (lines[i].charAt(0)=="#") continue;
    var fields = lines[i].split(" ");
    var x = +fields[opts.xcol]*opts.xsc;
    var y = +fields[opts.ycol]*opts.ysc;
    if (isFinite(x) && isFinite(y)){ data.push({x:x, y:y}); }
  }

  // push data to the global dataset, plot it
  plot.lines.push({
    xy:    data,
    color: opts.color,
    width: opts.width,
  });
  var lines = plot.line_layer.selectAll("path").data(plot.lines);
  lines.enter().append("path")
    .attr("class",        "line")
    .attr("fill",         "none")
    .attr("stroke",       function(d) { return d.color;})
    .attr("stroke-width", function(d) { return d.width;});
  lines.exit().remove();

  plot.ready+=1;
  plot_update(plot);
}

/***************************************************************/
// plot_add_point -- add point to the plot
function plot_add_point(plot, opts){
  plot_setdef(opts, plot_point_defs); // set default options

  plot.points.push({ x:opts.xpt, y:opts.ypt });

  var points = plot.point_layer.selectAll("circle").data(plot.points);
  points.enter().append("circle")
      .attr("class", opts.class)
      .attr("r",     opts.width/2.0)
      .attr("fill",  opts.color);
  points.exit().remove();

  plot.ready+=1;
  plot_update(plot);
}


/***************************************************************/
// plot_add_text -- add a text to the plot
function plot_add_text(plot, opts){
  plot_setdef(opts, plot_text_defs); // set default options

  plot.texts.push({ x:opts.x, y:opts.y });

  var texts = plot.text_layer.selectAll("text").data(plot.texts);
  texts.enter().append("text")
     .text(opts.text)
     .attr("class",       opts.class)
     .attr("fill",        opts.color)
     .attr("font-weight", opts.weight)
     .attr("text-anchor", opts.align);
  texts.exit().remove();

  plot.ready+=1;
  plot_update(plot);
}


/***************************************************************/
/***************************************************************/
// make a single plot

function make_plot(svg, opts){

  plot_setdef(opts, plot_defs); // set default options

  // Clear the svg canvas
  svg.text("");

  if (opts.data.length==0){ // no data - no plots!
    svg.attr("width", 0).attr("height", 0);
    return;
  }

  var plot = {
    total: opts.data.length, ready: 0, // total objects, processed objects
    lines: [], points: [], texts: []   // data for plotting objects
  };

  plot.svg = svg
    .attr("width", opts.width)
    .attr("height", opts.height)
    .append("g")
      .attr("transform",
        "translate(" + opts.lmarg + "," + opts.tmarg + ")");


  // Set the dimensions of the canvas / graph
  plot.width  = opts.width  - opts.lmarg - opts.rmarg,
  plot.height = opts.height - opts.tmarg - opts.bmarg;
  plot.opts = opts;

  // Set the ranges
  if (opts.xlog>0){
    plot.x = d3.scale.log().range([0, plot.width]).base(opts.ylog);
  } else {
    plot.x = d3.scale.linear().range([0, plot.width]);
  }
  if (opts.ylog>0){
    plot.y = d3.scale.log().range([plot.height,0]).base(opts.ylog);
  } else {
    plot.y = d3.scale.linear().range([plot.height, 0]);
  }

  // Define the axes
  plot.d3_xax = d3.svg.axis()
    .scale(plot.x)
    .orient("bottom")
    .tickSize(plot.height);
  plot.d3_yax = d3.svg.axis()
    .scale(plot.y)
    .orient("left")
    .tickSize(plot.width);

  if (opts.xlog>0){
    plot.d3_xax.ticks(opts.xticks, opts.xfmt);
  } else{
    plot.d3_xax.ticks(opts.xticks).tickFormat(d3.format(opts.xfmt));
  }
  if (opts.ylog>0){
    plot.d3_yax.ticks(opts.yticks, opts.yfmt);
  } else {
    plot.d3_yax.ticks(opts.yticks).tickFormat(d3.format(opts.yfmt));
  }

//  if (opts.yexp!=0){
//    yaxis.tickFormat(function(d) { return (d/Math.pow(10,opts.yexp)); })
//         .tickFormat(d3.format(',.2f/20'));
//    opts.ylabel = opts.ylabel + ", x1e" + opts.yexp;
//  }
//  if (opts.xexp!=0){
//    xaxis.tickFormat(function(d) { return (d/Math.pow(10,opts.xexp)); })
//         .tickFormat(d3.format(',.2f/20'));
//    opts.xlabel = opts.xlabel + ", x1e" + opts.xexp;
//  }

  // Define the line
  plot.valueline = d3.svg.line()
      .x(function(d) { return plot.x(d.x); })
      .y(function(d) { return plot.y(d.y); });

  // Add the X Axis
  plot.svg_xax = plot.svg.append("g").attr("class", "x axis")

  // Add the Y Axis
  plot.svg_yax = plot.svg.append("g").attr("class", "y axis")
    .attr("transform", "translate(" + plot.width + ",0)");

  // Layers for lines, points, texts
  plot.line_layer  = plot.svg.append("g").attr("class", "lines")
  plot.point_layer = plot.svg.append("g").attr("class", "points")
  plot.text_layer  = plot.svg.append("g").attr("class", "texts")

  // X label
  plot.svg_xlab = plot.svg.append("text")
    .text(opts.xlabel)
    .attr("class", "x label")
    .attr("x", plot.width/2)
    .attr("y", plot.height + opts.bmarg)
    .attr("text-anchor", "middle");

  // Y label
  plot.svg_ylab = plot.svg.append("text")
    .text(opts.ylabel)
    .attr("class", "y label")
    .attr("text-anchor", "middle")
    .attr("dominant-baseline", "text-before-edge")
    .attr("transform", "translate("+(-opts.lmarg)+","+(plot.height/2)+"),rotate(-90)");

  // make requests for getting and processing all elements
  for (var i=0; i<opts.data.length; i++){

    if (opts.data[i].url !== undefined){
      var f = function(error, text, i) {
        if (error) alert(error);
        var optsd = opts.data[i];
        plot_add_line(plot, text, optsd);
      }
      d3.text(opts.data[i].url, // I hate javascript!
        (function(i2) { return function(error, text){f(error, text, i2);}})(i)
      );
    }

    else if (opts.data[i].xpt !== undefined &&
             opts.data[i].ypt !== undefined){
        plot_add_point(plot, opts.data[i]);
    }

    else if (opts.data[i].text !== undefined){
        plot_add_text(plot, opts.data[i]);
    }

  }

  return plot;
}

/***************************************************************/
// make_mplot -- plot selector
// Usage:
//  make_mplot([
//    {title: 'Select a plot...'}
//
//    {title: 'Fermi-liquid parameter F0a vs pressure',
//     xlabel: "P, bar", ylabel: "F0a",
//     data: [{ url: "data/he3_f0a&0:0.2:34"}] },
//
//    {title: 'Fermi-liquid parameter F0s vs pressure',
//     xlabel: "P, bar", ylabel: "F0s",
//     data: [{ url: "data/he3_f0s&0:0.2:34"}] },
//  ]);

function make_mplot(data){

  if (data.length<1) return;

  var c   = d3.select("body").append("div"); // centered block
  var sel = c.append("form").append("select");  // the selector
  var svg = c.append("svg").attr("width", 0).attr("height", 0) // SVG

  // make selector options
  sel.selectAll("option").data(data)
    .enter().append("option")
      .text(function(d) { return d.title; })

  // callback with make_plot() on the selector change
  sel.on("change", function() {
    make_plot(svg, data[this.selectedIndex]);
  });

  // make default plot
  make_plot(svg, data[0]);
}
