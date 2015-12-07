function make_plot(opts){
  // setting default values:
  function set_def(par, def){
    if (typeof opts[par] === 'undefined') opts[par] = def;
  }

  var plot={};

  // global height/width
  set_def("width",  640);
  set_def("height", 480);
  // margins
  set_def("lmarg", 50);
  set_def("rmarg", 10);
  set_def("tmarg", 10);
  set_def("bmarg", 30);
  // x/y labels
  set_def("xlabel", "X");
  set_def("ylabel", "Y");
  // number of ticks
  set_def("xticks", 10);
  set_def("yticks", 10);
  // number of ticks
  set_def("xexp", 0);
  set_def("yexp", 0);
  // range
  set_def("xmin", "auto");
  set_def("xmax", "auto");
  set_def("ymin", "auto");
  set_def("ymax", "auto");

  // Set the dimensions of the canvas / graph
  plot.width  = opts.width  - opts.lmarg - opts.rmarg,
  plot.height = opts.height - opts.tmarg - opts.bmarg;

  // Set the ranges
  var x = d3.scale.linear().range([0, plot.width]);
  var y = d3.scale.linear().range([plot.height, 0]);

  // Define the axes
  var xaxis = d3.svg.axis()
    .scale(x)
    .orient("bottom")
    .ticks(opts.xticks)
    .tickSize(plot.height);
  var yaxis = d3.svg.axis()
    .scale(y)
    .orient("left")
    .ticks(opts.yticks)
    .tickSize(plot.width);

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
  var valueline = d3.svg.line()
      .x(function(d) { return x(d.X); })
      .y(function(d) { return y(d.Y); });

  // Adds the svg canvas
  plot.svg = d3.select("body")
    .append("svg")
      .attr("width", opts.width)
      .attr("height", opts.height)
    .append("g")
      .attr("transform",
            "translate(" + opts.lmarg + "," + opts.tmarg + ")");

  // Add the X Axis
  plot.svg.append("g")
    .attr("class", "x axis")
    .attr("id","xaxis")
    .call(xaxis);

  // Add the Y Axis
  plot.svg.append("g")
    .attr("class", "y axis")
    .attr("id","yaxis")
    .attr("transform", "translate(" + plot.width + ",0)")
    .call(yaxis);

  // X and Y labels
  plot.svg.append("text")
    .attr("class", "x label")
    .attr("id", "xlabel")
    .attr("x", plot.width/2)
    .attr("y", plot.height + opts.bmarg)
    .text(opts.xlabel)
    .attr("text-anchor", "middle");
  plot.svg.append("text")
    .attr("class", "y label")
    .attr("id", "ylabel")
    .text(opts.ylabel)
    .attr("text-anchor", "middle")
    .attr("transform", "translate("+(-opts.lmarg+10)+","+(plot.height/2)+"),rotate(-90)");
  plot.data=[];

  function update(){
    // update plots
    var lines = plot.svg.selectAll(".line").data(plot.data);
    lines.attr("d", function(d) { return valueline(d.data); });
    lines.enter().append("path").attr("class", "line")
      .attr("d", function(d) { return valueline(d.data); });
    lines.exit().remove();
  }


  // Get and plot the data
  plot.add = function(URL, opts) {
    function set_def(par, def){
      if (typeof opts[par] === 'undefined') opts[par] = def; }
    // x/y columns in the data
    set_def("xcol", 0);
    set_def("ycol", 1);
    set_def("xsc",  1);
    set_def("ysc",  1);
    // x/y labels
    set_def("xlabel", "");
    set_def("ylabel", "");
    set_def("color", "steelblue");

    d3.text(URL, function(error, text) {
      if (error) alert(error);
      var lines = text.split("\n");
      var data = [];
      for (i=0; i<lines.length; i++){ //>
        if (lines[i].charAt(0)=="#") continue;
        var fields = lines[i].split(" ");
        data.push({X: +fields[opts.xcol]*opts.xsc, Y: +fields[opts.ycol]*opts.ysc});
      }
      var xmin =  d3.min(data, function(d) { return d.X; });
      var xmax =  d3.max(data, function(d) { return d.X; });
      var ymin =  d3.min(data, function(d) { return d.Y; });
      var ymax =  d3.max(data, function(d) { return d.Y; });
      var ym = (ymax-ymin)*0.05;
      ymin = ymin-ym;
      ymax = ymax+ym;
      if (plot.data.length){
        xmin = Math.min(xmin, x.domain()[0]);
        xmax = Math.max(xmax, x.domain()[1]);
        ymin = Math.min(ymin, y.domain()[0]);
        ymax = Math.max(ymax, y.domain()[1]);
      }
      x.domain([xmin,xmax]);
      y.domain([ymin,ymax]);
      d3.select("#xaxis").call(xaxis);
      d3.select("#yaxis").call(yaxis);
      if (opts.xlabel) d3.select("#xlabel").text(opts.xlabel);
      if (opts.ylabel) d3.select("#ylabel").text(opts.ylabel);
      plot.data.push({data:data, xmin:xmin, xmax:xmax, ymin:ymin, ymax:ymax});
      update();
    });
  }

  plot.clear = function(){ plot.data=[]; update(); }

  return plot;
}
