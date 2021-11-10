import std.csv;
import std.conv;
import std.file;
import std.json;
import std.stdio;
import std.algorithm;
import std.typecons : Tuple;

import commandr;

import tree;
import array2d;
import parameters;

// todo: update the json out writer

void main(string[] args)
{
  auto a = new Program("stem-gas-transport model, cli-interface", "0.4")
    .summary("Code for manuscript")
    .author("Jani Anttila <janivaltteri@posteo.net>")
    .add(new Option("i", "infile",    "parameter json infile").required)
    .add(new Option("o", "outfile",   "outfile prefix").defaultValue(""))
    .add(new Option("s", "steps",     "maximum number of steps").defaultValue("20"))
    .add(new Option("v", "velocities","input axial flows infile").defaultValue(""))
    .add(new Option("m", "method",    "numerical method")
	 .defaultValue("euler").validate(new EnumValidator(["euler","rk4","original"])))
    .add(new Option("p", "printmode", "print either final state or all time-series")
	 .defaultValue("final").validate(new EnumValidator(["final","all"])))
    .add(new Flag("j", null, "json output").name("jsonout"))
    .add(new Flag("a", null, "storage output").name("storageout"))
    .add(new Flag("f", null, "flux output").name("fluxout"))
    .add(new Flag("e", null, "stop if equilibrium reached").name("equilibrate"))
    .add(new Option("t", "tolerance", "equilibrium tolerance").defaultValue("0.00001"))
    .parse(args);

  // input argument handling

  // parameter infile
  Parameters par;
  bool parameter_read_ok = par.read_json(a.option("infile"));
  par.print();

  // number of steps
  int nsteps; // = to!int(a.option("steps"));
  try{
    nsteps = to!int(a.option("steps"));
  }catch (ConvException ce){
    writeln("error parsing suggested steps value, using default 20");
    nsteps = 20;
  }

  // temporally variable input velocity
  bool velo_input_p = false;
  bool velo_read_ok = true;
  string velo_file = a.option("velocities");
  double[] velocities;
  if(velo_file != ""){
    velo_input_p = true;
    if(velo_file.exists()){
      velocities.length = nsteps;
      int counter = 0;
      auto vfile = File(velo_file, "r");
      foreach(record; vfile.byLine.joiner("\n").csvReader!(Tuple!(double))){
	if(counter >= nsteps){
	  writeln("error: ",velo_file," has more rows than steps");
	  velo_read_ok = false;
	  break;
	}else{
	  velocities[counter++] = record[0];
	}
      }
      if(counter < nsteps){
	writeln("error: ",velo_file," has fewer rows than steps");
	velo_read_ok = false;
      }
      vfile.close();
    }else{
      writeln("velocity file ",velo_file," does not exist!");
      velo_read_ok = false;
    }
  }

  // outfile names
  bool file_output;
  string outfile;
  if(a.option("outfile") == ""){
    file_output = false;
  }else{
    file_output = true;
    if(a.option("printmode") == "all"){
      outfile = a.option("outfile") ~ "-all.csv";
    }else if(a.option("printmode") == "final"){
      outfile = a.option("outfile") ~ "-final.csv";
    }else{
      outfile = a.option("outfile") ~ ".csv";
    }
    debug(2) writeln("outfile: ",outfile);
  }

  // numerical method
  int method;
  if(a.option("method") == "euler"){
    method = 0;
  }else if(a.option("method") == "rk4"){
    method = 1;
  }else if(a.option("method") == "original"){
    method = 2;
  }

  // json output
  bool json_output;
  string jsonfile;
  int json_output_flag = a.occurencesOf("jsonout");
  if(json_output_flag > 0){
    json_output = true;
    if(file_output){
      jsonfile = a.option("outfile") ~ ".json";
    }else{
      jsonfile = "out.json";
    }
    debug(2) writeln("json output to file ",jsonfile);
  }

  // full storage output
  bool storage_output;
  string storagefile;
  int storage_output_flag = a.occurencesOf("storageout");
  if(storage_output_flag > 0){
    storage_output = true;
    if(file_output){
      storagefile = a.option("outfile") ~ "-storage.csv";
    }else{
      storagefile = "out-storage.csv";
    }
    debug(2) writeln("storage output to file ",storagefile);
  }

  // elementwise flux output
  bool flux_output;
  string fluxfile;
  int flux_output_flag = a.occurencesOf("fluxout");
  if(flux_output_flag > 0){
    flux_output = true;
    if(file_output){
      fluxfile = a.option("outfile") ~ "-flux.csv";
    }else{
      fluxfile = "out-flux.csv";
    }
    debug(2) writeln("flux output to file ",fluxfile);
  }

  // equilibrium search
  bool equilibrium = false;
  // double[4] prev_summary = [ 0.0, 0.0, 0.0, 0.0 ]; // can't remember why I had this here
  double[4] summary = [ 0.0, 0.0, 0.0, 0.0 ];
  double tolerance;
  try{
    tolerance = to!double(a.option("tolerance"));
  }catch (ConvException ce){
    writeln("error parsing your suggested tolerance value, using default 0.00001");
    tolerance = 0.00001;
  }
  int equilibrium_search = a.occurencesOf("equilibrate");
  if(equilibrium_search > 0){
    writeln("stopping at equilibrium within tolerance ",tolerance);
  }

  // initialisation
  if(!parameter_read_ok){
    writeln("parameter reading failed, stopping");
  }else if(!velo_read_ok){
    writeln("velocity file reading failed, stopping");
  }else{
    Tree tree = new Tree(par);

    tree.initialise(par);
    
    File of; // regular output file
    File sf; // full storage output file
    File ff; // flux output file

    debug(2) writeln("simulating");

    // write headers for specified output files
    if(file_output){
      of = File(outfile,"w");
      of.writeln("time,storage,f_in,f_side,f_dtop,f_atop,c_in,c_side,c_dtop,c_atop");
    }
    if(storage_output){
      sf = File(storagefile,"w");
      sf.writeln("time,a,r,cair,cwat,nair,nwat,towt");
    }
    if(flux_output){
      ff = File(fluxfile,"w");
      ff.writeln("time,elem,type,value");
    }

    bool num_safe = true;

    // the actual simulation
    for(auto step = 0; step < nsteps; step++){

      // advance
      double velo;
      if(velo_input_p){
	velo = velocities[step];
      }else{
	velo = par.vel;
      }

      tree.euler_step(par, velo, method);

      // summary values are 0: flux in 1: flux side 2: flux top d 3: flux top a
      tree.fill_summary(summary);

      // check that these are positive and finite
      foreach(e; summary){
	if((e == double.nan) | (e < 0.0)){
	  num_safe = false;
	}
      }

      // write state
      if(a.option("printmode") == "all"){
	output(of,tree,file_output,false);
	// write storage
	if(storage_output){
	  storagewrite(sf,tree,par);
	}
	if(flux_output){
	  fluxwrite(ff,tree,par);
	}
      }

      if(num_safe){

	if(equilibrium_search > 0){
	  // test for equilibrium
	  double flux_out_sum = summary[1] + summary[2] + summary[3];
	  if((1.0 - (flux_out_sum / summary[0])) > tolerance){
	    equilibrium = false;
	  }else{
	    equilibrium = true;
	    writeln("equilibrium reached at time ",tree.time);
	    break;
	  }
	}

      }else{

	writeln("oops! numerical instability, stopping at step ",step);
	break;

      }
    }

    if(equilibrium_search > 0){
      if(!equilibrium){
	writeln("equilibrium not reached");
      }
    }

    // write final state
    if(a.option("printmode") == "final"){
      output(of,tree,file_output,false);
      if(storage_output){
	storagewrite(sf,tree,par);
      }
      if(flux_output){
	fluxwrite(ff,tree,par);
      }
    }

    // close files
    if(file_output) of.close();
    if(storage_output) sf.close();
    if(flux_output) ff.close();

    // json end state write
    if(json_output){
      json_write(jsonfile, tree, par, equilibrium, num_safe);
    }

  }

  debug(1) writeln("done!");
}

/**
 *  called if printmode == "all" or "final"
 */
void output(File f, ref Tree t,
	    bool file, bool with_influx)
{
  if(file){
    // write to file
    f.writeln(t.time,",",
	      t.storage_sum(),",",
	      t.fluxin_sum(),",",
	      t.fluxside_sum(),",",
	      t.fluxtop_d_sum(),",",
	      t.fluxtop_a_sum(),",",
	      t.fluxin_c,",",
	      t.fluxside_c,",",
	      t.fluxtop_d_c,",",
	      t.fluxtop_a_c);
  }else{
    // print to screen
    writeln("time: ",     t.time,
	    " storage: ", t.storage_sum(),
	    " fin: ",     t.fluxin_sum(),
	    " fside: ",   t.fluxside_sum(),
	    " fdtop: ",   t.fluxtop_d_sum(),
	    " fatop: ",   t.fluxtop_a_sum());
  }
}

/**
 *  called if -f is given
 */
void fluxwrite(File f, ref Tree t, ref Parameters p)
{
  for(auto i = 0; i < p.ny; i++){
    f.writeln(t.time,",",i,",side,",t.state.fluxside[i]);
  }
  for(auto j = 0; j < p.nr; j++){
    f.writeln(t.time,",",j,",top_d,",t.state.fluxtop_d[j]);
  }
  for(auto j = 0; j < p.nr; j++){
    f.writeln(t.time,",",j,",top_a,",t.state.fluxtop_a[j]);
  }
}

/**
 *  called if -a is given
 */
void storagewrite(File f, ref Tree t, ref Parameters p)
{
  for(auto i = 0; i < p.ny; i++){
    for(auto j = 0; j < p.nr; j++){
      f.writeln(t.time,",",i,",",j,",",
		t.state.cair[i,j],",",
		t.state.cwat[i,j],",",
		t.state.nair[i,j],",",
		t.state.nwat[i,j],",",
		t.derivative.towater[i,j]);
    }
  }
}

/**
 *  called if -j is given
 */
void json_write(in ref string jfile, ref Tree t, ref Parameters p,
		bool eql, bool num_safe){
  File jw = File(jfile,"w");
  JSONValue j = [ "time": t.time ];
  // parameters
  j.object["dt"]  = JSONValue(p.dt);
  j.object["dy"]  = JSONValue(p.dy);
  j.object["dr"]  = JSONValue(p.dr);
  j.object["nr"]  = JSONValue(p.nr);
  j.object["ny"]  = JSONValue(p.ny);
  j.object["nrs"] = JSONValue(p.nrs);
  j.object["res"] = JSONValue(p.res);
  j.object["vel"] = JSONValue(p.vel);
  j.object["gamma"] = JSONValue(p.gamma);
  j.object["diff_r"] = JSONValue(p.diff_r);
  j.object["diff_a"] = JSONValue(p.diff_a);
  j.object["diff_b"] = JSONValue(p.diff_b);
  j.object["henry"] = JSONValue(p.henry);
  j.object["height"] = JSONValue(p.height);
  j.object["radius"] = JSONValue(p.radius);
  j.object["camb"] = JSONValue(p.camb);
  j.object["csoil"] = JSONValue(p.csoil);
  j.object["fractwat"] = JSONValue(p.fractwat);
  j.object["fractair"] = JSONValue(p.fractair);
  j.object["fractwood"] = JSONValue(p.fractwood);
  if(eql){
    j.object["equilibrium"] = JSONValue(true);
  }else{
    j.object["equilibrium"] = JSONValue(false);
  }
  if(num_safe){
    // end state
    j.object["flux_in"]   = JSONValue(t.fluxin_sum());
    j.object["flux_side"] = JSONValue(t.fluxside_sum());
    j.object["flux_dtop"] = JSONValue(t.fluxtop_d_sum());
    j.object["flux_atop"] = JSONValue(t.fluxtop_a_sum());
    j.object["flux_top"]  = JSONValue(t.fluxtop_d_sum() + t.fluxtop_a_sum());
    j.object["cumu_in"]   = JSONValue(t.fluxin_c);
    j.object["cumu_side"] = JSONValue(t.fluxside_c);
    j.object["cumu_dtop"] = JSONValue(t.fluxtop_d_c);
    j.object["cumu_atop"] = JSONValue(t.fluxtop_a_c);
    j.object["cumu_top"]  = JSONValue(t.fluxtop_d_c + t.fluxtop_d_c);
    j.object["storage"]   = JSONValue(t.storage_sum());
    j.object["num_safe"]  = JSONValue(true);
  }else{
    j.object["num_safe"]  = JSONValue(false);
  }
  // write
  jw.writeln(j.toPrettyString);
  jw.close();
}
