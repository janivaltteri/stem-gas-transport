module parameters;

import std.file;
import std.json;
import std.stdio;
import std.conv;

struct Parameters
{
  bool initialised;
  
  int nr;
  int ny;
  int nrs;

  int res;

  double dy;
  double dr;
  double dt;
  
  double vel;
  double gamma;
  double tort;
  double tortax;
  double henry;
  double height;
  double width;
  double radius;
  double camb;
  double csoil;
  double fractair;
  double fractwat;
  double fractwood;

  double diff_r;
  double diff_a;
  double diff_b;

  void set_defaults(){

    this.nr = 6;
    this.ny = 7;
    this.nrs = 0;
    this.res = 1;
    
    this.vel = 0.75e-7; // 0.75e-7;
    this.gamma = 0.1;
    this.tort = 1.0e-5;
    this.tortax = 1.0e-2;
    this.henry = 0.000355;
    this.height = 0.48;
    this.radius = 0.005;
    this.camb = 0.0;
    this.csoil = 0.15;
    this.fractwat = 0.27;
    this.fractair = 0.44;
    this.fractwood = 1.0 - (this.fractair + this.fractwat);

    this.dt = 0.1;
    this.dy = this.height / to!double(this.ny);

    this.dr = this.radius / to!double(this.nr);

    this.diff_r = this.tort * 2.0 * 1.0e-5;
    this.diff_a = this.tortax * 2.0 * 1.0e-5;
    this.diff_b = this.diff_r;

    initialised = true;
  }

  bool read_json(in string file){
    debug(2) writeln("reading parameters from ",file);
    bool ok = true;
    if(file.exists()){
      string read = readText(file);
      JSONValue j = parseJSON(read);
      if("nr" in j){  nr  = to!int(j["nr"].integer);  }else{ writeln("no nr");  ok = false; }
      if("ny" in j){  ny  = to!int(j["ny"].integer);  }else{ writeln("no ny");  ok = false; }
      if("nrs" in j){ nrs = to!int(j["nrs"].integer); }else{ writeln("no nrs"); ok = false; }
      if("res" in j){ res = to!int(j["res"].integer); }else{ writeln("no res"); ok = false; }
      if("vel" in j){   vel   = j["vel"].floating;   }else{ writeln("no vel");   ok = false; }
      if("gamma" in j){ gamma = j["gamma"].floating; }else{ writeln("no gamma"); ok = false; }
      if("henry" in j){ henry = j["henry"].floating; }else{ writeln("no henry"); ok = false; }
      if("csoil" in j){ csoil = j["csoil"].floating; }else{ writeln("no csoil"); ok = false; }
      if("dt" in j){ dt = j["dt"].floating; }else{ writeln("no dt"); ok = false; }

      // height
      if("height" in j){
	height = j["height"].floating;
      }else{
	writeln("no height"); ok = false;
      }

      // radius
      if("radius" in j){
	radius = j["radius"].floating;
      }else{
	writeln("no radius"); ok = false;
      }

      // dr
      /*
      if("dr" in j){
	JSONValue dr_array = j["dr"].array; //get!(JSONValue[]); //j["dr"].floating;
	dr.length = ny;
	// todo: tarkista oikea määrä
	for(auto i = 0; i < ny; i++){
	  dr[i] = dr_array[i].get!float;
	}
      }else{
	writeln("no dr"); ok = false;
      }
      */

      // fractions
      if("fractair" in j){
	fractair = j["fractair"].floating;
      }else{
	writeln("no fractair"); ok = false;
      }
      if("fractwat" in j){
	fractwat = j["fractwat"].floating;
      }else{
	writeln("no fractwat"); ok = false;
      }
      if((fractair + fractwat) > 1.0){
	writeln("fractions of air and water larger than 1.0"); ok = false;
      }else{
	fractwood = 1.0 - (fractair + fractwat);
      }
      
      // either tort and tortax given or diffk and diffkax given
      if("diff_r" in j){
	if("diff_a" in j){
	  diff_r = j["diff_r"].floating;
	  diff_a = j["diff_a"].floating;
	  tort = diff_r / (2.0 * 1.0e-5);
	  tortax = diff_a / (2.0 * 1.0e-5);
	}else{
	  writeln("diff_r without diff_a"); ok = false;
	}
      }else if("tort" in j){
	if("tortax" in j){
	  tort = j["tort"].floating;
	  tortax = j["tortax"].floating;
	  diff_r = tort * 2.0 * 1.0e-5;
	  diff_a = tortax * 2.0 * 1.0e-5;
	}else{
	  writeln("tort without tortax"); ok = false;
	}
      }else{
	writeln("neither diff_r not tort given"); ok = false;
      }

      if("camb" in j){
	camb = j["camb"].floating;
      }else{
	camb = 0.0;
      }

      if("diff_b" in j){
	diff_b = j["diff_b"].floating;
	//debug(2) writeln("diff_b set to ",diff_b);
      }else{
	diff_b = diff_r;
	//debug(2) writeln("diff_b not set, using diff_r");
      }

      dy = height / to!double(ny);
      dr = radius / to!double(nr);
      
    }else{
      writeln("file ",file," not found");
      ok = false;
    }

    if(!ok) writeln("error in reading parameters");
    
    return ok;
  }

  void print(){
    writeln("parameters:");
    writefln("nr %s ny %s nrs %s", nr, ny, nrs);
    writefln("dy %s dr %s", dy, dr);
    writefln("diff_r %s diff_a %s diff_b %s gamma %s",
	     diff_r, diff_a, diff_b, gamma);
  }

}
