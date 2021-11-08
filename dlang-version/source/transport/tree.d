module tree;

import std.algorithm.iteration;
import std.stdio;
import std.range;
import std.conv;
import std.math;

import parameters;
import array2d;
import state;
import derivative;

immutable double pii = 3.1415926535;

class Tree
{
 public:

  State state;
  Derivative derivative;

  Matrix volume;
  Matrix volumeair;
  Matrix volumewat;
  Matrix crossect;
  Matrix ker1;
  Matrix ker2;

  int nr;
  int ny;

  double fluxin_c;
  double fluxside_c;
  double fluxtop_a_c;
  double fluxtop_d_c;

  bool initialised;

  @property double time() { return state.time; }

  this(in ref Parameters p){

    // integers
    this.nr = p.nr;
    this.ny = p.ny;

    // doubles
    this.fluxin_c = 0.0;
    this.fluxside_c = 0.0;
    this.fluxtop_a_c = 0.0;
    this.fluxtop_d_c = 0.0;

    // matrices
    this.ker1.initialise(p.ny,p.nr);
    this.ker2.initialise(p.ny,p.nr);
    this.volume.initialise(p.ny,p.nr);
    this.crossect.initialise(p.ny,p.nr);
    this.volumeair.initialise(p.ny,p.nr);
    this.volumewat.initialise(p.ny,p.nr);

    // state and derivative
    this.state = State(p.ny,p.nr);
    this.derivative = Derivative(p.ny,p.nr);

    this.initialised = false;
  }

  /**
   *  called from main if radial increments file is supplied
   *  needs more testing...
   */
  void initialise_dr_array(in ref Parameters p, const double[] dr_values){

    for(auto i = 0; i < p.ny; i++){
      for(auto j = 0; j < p.nr; j++){
	crossect[i,j] = crossect_area(j,dr_values[i]);
      }
    }

    for(auto i = 0; i < p.ny; i++){
      volume[i,0]    = p.dy * crossect[i,0];
      volumeair[i,0] = p.fractair * volume[i,0];
      volumewat[i,0] = p.fractwat * volume[i,0];
      ker1[i,0]      = double.max;
      ker2[i,0]      = log(2.0 * dr_values[i] / dr_values[i]);
      for(auto j = 1; j < p.nr; j++){
	volume[i,j]    = p.dy * crossect[i,j];
	volumeair[i,j] = p.fractair * volume[i,j];
	volumewat[i,j] = p.fractwat * volume[i,j];
	// old way
	/*
	ker1[i,j] = log((dr_values[i] * to!double(j + 1)) /
			(dr_values[i] * to!double(j)));
	ker2[i,j] = log((dr_values[i] * to!double(j + 2)) /
			(dr_values[i] * to!double(j + 1)));
	*/
	// new way
	ker1[i,j] = ((dr_values[i] * to!double(j + 1)) -
		     (dr_values[i] * to!double(j)));
	ker2[i,j] = ((dr_values[i] * to!double(j + 2)) -
		     (dr_values[i] * to!double(j + 1)));
      }
    }

    zero_state(p);

    debug(1){
      write("dr:");
      for(auto i = 0; i < p.ny; i++){ write(" ",dr_values[i]); }
      writeln("");
    }

    this.initialised = true;
  }

  /**
   *  called from main if a fixed dr value is read from par
   *  this should be the only function which uses p.dr
   */
  void initialise_dr_fixed(in ref Parameters p){

    for(auto i = 0; i < p.ny; i++){
      for(auto j = 0; j < p.nr; j++){
	crossect[i,j] = crossect_area(j,p.dr);
      }
    }

    for(auto i = 0; i < p.ny; i++){
      volume[i,0]    = p.dy * crossect[i,0];
      volumeair[i,0] = p.fractair * volume[i,0];
      volumewat[i,0] = p.fractwat * volume[i,0];
      ker1[i,0]      = double.max;
      ker2[i,0]      = log(2.0 * p.dr / p.dr);
      for(auto j = 1; j < p.nr; j++){
	volume[i,j]    = p.dy * crossect[i,j];
	volumeair[i,j] = p.fractair * volume[i,j];
	volumewat[i,j] = p.fractwat * volume[i,j];
	// old way
	/*
	ker1[i,j] = log((p.dr * to!double(j + 1)) /
			(p.dr * to!double(j)));
	ker2[i,j] = log((p.dr * to!double(j + 2)) /
			(p.dr * to!double(j + 1)));
	*/
	// new way
	ker1[i,j] = ((p.dr * to!double(j + 1)) -
		     (p.dr * to!double(j)));
	ker2[i,j] = ((p.dr * to!double(j + 2)) -
		     (p.dr * to!double(j + 1)));
      }
    }

    zero_state(p);

    debug(1){
      write("initialise_dr_fixed: crossects ");
      for(auto j = 0; j < p.nr; j++){ write(" ",crossect[0,j]); }
      writeln("");
    }

    this.initialised = true;
  }

  /**
   *  called from initialise_dr_() functions
   */
  void zero_state(in ref Parameters p)
  {
    for(auto i = 0; i < p.ny; i++){
      for(auto j = 0; j < p.nr; j++){
        state.cair[i,j] = p.camb;
        state.cwat[i,j] = state.cair[i,j]*p.henry;
        state.nair[i,j] = state.cair[i,j]*volumeair[i,j];
        state.nwat[i,j] = state.cwat[i,j]*volumewat[i,j];
      }
    }
  }

  /**
   *  called for setting cross section areas in the constructor
   */
  pure double crossect_area(in int idx, in double dr){
    if(idx == 0){
      return pii*dr*dr;
    }else{
      return pii*pow(to!double(idx+1)*dr,2) - pii*pow(to!double(idx)*dr,2);
    }
  }

  void euler_step(in ref Parameters p, in double ax_vel, in int method){

    for(auto step = 0; step < p.res; step++){

      if(method == 0){

	// euler first order method. default.

	radial_diffusion(derivative,state,p);
	axial_advection(derivative,state,p,ax_vel);
	axial_diffusion(derivative,state,p);
	phase_equilibration(derivative,state,p);
	update_cumulants(p);
	for(auto i = 0; i < p.ny; i++){
	  for(auto j = 0; j < p.nr; j++){
	    state.nair[i,j] = state.nair[i,j] + p.dt*(derivative.qair[i,j] +
						      derivative.axdiff[i,j] -
						      derivative.towater[i,j]);
	    state.nwat[i,j] = state.nwat[i,j] + p.dt*(derivative.qadv[i,j] +
						      derivative.towater[i,j]);
	    state.cair[i,j] = state.nair[i,j]/volumeair[i,j];
	    state.cwat[i,j] = state.nwat[i,j]/volumewat[i,j];
	  }
	}

      }else if(method == 1){

	// runge-kutta 4th order method. needs testing.

	State s1 = State(p.ny,p.nr);
	State s2 = State(p.ny,p.nr);
	State s3 = State(p.ny,p.nr);
	State s4 = State(p.ny,p.nr);
	State s1i = State(p.ny,p.nr);
	State s2i = State(p.ny,p.nr);
	State s3i = State(p.ny,p.nr);
	Derivative d2 = Derivative(p.ny,p.nr);
	Derivative d3 = Derivative(p.ny,p.nr);
	Derivative d4 = Derivative(p.ny,p.nr);

	radial_diffusion(derivative,state,p);
	axial_advection(derivative,state,p,ax_vel);
	axial_diffusion(derivative,state,p);
	phase_equilibration(derivative,state,p);

	update_cumulants(p);

	for(auto i = 0; i < p.ny; i++){
	  for(auto j = 0; j < p.nr; j++){
	    s1.nair[i,j] = p.dt*(derivative.qair[i,j] +
				 derivative.axdiff[i,j] -
				 derivative.towater[i,j]);
	    s1.nwat[i,j] = p.dt*(derivative.qadv[i,j] +
				 derivative.towater[i,j]);
	    //s1.cair[i,j] = s1.nair[i,j]/volumeair[i,j];
	    //s1.cwat[i,j] = s1.nwat[i,j]/volumewat[i,j];
	    s1i.nair[i,j] = state.nair[i,j] + 0.5*s1.nair[i,j];
	    s1i.nwat[i,j] = state.nwat[i,j] + 0.5*s1.nwat[i,j];
	    s1i.cair[i,j] = s1i.nair[i,j]/volumeair[i,j];
	    s1i.cwat[i,j] = s1i.nwat[i,j]/volumewat[i,j];
	  }
	}

	radial_diffusion(d2,s1i,p);
	axial_advection(d2,s1i,p,ax_vel);
	axial_diffusion(d2,s1i,p);
	phase_equilibration(d2,s1i,p);

	for(auto i = 0; i < p.ny; i++){
	  for(auto j = 0; j < p.nr; j++){
	    s2.nair[i,j] = p.dt*(d2.qair[i,j] +
				 d2.axdiff[i,j] -
				 d2.towater[i,j]);
	    s2.nwat[i,j] = p.dt*(d2.qadv[i,j] +
				 d2.towater[i,j]);
	    //s2.cair[i,j] = s2.nair[i,j]/volumeair[i,j];
	    //s2.cwat[i,j] = s2.nwat[i,j]/volumewat[i,j];
	    s2i.nair[i,j] = state.nair[i,j] + 0.5*s2.nair[i,j];
	    s2i.nwat[i,j] = state.nwat[i,j] + 0.5*s2.nwat[i,j];
	    s2i.cair[i,j] = s2i.nair[i,j]/volumeair[i,j];
	    s2i.cwat[i,j] = s2i.nwat[i,j]/volumewat[i,j];
	  }
	}

	radial_diffusion(d3,s2i,p);
	axial_advection(d3,s2i,p,ax_vel);
	axial_diffusion(d3,s2i,p);
	phase_equilibration(d3,s2i,p);

	for(auto i = 0; i < p.ny; i++){
	  for(auto j = 0; j < p.nr; j++){
	    s3.nair[i,j] = p.dt*(d3.qair[i,j] +
				 d3.axdiff[i,j] -
				 d3.towater[i,j]);
	    s3.nwat[i,j] = p.dt*(d3.qadv[i,j] +
				 d3.towater[i,j]);
	    //s3.cair[i,j] = s3.nair[i,j]/volumeair[i,j];
	    //s3.cwat[i,j] = s3.nwat[i,j]/volumewat[i,j];
	    s3i.nair[i,j] = state.nair[i,j] + s3.nair[i,j];
	    s3i.nwat[i,j] = state.nwat[i,j] + s3.nwat[i,j];
	    s3i.cair[i,j] = s3i.nair[i,j]/volumeair[i,j];
	    s3i.cwat[i,j] = s3i.nwat[i,j]/volumewat[i,j];
	  }
	}

	radial_diffusion(d4,s3i,p);
	axial_advection(d4,s3i,p,ax_vel);
	axial_diffusion(d4,s3i,p);
	phase_equilibration(d4,s3i,p);

	for(auto i = 0; i < p.ny; i++){
	  for(auto j = 0; j < p.nr; j++){
	    s4.nair[i,j] = p.dt*(d4.qair[i,j] +
				 d4.axdiff[i,j] -
				 d4.towater[i,j]);
	    s4.nwat[i,j] = p.dt*(d4.qadv[i,j] +
				 d4.towater[i,j]);
	    //s4.cair[i,j] = s4.nair[i,j]/volumeair[i,j];
	    //s4.cwat[i,j] = s4.nwat[i,j]/volumewat[i,j];
	  }
	}

	for(auto i = 0; i < p.ny; i++){
	  for(auto j = 0; j < p.nr; j++){
	    state.nair[i,j] = state.nair[i,j] + (s1.nair[i,j] +
						 s2.nair[i,j]*2.0 +
						 s3.nair[i,j]*2.0 +
						 s4.nair[i,j])/6.0;
	    state.nwat[i,j] = state.nwat[i,j] + (s1.nwat[i,j] +
						 s2.nwat[i,j]*2.0 +
						 s3.nwat[i,j]*2.0 +
						 s4.nwat[i,j])/6.0;
	    state.cair[i,j] = state.nair[i,j]/volumeair[i,j];
	    state.cwat[i,j] = state.nwat[i,j]/volumewat[i,j];
	  }
	}

      }else if(method == 2){

	// original two-step method used in Hölttä & Kolari publication
	// currently this is not operational...

	/*
	radial_diffusion(derivative,state,p);
	axial_advection(derivative,state,p,ax_vel);
	axial_diffusion(derivative,state,p);
	phase_equilibration(derivative,state,p);
	update_cumulants(p);
	
	for(auto i = 0; i < p.ny; i++){
	  for(auto j = 0; j < p.nr; j++){
	    state.nair[i,j] = state.nair[i,j] + p.dt*(state.qair[i,j] + axdiff[i,j]);
	    state.nwat[i,j] = state.nwat[i,j] + p.dt*state.qadv[i,j];
	    state.cair[i,j] = state.nair[i,j]/volumeair[i,j];
	    state.cwat[i,j] = state.nwat[i,j]/volumewat[i,j];
	  }
	}
	for(auto i = 0; i < p.ny; i++){
	  for(auto j = 0; j < p.nr; j++){
	    eqwat[i,j] = p.henry*s.cair[i,j];
	    towater[i,j] = (eqwat[i,j] - s.cwat[i,j])*volumewat[i,j]*p.gamma;
	    s.nair[i,j] = s.nair[i,j] - p.dt*towater[i,j];
	    s.nwat[i,j] = s.nwat[i,j] + p.dt*towater[i,j];
	    s.cair[i,j] = s.nair[i,j]/volumeair[i,j];
	    s.cwat[i,j] = s.nwat[i,j]/volumewat[i,j];
	  }
	}
	*/

      }

      state.time += p.dt;

    }
  }

  void radial_diffusion(ref Derivative dv, ref State st, in ref Parameters p){
    for(auto i = 0; i < p.ny; i++){
      for(auto j = 0; j < p.nr; j++){
	if(j == 0){
	  // innermost radial element
	  dv.qair[i,j] = -2.0*pii*p.diff_r*(st.cair[i,j] - st.cair[i,j+1]) / ker2[i,j];
	}else if(j == (nr - 1)){
	  // outermost radial element
	  dv.qair[i,j] = ((2.0*pii*p.diff_r*(st.cair[i,j-1] - st.cair[i,j]) / ker1[i,j]) -
			  (2.0*pii*p.diff_b*(st.cair[i,j] - p.camb) / ker2[i,j]));
	}else{
	  dv.qair[i,j] = ((2.0*pii*p.diff_r*(st.cair[i,j-1] - st.cair[i,j]) / ker1[i,j]) -
			  (2.0*pii*p.diff_r*(st.cair[i,j] - st.cair[i,j+1]) / ker2[i,j]));
	}
	// set radial fluxdiff from side
	st.fluxside[i] = 0.0;
	if(j == (nr - 1)){
	  st.fluxside[i] = 2.0*pii*p.diff_b*(st.cair[i,j] - p.camb) / ker2[i,j];
	}
      }
    }
  }

  void axial_advection(ref Derivative dv, ref State st,
		       in ref Parameters p, in double ax_vel){
    for(auto i = 0; i < p.ny; i++){
      for(auto j = 0; j < p.nr; j++){
	// modified: multiply every term with crossect
	if(j < p.nrs){
	  // zero inside the duramen
	  dv.qadv[i,j] = 0.0;
	}else{
	  if(i == 0){
	    // bottom element layer receives CH4 in water
	    dv.qadv[i,j] = ((ax_vel * p.fractwat * p.csoil * crossect[i,j]) -
			    (ax_vel * st.cwat[i,j] * crossect[i,j]));
	    st.fluxin[j] = ax_vel * p.fractwat * p.csoil * crossect[i,j];
	  }else if(i == (ny - 1)){
	    // the topmost elements simply lose CH4 in water
	    dv.qadv[i,j] = ((ax_vel*st.cwat[i-1,j]*crossect[i-1,j]) -
			    (ax_vel*st.cwat[i,j]*crossect[i,j]));
	  }else{
	    dv.qadv[i,j] = ((ax_vel*st.cwat[i-1,j]*crossect[i-1,j]) -
			    (ax_vel*st.cwat[i,j]*crossect[i,j]));
	  }
	}
	// set fluxadv
	if(j < p.nrs){
	  st.fluxtop_a[j] = 0.0;
	}else{
	  st.fluxtop_a[j] = ax_vel * st.cwat[ny-1,j] * crossect[i,j];
	}
      }
    }
  }

  void axial_diffusion(ref Derivative dv, ref State st, in ref Parameters p){
    for(auto i = 0; i < p.ny; i++){
      for(auto j = 0; j < p.nr; j++){
	if(i == 0){
	  // without axial diffusion in
	  dv.axdiff[i,j] = -(p.diff_a * crossect[i,j] *
			     (st.cair[i,j] - st.cair[i+1,j]) / p.dy);
	  // with axial diffusion in
	  /*
	    axdiff[i,j] = (p.diff_a*crossect[j]*(p.csoil - s.cair[i,j]) / p.dy) -
	    (p.diff_a*crossect[j]*(s.cair[i,j] - s.cair[i+1,j]) / p.dy);
	  */
	}else if(i == (ny - 1)){
	  dv.axdiff[i,j] = ((p.diff_a * crossect[i,j] *
			     (st.cair[i-1,j] - st.cair[i,j]) / p.dy) -
			    (p.diff_a * crossect[i,j] *
			     (st.cair[i,j] - p.camb) / p.dy));
	  st.fluxtop_d[j] = (p.diff_a * crossect[i,j] * (st.cair[i,j] - p.camb) / p.dy);
	}else{
	  dv.axdiff[i,j] = ((p.diff_a * crossect[i,j] *
			     (st.cair[i-1,j] - st.cair[i,j]) / p.dy) -
			    (p.diff_a * crossect[i,j] *
			     (st.cair[i,j] - st.cair[i+1,j]) / p.dy));
	}
      }
    }
  }

  void phase_equilibration(ref Derivative dv, ref State st, in ref Parameters p){
    for(auto i = 0; i < p.ny; i++){
      for(auto j = 0; j < p.nr; j++){
	dv.towater[i,j] = (p.henry*st.cair[i,j] - st.cwat[i,j])*volumewat[i,j]*p.gamma;
      }
    }
  }

  void update_cumulants(in ref Parameters p){
    fluxin_c += p.dt * fluxin_sum();
    fluxside_c += p.dt * fluxside_sum();
    fluxtop_d_c += p.dt * fluxtop_d_sum();
    fluxtop_a_c += p.dt * fluxtop_a_sum();
  }

  pure double storage_sum(){
    double sum = 0.0;
    for(auto i = 0; i < ny; i++){
      for(auto j = 0; j < nr; j++){
	sum += state.nwat[i,j] + state.nair[i,j];
      }
    }
    return sum;
  }

  pure double fluxin_sum(){
    return sum(state.fluxin);
  }

  pure double fluxside_sum(){
    return sum(state.fluxside);
  }

  pure double fluxtop_a_sum(){
    return sum(state.fluxtop_a);
  }

  pure double fluxtop_d_sum(){
    return sum(state.fluxtop_d);
  }

  /**
   *  fills the 4 value array used for checking equilibrium
   */
  void fill_summary(ref double[4] s){
    // previously: s[0] = storage_sum();
    s[0] = fluxin_sum();
    s[1] = fluxside_sum();
    s[2] = fluxtop_d_sum();
    s[3] = fluxtop_a_sum();
  }
  
};

unittest
{
  Parameters p;
  p.set_defaults();
  Tree t = new Tree(p);
  assert(t.state.time < 0.1);
}
