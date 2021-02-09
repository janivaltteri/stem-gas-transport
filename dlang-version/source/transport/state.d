module state;

import array2d;

struct State
{
public:

  double time;

  double[] fluxin;
  double[] fluxside;
  double[] fluxtop_d;
  double[] fluxtop_a;
  
  Matrix cair;
  Matrix cwat;
  Matrix nair;
  Matrix nwat;

  this(in int in_ny, in int in_nr){

    this.time = 0.0;

    this.fluxin.length = in_nr;
    this.fluxside.length = in_ny;
    this.fluxtop_d.length = in_nr;
    this.fluxtop_a.length = in_nr;
    
    this.cair.initialise(in_ny,in_nr);
    this.cwat.initialise(in_ny,in_nr);
    this.nair.initialise(in_ny,in_nr);
    this.nwat.initialise(in_ny,in_nr);

    foreach(ref e; fluxin) e = 0.0;
    foreach(ref e; fluxside) e = 0.0;
    foreach(ref e; fluxtop_d) e = 0.0;
    foreach(ref e; fluxtop_a) e = 0.0;

  }

}

unittest
{
  State s = State(5,5);

  assert(s.time < 1.0);
  assert(s.cair[0,0] < 0.1);
  assert(s.cair[0,0] > -0.1);
  assert(s.fluxside[0] < 0.1);
  assert(s.fluxside[0] > -0.1);
  assert(s.fluxtop_d[0] < 0.1);
  assert(s.fluxtop_d[0] > -0.1);
  assert(s.fluxtop_a[0] < 0.1);
  assert(s.fluxtop_a[0] > -0.1);
}
