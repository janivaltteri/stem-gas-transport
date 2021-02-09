module derivative;

import array2d;

struct Derivative
{
public:

  Matrix qair;
  Matrix qadv;
  Matrix axdiff;
  Matrix towater;

  this(in int in_ny, in int in_nr){
    this.qair.initialise(in_ny,in_nr);
    this.qadv.initialise(in_ny,in_nr);
    this.axdiff.initialise(in_ny,in_nr);
    this.towater.initialise(in_ny,in_nr);
  }

  
}
