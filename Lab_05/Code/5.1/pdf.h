#ifndef _pdf_h_
#define _pdf_h_

#include "position.h"

#include <cmath>

//double _a0 = 0.0529e-9;

class pdf{
  public:
    virtual double eval(Position *pos)=0;
};


class hydro_100 : public pdf {
  public:
		double eval(Position *pos);
};


class hydro_210 : public pdf {
  public:
    double eval(Position *pos);
};

#endif
