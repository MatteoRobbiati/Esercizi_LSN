#ifndef _functions_h_
#define _functions_h_

#include <cmath>

class BaseFunction{
  public:
    virtual double eval(double x)=0;
};


class MyCos : public BaseFunction {
  public:
		double eval(double x);
};


class MyGeneral : public BaseFunction {
  public:
    double eval(double x);
};

#endif
