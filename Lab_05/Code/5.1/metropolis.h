#ifndef _metropolis_h_
#define _metropolis_h_

#include "random.h"
#include "position.h"
#include "pdf.h"

using namespace std;

class Metropolis : public Measure, public pdf{
  private:
    int _steps, _accepted;                            // number of steps, counter of accepted steps
    pdf *p_function;                              // density function
    Position *_xold;                              // x_old
    Position *_x;                                 // x_new
    Random _rnd;

  public:
    Metropolis(int N, pdf *mypdf, Position *start, Random rnd);
    ~Metropolis();
    double T_uni();
    //double T_gauss();
    double A_eval();
    void make_step();
    void walker_at_home();
};

#endif
