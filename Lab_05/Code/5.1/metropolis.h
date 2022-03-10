#ifndef _metropolis_h_
#define _metropolis_h_

#include "random.h"
#include "position.h"
#include "pdf.h"
#include "measure.h"
#include <string>
#include <vector>

using namespace std;

class Metropolis : public Measure{
  private:
    int _steps, _accepted;                        // number of steps, counter of accepted steps
    double _stepsize;
    pdf *p_function;                              // density function
    Position *_x;                                 // x_new
    Random* _rnd;
    int _attempted;
    string _method;

  public:
    Metropolis(int N, pdf *mypdf, Position *start, Random* rnd, double stepsize, string method);
    ~Metropolis();

    void Equilibrate(int Nequi);
    void try_step();
    void run(string filename);
    double rate_of_acceptance();
    vector<double> get_measure();
    int get_dimension();
};

#endif
