#ifndef _measure_h_
#define _measure_h_

#include <vector>
using namespace std;

// classe astratta Measure
// in ogni problema in cui la usiamo in `blocking` bisogna concretizzarla
// serve implementare un metodo di misura vera e propria che fornisca una per una
// le componenti della misura sotto forma di vector<double>
// e un metodo che fornisca il numero delle misure effettuate

class Measure{
  public:
    virtual vector<double> get_measure()=0;
    virtual int get_dimension()=0;
};

#endif
