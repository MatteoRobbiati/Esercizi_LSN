#ifndef _randomwalk_h_
#define _randomwalk_h_

#include "random.h"
#include "measure.h"
#include "statistic.h"
#include <vector>
#include <string>
#include <cmath>

using namespace std;

class RandomWalk : public Random, public Measure{
  private:
    double _d;                      // lattice constant
    int _N, _dim;                   // number of steps and dimension of the hypervolume
                                    // _dim = 3 in our (x,y,z) case
    vector<double> _pos;            // position of the walker
    Random _rnd;
    string _metric;                 // could be discrete or continuum

  public:
    RandomWalk(double d, int N, int dim, string metric, Random rnd);
    ~RandomWalk();
    vector<double> get_position();                        // get coordinates
    double get_squared_mod();                             // get squared mod by position
    void walker_at_home();                                // bring walker to the origin
    void make_discrete_step();                            // one discrete step
    void make_continuum_step();                           // one continuum step
    void run_walk(string filename);                       // one discrete random walk
    vector<double> get_measure();
    int get_dimension();
};


#endif
