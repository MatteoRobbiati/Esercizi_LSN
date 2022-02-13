#ifndef _position_h_
#define _position_h_

#include <vector>
#include <iostream>
#include "random.h"

class Position{
  private:
    double _x, _y, _z;
    Random* _rnd;
  public:
    Position(Random* rnd);
    Position(double x, double y, double z, Random* rnd);
    ~Position();
    vector<double> get_coordinates();
    void set_coordinates(vector<double> coord);
    void to_zero();
    double get_theta();
    double get_radius();
    void uniform_step(double stepsize);
    void gaussian_step(double stepsize);
    void show_position();
};

#endif
