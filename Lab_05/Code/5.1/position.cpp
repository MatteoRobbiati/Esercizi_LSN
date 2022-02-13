#include "position.h"
#include "random.h"
#include <cmath>

using namespace std;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ constructor and distructor ~~~~~~~~~~~~~~~~~~

Position::Position(Random* rnd){
  _x = 0.;
  _y = 0.;
  _z = 0.;
  _rnd = rnd;
}

Position::Position(double x, double y, double z, Random* rnd){
  _x = x;
  _y = y;
  _z = z;
  _rnd = rnd;
}

Position::~Position(){
  _rnd->SaveSeed();
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ playing with coordinates ~~~~~~~~~~~~~~~~~~~~~

void Position::set_coordinates(vector<double> coord){
  _x = coord.at(0);
  _y = coord.at(1);
  _z = coord.at(2);
  return;
}

vector<double> Position::get_coordinates(){
  vector<double> coor(3,0);
  coor.at(0)=_x;
  coor.at(1)=_y;
  coor.at(2)=_z;
  return coor;
}

void Position::to_zero(){
  _x = 0.;
  _y = 0.;
  _z = 0.;
}

double Position::get_radius(){
  return sqrt(_x*_x+_y*_y+_z*_z);
}

double Position::get_theta(){
  return atan(_x*_x+_y*_y)/_z;
}

void Position::uniform_step(double stepsize){
  double X = _rnd->Rannyu(-stepsize, stepsize);
  double Y = _rnd->Rannyu(-stepsize, stepsize);
  double Z = _rnd->Rannyu(-stepsize, stepsize);
  _x+=X;
  _y+=Y;
  _z+=Z;
}

void Position::gaussian_step(double stepsize){
  double X = _rnd->Gauss(0, stepsize);
  double Y = _rnd->Gauss(0, stepsize);
  double Z = _rnd->Gauss(0, stepsize);
  _x+=X;
  _y+=Y;
  _z+=Z;
}


void Position::show_position(){
  cout << "Position: (" << _x << ", " << _y << ", " << _z << ")" << endl;
  return;
}
