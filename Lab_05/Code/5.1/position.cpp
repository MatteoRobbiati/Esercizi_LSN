#include "position.h"

using namespace std;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ constructor and distructor ~~~~~~~~~~~~~~~~~~

Position::Position(){
  _r = 0.;
  _theta = 0.;
  _phi = 0.;
}

Position::~Position();

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ playing with coordinates ~~~~~~~~~~~~~~~~~~~~~

Position Position::get_position(){
  Position *temPos = new Position;
  temPos->_r =_r;
  temPos->_theta = _theta;
  temPos->_phi = _phi;
  return tempPos;
}

void Position::set_position(Position *Pos){
  _r = Pos->_r;
  _theta = Pos->_theta;
  _phi = Pos->_phi;
  return;
}


void Position::to_zero(){
  _r = 0.;
  _theta = 0.;
  _phi = 0.;
}
