#ifndef _position_h_
#define _position_h_

#include <vector>

class Position{
  private:
    double _r, _theta, _phi;
  public:
    Position();
    ~Position();
    Position get_position();
    void set_position(Position);
    void to_zero();
};

#endif
