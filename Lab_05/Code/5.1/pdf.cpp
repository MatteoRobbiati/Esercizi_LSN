#include "pdf.h"
#include <cmath>

using namespace std;

double Hydro_100::eval(Position *pos){
	return (1./(sqrt(M_PI)))*exp(-(pos->_r));
}

double hydro_210::eval(Position *pos){
	return 1./8.*(sqrt(2./M_PI))*(pos->_r)*exp(-(pos->_r)/2.)*cos(pos->_theta);
}
