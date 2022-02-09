#include "pdf.h"
#include <cmath>

using namespace std;

double hydro_100::eval(Position *pos){
	return (1./M_PI)*exp(-2.*(pos->get_radius()));
}

double hydro_210::eval(Position *pos){
	return 1./64.*(sqrt(2./M_PI))*(pos->get_radius())*exp(-(pos->get_radius())/2.)*cos(pos->get_theta());
}
