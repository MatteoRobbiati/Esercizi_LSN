#include "pdf.h"
#include <cmath>

using namespace std;

double hydro_100::eval(Position *pos){
	return (1./M_PI)*exp(-2.*(pos->get_radius()));
}

double hydro_210::eval(Position *pos){
	vector<double> coord = pos->get_coordinates();
	return (1./(32.*M_PI))*(coord.at(2)*coord.at(2))*exp(-(pos->get_radius()));
}
