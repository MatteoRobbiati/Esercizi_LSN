#include "functions.h"
#include <cmath>

using namespace std;

double MyCos::eval(double x){
	return (0.5*M_PI)*cos(0.5*x*M_PI);
}

double MyGeneral::eval(double x){
	return ((0.5*M_PI)*cos(0.5*x*M_PI))/(2*(1-x));
}
