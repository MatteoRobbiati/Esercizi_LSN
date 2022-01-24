#ifndef __myfunction_h__
#define __myfunction_h__

#include "basefunction.h"

class MyCos : public BaseFunction {
	public:
		MyCos();
		~MyCos();
		double eval(double x);
};


class MyGeneral : public BaseFuction {
	public:
		MyGeneral();
		~MyGeneral();
		double eval(double x);
}

#endif
