#pragma once
#define GNUPLOT_ENABLE_PTY
#include "gnuplot-iostream.h"
#include <iostream>
#include <complex>
#include <chrono>
#include <map>
#include <vector>
#include <cmath>
#include <stdio.h>
#include <boost/tuple/tuple.hpp>

#define pi 3.1415962
#define Field std::complex<float>

#define dfunc(W) Field (*W)(Field*,int,float,int,int,int,int)

#define stepf(X) void (*X)(Field,Field*,Field*,float,float,float,float,int,int,int,int,int, dfunc(),dfunc(),dfunc())

#define stepf2(X) void (*X)(Field*,Field*,int,SysParams*,dfunc(),dfunc(),dfunc());

#define stepC(X) Field (*X)(Field*,int,float,float,float,int,int,int,int)

#define stepC2(X) Field (*X)(Field*,int,SysParams*)


class SysParams {
	public:
	int NF;
	int NX;
	int NY;
	int NZ;

	int NT;

	float dt;
	float dx;
	float dy;
	float dz;

	float* X;
	float* Y;
	float* Z;

	float tmin;
	float tmax;

	Field im = std::complex<float>(0,1);

	float cfl = 0.1f;

	float diss = -0.1f;

	stepf2(Step);

	stepC2(Constraint);	

	dfunc(DX);
	dfunc(DY);
	dfunc(DZ);

	dfunc(DX0);
	dfunc(DY0);
	dfunc(DZ0);

	dfunc(DXN);
	dfunc(DYN);
	dfunc(DZN);

	dfunc(DISS_DX);
	dfunc(DISS_DY);
	dfunc(DISS_DZ);

	// Plotting Variables

	Gnuplot* gp;

	bool enable_plots = true;

	bool evolve_exterior = false;

	SysParams() {
			NF = 0;
			NX = 0;
			NY = 0;
			NZ = 0;
			NT = 0;
			dt = 0.0f;
			dx = 0.0f;
			dy = 0.0f;
			dz = 0.0f;

			tmin = 0.0f;
			tmax = 0.0f;
	}
};

#include "Derivatives.h"
#include "Evolutions.h"
