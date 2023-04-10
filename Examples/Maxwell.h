#pragma once
#include "Definitions.h"

void Maxwell_Org_Evolution(Field* u, Field* du,int I,SysParams* sys, dfunc(DX), dfunc(DY), dfunc(DZ)) {
	du[0+I] = ( sys->im * DY(u,0+3-1+I,sys->dy,sys->NF,sys->NX,sys->NY,sys->NZ) - sys->im * DZ(u,0+2-1+I,sys->dz,sys->NF,sys->NX,sys->NY,sys->NZ));
	du[1+I] = ( - sys->im * DX(u,0+3-1+I,sys->dx,sys->NF,sys->NX,sys->NY,sys->NZ) + sys->im * DZ(u,0+1-1+I,sys->dz,sys->NF,sys->NX,sys->NY,sys->NZ));
	du[2+I] = ( sys->im * DX(u,0+2-1+I,sys->dx,sys->NF,sys->NX,sys->NY,sys->NZ) - sys->im * DY(u,0+1-1+I,sys->dy,sys->NF,sys->NX,sys->NY,sys->NZ));
}

void Maxwell_Mod_Evolution(Field* u, Field* du,int I,SysParams* sys, dfunc(DX), dfunc(DY), dfunc(DZ)) {
	float b = 0.0f;
	float a = 1.0f;
	du[0+I] = ( sys->im * DY(u,0+3-1+I,sys->dy,sys->NF,sys->NX,sys->NY,sys->NZ) - sys->im * DZ(u,0+2-1+I,sys->dz,sys->NF,sys->NX,sys->NY,sys->NZ)) + a * DX(u,4-1+I,sys->dx,sys->NF,sys->NX,sys->NY,sys->NZ);
	du[1+I] = ( - sys->im * DX(u,0+3-1+I,sys->dx,sys->NF,sys->NX,sys->NY,sys->NZ) + sys->im * DZ(u,0+1-1+I,sys->dz,sys->NF,sys->NX,sys->NY,sys->NZ)) + a * DY(u,4-1+I,sys->dy,sys->NF,sys->NX,sys->NY,sys->NZ);
	du[2+I] = ( sys->im * DX(u,0+2-1+I,sys->dx,sys->NF,sys->NX,sys->NY,sys->NZ) - sys->im * DY(u,0+1-1+I,sys->dy,sys->NF,sys->NX,sys->NY,sys->NZ)) + a * DZ(u,4-1+I,sys->dz,sys->NF,sys->NX,sys->NY,sys->NZ);
	du[3+I] = (1/a)*(DX(u,I,sys->dx,sys->NF,sys->NX,sys->NY,sys->NZ) + DY(u,1+I,sys->dy,sys->NF,sys->NX,sys->NY,sys->NZ) + DZ(u,2+I,sys->dz,sys->NF,sys->NX,sys->NY,sys->NZ)) + b * u[3+I];	
}

Field Maxwell_Org_Constraint(Field* u, int I, SysParams* sys){
	return sys->DX(u,0+I,sys->dx,sys->NF,sys->NX,sys->NY,sys->NZ) + sys->DY(u,1+I,sys->dy,sys->NF,sys->NX,sys->NY,sys->NZ) + sys->DZ(u,2+I,sys->dz,sys->NF,sys->NX,sys->NY,sys->NZ);
}

SysParams Init_MOPE(Field*& u, int NX, int NY, int NZ) {
	SysParams sys;
	sys.NF = 3;
	sys.NX = NX;
	sys.NY = NY;
	sys.NZ = NZ;
	u = new Field[sys.NF*NX*NY*NZ];
	// float s3 = 1.0f;
	float s3 = std::sqrt(3);
	float L = 4*pi*s3;
	sys.X = new float[NX];
	sys.Y = new float[NY];
	sys.Z = new float[NZ];
	sys.dx = L/NX;
	sys.dy = L/NY;
	sys.dz = L/NZ;

	sys.dt = sys.cfl*std::min(sys.dx,std::min(sys.dy,sys.dz));
	int I = 0;

	sys.tmin = 0.0f;
	sys.tmax = 100.0f;

	float t_temp = sys.tmin;
	while (t_temp < sys.tmax) {
		t_temp+=sys.dt;
		t_temp+=sys.dt;
		sys.NT += 2;
	}
	for(int i = 0; i < NX; i++) for(int j = 0; j < NY; j++) for(int k = 0; k < NZ; k++) {
		I = i*sys.NF + j*sys.NF*NX + k*sys.NF*NX*NY;

		float x = -L/2 + i*sys.dx;
		float y = -L/2 + j*sys.dy;
		float z = -L/2 + k*sys.dz;
		sys.X[i] = x;
		sys.Y[j] = y;
		sys.Z[k] = z;

		for(int n = 0; n < sys.NF; n++) {
			u[n + I] = 0.0;
		}

		// kx = ky = kz = 1/sqrt(3);
		u[0+I] = std::exp(sys.im*(x+y+z)/s3);
		u[1+I] = -(std::complex<float>(1,0)+sys.im*s3)*float(0.5) * std::exp(sys.im*(x+y+z)/s3);
		u[2+I] = (-std::complex<float>(1,0)+sys.im*s3)*float(0.5) * std::exp(sys.im*(x+y+z)/s3);

		// kx = 1, ky = 0, kz = 0;
		// u[0+I] = 0.0f;
		// u[1+I] = std::exp(sys.im*x);
		// u[2+I] = -sys.im * std::exp(sys.im*x);

		
		// float r = std::sqrt(x*x+y*y+z*z);
		// u[0+I] = std::exp(-r*r);
	}
	sys.DX = DerX;
	sys.DY = DerY;
	sys.DZ = DerZ;

	sys.DX0 = DerX_0_Periodic;
	sys.DY0 = DerY_0_Periodic;
	sys.DZ0 = DerZ_0_Periodic;

	sys.DXN = DerX_N_Periodic;
	sys.DYN = DerY_N_Periodic;
	sys.DZN = DerZ_N_Periodic;
	
	sys.DISS_DX = KO4_Diss_DX;
	sys.DISS_DY = KO4_Diss_DY;
	sys.DISS_DZ = KO4_Diss_DZ;

	sys.Step = Maxwell_Org_Evolution;
	sys.Constraint = Maxwell_Org_Constraint;
	sys.evolve_exterior = true;
	
	return sys;	
}

SysParams Init_MORE(Field*& u, int NX, int NY, int NZ) {
	SysParams sys;
	sys.NF = 3;
	sys.NX = NX;
	sys.NY = NY;
	sys.NZ = NZ;
	u = new Field[sys.NF*NX*NY*NZ];
	float L = 5.0f;
	sys.X = new float[NX];
	sys.Y = new float[NY];
	sys.Z = new float[NZ];
	sys.dx = L/NX;
	sys.dy = L/NY;
	sys.dz = L/NZ;

	sys.dt = sys.cfl*std::min(sys.dx,std::min(sys.dy,sys.dz));
	int I = 0;

	sys.tmin = 0.0f;
	sys.tmax = 1000.0f;

	float t_temp = sys.tmin;
	while (t_temp < sys.tmax) {
		t_temp+=sys.dt;
		t_temp+=sys.dt;
		sys.NT += 2;
	}
	for(int i = 0; i < NX; i++) for(int j = 0; j < NY; j++) for(int k = 0; k < NZ; k++) {
		I = i*sys.NF + j*sys.NF*NX + k*sys.NF*NX*NY;

		float x = -L/2 + i*sys.dx;
		float y = -L/2 + j*sys.dy;
		float z = -L/2 + k*sys.dz;
		sys.X[i] = x;
		sys.Y[j] = y;
		sys.Z[k] = z;
		float r = std::sqrt(x*x+y*y+z*z);

		for(int n = 0; n < sys.NF; n++) {
			u[n + I] = 0.0;
		}

		float r3 = (r*r*r + 0.001f);
		u[0 + I] = x/r3;
		u[1 + I] = y/r3;
		u[2 + I] = z/r3;
	}
	sys.DX = DerX;
	sys.DY = DerY;
	sys.DZ = DerZ;

	sys.DX0 = DerX_0_Evolve;
	sys.DY0 = DerY_0_Evolve;
	sys.DZ0 = DerZ_0_Evolve;

	sys.DXN = DerX_N_Evolve;
	sys.DYN = DerY_N_Evolve;
	sys.DZN = DerZ_N_Evolve;

	sys.DISS_DX = KO4_Diss_DX;
	sys.DISS_DY = KO4_Diss_DY;
	sys.DISS_DZ = KO4_Diss_DZ;

	sys.Step = Maxwell_Org_Evolution;
	sys.Constraint = Maxwell_Org_Constraint;
	sys.evolve_exterior = false;
	return sys;
}

SysParams Init_MMRE(Field*& u, int NX, int NY, int NZ) {
	SysParams sys;
	sys.NF = 4;
	sys.NX = NX;
	sys.NY = NY;
	sys.NZ = NZ;
	u = new Field[sys.NF*NX*NY*NZ];
	float L = 5.0f;
	sys.X = new float[NX];
	sys.Y = new float[NY];
	sys.Z = new float[NZ];
	sys.dx = L/NX;
	sys.dy = L/NY;
	sys.dz = L/NZ;

	sys.dt = sys.cfl*std::min(sys.dx,std::min(sys.dy,sys.dz));
	int I = 0;

	sys.tmin = 0.0f;
	sys.tmax = 10.0f;

	float t_temp = sys.tmin;
	while (t_temp < sys.tmax) {
		t_temp+=sys.dt;
		t_temp+=sys.dt;
		sys.NT += 2;
	}
	for(int i = 0; i < NX; i++) for(int j = 0; j < NY; j++) for(int k = 0; k < NZ; k++) {
		I = i*sys.NF + j*sys.NF*NX + k*sys.NF*NX*NY;

		float x = -L/2 + i*sys.dx;
		float y = -L/2 + j*sys.dy;
		float z = -L/2 + k*sys.dz;
		sys.X[i] = x;
		sys.Y[j] = y;
		sys.Z[k] = z;
		float r = std::sqrt(x*x+y*y+z*z);

		for(int n = 0; n < sys.NF; n++) {
			u[n + I] = 0.0;
		}

		float r3 = (r*r*r + 0.001f);
		u[0 + I] = x/r3;
		u[1 + I] = y/r3;
		u[2 + I] = z/r3;
	}
	sys.DX = DerX;
	sys.DY = DerY;
	sys.DZ = DerZ;

	sys.DX0 = DerX_0_Evolve;
	sys.DY0 = DerY_0_Evolve;
	sys.DZ0 = DerZ_0_Evolve;

	sys.DXN = DerX_N_Evolve;
	sys.DYN = DerY_N_Evolve;
	sys.DZN = DerZ_N_Evolve;

	sys.DISS_DX = KO4_Diss_DX;
	sys.DISS_DY = KO4_Diss_DY;
	sys.DISS_DZ = KO4_Diss_DZ;

	sys.Step = Maxwell_Mod_Evolution;
	sys.Constraint = Maxwell_Org_Constraint;
	sys.evolve_exterior = false;
	return sys;
}

SysParams Init_MMPE(Field*& u, int NX, int NY, int NZ) {
	SysParams sys;
	sys.NF = 3;
	sys.NX = NX;
	sys.NY = NY;
	sys.NZ = NZ;
	u = new Field[sys.NF*NX*NY*NZ];
	// float s3 = 1.0f;
	float s3 = std::sqrt(3);
	float L = 4*pi*s3;
	sys.X = new float[NX];
	sys.Y = new float[NY];
	sys.Z = new float[NZ];
	sys.dx = L/NX;
	sys.dy = L/NY;
	sys.dz = L/NZ;

	sys.dt = sys.cfl*std::min(sys.dx,std::min(sys.dy,sys.dz));
	int I = 0;

	sys.tmin = 0.0f;
	sys.tmax = 100.0f;

	float t_temp = sys.tmin;
	while (t_temp < sys.tmax) {
		t_temp+=sys.dt;
		t_temp+=sys.dt;
		sys.NT += 2;
	}
	for(int i = 0; i < NX; i++) for(int j = 0; j < NY; j++) for(int k = 0; k < NZ; k++) {
		I = i*sys.NF + j*sys.NF*NX + k*sys.NF*NX*NY;

		float x = -L/2 + i*sys.dx;
		float y = -L/2 + j*sys.dy;
		float z = -L/2 + k*sys.dz;
		sys.X[i] = x;
		sys.Y[j] = y;
		sys.Z[k] = z;

		for(int n = 0; n < sys.NF; n++) {
			u[n + I] = 0.0;
		}

		// kx = ky = kz = 1/sqrt(3);
		u[0+I] = std::exp(sys.im*(x+y+z)/s3);
		u[1+I] = -(std::complex<float>(1,0)+sys.im*s3)*float(0.5) * std::exp(sys.im*(x+y+z)/s3);
		u[2+I] = (-std::complex<float>(1,0)+sys.im*s3)*float(0.5) * std::exp(sys.im*(x+y+z)/s3);

		// kx = 1, ky = 0, kz = 0;
		// u[0+I] = 0.0f;
		// u[1+I] = std::exp(sys.im*x);
		// u[2+I] = -sys.im * std::exp(sys.im*x);

		
		// float r = std::sqrt(x*x+y*y+z*z);
		// u[0+I] = std::exp(-r*r);
	}
	sys.DX = DerX;
	sys.DY = DerY;
	sys.DZ = DerZ;

	sys.DX0 = DerX_0_Periodic;
	sys.DY0 = DerY_0_Periodic;
	sys.DZ0 = DerZ_0_Periodic;

	sys.DXN = DerX_N_Periodic;
	sys.DYN = DerY_N_Periodic;
	sys.DZN = DerZ_N_Periodic;
	
	sys.DISS_DX = KO4_Diss_DX;
	sys.DISS_DY = KO4_Diss_DY;
	sys.DISS_DZ = KO4_Diss_DZ;

	sys.Step = Maxwell_Org_Evolution;
	sys.Constraint = Maxwell_Org_Constraint;
	sys.evolve_exterior = true;
	
	return sys;	
}
