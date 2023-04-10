#pragma once
#include "Definitions.h"
#include "Plotting.h"

float Constraint(Field* u, SysParams* sys) {
	float G = 0.0f;
	int I;
	for(int i = 1; i < sys->NX-1; i++) for(int j = 1; j < sys->NY-1; j++) for(int k = 1; k < sys->NZ-1; k++) {
		I = i*sys->NF + j*sys->NF*sys->NX + k*sys->NF*sys->NX*sys->NY;
		G += std::abs(std::pow(sys->Constraint(u,I,sys),2));
	}
	return std::log(std::sqrt(G/(sys->NX*sys->NY*sys->NZ)));
}

void Add_Dissapation(Field* u, Field* du, SysParams* sys) {
	int I;
	float d2;
	for(int i = 3; i < sys->NX-2; i++) for(int j = 3; j < sys->NY-2; j++) for(int k = 3; k < sys->NZ-2; k++) {
		for(int n = 0; n < sys->NF; n++) {
			I = n + i*sys->NF + j*sys->NF*sys->NX + k*sys->NF*sys->NX*sys->NY;
			d2 = sys->dx*sys->dx;
			du[I] += sys->diss * (d2*sys->dx) * sys->DISS_DX(u,I,sys->dx,sys->NF,sys->NX,sys->NY,sys->NZ);
			d2 = sys->dy*sys->dy;
			du[I] += sys->diss * (d2*sys->dy) * sys->DISS_DY(u,I,sys->dy,sys->NF,sys->NX,sys->NY,sys->NZ);
			d2 = sys->dz*sys->dz;
			du[I] += sys->diss * (d2*sys->dz)  * sys->DISS_DZ(u,I,sys->dz,sys->NF,sys->NX,sys->NY,sys->NZ);
		}
	}
}

void Step_Ins(Field* u, Field* du, SysParams* sys) {
	int NF = sys->NF;
	int NX = sys->NX;
	int NY = sys->NY;
	int NZ = sys->NZ;
	int I = 0;
	for(int i = 1; i < NX-1; i++) for(int j = 1; j < NY-1; j++) for(int k = 1; k < NZ-1; k++) {
		I = i*NF + j*NX*NF + k*NX*NY*NF;
		sys->Step(u,du,I,sys,sys->DX,sys->DY,sys->DZ);
	}
}

void Step_Ends(Field* u, Field* du,SysParams* sys) {
	int I = 0;
	int NF = sys->NF;
	int NX = sys->NX;
	int NY = sys->NY;
	int NZ = sys->NZ;
	Step_Ins(u,du,sys);

	if(!sys->evolve_exterior) {
		return;
	}

	for(int i = 1; i < NX-1; i++) for(int j = 1; j < NY-1; j++) {
		 I = i*NF + j*NF*NX;
		 sys->Step(u,du,I,sys,sys->DX,sys->DY,sys->DZ0);
		 I = i*NF + j*NF*NX + (NZ-1)*NF*NX*NY;
		 sys->Step(u,du,I,sys,sys->DX,sys->DY,sys->DZN);
	}

	// XZ Face
	for(int i = 1; i < NX-1; i++) for(int k = 1; k < NZ-1; k++) {
		I = i*NF + k*NF*NX*NY;
		sys->Step(u,du,I,sys,sys->DX,sys->DY0,sys->DZ);
		I = i*NF + (NY-1)*NF*NX + k*NF*NX*NY;
		sys->Step(u,du,I,sys,sys->DX,sys->DYN,sys->DZ);
	}
	
	// YZ Face
	for(int j = 1; j < NY-1; j++) for(int k = 1; k < NZ-1; k++) {
		I = j*NF*NX + k*NF*NX*NY;
		sys->Step(u,du,I,sys,sys->DX0,sys->DY,sys->DZ);
	
		I = (NX-1)*NF + j*NF*NX + k*NF*NX*NY;
		sys->Step(u,du,I,sys,sys->DXN,sys->DY,sys->DZ);
	}

	for(int i = 1; i < NX-1; i++) {
		I = i*NF;
		sys->Step(u,du,I,sys,sys->DX,sys->DY0,sys->DZ0);
		I = i*NF + (NZ-1)*NF*NX*NY;
		sys->Step(u,du,I,sys,sys->DX,sys->DY0,sys->DZN);
		I = i*NF + (NY-1)*NF*NX;
		sys->Step(u,du,I,sys,sys->DX,sys->DYN,sys->DZ0);
		I = i*NF + (NY-1)*NF*NX + (NZ-1)*NF*NX*NY;
		sys->Step(u,du,I,sys,sys->DX,sys->DYN,sys->DZN);
	
	}
	
	// Y Line
	for(int j = 1; j < NY-1; j++) {
		I = j*NF*NX;
		sys->Step(u,du,I,sys,sys->DX0,sys->DY,sys->DZ0);
		I = j*NF*NX + (NZ-1)*NF*NX*NY;
		sys->Step(u,du,I,sys,sys->DX0,sys->DY,sys->DZN);
		I = (NX-1)*NF + j*NF*NX;
		sys->Step(u,du,I,sys,sys->DXN,sys->DY,sys->DZ0);
		I = (NX-1)*NF + j*NF*NX + (NZ-1)*NF*NX*NY;
		sys->Step(u,du,I,sys,sys->DXN,sys->DY,sys->DZN);
	}
	
	// Z Lines
	for(int k = 1; k < NZ-1; k++) {
		I = k*NF*NX*NY;
		sys->Step(u,du,I,sys,sys->DX0,sys->DY0,sys->DZ);
		I = (NY-1)*NF*NX + k*NF*NX*NY;
		sys->Step(u,du,I,sys,sys->DX0,sys->DYN,sys->DZ);
		I = (NX-1)*NF + k*NF*NX*NY;
		sys->Step(u,du,I,sys,sys->DXN,sys->DY0,sys->DZ);
		I = (NX-1)*NF + (NY-1)*NF*NX + k*NF*NX*NY;
		sys->Step(u,du,I,sys,sys->DXN,sys->DYN,sys->DZ);
	
	}


		I = 0;
		sys->Step(u,du,I,sys,sys->DX0,sys->DY0,sys->DZ0);
		I = (NX-1)*NF;
		sys->Step(u,du,I,sys,sys->DXN,sys->DY0,sys->DZ0);
		I = (NY-1)*NF*NX;
		sys->Step(u,du,I,sys,sys->DX0,sys->DYN,sys->DZ0);
		I = (NZ-1)*NF*NX*NY;
		sys->Step(u,du,I,sys,sys->DX0,sys->DY0,sys->DZN);
		I = (NX-1)*NF + (NY-1)*NF*NX;
		sys->Step(u,du,I,sys,sys->DXN,sys->DYN,sys->DZ0);
		I = (NY-1)*NF*NX + (NZ-1)*NF*NX*NY;
		sys->Step(u,du,I,sys,sys->DX0,sys->DYN,sys->DZN);
		I = (NX-1)*NF + (NZ-1)*NF*NX*NY;
		sys->Step(u,du,I,sys,sys->DXN,sys->DY0,sys->DZN);
		I = (NX-1)*NF + (NY-1)*NF*NX + (NZ-1)*NF*NX*NY;
		sys->Step(u,du,I,sys,sys->DXN,sys->DYN,sys->DZN);
}

void ICN_Evolve(Field* u, Field* u_new,Field* du, Field* u_half, SysParams* sys) {
	Step_Ends(u,du,sys);
	// Add_Dissapation(u, du,sys);
	int I = 0;
	for(I = 0; I < sys->NF*sys->NX*sys->NY*sys->NZ; I++) {
		u_new[I] = u[I] + sys->dt * du[I];
		u_half[I] = (u[I] + u_new[I])/2.0f;
		du[I] = 0.0f;
	}
	
	Step_Ends(u_half,du,sys);
	// Add_Dissapation(u_half,du,sys);
	for(I = 0; I < sys->NF*sys->NX*sys->NY*sys->NZ; I++) {
		u_new[I] = u[I] + sys->dt * du[I];
		u_half[I] = (u[I] + u_new[I])/2.0f;
		du[I] = 0.0f;
	}

	Step_Ends(u_half,du,sys);
	// Add_Dissapation(u_half,du,sys);
	for(I = 0; I < sys->NF*sys->NX*sys->NY*sys->NZ; I++) {
		u_new[I] = u[I] + sys->dt * du[I];
		du[I] = 0.0f;
	}
}

void Integrate(Field* u, Field* u_new, SysParams* sys) {

	int NN = sys->NF*sys->NX*sys->NY*sys->NZ;
	Field* u_old  = new Field[NN];
	Field* u_half = new Field[NN];
	Field* du     = new Field[NN];
	u_new = new Field[NN];
	for(int i = 0; i < NN; i++) {
			u_old[i] = u[i];
			u_half[i] = 0.0f;
			du[i] = 0.0f;
		}
	float t      = sys->tmin;
	float t_frac = 0.0f;
	int t_count  = 100;

	float* G = new float[sys->NT];
	
	for(int i = 0; i < t_count; i++) {
		std::cout << "#";
	}
	std::cout<< "#" << std::endl;

	int I = 0;
	while(t  < sys->tmax) {
		ICN_Evolve(u,u_new,du,u_half,sys);
		t += sys->dt;
		G[I] = Constraint(u_new,sys);

		ICN_Evolve(u_new,u,du,u_half,sys);
		G[I+1] = Constraint(u_new,sys);
		t += sys->dt;
		I+=2;
		while((t-sys->tmin)/sys->tmax > t_frac) {
			t_frac += 1.0f/float(t_count);
			std::cout << "#" << std::flush;
			if(sys->enable_plots) {
				Plot2D(u,1,1,sys->NY/2,sys);
			}
		}
	}
	std::cout << std::endl;
	PlotConstraint(G,sys);
	std::cin.get();
	for(int i = 0; i < NN; i++) {
		u_new[i] = u[i];
		u[i] = u_old[i];
	}
}
