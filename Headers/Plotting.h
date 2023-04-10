#pragma once
#include "Definitions.h"

void Plot1D(Field* u,int comp, int dir, int pos1, int  pos2,SysParams* sys) {	
	
	std::vector<std::pair<float,float>> data;
	if(dir == 0) {
		for(int i = 0; i < sys->NX; i++) {
			data.push_back(std::pair<float,float>(sys->X[i], u[comp + i*sys->NF + pos1*sys->NF*sys->NX + pos2*sys->NF*sys->NX*sys->NY].real()));
		}
	} else if(dir == 1) {
		for(int i = 0; i < sys->NY; i++) {
			data.push_back(std::pair<float,float>(sys->Y[i], u[comp + pos1*sys->NF + i*sys->NF*sys->NX + pos2*sys->NF*sys->NX*sys->NY].real()));
		}
	} else if(dir == 2) {
		for(int i = 0; i < sys->NZ; i++) {
			data.push_back(std::pair<float,float>(sys->Z[i],  u[comp + pos1*sys->NF + pos2*sys->NF*sys->NX + i*sys->NF*sys->NX*sys->NY].real()));
		}
	} else {
		return;
	}
	*(sys->gp) << "set yrange [-2:2]\n";
	*(sys->gp) << "plot '-' with lines notitle\n";
	*(sys->gp)->send1d(data);
		
}

void Plot2D(Field* u, int comp, int dir, int pos,SysParams* sys) {
	if (dir == 0) {
		std::vector<std::vector<std::tuple<float,float,float>>> data(sys->NY);
		int I;
		for(int j = 0; j < sys->NY; j++) {
			data[j].resize(sys->NZ);
			for(int k = 0; k < sys->NZ; k++) {
				I = comp + pos*sys->NF + j*sys->NF*sys->NX + k*sys->NF*sys->NX*sys->NY;	
				data[j][k] = std::make_tuple(sys->Y[j],sys->Z[k],u[I].real());
			}
		}
		*(sys->gp) << "splot '-' with pm3d notitle\n";
		*(sys->gp)->send2d(data);
	} else if(dir == 1) {
		std::vector<std::vector<std::tuple<float,float,float>>> data(sys->NX);
		int I;
		for(int i = 0; i < sys->NX; i++) {
		data[i].resize(sys->NZ);
			for(int k = 0; k < sys->NZ; k++) {
				I = comp + i*sys->NF + pos*sys->NF*sys->NX + k*sys->NF*sys->NX*sys->NY;	
				data[i][k] = std::make_tuple(sys->X[i],sys->Z[k],u[I].real());
			}
		}
		*(sys->gp) << "splot '-' with pm3d notitle\n";
		*(sys->gp)->send2d(data);
	} else if (dir == 2) {
		std::vector<std::vector<std::tuple<float,float,float>>> data(sys->NX);
		int I;
		for(int i = 0; i < sys->NX; i ++) {
			data[i].resize(sys->NY);
			for(int j = 0; j < sys->NY; j++) {
				I = comp + i*sys->NF + j*sys->NF*sys->NX + pos*sys->NF*sys->NX*sys->NY;
				data[i][j] = std::make_tuple(sys->X[i],sys->Y[j],u[I].real());
			}
		}
		*(sys->gp) << "splot '-' with pm3d notitle\n";
		*(sys->gp)->send2d(data);
	}
	
}

void PlotConstraint(float* G, SysParams* sys) {
	std::vector<std::pair<float,float>> data(sys->NT);
	float t = sys->tmin;
	for(int i =0; i < sys->NT; i++) {
		t += float(sys->dt);
		data[i] = std::pair<float,float>(t,float(G[i]));

	}
	*(sys->gp) << "plot '-' with lines notitle\n";
	*(sys->gp)->send1d(data);
}
