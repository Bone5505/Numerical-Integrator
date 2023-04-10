#pragma once
#include "Definitions.h"

Field DerX(Field* u, int I, float dx, int NF, int NX,int NY, int NZ) {
	return (u[I+NF]-u[I-NF])/(2*dx);
}

Field DerY(Field* u, int I, float dy, int NF, int NX, int NY, int NZ) {
	return (u[I+NF*NX] - u[I-NF*NX])/(2*dy);
}

Field DerZ(Field* u, int I, float dz, int NF, int NX, int NY, int NZ) {
	return (u[I+NF*NY*NZ]-u[I-NF*NX*NY])/(2*dz);
}

// Periodic Boundary Conditions. 


Field DerX_0_Periodic(Field* u, int I, float dx, int NF, int NX,int NY, int NZ) {
	return (u[I+NF]-u[I+(NX-1)*NF])/(2*dx);
}

Field DerY_0_Periodic(Field* u, int I, float dy, int NF, int NX, int NY, int NZ) {
	return (u[I+NF*NX]-u[I+(NY-1)*NF*NX])/(2*dy);
}

Field DerZ_0_Periodic(Field* u, int I, float dz, int NF, int NX, int NY, int NZ) {
	return (u[I+NF*NX*NY]-u[I+NF*NX*NY*(NZ-1)])/(2*dz);
}

Field DerX_N_Periodic(Field* u, int I, float dx, int NF, int NX,int NY, int NZ) {
	return (u[I-(NX-1)*NF]-u[I-NF])/(2*dx);
}

Field DerY_N_Periodic(Field* u, int I, float dy, int NF, int NX, int NY, int NZ) {
	return (u[I-(NY-1)*NF*NX]-u[I-NF*NX])/(2*dy);
}

Field DerZ_N_Periodic(Field* u, int I, float dz, int NF, int NX, int NY, int NZ) {
	return (u[I-(NZ-1)*NF*NX*NY]-u[I-NF*NX*NY])/(2*dz);
}


// Evolve the ends using up/down-wind scheme


Field DerX_0_Evolve(Field* u, int I, float dx, int NF, int NX,int NY, int NZ) {
	return (-3.0f*u[I] +4.0f*u[I+NF]-u[I+2*NF])/(2*dx);
	// return (-u[I]+u[I+NF])/dx;
}

Field DerY_0_Evolve(Field* u, int I, float dy, int NF, int NX, int NY, int NZ) {
	return (-3.0f*u[I] +4.0f*u[I+NF*NX]-u[I+2*NF*NX])/(2*dy);
	// return (-u[I]+u[I+NF*NX])/dy;
}

Field DerZ_0_Evolve(Field* u, int I, float dz, int NF, int NX, int NY, int NZ) {
	return (-3.0f*u[I] +4.0f*u[I+NF*NX*NY]-u[I+2*NF*NX*NY])/(2*dz);
	// return (-u[I]+u[I+NF*NX*NY])/dz;
}


Field DerX_N_Evolve(Field* u, int I, float dx, int NF, int NX,int NY, int NZ) {
	return (u[I-2*NF] - 4.0f*u[I-NF]+3.0f*u[I])/(2*dx);
	// return (-u[I-NF]+u[I])/dx;
}

Field DerY_N_Evolve(Field* u, int I, float dy, int NF, int NX, int NY, int NZ) {
	return (u[I-2*NF*NX] - 4.0f*u[I-NF*NX]+3.0f*u[I])/(2*dy);
	// return (-u[I-NF*NX]+u[I])/dy;
}

Field DerZ_N_Evolve(Field* u, int I, float dz, int NF, int NX, int NY, int NZ) {
	return (u[I-2*NF*NX*NY] - 4.0f*u[I-NF*NX*NY]+3.0f*u[I])/(2*dz);
	// return (-u[I-NF*NX*NY]+u[I])/dz;
}


Field KO4_Diss_DX(Field* u,int I, float dx, int NF, int NX, int NY, int NZ) {
	float dx2 = dx*dx;
	return (u[I+2*NF]-4.0f*u[I+NF]+6.0f*u[I]-4.0f*u[I-NF]+u[I-2*NF])/(dx2*dx2);
}

Field KO4_Diss_DY(Field* u,int I, float dy, int NF, int NX, int NY, int NZ) {
	float dy2 = dy*dy;
	return (u[I+2*NF*NX]-4.0f*u[I+NF*NX]+6.0f*u[I]-4.0f*u[I-NF*NX]+u[I-2*NF*NX])/(dy2*dy2);
}

Field KO4_Diss_DZ(Field* u,int I, float dz, int NF, int NX, int NY, int NZ) {
	float dz2 = dz*dz;
	return (u[I+2*NF*NX*NY]-4.0f*u[I+NF*NX*NY]+6.0f*u[I]-4.0f*u[I-NF*NX*NY]+u[I-2*NF*NX*NY])/(dz2*dz2);
}
