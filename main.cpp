#include "Definitions.h"
#include "Maxwell.h"

int main() {
	int NP = 32;
	Field* u;
	Field* u_new;
	// SysParams sys =  Init_MORE(u,NP,NP,NP);
	SysParams sys =  Init_MMPE(u,NP,NP,NP);
	// SysParams sys = Init_MMRE(u,NP,NP,NP);
	Gnuplot gp;
	sys.gp = &gp;


	printf("dx: %f, dy: %f, dz: %f\n",sys.dx,sys.dy,sys.dz);
	Plot2D(u,1,1,sys.NY/2,&sys);
	
	Integrate(u,u_new,&sys);

	printf("Finished!\n");
}
