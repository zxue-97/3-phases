
#include "constants.h"

struct CFDValues cfdbv;
clock_t simulation_str_time, simulation_end_time;
double simulation_time_total;

int main(int argc, char** argv)
{
	simulation_str_time = clock();
	;
	numericalmainvalues(argv, argc, &cfdbv);
	size(cfdbv.domainsize);
#if AXI
	;
#else
	origin(0, -cfdbv.domainsize / 2., -cfdbv.domainsize / 2.);
#endif
	int initialgrid = pow(2, LEVELmin);
	init_grid(initialgrid);
	;

int Viscous_velocity = 1; // 1 for capillary-viscous velocity scaling and 
							// 0 for capillary-inertial veloity scaling 
	f1.sigma = cfdbv.sigma_1;
	f2.sigma = cfdbv.sigma_2;
	f3.sigma = cfdbv.sigma_3;
	rho1 = cfdbv.rho1;
	rho2 = cfdbv.rho2;
	rho3 = cfdbv.rho3;
	mu1 = cfdbv.mu_1;
	mu2 = cfdbv.mu_2;
	mu3 = cfdbv.mu_3;
	;
	printf("initialgrid: %d\n", initialgrid);
	printf("cfdbv.domainsize: %f\n", cfdbv.domainsize);
	printf("cfdbv.diameter: %f\n", cfdbv.diameter);
	printf("cfdbv.refinegap: %f\n", cfdbv.refinegap);
	printf("cfdbv.Laplas: %f\n", cfdbv.Laplace);
	printf("cfdbv.sigma_12: %f\n", cfdbv.sigma_12);
	printf("cfdbv.sigma_23: %f\n", cfdbv.sigma_23);
	printf("cfdbv.sigma_31: %f\n", cfdbv.sigma_31);
	printf("cfdbv.sigma_23_12: %f\n", cfdbv.sigma_23_12);
	printf("cfdbv.sigma_31_12: %f\n", cfdbv.sigma_31_12);
	printf("LEVELmin: %d\n", LEVELmin);
	printf("LEVELmax: %d\n", LEVELmax);
	;
	TOLERANCE = 1e-6;
	run();
	;
	return 1;
}

event init(t = 0)
{
	double r0 = cfdbv.diameter;
	double x0 = cfdbv.pooldepth;
	double pr = cfdbv.refinegap;
	//
	refine(sq(x - x0) + sq(y) + sq(z) < sq((1.0 + pr) * r0) && sq(x - x0) + sq(y) + sq(z) > sq((1.0 - pr) * r0) && level < LEVELmax);
	refine((x < (1.0 + pr)* x0&& x >(1.0 - pr) * x0 && y > r0) && level < LEVELmax);
	//
	foreach()
	{
		if (sq(x - x0) + sq(y) + sq(z) < sq(r0))
		{
			f1[] = 0.0;
			f2[] = 0.0;
			f3[] = 1.0;
		}
		else if (x < x0)
		{
			f1[] = 1.0;
			f2[] = 0.0;
			f3[] = 0.0;
		}
		else
		{
			f1[] = 0.0;
			f2[] = 1.0;
			f3[] = 0.0;
		}
	}
}

event adapt(i++)
{
	adapt_wavelet(REFINE_VAR, (double[]) REFINE_VAL, maxlevel = LEVELmax, minlevel = LEVELmin);
}

event writedata(t += SAVE_FILE_EVERY)
{
	printf("time: %g\n", t);
	;
	char name[100];
	sprintf(name, "dump-%06d", (int)round(t * 1000));
	dump(name);
	;
	dump(file = FILENAME_LASTFILE);
}

event pngfiles(t += SAVE_FILE_EVERY)
{
	char name[100];
	scalar fall[];
	foreach()
		fall[] = f1[] + 2.0 * f2[] + 3.0 * f3[];
	sprintf(name, "snapshot-%06d.png", (int)round(t * 1000));
	output_ppm(fall, file = name, n = 2000, min = 0, max = 2);
}

event triplepoint(t += SAVE_FILE_EVERY)
{
	char name[100];
	sprintf(name, "triple-point-%06d.txt", (int)round(t * 1000));
	find_triple_point(name, f1, f2, f3);
}

event vtkfiles(t += SAVE_FILE_EVERY)
{
	char name[100];
	const int angle = 90.0;
	sprintf(name, "contour-%06d.vtk", (int)round(t * 1000));
	output_vtk_contour(name, t, angle);
	//
	sprintf(name, "interface-f1-%06d.vtk", (int)round(t * 1000));
	output_vtk_interface(name, f1, t, angle);
	sprintf(name, "interface-f2-%06d.vtk", (int)round(t * 1000));
	output_vtk_interface(name, f2, t, angle);
	sprintf(name, "interface-f3-%06d.vtk", (int)round(t * 1000));
	output_vtk_interface(name, f3, t, angle);
}

event pltfiles(t += SAVE_FILE_EVERY)
{
	char name[100];
	const int angle = 90.0;
	scalar omega[];
	;
	vorticity(u, omega);
	sprintf(name, "contour-%06d.plt", (int)round(t * 1000));
	output_tecplot2D_CC(name, t, {u.x, u.y, p, omega}, angle);
	//
	sprintf(name, "interface-f1-%06d.plt", (int)round(t * 1000));
	output_tecplot2D_FI(name, f1, t, angle);
	sprintf(name, "interface-f2-%06d.plt", (int)round(t * 1000));
	output_tecplot2D_FI(name, f2, t, angle);
	sprintf(name, "interface-f3-%06d.plt", (int)round(t * 1000));
	output_tecplot2D_FI(name, f3, t, angle);
}




event final(t = MAX_TIME)
{
	printf("final iteration");
}

event running_time_estimation(t += SAVE_FILE_EVERY)
{
	double estimatetimeleft;
	char TDc[100], ETLc[100];
	;
	simulation_end_time = clock();
	;
	simulation_time_total += (double)(simulation_end_time - simulation_str_time) / CLOCKS_PER_SEC;
	if (t == 0.0)
		estimatetimeleft = 0.0;
	else
		estimatetimeleft = simulation_time_total * (MAX_TIME) / t - simulation_time_total;
	timecalculation(simulation_time_total, TDc);
	timecalculation(estimatetimeleft, ETLc);
	simulation_str_time = clock();
	;
	switch (pid())
	{
	case 0:
	{
		printf("\r\nDuration until now: %s\t\tTime left: %s\r\n", TDc, ETLc);
		break;
	}
	}
}


event movies(i += 10)
{
	clear();
	view(fov = 6.989, quat = {0,0,-0.707,0.707},
		tx = 1e-6, ty = -0.5, width = 780, height = 382);
	draw_vof("f1", lw = 2);
	draw_vof("f2", lw = 2);
	draw_vof("f3", lw = 2);
	mirror({1}) {
		draw_vof("f1", lw = 2);
		draw_vof("f2", lw = 2);
		draw_vof("f3", lw = 2);
	}
	//  squares ("u.x", spread = -1);
	//  squares ("p");
	cells();
	save("movie.mp4");

	clear();
	view(fov = 0.593863, tx = 0.122977, ty = -0.489743);
	draw_vof("f1", lw = 2, lc = {1,0,0});
	draw_vof("f2", lw = 2, lc = {0,1,0});
	draw_vof("f3", lw = 2, lc = {0,1,1});
	cells();
	save("zoom.mp4");
}