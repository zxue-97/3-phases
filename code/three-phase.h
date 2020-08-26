/**
This file helps setup simulations for flows of three fluids separated by
three interfaces (i.e. immiscible fluids). It is typically used in
combination with a [Navier--Stokes solver](navier-stokes/centered.h).

The interfaces between the fluids is tracked with a Volume-Of-Fluid
method. The volume fraction in fluid 1 is $f1 = 1, f2 = 0, f3 = 0$ and in fluid 2 will be $f2 = 1, f3 = 0, f1 = 0$ and in fluid 3 is $f3 = 1, f1 = 0, f2 = 0$.
The densities and dynamic viscosities for fluid 1, 2 and 3 are *rho1*, *mu1*, *rho2*, *mu2*, *rho3*, *mu3*, respectively. */

scalar f1[], f2[], f3[], * interfaces = { f1, f2, f3 };

/**
Please be noted that we are using a different vof.h, since the interface calculation is modified
*/

#include "vof-3p.h"

double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0., rho3 = 1., mu3 = 0.;

face vector alphav[];
scalar rhov[];

event defaults(i = 0)
{
	alpha = alphav;
	rho = rhov;
	if (mu1 || mu2)
		mu = new face vector;
}

/**
The density and viscosity are defined using arithmetic averages by
default. The user can overload these definitions to use other types of
averages (i.e. harmonic). */

#ifndef rho
#define rho(f1, f2, f3) (clamp(f1, 0., 1.) * rho1 + clamp(f2, 0., 1.) * rho2 + clamp(f3, 0., 1.) * rho3)
#endif
#ifndef mu
#define mu(f1, f2, f3) (clamp(f1, 0., 1.) * mu1 + clamp(f2, 0., 1.) * mu2 + clamp(f3, 0., 1.) * mu3)
#endif

//#if R_FILTERED_FRACTION_Y_N == 'y'
//scalar sf1[], sf2[], sf3[];
//#else
//# define sf1 f1
//# define sf2 f2
//# define sf3 f3
//#endif

event properties(i++)
{
	if (R_CORRECT_FRACTIONS_Y_N == 'y')
	{
		foreach()
		{
			double ft;
			ft = f1[] + f2[] + f3[];
			f1[] /= ft;
			f2[] /= ft;
			f3[] /= ft;
		}
	}
	if (R_FILTERED_FRACTION_Y_N == 'y')
	{
		scalar sf1[], sf2[], sf3[];
#if dimension <= 2
		foreach()
		{
			sf1[] = (4.0 * f1[] + 2.0 * (f1[0, 1] + f1[0, -1] + f1[1, 0] + f1[-1, 0]) +
				f1[-1, -1] + f1[1, -1] + f1[1, 1] + f1[-1, 1]) / 16.0;
			sf2[] = (4.0 * f2[] + 2.0 * (f2[0, 1] + f2[0, -1] + f2[1, 0] + f2[-1, 0]) +
				f2[-1, -1] + f2[1, -1] + f2[1, 1] + f2[-1, 1]) / 16.0;
			sf3[] = (4.0 * f3[] + 2.0 * (f3[0, 1] + f3[0, -1] + f3[1, 0] + f3[-1, 0]) +
				f3[-1, -1] + f3[1, -1] + f3[1, 1] + f3[-1, 1]) / 16.0;
		}
#else // dimension == 3
		foreach()
		{
			sf1[] = (8. * f1[] +
				4. * (f1[-1] + f1[1] + f1[0, 1] + f1[0, -1] + f1[0, 0, 1] + f1[0, 0, -1]) +
				2. * (f1[-1, 1] + f1[-1, 0, 1] + f1[-1, 0, -1] + f1[-1, -1] +
					f1[0, 1, 1] + f1[0, 1, -1] + f1[0, -1, 1] + f1[0, -1, -1] +
					f1[1, 1] + f1[1, 0, 1] + f1[1, -1] + f1[1, 0, -1]) +
				f1[1, -1, 1] + f1[-1, 1, 1] + f1[-1, 1, -1] + f1[1, 1, 1] +
				f1[1, 1, -1] + f1[-1, -1, -1] + f1[1, -1, -1] + f1[-1, -1, 1]) / 64.;
			sf2[] = (8. * f2[] +
				4. * (f2[-1] + f2[1] + f2[0, 1] + f2[0, -1] + f2[0, 0, 1] + f2[0, 0, -1]) +
				2. * (f2[-1, 1] + f2[-1, 0, 1] + f2[-1, 0, -1] + f2[-1, -1] +
					f2[0, 1, 1] + f2[0, 1, -1] + f2[0, -1, 1] + f2[0, -1, -1] +
					f2[1, 1] + f2[1, 0, 1] + f2[1, -1] + f2[1, 0, -1]) +
				f2[1, -1, 1] + f2[-1, 1, 1] + f2[-1, 1, -1] + f2[1, 1, 1] +
				f2[1, 1, -1] + f2[-1, -1, -1] + f2[1, -1, -1] + f2[-1, -1, 1]) / 64.;
			sf3[] = (8. * f3[] +
				4. * (f3[-1] + f3[1] + f3[0, 1] + f3[0, -1] + f3[0, 0, 1] + f3[0, 0, -1]) +
				2. * (f3[-1, 1] + f3[-1, 0, 1] + f3[-1, 0, -1] + f3[-1, -1] +
					f3[0, 1, 1] + f3[0, 1, -1] + f3[0, -1, 1] + f3[0, -1, -1] +
					f3[1, 1] + f3[1, 0, 1] + f3[1, -1] + f3[1, 0, -1]) +
				f3[1, -1, 1] + f3[-1, 1, 1] + f3[-1, 1, -1] + f3[1, 1, 1] +
				f3[1, 1, -1] + f3[-1, -1, -1] + f3[1, -1, -1] + f3[-1, -1, 1]) / 64.;
		}
#endif
#if TREE
		sf1.prolongation = refine_bilinear;
		sf2.prolongation = refine_bilinear;
		sf3.prolongation = refine_bilinear;
		boundary({ sf1, sf2, sf3 });
#endif
		foreach_face()
		{
			double ff1 = (sf1[] + sf1[-1]) / 2.;
			double ff2 = (sf2[] + sf2[-1]) / 2.;
			double ff3 = (sf3[] + sf3[-1]) / 2.;
			alphav.x[] = fm.x[] / rho(ff1, ff2, ff3);
			if (mu1 || mu2)
			{
				face vector muv = mu;
				muv.x[] = fm.x[] * mu(ff1, ff2, ff3);
			}
	}
		foreach()
			rhov[] = cm[] * rho(sf1[], sf2[], sf3[]);

#if TREE
		sf1.prolongation = fraction_refine;
		sf2.prolongation = fraction_refine;
		sf3.prolongation = fraction_refine;
		boundary({ sf1, sf2, sf3 });
#endif
	}
	else
	{
#if TREE
		f1.prolongation = refine_bilinear;
		f2.prolongation = refine_bilinear;
		f3.prolongation = refine_bilinear;
		boundary({ f1, f2, f3 });
#endif
		foreach_face()
		{
			double ff1 = (f1[] + f1[-1]) / 2.;
			double ff2 = (f2[] + f2[-1]) / 2.;
			double ff3 = (f3[] + f3[-1]) / 2.;
			alphav.x[] = fm.x[] / rho(ff1, ff2, ff3);
			if (mu1 || mu2)
			{
				face vector muv = mu;
				muv.x[] = fm.x[] * mu(ff1, ff2, ff3);
			}
		}
		foreach()
			rhov[] = cm[] * rho(f1[], f2[], f3[]);

#if TREE
		f1.prolongation = fraction_refine;
		f2.prolongation = fraction_refine;
		f3.prolongation = fraction_refine;
		boundary({ f1, f2, f3 });
#endif
	}
}