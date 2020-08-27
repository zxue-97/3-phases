
char R_CORRECT_INTERFACE_Y_N = 'n'; // y: yes *** n: no --- in command line use ci to enable this
char R_FILTERED_FRACTION_Y_N = 'n'; // y: yes *** n: no --- in command line use ff to enable this
char R_CORRECT_FRACTIONS_Y_N = 'n'; // y: yes *** n: no --- in command line use cf to enable this
#define R_VOFLIMIT					1.0e-9
#define R_PI						3.1415926535897932384626433832795

#include "axi.h"
#include "navier-stokes/centered.h"
#include "three-phase.h"
#include "tension.h"
#include "tag.h"
#include "triple-point-analysis.h"
#include "tecplot-3p.h"
#include "view-3p.h"

#define DIM_NONDIM_EXP				'n' // d: dimension; n: nondimension;

#if DIM_NONDIM_EXP == 'd' || DIM_NONDIM_EXP == 'D'

#define DROP_DIAMETER				1.0
#define RHO_1						1000.0
#define RHO_2						1000.0
#define RHO_3						1000.0
#define MU_1						1.0
#define MU_2						1.0
#define MU_3						1.0
#define SIGMA_12					1.0
#define SIGMA_23					1.0
#define SIGMA_31					1.0

#define LAPLACE						0.0
#define RHO_2_OVER_RHO_1			0.0
#define RHO_3_OVER_RHO_1			0.0
#define MU_2_OVER_MU_1				0.0
#define MU_3_OVER_MU_1				0.0
#define SIGMA_23_OVER_SIGMA_12		0.0
#define SIGMA_31_OVER_SIGMA_12		0.0

#elif DIM_NONDIM_EXP == 'n' || DIM_NONDIM_EXP == 'N'
// index 1 is the pool, index 2 is the drop, index 3 is above pool
#define LAPLACE						10.0
#define RHO_2_OVER_RHO_1			1.0
#define RHO_3_OVER_RHO_1			1.0
#define MU_2_OVER_MU_1				1.0
#define MU_3_OVER_MU_1				1.0
#define SIGMA_23_OVER_SIGMA_12		1.250 // in the command line call it by sj
#define SIGMA_31_OVER_SIGMA_12		1.750 // in the command line call it by sk

#define DROP_DIAMETER				0.0
#define RHO_1						0.0
#define RHO_2						0.0
#define RHO_3						0.0
#define MU_1						0.0
#define MU_2						0.0
#define MU_3						0.0
#define SIGMA_12					0.0
#define SIGMA_23					0.0
#define SIGMA_31					0.0

#endif

#define INITAL_GRID_LEVEL			6
#define MAX_GRID_LEVEL				7
#define DOMAIN_WIDTH				5.00
#define POOL_DEPTH					2.50
#define REFINE_GAP					0.10
#define MAX_TIME					30.0
#define SAVE_FILE_EVERY				0.10

#define REFINE_VAR					{f1, f2, f3, u.x, u.y}
#define REFINE_VAL					{1.0e-9, 1.0e-9, 1.0e-9, 1.0e-3, 1.0e-3}

#define FILENAME_DATA				"data"
#define FILENAME_DURATION			"duration"
#define FILENAME_PARAMETERS			"parameters.txt"
#define FILENAME_ENDOFRUN			"endofrun"
#define FILENAME_LASTFILE			"lastfile"
#define FILENAME_CONSTANT			"constant.txt"

int LEVELmin = INITAL_GRID_LEVEL, LEVELmax = MAX_GRID_LEVEL;
double maxruntime = HUGE;

struct CFDValues {
	double diameter, rho_1, rho_2, rho_3, mu_1, mu_2, mu_3, sigma_12, sigma_23, sigma_31, sigma_1, sigma_2, sigma_3;
	double Laplace;
	double domainsize, refinegap, pooldepth;
	double timecontact, timeend, timestep;
	double sigma_31_12, sigma_23_12;
};

void readfromarg(char** argv, int argc, struct CFDValues* bvalues);
void output_vtk_contour(char* name, double time, double rotationangle);
void output_vtk_interface(char* name, scalar intrfc, double time, double rotationangle);

void output_vtk_contour(char* name, double time, double rotationangle)
{
	int cno, c;
	double xtmp, ytmp;
	const double angle = rotationangle * R_PI / 180.0;
	FILE* fp;
	;
	cno = 0;
	foreach()
	{
		cno++;
	}
	fp = fopen(name, "w");
	;
	fprintf(fp, "# vtk DataFile Version 2.0\n");
	fprintf(fp, "Example of Unstructured Grid\n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
	;
	fprintf(fp, "POINTS %d float\n", 4 * cno);
	foreach()
	{
		xtmp = x - 0.50 * Delta;
		ytmp = y - 0.50 * Delta;
		fprintf(fp, "%g %g %g ", xtmp * cos(angle) - ytmp * sin(angle), ytmp * cos(angle) + xtmp * sin(angle), 0.0);
		xtmp = x - 0.50 * Delta;
		ytmp = y + 0.50 * Delta;
		fprintf(fp, "%g %g %g ", xtmp * cos(angle) - ytmp * sin(angle), ytmp * cos(angle) + xtmp * sin(angle), 0.0);
		xtmp = x + 0.50 * Delta;
		ytmp = y + 0.50 * Delta;
		fprintf(fp, "%g %g %g ", xtmp * cos(angle) - ytmp * sin(angle), ytmp * cos(angle) + xtmp * sin(angle), 0.0);
		xtmp = x + 0.50 * Delta;
		ytmp = y - 0.50 * Delta;
		fprintf(fp, "%g %g %g ", xtmp * cos(angle) - ytmp * sin(angle), ytmp * cos(angle) + xtmp * sin(angle), 0.0);
		fprintf(fp, "\n");
	}
	fprintf(fp, "CELLS %d %d\n", cno, cno + 4 * cno);
	for (c = 0; c < cno; c++)
		fprintf(fp, "%d %d %d %d %d ", 4, 4 * c + 0, 4 * c + 1, 4 * c + 2, 4 * c + 3);
	fprintf(fp, "\n");
	;
	fprintf(fp, "CELL_TYPES %d\n", cno);
	for (c = 0; c < cno; c++)
		fprintf(fp, "9 ");
	fprintf(fp, "\n");
	;
	fprintf(fp, "CELL_DATA %d\n", cno);
	;
	fprintf(fp, "SCALARS fall float\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	foreach()
	{
		fprintf(fp, "%g ", f1[] + 2.0 * f2[] + 3.0 * f3[]);
	}
	fprintf(fp, "\n");
	;
	fprintf(fp, "SCALARS f1 float\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	foreach()
	{
		fprintf(fp, "%g ", f1[]);
	}
	fprintf(fp, "\n");
	;
	fprintf(fp, "SCALARS f2 float\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	foreach()
	{
		fprintf(fp, "%g ", f2[]);
	}
	fprintf(fp, "\n");
	;
	fprintf(fp, "SCALARS f3 float\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	foreach()
	{
		fprintf(fp, "%g ", f3[]);
	}
	fprintf(fp, "\n");
	;
	fprintf(fp, "SCALARS P float\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	foreach()
	{
		fprintf(fp, "%g ", p[]);
	}
	fprintf(fp, "\n");
	;
	fclose(fp);
}

void output_vtk_interface(char* name, scalar intrfc, double time, double rotationangle)
{
	int interfacepoints, cellnumber, iii, interfacecells;
	char typei[2], * typeintersect;
	double xi[2], yi[2], * xintersect, * yintersect;
	FILE* fp;
	const double angle = rotationangle * R_PI / 180.0;
	scalar alpha[];
	vector n[];
	;
	cellnumber = 0;
	foreach()
		cellnumber++;
	xintersect = (double*)calloc(sizeof(double), cellnumber);
	yintersect = (double*)calloc(sizeof(double), cellnumber);
	typeintersect = (char*)calloc(sizeof(char), cellnumber);
	;
	reconstruction(intrfc, n, alpha);
	interfacepoints = 0;
	foreach_leaf()
	{
		if (intrfc[] > R_VOFLIMIT && intrfc[] < 1.0 - R_VOFLIMIT)
		{
			findintersectionpoints(n.x[], n.y[], alpha[], x, y, Delta, xi, yi, typei);
			xintersect[interfacepoints] = xi[0] * cos(angle) - yi[0] * sin(angle);
			yintersect[interfacepoints] = yi[0] * cos(angle) + xi[0] * sin(angle);
			typeintersect[interfacepoints] = typei[0];
			interfacepoints++;
			xintersect[interfacepoints] = xi[1] * cos(angle) - yi[1] * sin(angle);
			yintersect[interfacepoints] = yi[1] * cos(angle) + xi[1] * sin(angle);
			typeintersect[interfacepoints] = typei[1];
			interfacepoints++;
		}
	}
	interfacecells = interfacepoints / 2;
	fp = fopen(name, "w");
	fprintf(fp, "# vtk DataFile Version 2.0\n");
	fprintf(fp, "Example of Unstructured Grid\n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
	switch (interfacepoints)
	{
	case 0:
	{
		fprintf(fp, "POINTS %d float\n", 3);
		fprintf(fp, "%g %g %g ", 0.0, 0.0, 0.0);
		fprintf(fp, "%g %g %g ", 0.0, 0.0, 0.0);
		fprintf(fp, "%g %g %g ", 0.0, 0.0, 0.0);
		fprintf(fp, "\n");
		fprintf(fp, "CELLS %d %d\n", 1, 4);
		fprintf(fp, "%d %d %d %d ", 3, 0, 1, 2);
		fprintf(fp, "\n");
		fprintf(fp, "CELL_TYPES %d\n", 1);
		fprintf(fp, "5 ");
		fprintf(fp, "\n");
		fprintf(fp, "CELL_DATA %d\n", 1);
		fprintf(fp, "SCALARS fall float\n");
		fprintf(fp, "LOOKUP_TABLE default\n");
		fprintf(fp, "%g ", 1.0);
		fprintf(fp, "\n");
		break;
	}
	default:
	{
		fprintf(fp, "POINTS %d float\n", 3 * interfacecells);
		for (iii = 0; iii < interfacepoints - 1; iii += 2)
		{
			fprintf(fp, "%g %g %g ", xintersect[iii], yintersect[iii], 0.0);
			fprintf(fp, "%g %g %g ", xintersect[iii + 1], yintersect[iii + 1], 0.0);
			fprintf(fp, "%g %g %g ", xintersect[iii], yintersect[iii], 0.0);
		}
		fprintf(fp, "\n");
		;
		fprintf(fp, "CELLS %d %d\n", interfacecells, interfacecells + 3 * interfacecells);
		for (iii = 0; iii < interfacecells; iii++)
			fprintf(fp, "%d %d %d %d ", 3, 3 * iii + 0, 3 * iii + 1, 3 * iii + 2);
		fprintf(fp, "\n");
		;
		fprintf(fp, "CELL_TYPES %d\n", interfacecells);
		for (iii = 0; iii < interfacecells; iii++)
			fprintf(fp, "5 ");
		fprintf(fp, "\n");
		;
		fprintf(fp, "CELL_DATA %d\n", interfacecells);
		;
		fprintf(fp, "SCALARS fall float\n");
		fprintf(fp, "LOOKUP_TABLE default\n");
		for (iii = 0; iii < interfacecells; iii++)
		{
			fprintf(fp, "%g ", 1.0);
		}
		fprintf(fp, "\n");
		break;
	}
	}
	fclose(fp);
	;
	free(xintersect);
	free(yintersect);
}

int numericalmainvalues(char** argv, int argc, struct CFDValues* bvalues)
{
	bvalues->sigma_12 	= 1.0;
	bvalues->mu_1 		= 1.0;
	bvalues->diameter 	= 1.0;
	;
	bvalues->Laplace = -1.0;
	bvalues->pooldepth = -1.0;
	bvalues->timeend = -1.0;
	bvalues->timestep = -1.0;
	;
	bvalues->sigma_23_12 = -1.0;
	bvalues->sigma_31_12 = -1.0;
	;
	readfromarg(argv, argc, bvalues);
	switch (DIM_NONDIM_EXP)
	{
	case 'd':
	case 'D':
	{
		bvalues->rho_1 = RHO_1;
		bvalues->rho_2 = RHO_2;
		bvalues->rho_3 = RHO_3;
		bvalues->mu_1 = MU_1;
		bvalues->mu_2 = MU_2;
		bvalues->mu_3 = MU_3;
		bvalues->sigma_12 = SIGMA_12;
		bvalues->sigma_23 = SIGMA_23;
		bvalues->sigma_31 = SIGMA_31;
		bvalues->diameter = DROP_DIAMETER;
		bvalues->Laplace = bvalues->sigma_12 * bvalues->rho_1 * bvalues->diameter / (bvalues->mu_1 * bvalues->mu_1);
		break;
	}
	case 'n':
	case 'N':
	{
		if (bvalues->Laplace < 0.0)
			bvalues->Laplace = LAPLACE;
		bvalues->rho_1 = bvalues->Laplace * bvalues->mu_1 * bvalues->mu_1 / (bvalues->sigma_12 * bvalues->diameter);
		bvalues->rho_2 = RHO_2_OVER_RHO_1 * bvalues->rho_1;
		bvalues->rho_3 = RHO_3_OVER_RHO_1 * bvalues->rho_1;
		bvalues->mu_2 = MU_2_OVER_MU_1 * bvalues->mu_1;
		bvalues->mu_3 = MU_3_OVER_MU_1 * bvalues->mu_1;
		if (bvalues->sigma_23_12 < 0.0)
			bvalues->sigma_23_12 = SIGMA_23_OVER_SIGMA_12;
		bvalues->sigma_23 = bvalues->sigma_23_12 * bvalues->sigma_12;
		if (bvalues->sigma_31_12 < 0.0)
			bvalues->sigma_31_12 = SIGMA_31_OVER_SIGMA_12;
		bvalues->sigma_31 = bvalues->sigma_31_12 * bvalues->sigma_12;
		break;
	}
	}
	bvalues->domainsize = DOMAIN_WIDTH * bvalues->diameter;
	if (bvalues->pooldepth < 0.0)
		bvalues->pooldepth = POOL_DEPTH * bvalues->diameter;
	bvalues->refinegap = REFINE_GAP * bvalues->diameter;
	;
	if (bvalues->timeend < 0.0)
		bvalues->timeend = MAX_TIME;
	if (bvalues->timestep < 0.0)
		bvalues->timestep = SAVE_FILE_EVERY;
	;
	bvalues->sigma_1 = 0.50 * (bvalues->sigma_12 + bvalues->sigma_31 - bvalues->sigma_23);
	bvalues->sigma_2 = 0.50 * (bvalues->sigma_23 + bvalues->sigma_12 - bvalues->sigma_31);
	bvalues->sigma_3 = 0.50 * (bvalues->sigma_31 + bvalues->sigma_23 - bvalues->sigma_12);
	return 1;
}

void readfromarg(char** argv, int argc, struct CFDValues* bvalues)
{
	int i, j;
	char tmp[100];
	if (argc < 2)
		return;
	for (i = 1; i < argc; i++)
	{
		switch (argv[i][0])
		{
		case 'l':
		case 'L':
		{
			for (j = 1; j < (int)strlen(argv[i]); j++)
				tmp[j - 1] = argv[i][j];
			tmp[j - 1] = '\0';
			bvalues->Laplace = atof(tmp);
			break;
		}
		case 'c':
		case 'C':
		{
			switch (argv[i][1])
			{
			case 'i':
			case 'I':
			{
				R_CORRECT_INTERFACE_Y_N = 'y';
				break;
			}
			case 'f':
			case 'F':
			{
				R_CORRECT_FRACTIONS_Y_N = 'y';
				break;
			}
			}
			break;
		}
		case 'f':
		case 'F':
		{
			switch (argv[i][1])
			{
			case 'f':
			case 'F':
			{
				R_FILTERED_FRACTION_Y_N = 'y';
				break;
			}
			}
			break;
		}
		case 's':
		case 'S':
		{
			switch (argv[i][1])
			{
			case 'j':
			case 'J':
			{
				for (j = 2; j < (int)strlen(argv[i]); j++)
					tmp[j - 2] = argv[i][j];
				tmp[j - 2] = '\0';
				bvalues->sigma_23_12 = atof(tmp);
				break;
			}
			case 'k':
			case 'K':
			{
				for (j = 2; j < (int)strlen(argv[i]); j++)
					tmp[j - 2] = argv[i][j];
				tmp[j - 2] = '\0';
				bvalues->sigma_31_12 = atof(tmp);
				break;
			}
			}
			break;
		}
		case 'x':
		case 'X':
		{
			for (j = 1; j < (int)strlen(argv[i]); j++)
				tmp[j - 1] = argv[i][j];
			tmp[j - 1] = '\0';
			LEVELmax = atoi(tmp);
			break;
		}
		case 'n':
		case 'N':
		{
			for (j = 1; j < (int)strlen(argv[i]); j++)
				tmp[j - 1] = argv[i][j];
			tmp[j - 1] = '\0';
			LEVELmin = atoi(tmp);
			break;
		}
		case 't':
		case 'T':
		{
			switch (argv[i][1])
			{
			case 's':
			case 'S':
			{
				for (j = 2; j < (int)strlen(argv[i]); j++)
					tmp[j - 2] = argv[i][j];
				tmp[j - 2] = '\0';
				bvalues->timestep = atof(tmp);
				break;
			}
			case 'e':
			case 'E':
			{
				for (j = 2; j < (int)strlen(argv[i]); j++)
					tmp[j - 2] = argv[i][j];
				tmp[j - 2] = '\0';
				bvalues->timeend = atof(tmp);
				break;
			}
			}
			break;
		}
		}
	}
}

int findintersectionpoints(double nx, double ny, double alpha, double xc, double yc, double delta, double xi[2], double yi[2], char typei[2])
{
	int ppp = 0;
	double xtmp[2], ytmp[2], underflow = 1.0e-6;
	if (fabs(nx) < underflow)
	{
		ytmp[0] = (alpha - nx * (-0.50)) / ny;
		ytmp[1] = (alpha - nx * (+0.50)) / ny;
		xi[ppp] = xc + (-0.50) * delta;
		yi[ppp] = yc + (ytmp[0]) * delta;
		typei[ppp] = 'l';
		(ppp)++;
		xi[ppp] = xc + (+0.50) * delta;
		yi[ppp] = yc + (ytmp[1]) * delta;
		typei[ppp] = 'r';
		(ppp)++;
	}
	else if (fabs(ny) < underflow)
	{
		xtmp[0] = (alpha - ny * (-0.50)) / nx;
		xtmp[1] = (alpha - ny * (+0.50)) / nx;
		xi[ppp] = xc + (xtmp[0]) * delta;
		yi[ppp] = yc + (-0.50) * delta;
		typei[ppp] = 'b';
		(ppp)++;
		xi[ppp] = xc + (xtmp[1]) * delta;
		yi[ppp] = yc + (+0.50) * delta;
		typei[ppp] = 't';
		(ppp)++;
	}
	else
	{
		xtmp[0] = (alpha - ny * (-0.50)) / nx;
		xtmp[1] = (alpha - ny * (+0.50)) / nx;
		ytmp[0] = (alpha - nx * (-0.50)) / ny;
		ytmp[1] = (alpha - nx * (+0.50)) / ny;

		if (-0.50 <= ytmp[0] && ytmp[0] <= +0.50)
		{
			xi[ppp] = xc + (-0.50) * delta;
			yi[ppp] = yc + (ytmp[0]) * delta;
			typei[ppp] = 'l';
			(ppp)++;
		}
		if (-0.50 <= xtmp[0] && xtmp[0] <= +0.50)
		{
			xi[ppp] = xc + (xtmp[0]) * delta;
			yi[ppp] = yc + (-0.50) * delta;
			typei[ppp] = 'b';
			(ppp)++;
		}
		if (-0.50 <= ytmp[1] && ytmp[1] <= +0.50)
		{
			xi[ppp] = xc + (+0.50) * delta;
			yi[ppp] = yc + (ytmp[1]) * delta;
			typei[ppp] = 'r';
			(ppp)++;
		}
		if (-0.50 <= xtmp[1] && xtmp[1] <= +0.50)
		{
			xi[ppp] = xc + (xtmp[1]) * delta;
			yi[ppp] = yc + (+0.50) * delta;
			typei[ppp] = 't';
			(ppp)++;
		}
	}
	return true;
}

int timecalculation(double t, char* chartime)
{
	int d, h, m, s;
	if (t < 60.0)
	{
		d = 0;
		h = 0;
		m = 0;
		s = (int)t;
	}
	else if (t < 3600.0)
	{
		d = 0;
		h = 0;
		m = (int)(t / 60.0);
		s = (int)(t - m * 60.0);
	}
	else if (t < 3600.0 * 24.0)
	{
		d = 0;
		h = (int)(t / 3600.0);
		m = (int)((t - h * 3600.0) / 60.0);
		s = (int)(t - h * 3600.0 - m * 60.0);
	}
	else
	{
		d = (int)(t / 3600.0 / 24.0);
		h = (int)((t - d * 3600.0 * 24.0) / 3600.0);
		m = (int)((t - d * 3600.0 * 24.0 - h * 3600.0) / 60.0);
		s = (int)(t - d * 3600.0 * 24.0 - h * 3600.0 - m * 60.0);
	}
	sprintf(chartime, "%d:%02d:%02d.%02d (d:hh:mm.ss)", d, h, m, s);
	return 1;
}
