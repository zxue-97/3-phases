
struct CellValues {
	double xc, yc, delta;
	double xtfw, ytfw, vvfw;
	double xtar, ytar, vvar;
	double xf[2][3], yf[2][3], nxf[3], nyf[3], f[3];
	double theta[4]; // 0: TR, 1: TL, 2: BL, 3: BR
	double xtp[3][3], ytp[3][3]; // two interfaces, unique line
	double area_integrated[3];
};

void find_void_area(struct CellValues* cv, int res);
void find_four_corner_angles(struct CellValues* cv, double xbase, double ybase);
void find_all_interface_lines_f2f1f0(struct CellValues* cv, double xbase, double ybase);
void find_all_interface_lines_f0f1f2(struct CellValues* cv, double xbase, double ybase);
bool point_in_which_side_line(double xl, double yl, double nx, double ny, double xp, double yp);
int check_if_p2_is_CW_cell_coordinates(double xp1, double yp1, double xp2, double yp2);
double triangle_area(double x[3], double y[3]);
double find_area_pt_and_two_border_points_cell_coordinates(double xpt, double ypt, double xb[2], double yb[2]);
void find_next_CCW_point_cell_coordinates(double* xxx, double* yyy, int res);
void integration_by_rotating_around_triple_point(double fraction, double xtcc, double ytcc, double xtmpcc[2], double ytmpcc[2], double* area);
int findintersectionpoints(double nx, double ny, double alpha, double xc, double yc, double delta, double xi[2], double yi[2], char typei[2]);

void find_triple_point(char* name, scalar f1, scalar f2, scalar f3)
{
	int ci, c;
	const double limit = 0.01;
	char typei[2];
	double xi[2], yi[2];
	struct CellValues cv[100];
	scalar a1[], a2[], a3[];
	vector n1[], n2[], n3[];
	FILE* fp;
	;
	reconstruction(f1, n1, a1);
	reconstruction(f2, n2, a2);
	reconstruction(f3, n3, a3);
	;
	ci = 0;
	foreach()
	{
		if (f1[] > limit && f1[] < 1.0 - limit)
		{
			if (f2[] > limit && f2[] < 1.0 - limit)
			{
				if (f3[] > limit && f3[] < 1.0 - limit)
				{
					findintersectionpoints(n1.x[], n1.y[], a1[], x, y, Delta, xi, yi, typei);
					cv[ci].xf[0][0] = xi[0];
					cv[ci].yf[0][0] = yi[0];
					cv[ci].xf[1][0] = xi[1];
					cv[ci].yf[1][0] = yi[1];
					;
					findintersectionpoints(n2.x[], n2.y[], a2[], x, y, Delta, xi, yi, typei);
					cv[ci].xf[0][1] = xi[0];
					cv[ci].yf[0][1] = yi[0];
					cv[ci].xf[1][1] = xi[1];
					cv[ci].yf[1][1] = yi[1];
					;
					findintersectionpoints(n3.x[], n3.y[], a3[], x, y, Delta, xi, yi, typei);
					cv[ci].xf[0][2] = xi[0];
					cv[ci].yf[0][2] = yi[0];
					cv[ci].xf[1][2] = xi[1];
					cv[ci].yf[1][2] = yi[1];
					;
					cv[ci].xc = x;
					cv[ci].yc = y;
					cv[ci].delta = Delta;
					cv[ci].nxf[0] = n1.x[];
					cv[ci].nyf[0] = n1.y[];
					cv[ci].nxf[1] = n2.x[];
					cv[ci].nyf[1] = n2.y[];
					cv[ci].nxf[2] = n3.x[];
					cv[ci].nyf[2] = n3.y[];
					cv[ci].f[0] = f1[];
					cv[ci].f[1] = f2[];
					cv[ci].f[2] = f3[];
					;
					ci++;
				}
			}
		}
	}
	if (ci > 0)
	{
		fp = fopen(name, "w");
		for (c = 0; c < ci; c++)
		{
			fprintf(fp, "x\ty\tdelta\t%f\t%f\t%f\n", cv[c].xc, cv[c].yc, cv[c].delta);
			fprintf(fp, "f1\tx1\tx2\ty1\ty2\t%f\t%f\t%f\t%f\t%f\n", cv[c].f[0], cv[c].xf[0][0], cv[c].xf[1][0], cv[c].yf[0][0], cv[c].yf[1][0]);
			fprintf(fp, "f2\tx1\tx2\ty1\ty2\t%f\t%f\t%f\t%f\t%f\n", cv[c].f[1], cv[c].xf[0][1], cv[c].xf[1][1], cv[c].yf[0][1], cv[c].yf[1][1]);
			fprintf(fp, "f3\tx1\tx2\ty1\ty2\t%f\t%f\t%f\t%f\t%f\n", cv[c].f[2], cv[c].xf[0][2], cv[c].xf[1][2], cv[c].yf[0][2], cv[c].yf[1][2]);
			;
			find_void_area(&(cv[c]), 1000);
			;
			find_four_corner_angles(&(cv[c]), cv[c].xtar, cv[c].ytar);
			find_all_interface_lines_f0f1f2(&(cv[c]), cv[c].xtar, cv[c].ytar);
			fprintf(fp, "Area\n");
			fprintf(fp, "Void Area\tVoid Fraction\t%f\t%f\n", cv[c].vvar, cv[c].vvar / (cv[c].delta * cv[c].delta));
			fprintf(fp, "f12\tXT\tYT\t%f\t%f\n", cv[c].xtar, cv[c].ytar);
			fprintf(fp, "f12\tXT\tYT\t%f\t%f\n", cv[c].xtp[0][1], cv[c].ytp[0][1]);
			fprintf(fp, "f12\tXT\tYT\t%f\t%f\n", cv[c].xtp[1][0], cv[c].ytp[1][0]);
			fprintf(fp, "f23\tXT\tYT\t%f\t%f\n", cv[c].xtar, cv[c].ytar);
			fprintf(fp, "f23\tXT\tYT\t%f\t%f\n", cv[c].xtp[1][2], cv[c].ytp[1][2]);
			fprintf(fp, "f23\tXT\tYT\t%f\t%f\n", cv[c].xtp[2][1], cv[c].ytp[2][1]);
			fprintf(fp, "f31\tXT\tYT\t%f\t%f\n", cv[c].xtar, cv[c].ytar);
			fprintf(fp, "f31\tXT\tYT\t%f\t%f\n", cv[c].xtp[2][0], cv[c].ytp[2][0]);
			fprintf(fp, "f31\tXT\tYT\t%f\t%f\n", cv[c].xtp[0][2], cv[c].ytp[0][2]);
			fprintf(fp, "\n");
			find_all_interface_lines_f2f1f0(&(cv[c]), cv[c].xtar, cv[c].ytar);
			fprintf(fp, "f12\tXT\tYT\t%f\t%f\n", cv[c].xtar, cv[c].ytar);
			fprintf(fp, "f12\tXT\tYT\t%f\t%f\n", cv[c].xtp[0][1], cv[c].ytp[0][1]);
			fprintf(fp, "f12\tXT\tYT\t%f\t%f\n", cv[c].xtp[1][0], cv[c].ytp[1][0]);
			fprintf(fp, "f23\tXT\tYT\t%f\t%f\n", cv[c].xtar, cv[c].ytar);
			fprintf(fp, "f23\tXT\tYT\t%f\t%f\n", cv[c].xtp[1][2], cv[c].ytp[1][2]);
			fprintf(fp, "f23\tXT\tYT\t%f\t%f\n", cv[c].xtp[2][1], cv[c].ytp[2][1]);
			fprintf(fp, "f31\tXT\tYT\t%f\t%f\n", cv[c].xtar, cv[c].ytar);
			fprintf(fp, "f31\tXT\tYT\t%f\t%f\n", cv[c].xtp[2][0], cv[c].ytp[2][0]);
			fprintf(fp, "f31\tXT\tYT\t%f\t%f\n", cv[c].xtp[0][2], cv[c].ytp[0][2]);
			fprintf(fp, "\n");
			;
			find_four_corner_angles(&(cv[c]), cv[c].xtfw, cv[c].ytfw);
			find_all_interface_lines_f0f1f2(&(cv[c]), cv[c].xtfw, cv[c].ytfw);
			fprintf(fp, "Fraction Weighted\n");
			fprintf(fp, "Void Area (FW)\tVoid Fraction (FW)\t%f\t%f\n", cv[c].vvfw, cv[c].vvfw / (cv[c].delta * cv[c].delta));
			fprintf(fp, "f12\tXT\tYT\t%f\t%f\n", cv[c].xtfw, cv[c].ytfw);
			fprintf(fp, "f12\tXT\tYT\t%f\t%f\n", cv[c].xtp[0][1], cv[c].ytp[0][1]);
			fprintf(fp, "f12\tXT\tYT\t%f\t%f\n", cv[c].xtp[1][0], cv[c].ytp[1][0]);
			fprintf(fp, "f23\tXT\tYT\t%f\t%f\n", cv[c].xtfw, cv[c].ytfw);
			fprintf(fp, "f23\tXT\tYT\t%f\t%f\n", cv[c].xtp[1][2], cv[c].ytp[1][2]);
			fprintf(fp, "f23\tXT\tYT\t%f\t%f\n", cv[c].xtp[2][1], cv[c].ytp[2][1]);
			fprintf(fp, "f31\tXT\tYT\t%f\t%f\n", cv[c].xtfw, cv[c].ytfw);
			fprintf(fp, "f31\tXT\tYT\t%f\t%f\n", cv[c].xtp[2][0], cv[c].ytp[2][0]);
			fprintf(fp, "f31\tXT\tYT\t%f\t%f\n", cv[c].xtp[0][2], cv[c].ytp[0][2]);
			fprintf(fp, "\n");
			find_all_interface_lines_f2f1f0(&(cv[c]), cv[c].xtfw, cv[c].ytfw);
			fprintf(fp, "f12\tXT\tYT\t%f\t%f\n", cv[c].xtfw, cv[c].ytfw);
			fprintf(fp, "f12\tXT\tYT\t%f\t%f\n", cv[c].xtp[0][1], cv[c].ytp[0][1]);
			fprintf(fp, "f12\tXT\tYT\t%f\t%f\n", cv[c].xtp[1][0], cv[c].ytp[1][0]);
			fprintf(fp, "f23\tXT\tYT\t%f\t%f\n", cv[c].xtfw, cv[c].ytfw);
			fprintf(fp, "f23\tXT\tYT\t%f\t%f\n", cv[c].xtp[1][2], cv[c].ytp[1][2]);
			fprintf(fp, "f23\tXT\tYT\t%f\t%f\n", cv[c].xtp[2][1], cv[c].ytp[2][1]);
			fprintf(fp, "f31\tXT\tYT\t%f\t%f\n", cv[c].xtfw, cv[c].ytfw);
			fprintf(fp, "f31\tXT\tYT\t%f\t%f\n", cv[c].xtp[2][0], cv[c].ytp[2][0]);
			fprintf(fp, "f31\tXT\tYT\t%f\t%f\n", cv[c].xtp[0][2], cv[c].ytp[0][2]);
			fprintf(fp, "-------------------------------\n");
		}
		fclose(fp);
	}
}

void find_all_interface_lines_f2f1f0(struct CellValues* cv, double xbase, double ybase)
{
	int i, j, f1, f2, f3;
	double dx, dy, th, xtmp, ytmp, xtcc, ytcc, xtmpcc[2], ytmpcc[2], area_tmp;
	bool side_f_fn_point[3][3][2];
	int number_f_fn_point[3][3][2];
	;
	for (f1 = 0; f1 < 3; f1++)
	{
		for (f2 = 0; f2 < 3; f2++)
		{
			side_f_fn_point[f1][f2][0] = point_in_which_side_line(cv->xf[0][f1], cv->yf[0][f1], cv->nxf[f1], cv->nyf[f1], cv->xf[0][f2], cv->yf[0][f2]);
			side_f_fn_point[f1][f2][1] = point_in_which_side_line(cv->xf[0][f1], cv->yf[0][f1], cv->nxf[f1], cv->nyf[f1], cv->xf[1][f2], cv->yf[1][f2]);
			cv->xtp[f1][f2] = -1e6;
			cv->ytp[f1][f2] = -1e6;
			number_f_fn_point[f1][f2][0] = -1;
			number_f_fn_point[f1][f2][1] = -1;
		}
	}
	;
	for (f1 = 2; f1 > -1; f1--)
	{
		f2 = (f1 + 1) % 3;
		f3 = (f1 + 2) % 3;
		{
			for (i = 0; i < 2; i++)
			{
				for (j = 0; j < 2; j++)
				{
					if (side_f_fn_point[f1][f2][i] && side_f_fn_point[f2][f1][j] && !side_f_fn_point[f1][f2][1 - i] && !side_f_fn_point[f2][f1][1 - j])
					{
						cv->xtp[f1][f2] = 0.50 * (cv->xf[j][f1] + cv->xf[i][f2]);
						cv->ytp[f1][f2] = 0.50 * (cv->yf[j][f1] + cv->yf[i][f2]);
						number_f_fn_point[f1][f2][0] = i;
						number_f_fn_point[f2][f1][1] = j;
						;
						dx = cv->xtp[f1][f2] - xbase;
						dy = cv->ytp[f1][f2] - ybase;
						th = atan2(dy, dx);
						if (th < 0.0)
							th += 2.0 * R_PI;
						if (cv->theta[0] <= th && th <= cv->theta[1])
						{
							ytmp = cv->yc + 0.50 * cv->delta;
							xtmp = xbase + (dx / dy) * (ytmp - ybase);
						}
						else if (cv->theta[1] <= th && th <= cv->theta[2])
						{
							xtmp = cv->xc - 0.50 * cv->delta;
							ytmp = ybase + (dy / dx) * (xtmp - xbase);
						}
						else if (cv->theta[2] <= th && th <= cv->theta[3])
						{
							ytmp = cv->yc - 0.50 * cv->delta;
							xtmp = xbase + (dx / dy) * (ytmp - ybase);
						}
						else
						{
							xtmp = cv->xc + 0.50 * cv->delta;
							ytmp = ybase + (dy / dx) * (xtmp - xbase);
						}
						cv->xtp[f1][f2] = xtmp;
						cv->ytp[f1][f2] = ytmp;
						cv->xtp[f2][f1] = cv->xtp[f1][f2];
						cv->ytp[f2][f1] = cv->ytp[f1][f2];
						;
						i = 10;
						j = 10;
						;
						f1 = -10;
					}
				}
			}
		}
	}
	;
	for (f1 = 2; f1 > -1; f1--)
	{
		f2 = (f1 + 1) % 3;
		f3 = (f1 + 2) % 3;
		if (number_f_fn_point[f1][f2][0] != -1)
		{
			i = number_f_fn_point[f1][f2][0];
			j = number_f_fn_point[f2][f1][1];
			xtcc = (xbase - cv->xc) / cv->delta;
			ytcc = (ybase - cv->yc) / cv->delta;
			xtmpcc[0] = (cv->xtp[f1][f2] - cv->xc) / cv->delta;
			ytmpcc[0] = (cv->ytp[f1][f2] - cv->yc) / cv->delta;
			xtmpcc[1] = (cv->xf[j][f1] - cv->xc) / cv->delta;
			ytmpcc[1] = (cv->yf[j][f1] - cv->yc) / cv->delta;
			if (check_if_p2_is_CW_cell_coordinates(xtmpcc[0], ytmpcc[0], xtmpcc[1], ytmpcc[1]) == -1)
			{
				integration_by_rotating_around_triple_point(cv->f[f1], xtcc, ytcc, xtmpcc, ytmpcc, &area_tmp);
				;
				cv->area_integrated[f1] = area_tmp;
				cv->xtp[f1][f3] = cv->xc + xtmpcc[1] * cv->delta;
				cv->ytp[f1][f3] = cv->yc + ytmpcc[1] * cv->delta;
				cv->xtp[f3][f1] = cv->xtp[f1][f3];
				cv->ytp[f3][f1] = cv->ytp[f1][f3];
				;
				xtmpcc[0] = xtmpcc[1];
				ytmpcc[0] = ytmpcc[1];
				;
				integration_by_rotating_around_triple_point(cv->f[f3], xtcc, ytcc, xtmpcc, ytmpcc, &area_tmp);
				;
				cv->area_integrated[f3] = area_tmp;
				cv->xtp[f2][f3] = cv->xc + xtmpcc[1] * cv->delta;
				cv->ytp[f2][f3] = cv->yc + ytmpcc[1] * cv->delta;
				cv->xtp[f3][f2] = cv->xtp[f2][f3];
				cv->ytp[f3][f2] = cv->ytp[f2][f3];
				;
				cv->area_integrated[f2] = 1.0 - cv->area_integrated[f1] - cv->area_integrated[f3];
			}
			else
			{
				xtmpcc[0] = (cv->xtp[f1][f2] - cv->xc) / cv->delta;
				ytmpcc[0] = (cv->ytp[f1][f2] - cv->yc) / cv->delta;
				xtmpcc[1] = (cv->xf[i][f2] - cv->xc) / cv->delta;
				ytmpcc[1] = (cv->yf[i][f2] - cv->yc) / cv->delta;
				if (check_if_p2_is_CW_cell_coordinates(xtmpcc[0], ytmpcc[0], xtmpcc[1], ytmpcc[1]) == -1)
				{
					integration_by_rotating_around_triple_point(cv->f[f2], xtcc, ytcc, xtmpcc, ytmpcc, &area_tmp);
					;
					cv->area_integrated[f2] = area_tmp;
					cv->xtp[f2][f3] = cv->xc + xtmpcc[1] * cv->delta;
					cv->ytp[f2][f3] = cv->yc + ytmpcc[1] * cv->delta;
					cv->xtp[f3][f2] = cv->xtp[f2][f3];
					cv->ytp[f3][f2] = cv->ytp[f2][f3];
					;
					xtmpcc[0] = xtmpcc[1];
					ytmpcc[0] = ytmpcc[1];
					;
					integration_by_rotating_around_triple_point(cv->f[f3], xtcc, ytcc, xtmpcc, ytmpcc, &area_tmp);
					;
					cv->area_integrated[f3] = area_tmp;
					cv->xtp[f1][f3] = cv->xc + xtmpcc[1] * cv->delta;
					cv->ytp[f1][f3] = cv->yc + ytmpcc[1] * cv->delta;
					cv->xtp[f3][f1] = cv->xtp[f1][f3];
					cv->ytp[f3][f1] = cv->ytp[f1][f3];
					;
					cv->area_integrated[f1] = 1.0 - cv->area_integrated[f2] - cv->area_integrated[f3];
				}
				else
				{
					printf("ERRRRROOOORRRR\n");
				}
			}
			f1 = -10;
		}
	}
}

void find_all_interface_lines_f0f1f2(struct CellValues* cv, double xbase, double ybase)
{
	int i, j, f1, f2, f3;
	double dx, dy, th, xtmp, ytmp, xtcc, ytcc, xtmpcc[2], ytmpcc[2], area_tmp;
	bool side_f_fn_point[3][3][2];
	int number_f_fn_point[3][3][2];
	;
	for (f1 = 0; f1 < 3; f1++)
	{
		for (f2 = 0; f2 < 3; f2++)
		{
			side_f_fn_point[f1][f2][0] = point_in_which_side_line(cv->xf[0][f1], cv->yf[0][f1], cv->nxf[f1], cv->nyf[f1], cv->xf[0][f2], cv->yf[0][f2]);
			side_f_fn_point[f1][f2][1] = point_in_which_side_line(cv->xf[0][f1], cv->yf[0][f1], cv->nxf[f1], cv->nyf[f1], cv->xf[1][f2], cv->yf[1][f2]);
			cv->xtp[f1][f2] = -1e6;
			cv->ytp[f1][f2] = -1e6;
			number_f_fn_point[f1][f2][0] = -1;
			number_f_fn_point[f1][f2][1] = -1;
		}
	}
	;
	for (f1 = 0; f1 < 3; f1++)
	{
		f2 = (f1 + 1) % 3;
		f3 = (f1 + 2) % 3;
		{
			for (i = 0; i < 2; i++)
			{
				for (j = 0; j < 2; j++)
				{
					if (side_f_fn_point[f1][f2][i] && side_f_fn_point[f2][f1][j] && !side_f_fn_point[f1][f2][1 - i] && !side_f_fn_point[f2][f1][1 - j])
					{
						cv->xtp[f1][f2] = 0.50 * (cv->xf[j][f1] + cv->xf[i][f2]);
						cv->ytp[f1][f2] = 0.50 * (cv->yf[j][f1] + cv->yf[i][f2]);
						number_f_fn_point[f1][f2][0] = i;
						number_f_fn_point[f2][f1][1] = j;
						;
						dx = cv->xtp[f1][f2] - xbase;
						dy = cv->ytp[f1][f2] - ybase;
						th = atan2(dy, dx);
						if (th < 0.0)
							th += 2.0 * R_PI;
						if (cv->theta[0] <= th && th <= cv->theta[1])
						{
							ytmp = cv->yc + 0.50 * cv->delta;
							xtmp = xbase + (dx / dy) * (ytmp - ybase);
						}
						else if (cv->theta[1] <= th && th <= cv->theta[2])
						{
							xtmp = cv->xc - 0.50 * cv->delta;
							ytmp = ybase + (dy / dx) * (xtmp - xbase);
						}
						else if (cv->theta[2] <= th && th <= cv->theta[3])
						{
							ytmp = cv->yc - 0.50 * cv->delta;
							xtmp = xbase + (dx / dy) * (ytmp - ybase);
						}
						else
						{
							xtmp = cv->xc + 0.50 * cv->delta;
							ytmp = ybase + (dy / dx) * (xtmp - xbase);
						}
						cv->xtp[f1][f2] = xtmp;
						cv->ytp[f1][f2] = ytmp;
						cv->xtp[f2][f1] = cv->xtp[f1][f2];
						cv->ytp[f2][f1] = cv->ytp[f1][f2];
						;
						i = 10;
						j = 10;
						;
						f1 = 10;
					}
				}
			}
		}
	}
	;
	for (f1 = 0; f1 < 3; f1++)
	{
		f2 = (f1 + 1) % 3;
		f3 = (f1 + 2) % 3;
		if (number_f_fn_point[f1][f2][0] != -1)
		{
			i = number_f_fn_point[f1][f2][0];
			j = number_f_fn_point[f2][f1][1];
			xtcc = (xbase - cv->xc) / cv->delta;
			ytcc = (ybase - cv->yc) / cv->delta;
			xtmpcc[0] = (cv->xtp[f1][f2] - cv->xc) / cv->delta;
			ytmpcc[0] = (cv->ytp[f1][f2] - cv->yc) / cv->delta;
			xtmpcc[1] = (cv->xf[j][f1] - cv->xc) / cv->delta;
			ytmpcc[1] = (cv->yf[j][f1] - cv->yc) / cv->delta;
			if (check_if_p2_is_CW_cell_coordinates(xtmpcc[0], ytmpcc[0], xtmpcc[1], ytmpcc[1]) == -1)
			{
				integration_by_rotating_around_triple_point(cv->f[f1], xtcc, ytcc, xtmpcc, ytmpcc, &area_tmp);
				;
				cv->area_integrated[f1] = area_tmp;
				cv->xtp[f1][f3] = cv->xc + xtmpcc[1] * cv->delta;
				cv->ytp[f1][f3] = cv->yc + ytmpcc[1] * cv->delta;
				cv->xtp[f3][f1] = cv->xtp[f1][f3];
				cv->ytp[f3][f1] = cv->ytp[f1][f3];
				;
				xtmpcc[0] = xtmpcc[1];
				ytmpcc[0] = ytmpcc[1];
				;
				integration_by_rotating_around_triple_point(cv->f[f3], xtcc, ytcc, xtmpcc, ytmpcc, &area_tmp);
				;
				cv->area_integrated[f3] = area_tmp;
				cv->xtp[f2][f3] = cv->xc + xtmpcc[1] * cv->delta;
				cv->ytp[f2][f3] = cv->yc + ytmpcc[1] * cv->delta;
				cv->xtp[f3][f2] = cv->xtp[f2][f3];
				cv->ytp[f3][f2] = cv->ytp[f2][f3];
				;
				cv->area_integrated[f2] = 1.0 - cv->area_integrated[f1] - cv->area_integrated[f3];
			}
			else
			{
				xtmpcc[0] = (cv->xtp[f1][f2] - cv->xc) / cv->delta;
				ytmpcc[0] = (cv->ytp[f1][f2] - cv->yc) / cv->delta;
				xtmpcc[1] = (cv->xf[i][f2] - cv->xc) / cv->delta;
				ytmpcc[1] = (cv->yf[i][f2] - cv->yc) / cv->delta;
				if (check_if_p2_is_CW_cell_coordinates(xtmpcc[0], ytmpcc[0], xtmpcc[1], ytmpcc[1]) == -1)
				{
					integration_by_rotating_around_triple_point(cv->f[f2], xtcc, ytcc, xtmpcc, ytmpcc, &area_tmp);
					;
					cv->area_integrated[f2] = area_tmp;
					cv->xtp[f2][f3] = cv->xc + xtmpcc[1] * cv->delta;
					cv->ytp[f2][f3] = cv->yc + ytmpcc[1] * cv->delta;
					cv->xtp[f3][f2] = cv->xtp[f2][f3];
					cv->ytp[f3][f2] = cv->ytp[f2][f3];
					;
					xtmpcc[0] = xtmpcc[1];
					ytmpcc[0] = ytmpcc[1];
					;
					integration_by_rotating_around_triple_point(cv->f[f3], xtcc, ytcc, xtmpcc, ytmpcc, &area_tmp);
					;
					cv->area_integrated[f3] = area_tmp;
					cv->xtp[f1][f3] = cv->xc + xtmpcc[1] * cv->delta;
					cv->ytp[f1][f3] = cv->yc + ytmpcc[1] * cv->delta;
					cv->xtp[f3][f1] = cv->xtp[f1][f3];
					cv->ytp[f3][f1] = cv->ytp[f1][f3];
					;
					cv->area_integrated[f1] = 1.0 - cv->area_integrated[f2] - cv->area_integrated[f3];
				}
				else
				{
					printf("ERRRRROOOORRRR\n");
				}
			}
			f1 = 10;
		}
	}
}

void integration_by_rotating_around_triple_point(double fraction, double xtcc, double ytcc, double xtmpcc[2], double ytmpcc[2], double* area)
{
	const int INTEGRATION_LEVEL = 5;
	int l, L;
	double area_tmp, area_tmp_last;
	;
	area_tmp = 0.0;
	area_tmp_last = 0.0;
	L = 100;
	;
	for (l = 1; l <= INTEGRATION_LEVEL; l++)
	{
		L *= 10;
		area_tmp -= area_tmp_last;
		while (area_tmp < fraction)
		{
			xtmpcc[1] = xtmpcc[0];
			ytmpcc[1] = ytmpcc[0];
			find_next_CCW_point_cell_coordinates(&(xtmpcc[0]), &(ytmpcc[0]), L);
			area_tmp_last = find_area_pt_and_two_border_points_cell_coordinates(xtcc, ytcc, xtmpcc, ytmpcc);
			area_tmp += area_tmp_last;
		}
	}
	*area = area_tmp;
}

double find_area_pt_and_two_border_points_cell_coordinates(double xpt, double ypt, double xb[2], double yb[2])
{
	double x[3], y[3], xtmp, ytmp, areatmp;
	if ((fabs(xb[0] - xb[1]) < R_VOFLIMIT) || (fabs(yb[0] - yb[1]) < R_VOFLIMIT))
	{
		x[0] = xpt;
		y[0] = ypt;
		x[1] = xb[0];
		y[1] = yb[0];
		x[2] = xb[1];
		y[2] = yb[1];
		return triangle_area(x, y);
	}
	xtmp = 0.0;
	ytmp = 0.0;
	if ((fabs(0.50 - xb[0]) < R_VOFLIMIT && fabs(0.50 - yb[1]) < R_VOFLIMIT) || (fabs(0.50 - xb[1]) < R_VOFLIMIT && fabs(0.50 - yb[0]) < R_VOFLIMIT))
	{
		xtmp = 0.50;
		ytmp = 0.50;
	}
	else if ((fabs(0.50 + xb[0]) < R_VOFLIMIT && fabs(0.50 - yb[1]) < R_VOFLIMIT) || (fabs(0.50 + xb[1]) < R_VOFLIMIT && fabs(0.50 - yb[0]) < R_VOFLIMIT))
	{
		xtmp = -0.50;
		ytmp = 0.50;
	}
	else if ((fabs(0.50 + xb[0]) < R_VOFLIMIT && fabs(0.50 + yb[1]) < R_VOFLIMIT) || (fabs(0.50 + xb[1]) < R_VOFLIMIT && fabs(0.50 + yb[0]) < R_VOFLIMIT))
	{
		xtmp = -0.50;
		ytmp = -0.50;
	}
	else if ((fabs(0.50 - xb[0]) < R_VOFLIMIT && fabs(0.50 + yb[1]) < R_VOFLIMIT) || (fabs(0.50 - xb[1]) < R_VOFLIMIT && fabs(0.50 + yb[0]) < R_VOFLIMIT))
	{
		xtmp = 0.50;
		ytmp = -0.50;
	}
	if (xtmp == 0.0 || ytmp == 0.0)
	{
		x[0] = xpt;
		y[0] = ypt;
		x[1] = xb[0];
		y[1] = yb[0];
		x[2] = xb[1];
		y[2] = yb[1];
		printf("AREAAAAAA\n");
		return triangle_area(x, y);
	}
	areatmp = 0.0;
	x[0] = xpt;
	y[0] = ypt;
	x[1] = xb[0];
	y[1] = yb[0];
	x[2] = xtmp;
	y[2] = ytmp;
	areatmp += triangle_area(x, y);
	x[0] = xpt;
	y[0] = ypt;
	x[1] = xb[1];
	y[1] = yb[1];
	x[2] = xtmp;
	y[2] = ytmp;
	areatmp += triangle_area(x, y);
	return areatmp;
}

double triangle_area(double x[3], double y[3])
{
	return 0.50 * fabs((x[0] * y[1] - x[1] * y[0]) + (x[1] * y[2] - x[2] * y[1]) + (x[2] * y[0] - x[0] * y[2]));
}

void find_next_CCW_point_cell_coordinates(double* xxx, double* yyy, int res)
{
	double step = 1.0 / (1.0 * res);
	if (fabs(0.50 - *xxx) < R_VOFLIMIT)
	{
		if (*yyy + step <= 0.50)
		{
			*yyy += step;
		}
		else
		{
			*xxx -= ((*yyy + step) - (+0.50));
			*yyy = +0.50;
		}
	}
	else if (fabs(0.50 - *yyy) < R_VOFLIMIT)
	{
		if (*xxx - step >= -0.50)
		{
			*xxx -= step;
		}
		else
		{
			*yyy += ((*xxx - step) - (-0.50));
			*xxx = -0.50;
		}
	}
	else if (fabs(0.50 + *xxx) < R_VOFLIMIT)
	{
		if (*yyy - step >= -0.50)
		{
			*yyy -= step;
		}
		else
		{
			*xxx -= ((*yyy - step) - (-0.50));
			*yyy = -0.50;
		}
	}
	else if (fabs(0.50 + *yyy) < R_VOFLIMIT)
	{
		if (*xxx + step <= 0.50)
		{
			*xxx += step;
		}
		else
		{
			*yyy += ((*xxx + step) - (+0.50));
			*xxx = +0.50;
		}
	}
}

int check_if_p2_is_CW_cell_coordinates(double xp1, double yp1, double xp2, double yp2)
{
	if (fabs(0.50 - xp1) < R_VOFLIMIT)
	{
		if (fabs(0.50 + xp2) < R_VOFLIMIT)
		{
			return 0;
		}
		else
		{
			if (yp2 > yp1)
			{
				return 1;
			}
			else
			{
				return -1;
			}
		}
	}
	else if (fabs(0.50 + xp1) < R_VOFLIMIT)
	{
		if (fabs(0.50 - xp2) < R_VOFLIMIT)
		{
			return 0;
		}
		else
		{
			if (yp2 < yp1)
			{
				return 1;
			}
			else
			{
				return -1;
			}
		}
	}
	else if (fabs(0.50 - yp1) < R_VOFLIMIT)
	{
		if (fabs(0.50 + yp2) < R_VOFLIMIT)
		{
			return 0;
		}
		else
		{
			if (xp2 < xp1)
			{
				return 1;
			}
			else
			{
				return -1;
			}
		}
	}
	else if (fabs(0.50 + yp1) < R_VOFLIMIT)
	{
		if (fabs(0.50 - yp2) < R_VOFLIMIT)
		{
			return 0;
		}
		else
		{
			if (xp2 > xp1)
			{
				return 1;
			}
			else
			{
				return -1;
			}
		}
	}
	return 0;
}

bool point_in_which_side_line(double xl, double yl, double nx, double ny, double xp, double yp)
{
	double vx, vy;
	vx = xp - xl;
	vy = yp - yl;
	if (vx * nx + vy * ny >= 0.0)
		return false;
	return true;
}

void find_four_corner_angles(struct CellValues* cv, double xbase, double ybase)
{
	double dx, dy;
	dx = cv->xc + 0.50 * cv->delta - xbase;
	dy = cv->yc + 0.50 * cv->delta - ybase;
	cv->theta[0] = atan2(dy, dx);
	if (cv->theta[0] < 0.0)
		cv->theta[0] += 2.0 * R_PI;
	dx = cv->xc - 0.50 * cv->delta - xbase;
	dy = cv->yc + 0.50 * cv->delta - ybase;
	cv->theta[1] = atan2(dy, dx);
	if (cv->theta[1] < 0.0)
		cv->theta[1] += 2.0 * R_PI;
	dx = cv->xc - 0.50 * cv->delta - xbase;
	dy = cv->yc - 0.50 * cv->delta - ybase;
	cv->theta[2] = atan2(dy, dx);
	if (cv->theta[2] < 0.0)
		cv->theta[2] += 2.0 * R_PI;
	dx = cv->xc + 0.50 * cv->delta - xbase;
	dy = cv->yc - 0.50 * cv->delta - ybase;
	cv->theta[3] = atan2(dy, dx);
	if (cv->theta[3] < 0.0)
		cv->theta[3] += 2.0 * R_PI;
}

void find_void_area(struct CellValues* cv, int res)
{
	int f;
	double x, y, ymax, ymin, cymax, cymin;
	double dx, dv, de, dxc, dyc, dvf, def, dxcf, dycf, wfmin, wfmax;
	;
	dv = 0.0;
	dxc = 0.0;
	dyc = 0.0;
	dx = cv->delta / (1.0 * res);
	ymin = cv->yc - 0.50 * cv->delta;
	ymax = cv->yc + 0.50 * cv->delta;
	cymin = ymin;
	cymax = ymax;
	;
	dvf = 0.0;
	dxcf = 0.0;
	dycf = 0.0;
	;
	for (x = cv->xc - 0.50 * cv->delta + 0.50 * dx; x < cv->xc + 0.50 * cv->delta; x += dx)
	{
		ymin = cv->yc - 0.50 * cv->delta;
		ymax = cv->yc + 0.50 * cv->delta;
		//wfmin = 3.0;
		//wfmax = 3.0;
		wfmin = 2.0 / 3.0;
		wfmax = 2.0 / 3.0;
		for (f = 0; f < 3; f++)
		{
			if (cv->nyf[f] < -R_VOFLIMIT)
			{
				y = cv->yf[0][f] + (cv->yf[1][f] - cv->yf[0][f]) / (cv->xf[1][f] - cv->xf[0][f]) * (x - cv->xf[0][f]);
				if (ymax > y)
				{
					ymax = y;
					//wfmax = 1.0 / cv->f[f];
					wfmax = 1.0 - cv->f[f];
					//wfmax = cv->f[f];
				}
			}
			else if (cv->nyf[f] > R_VOFLIMIT)
			{
				y = cv->yf[0][f] + (cv->yf[1][f] - cv->yf[0][f]) / (cv->xf[1][f] - cv->xf[0][f]) * (x - cv->xf[0][f]);
				if (ymin < y)
				{
					ymin = y;
					//wfmin = 1.0 / cv->f[f];
					wfmin = 1.0 - cv->f[f];
					//wfmin = cv->f[f];
				}
			}
			else
			{
				if ((x < 0.50 * (cv->xf[1][f] + cv->xf[0][f]) && cv->nxf[f] > R_VOFLIMIT)
					|| (x > 0.50 * (cv->xf[1][f] + cv->xf[0][f]) && cv->nxf[f] < -R_VOFLIMIT))
				{
					ymin = cymax;
					ymax = cymin;
				}
			}
		}
		if (ymax > ymin)
		{
			de = dx * (ymax - ymin);
			dv += de;
			dxc += x * de;
			dyc += 0.50 * (ymax + ymin) * de;
			;
			def = de * 0.50 * (wfmin + wfmax);
			dvf += def;
			dxcf += x * def;
			dycf += (ymax * wfmax + ymin * wfmin) / (wfmax + wfmin) * def;
		}
	}
	cv->vvar = dv;
	if (dv > R_VOFLIMIT || dv < -R_VOFLIMIT)
	{
		cv->xtar = dxc / dv;
		cv->ytar = dyc / dv;
	}
	;
	cv->vvfw = dvf;
	if (dvf > R_VOFLIMIT || dvf < -R_VOFLIMIT)
	{
		cv->xtfw = dxcf / dvf;
		cv->ytfw = dycf / dvf;
	}
}
