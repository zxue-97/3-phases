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


void output_tecplot2D_FI(char *name, scalar intrfc, double time, double rotationangle)
{
	int interfacepoints, cellnumber, iii;
	char typei[2], *typeintersect;
	double xi[2], yi[2], *xintersect, *yintersect;
	FILE *fp;
	const double angle = rotationangle * R_PI / 180.0;
	scalar alpha[];
	vector n[];
	;
	cellnumber = 0;
	foreach ()
		cellnumber++;
	xintersect = (double *)calloc(sizeof(double), cellnumber);
	yintersect = (double *)calloc(sizeof(double), cellnumber);
	typeintersect = (char *)calloc(sizeof(char), cellnumber);
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
	fp = fopen(name, "w");
	switch (interfacepoints)
	{
	case 0:
	{
		fprintf(fp,
				"\r\nvariables = x, y"
				"\r\nZONE T=InterfaceData"
				"\r\nN = 3, E = 1, DATAPACKING = BLOCK, ZONETYPE = FETRIANGLE"
				"\r\nSOLUTIONTIME = %.4e",
				time);
		foreach ()
		{
			fprintf(fp, "\r\n%e %e %e", x, x, x);
			fprintf(fp, "\r\n%e %e %e", y, y, y);
			fprintf(fp, "\r\n1 2 3");
			break;
		}
		break;
	}
	default:
	{
		fprintf(fp,
				"\r\nvariables = x, y"
				"\r\nZONE T=InterfaceData"
				"\r\nN = %d, E = %d, DATAPACKING = BLOCK, ZONETYPE = FETRIANGLE"
				"\r\nSOLUTIONTIME = %.4e",
				3 * interfacepoints / 2, interfacepoints / 2, time);
		for (iii = 0; iii < interfacepoints - 1; iii += 2)
			fprintf(fp, "\r\n%e %e %e", xintersect[iii], xintersect[iii + 1], xintersect[iii]);
		for (iii = 0; iii < interfacepoints - 1; iii += 2)
			fprintf(fp, "\r\n%e %e %e", yintersect[iii], yintersect[iii + 1], yintersect[iii]);
		for (iii = 1; iii < 3 * interfacepoints / 2; iii += 3)
			fprintf(fp, "\r\n%d %d %d", iii, iii + 1, iii + 2);
		break;
	}
	}
	fclose(fp);
	;
	free(xintersect);
	free(yintersect);
}