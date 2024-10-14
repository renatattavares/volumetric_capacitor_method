#include "..\Headers\Include.h"

void map(int type, int o, int Nx, int Ny, double dx, double dy, double D, double t, double l, double t_fin, double t_pcm, double kx_bat, double ky_bat, double k_pcm, double k_alu, double k_cpcm, double kx_gra, double ky_gra, double rho_bat, double rho_pcm, double rho_alu,	double rho_cpcm, double rho_gra, double cp_bat, double cp_pcm, double cp_alu, double cp_cpcm,	double cp_gra, int* pp, int* R, double* kx, double* ky, double* rho, double* cp)
{
	// Simulation timer
	clock_t start, end;
	start = clock();		// Start timer!

	if (type == 1) // Cylindrical cell 
	{
		int a, i, l, m;
		double x_center, y_center, x, y, radius;
		double cell_radius = D / 2;
		double x_cell = Nx * dx / 2;
		double y_cell = Ny * dy / 2;

		for (a = 0; a < o; a++) 
		{
			m = pp[a];
			i = m / (Ny + 2);
			l = m - i * (Ny + 2);

			// Mesh Regions and Properties
			x_center = (dx * 0.5) + (i-1) * dx;
			y_center = (dy * 0.5) + (l-1) * dy;
			x = fabs(x_center - x_cell);
			y = fabs(y_center - y_cell);
			radius = sqrt(pow(x, 2) + pow(y, 2));
			
			if (radius <= cell_radius)
			{
				R[m] = 1; // Battery Region
				kx[m] = kx_bat;
				ky[m] = ky_bat;
				rho[m] = rho_bat;
				cp[m] = cp_bat;
			}
			else 
			{
				R[m] = 0; // PCM Region
				kx[m] = k_pcm;
				ky[m] = k_pcm;
				rho[m] = rho_pcm;
				cp[m] = cp_pcm;
			}
			//printf("CV = %4i\t x_center = %4.5f\t y_center = %4.5f\t radius = %4.5f\t R = %4i\n", m, x_center, y_center, radius, R[m]);
		}
	}

	else if(type == 2) // Pouch cell - Baseline design 
	{
		int i, l, a, m, b;
		double x_center;
		double x_batt[16];
		int batts = 8;
		int Nx_ghosts = Nx + 2;
		int	Ny_ghosts = Ny + 2;

		// Battery regions
		for (int i = 0; i < batts; i++)
		{
			x_batt[2 * i] = (t_fin) + i * (t_fin + t);
			x_batt[(2 * i) + 1] = (t_fin + t) + i * (t_fin + t);
		}

		for (b = 0; b < o; b++) 
		{
			m = pp[b];
			i = (int)m / Ny_ghosts;
			l = (int)m - i * Ny_ghosts;

			// Mesh Regions and Properties
			x_center = (dx * 0.5) + (i-1) * dx;

			for (a = 0; a < batts; a++)
			{
				if ((x_center > x_batt[2 * a]) && (x_center < x_batt[(2 * a) + 1]))
				{
					R[m] = 1; // Battery Region
					kx[m] = kx_bat;
					ky[m] = ky_bat;
					rho[m] = rho_bat;
					cp[m] = cp_bat;
					break;
				}
				else {
					R[m] = 3; // Aluminum Region
					kx[m] = k_alu;
					ky[m] = k_alu;
					rho[m] = rho_alu;
					cp[m] = cp_alu;
				}
			}
		}
	}

	else if (type == 3) 
	{
		int batts = 8;
		int fins = 9;
		int pcms = 16;
		int i, l, a, b, c, m, z;
		double x_center, y_center;
		double x_batt[16];
		double x_fin[18];
		double x_pcm1[32];
		double x_pcm2[32];
		int Nx_ghosts = Nx + 2;
		int	Ny_ghosts = Ny + 2;

		// Battery, fins and PCM regions
		for (int i = 0; i < batts; i++)
		{	
			x_batt[2 * i] = (t_fin + t_pcm) + i*(t_fin + 2 * t_pcm + t);
			x_batt[(2 * i) + 1] = (t_fin + t_pcm + t) + i * (t_fin + 2 * t_pcm + t);
		}
		for (int b = 0; b < fins; b++)
		{
			x_fin[2 * b] = b * (t_fin + 2 * t_pcm + t);
			x_fin[(2 * b) + 1] = (t_fin) + b * (t_fin + 2 * t_pcm + t);
		}
		for (int c = 0; c < pcms; c++)
		{
			
			x_pcm1[2 * c] = (t_fin) + c * (t + t_fin + 2 * t_pcm);
			x_pcm1[(2 * c) + 1] = (t_fin + t_pcm) + c * (t + t_fin + 2 * t_pcm);
			x_pcm2[2 * c] = (t_fin + t + t_pcm) + c * (t + t_fin + 2 * t_pcm);
			x_pcm2[(2 * c) + 1] = (t_fin + t + 2 * t_pcm) + c * (t + t_fin +  2 * t_pcm);
			
		}
		
		for (z = 0; z < o; z++)
		{
			m = pp[z];
			i = (int)m / Ny_ghosts;
			l = (int)m - i * Ny_ghosts;

			// Mesh Regions and Properties

			x_center = (dx * 0.5) + (i-1) * dx;

			for (a = 0; a < batts; a++)
			{
				if ((x_center > x_batt[2 * a]) && (x_center < x_batt[(2 * a) + 1]))
				{
					R[m] = 1; // Battery Region
					kx[m] = kx_bat;
					ky[m] = ky_bat;
					rho[m] = rho_bat;
					cp[m] = cp_bat;
					break;
				}
				else
				{
					for (b = 0; b < fins; b++)
					{
						if ((x_center > x_fin[2 * b]) && (x_center < x_fin[(2 * b) + 1]))
						{

							R[m] = 6; // Graphite Region
							kx[m] = kx_gra;
							ky[m] = ky_gra;
							rho[m] = rho_gra;
							cp[m] = cp_gra;
							break;
						}
						else
						{
							for (c = 0; c < pcms; c++)
							{
								if (((x_center > x_pcm1[2 * c]) && (x_center < x_pcm1[(2 * c) + 1])) || 
									((x_center > x_pcm2[2 * c]) && (x_center < x_pcm2[(2 * c) + 1])))
								{
									y_center = (dy * 0.5) + (l-1) * dy;

									if (y_center >= l/2) {
										R[m] = 5; // CPCM Region
										kx[m] = k_cpcm;
										ky[m] = k_cpcm;
										rho[m] = rho_cpcm;
										cp[m] = cp_cpcm;
										break;	
									}
									else
									{
										R[m] = 0; // PCM Region
										kx[m] = k_pcm;
										ky[m] = k_pcm;
										rho[m] = rho_pcm;
										cp[m] = cp_pcm;
										break;
									}
								}
							}
						}
					}
				}
			}
		}
	}

	end = clock(); // Stop timer!
	double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
	cout << "\nMap time: " << time_taken << " s " << endl;
}


void assembly(int o, int* pp, int type, int Nx, int Ny, double dx, double dy, int* ww, int* ee, int* nn, int* ss, double* kx, double* ky, double* ap, double* aw, double* ae, double* an, double* as, double* su, double* sp, double* b, double Tw, double Te, double Tn, double Ts, double qw, double qe, double qn, double qs, double* T, double* Ti, double* rho, double* cp,	double L_pcm, double L_cpcm, double* f, double* fi, int* R, double Q, double dt, double h_cp, double T_cp, double h_air, double T_air, double w, int active) 
{

	double L = 0.0;
	int i, j, m, l;
	
	for (m = 0; m < o; m++) 
	{
		L = 0.0;

		l = pp[m];
		i = l / (Ny + 2);
		j = l - i * (Ny + 2);

		su[l] = 0.0;
		sp[l] = 0.0;

		aw[l] = (2 * kx[ww[l]] * kx[l]) / (dx * (kx[ww[l]] * dx + kx[l] * dx));
		ae[l] = (2 * kx[ee[l]] * kx[l]) / (dx * (kx[ee[l]] * dx + kx[l] * dx));
		an[l] = (2 * ky[nn[l]] * ky[l]) / (dy * (ky[nn[l]] * dy + ky[l] * dy));
		as[l] = (2 * ky[ss[l]] * ky[l]) / (dy * (ky[ss[l]] * dy + ky[l] * dy));

		// Boundary conditions - West
		if (i == 1)
		{
			aw[l] = 0;
			su[l] = su[l] + (2 * h_air * T_air) / dx;
			sp[l] = sp[l] - (2 * h_air) / dx;
		}

		// Boundary conditions - East
		if (i == Nx)
		{
			ae[l] = 0;
			su[l] = su[l] + (2 * h_air * T_air) / dx;
			sp[l] = sp[l] - (2 * h_air) / dx;
		}

		// Boundary conditions - North
		if (j == 1)
		{
			an[l] = 0;
			su[l] = su[l] + (2 * h_air * T_air) / dy;
			sp[l] = sp[l] - (2 * h_air) / dy;

		}

		// Boundary conditions - South
		if (j == Ny)
		{
			as[l] = 0;
			if (type == 3 && active == 1)
			{
				su[l] = su[l] + (2 * h_cp * T_cp) / dy;
				sp[l] = sp[l] - (2 * h_cp) / dy;
			}
			else if (type == 2)
			{
				su[l] = su[l] + (2 * h_cp * T_cp) / dy;
				sp[l] = sp[l] - (2 * h_cp) / dy;
			}
		}

		//Boundary conditions - Source term
		if (R[l] == 1)
		{
			su[l] = su[l] + Q;
		}

		// Latent Heat of Fusion
		if (R[l] == 0)
		{
			L = L_pcm;
		}
		else if (R[l] == 5)
		{
			L = L_cpcm;
		}

		//Calculo do ap e do b
		ap[l] = aw[l] + ae[l] + an[l] + as[l] + ((rho[l] * cp[l]) / dt) - sp[l];
		b[l] = su[l] + ((rho[l] * cp[l] * Ti[l]) / dt) - ((rho[l] * L * (f[l] - fi[l])) / dt);
		//printf("pp = %i\t ww = %5.1i\t ee = %5.1i\t nn = %5.1i\t ss = %5.1i\n", l, ww[l], ee[l], nn[l], ss[l]);
		//printf("pp = %i\t aw = %5.1E\t ae = %5.1E\t an = %5.1E\t as = %5.1E\t i = %i\t j = %i\n", l, aw[l], ae[l], an[l], as[l], i, j);
		//printf("pp = %i\t aw = %5.1E\t ae = %5.1E\t an = %5.1E\t as = %5.1E\t ap = %5.1E\t b = %5.1E\n", l, aw[l], ae[l], an[l], as[l], ap[l], b[l]);
	}
}


void set_mesh_problem(int type, int& N, int& Nx, int& Ny, double& Lx, double& Ly,
	double dx, double dy, double t, double l, double t_fin, double t_pcm) {

	if (type == 1)
	{
		Lx = 0.05;			// Mesh size in x direcition [m] 
		Ly = 0.05;			// Mesh size in y direcition [m]
	}
	else if (type == 2) 
	{	
		int batts = 8;
		int fins = 9;
		Lx = batts * t + fins * t_fin;			// Mesh size in x direcition [m] 
		Ly = l;									// Mesh size in y direcition [m]
	}
	else if (type == 3)
	{
		int batts = 8;
		int fins = 9;
		int pcms = 16;
		Lx = batts * t + fins * t_fin + pcms * t_pcm;	// Mesh size in x direcition [m] 
		Ly = l;											// Mesh size in y direcition [m]
	}
	Nx = (Lx / dx);				// Volumes in x direction
	Ny = (Ly / dy);				// Volumes in y direction 
	N = (Nx + 2) * (Ny + 2);	// Total number of CVs. Include ghost elements in boundaries.

	cout << "Mesh length in x direction [m]: " << Lx << endl;
	cout << "Mesh length in y direction [m]: " << Ly << endl;
	cout << "CV dimension in x direction [m]: " << dx << endl;
	cout << "CV dimension in y direction [m]: " << dy << endl;
	cout << "Number of cells in x direction: " << Nx << endl;
	cout << "Number of cells in y direction: " << Ny << endl;
	cout << "Total number of true cells: " << Nx * Ny << endl;
	cout << "Total number of cells: " << N << endl;
}


int path(int Nx, int Ny, int* pp, int* ee, int* ww, int* nn, int* ss) {

	int i, j, l, o, m;

	o = 0;

	for (i = 1; i <= Nx; i++)
	{
		for (j = 1; j <= Ny; j++)
		{
			pp[o] = i * (Ny + 2) + j;
			//printf("o = %i\t pp = %i\t i = %i\t j = %i\n", o, pp[o], i, j);
			o++;
		}
	}

	for (m = 0; m < o; m++)
	{
		l = pp[m];
		i = l / (Ny + 2);
		j = l - i * (Ny + 2);
		ee[l] = (i + 1) * (Ny + 2) + j;
		ww[l] = (i - 1) * (Ny + 2) + j;
		nn[l] = i * (Ny + 2) + (j - 1);
		ss[l] = i * (Ny + 2) + (j + 1);
		//printf("pp = %i\t ww = %5.1i\t ee = %5.1i\t nn = %5.1i\t ss = %5.1i\n", l, ww[l], ee[l], nn[l], ss[l]);
	}
	return o;
}
