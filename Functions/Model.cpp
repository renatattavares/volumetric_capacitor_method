#include "..\Headers\Include.h"

void map(int type, int o, int Nx, int Ny, double dx, double dy, double D, double t, double l,
	double t_fin, double t_pcm, double kx_bat, double ky_bat, double k_pcm, double k_alu,	
	double k_cpcm, double kx_gra, double ky_gra, double rho_bat, double rho_pcm, double rho_alu,
	double rho_cpcm, double rho_gra, double cp_bat, double cp_pcm, double cp_alu, double cp_cpcm,
	double cp_gra, int* pp, int* R, double* kx, double* ky, double* rho, double* cp)
{
	// Simulation timer
	clock_t start, end;
	start = clock();		// Start timer!

	if (type == 1) // Cylindrical cell 
	{
		int a, i, j, m;
		double x_center, y_center, x, y, radius;
		double cell_radius = D / 2;
		double x_cell = Nx * dx / 2;
		double y_cell = Ny * dy / 2;
		int Nx_ghosts = Nx + 2;
		int	Ny_ghosts = Ny + 2;

		for (a = 0; a < o; a++) 
		{
			m = pp[a];
			i = (int)m / Ny_ghosts;
			j = (int)m - i * Ny_ghosts;

			// Mesh Regions and Properties
			x_center = (dx * 0.5) + (i-1) * dx;
			y_center = (dy * 0.5) + (j-1) * dy;
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
		int i, j, a, m, b;
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
			j = (int)m - i * Ny_ghosts;

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
		int i, j, a, b, c, m, z;
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
			j = (int)m - i * Ny_ghosts;

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
									y_center = (dy * 0.5) + (j-1) * dy;

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


void assembly(int o, int* pp, int type, int Nx, int Ny, double dx, double dy, int* ww, int* ee,
	int* nn, int* ss, double* kx, double* ky, double* ap, double* aw, double* ae, double* an,
	double* as, double* su, double* sp, double* b, double Tw, double Te, double Tn, double Ts, 
	double qw, double qe, double qn, double qs, double* T, double* Ti, double* rho, double* cp,
	double L_pcm, double L_cpcm, double* f, double* fi, int* R, double Q, double dt, double h_cp,
	double T_cp, double h_air, double T_air, double w, int active) 
{

	double L = 0.0;
	int i, j, x, y;
	int Nx_ghosts = Nx + 2;
	int	Ny_ghosts = Ny + 2;
	
	for (i = 0; i < o; i++) 
	{
		if (pp[i] != 0)
		{
			L = 0.0;

			j = pp[i];
			x = (int)j / Ny_ghosts;
			y = (int)j - x * Ny_ghosts;

			su[j] = 0.0;
			sp[j] = 0.0;

			aw[j] = (2 * kx[ww[j]] * kx[j]) / (dx * (kx[ww[j]] * dx + kx[j] * dx));
			ae[j] = (2 * kx[ee[j]] * kx[j]) / (dx * (kx[ee[j]] * dx + kx[j] * dx));
			an[j] = (2 * ky[nn[j]] * ky[j]) / (dy * (ky[nn[j]] * dy + ky[j] * dy));
			as[j] = (2 * ky[ss[j]] * ky[j]) / (dy * (ky[ss[j]] * dy + ky[j] * dy));

			// Boundary conditions - West
			if (x == 1)
			{
				aw[j] = 0;
				su[j] = su[j] + (2 * h_air * T_air) / dx;
				sp[j] = sp[j] - (2 * h_air) / dx;
			}

			// Boundary conditions - East
			if (x == Nx)
			{
				ae[j] = 0;
				su[j] = su[j] + (2 * h_air * T_air) / dx;
				sp[j] = sp[j] - (2 * h_air) / dx;
			}

			// Boundary conditions - North
			if (y == Ny)
			{
				an[j] = 0;
				su[j] = su[j] + (2 * h_air * T_air) / dy;
				sp[j] = sp[j] - (2 * h_air) / dy;

			}

			// Boundary conditions - South
			if (y == 1)
			{
				as[j] = 0;
				if (type == 3 && active == 1)
				{
					su[j] = su[j] + (2 * h_cp * T_cp) / dy;
					sp[j] = sp[j] - (2 * h_cp) / dy;
				}
				else if (type == 2)
				{
					su[j] = su[j] + (2 * h_cp * T_cp) / dy;
					sp[j] = sp[j] - (2 * h_cp) / dy;
				}
			}

			//Boundary conditions - Source term
			if (R[j] == 1)
			{
				su[j] = su[j] + Q;
			}

			// Latent Heat of Fusion
			if (R[j] == 0)
			{
				L = L_pcm;
			}
			else if (R[j] == 5)
			{
				L = L_cpcm;
			}

			//Calculo do ap e do b
			ap[j] = aw[j] + ae[j] + an[j] + as[j] + ((rho[j] * cp[j]) / dt) - sp[j];
			b[j] = su[j] + ((rho[j] * cp[j] * Ti[j]) / dt) - ((rho[j] * L * (f[j] - fi[j])) / dt);
			//printf("\nCV = %i\t ww = %5.1i\t ee = %5.1i\t nn = %5.1i\t ss = %5.1i\n", i, ww[i], ee[i], nn[i], ss[i]);
			//printf("CV = %i\t aw = %5.1E\t ae = %5.1E\t an = %5.1E\t as = %5.1E\t ap = %5.1E\t b = %5.1E\n", j, aw[j], ae[j], an[j], as[j], ap[j], b[j]);
		}
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
	int Nx_ghosts = Nx + 2;
	int	Ny_ghosts = Ny + 2;

	for (i = 1; i <= Nx; i++) 
	{
		for (j = 1; j <= Ny; j++)
		{
			pp[o] = i * Ny_ghosts + j;
			o++;
		}
	}

	for (m = 0; m < o; m++)
	{
		if (pp[m] != 0)
		{
			l = pp[m];
			i = (int)l / Ny_ghosts;
			j = (int)l - i * Ny_ghosts;

			ee[l] = (i + 1) * Ny_ghosts + j;
			ww[l] = (i - 1) * Ny_ghosts + j;
			nn[l] = i * Ny_ghosts + (j + 1);
			ss[l] = i * Ny_ghosts + (j - 1);

		}
	}
	return o;
}
