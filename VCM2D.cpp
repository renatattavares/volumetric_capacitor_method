 #include "Headers/Include.h"

int main() {

	// ---------------------- SIMULATION SETUP ------------------------- //
	// Problem type
	const int type = 1;				// Cylindrical = 1 / Pouch Baseline = 2 / Pouch proposed = 3
	const int solver = 2;

	// Time Setup
	double TotalTime = 1;			// Total simulation time [s]
	double dt = 1;				// Time step [s]

	// -------------------------- MESH SETUP --------------------------- //
	// Control Volumes
	double dx = 0.001;				// CV dimension in x direction[m]
	double dy = 0.001;				// CV dimension in y direction[m]
	int N, Nx, Ny;
	double Lx, Ly;

	// Battery Cell
	// Cylindrical cell (Top view) - Type 1
	double D = 0.018;				// Cell diameter [m]
	double h = 0.065;				// Cell height [m]
	 
	// Pouch cell - Baseline Design (Side view) - Type 2
	double t = 0.008;				// Cell thickness [m]
	double l = 0.188;				// Cell length [m]
	double w = 0.100;				// Cell width [m]
	double t_fin = 0.001;			// Aluminum fin thickness [m]

	// Pouch cell - Proposed Design (Side view) - Type 3
	double t_pcm = 0.001;			// PCM/CPCM thickness [m]

	set_mesh_problem(type, N, Nx, Ny, Lx, Ly, dx, dy, t, l, t_fin, t_pcm);

	// ----------------- BOUNDARY AND INITIAL CONDITIONS --------------- //
	// Initial Condition
	double Tinitial = 25;			// Temperature [°C] in time = 0s
	double Q = 200000;				// W/m³

	// Temperatures for boundary conditions
	double Tn = 0.0;				// North
	double Ts = 0.0;				// South 
	double Te = 0.0;				// East
	double Tw = 0.0;				// West

	// Heat flux for boundary conditions
	double qn = 0.0;				// North
	double qs = 0.0;				// South 
	double qe = 0.0;				// East
	double qw = 0.0;				// West

	// ------------------------ PHYSICAL PROPERTIES --------------------- //
	// Option 1: RT 18 HC [T_melt = 18°C, Latent heat of fusion = 250 kJ/kg, Cp = 2 kJ/kg. K, k = 0.2 W/kg.K, rho = 880 <-> 770 kg/m³ (solid/liquid)]
	// Option 2: RT 21 HC [T_melt = 21°C, Latent heat of fusion = 250 kJ/kg, Cp = 2 kJ/kg. K, k = 0.2 W/kg.K, rho = 880 <-> 770 kg/m³ (solid/liquid)]

	// Battery
	double rho_bat = 2250.0;
	double kx_bat = 0.97;
	double ky_bat = 26.57;
	double cp_bat = 1230.0;
	
	// PCM
	double rho_pcm = 776.0;
	double k_pcm = 0.425;
	double cp_pcm = 2150.0;
	double L_pcm = 247050.0;		
	double Tmelt = 36.1;			

	// Aluminum
	double rho_alu = 2719.0;	
	double k_alu = 202.4;
	double cp_alu = 871.0;

	// Graphite
	double rho_gra = 2235.0;
	double kx_gra = 6.5;
	double ky_gra = 1659.3;
	double cp_gra = 762.5;

	// CPCM - Copper Foam
	double porosity = 0.92;			// Copper foam
	double rho_cop = 8978.0;
	double k_cop = 387.6;
	double cp_cop = 381.0;
	
	double rho_cpcm = porosity * rho_pcm + (1 - porosity) * rho_cop;
	double k_cpcm = porosity * k_pcm + (1 - porosity) * k_cop;
	double cp_cpcm = porosity * cp_pcm + (1 - porosity) * cp_cop;
	double L_cpcm = porosity * L_pcm;

	// Cold Plate
	double h_cp = 200;
	double T_cp = 25;

	// Natural Convection
	double h_air = 2;
	double T_air = 25;

	// --------------------------- RESULTS DATA ------------------------- //
	string results = create_results_folder();

	// -------------- MEMORY ALLOCATION AND INITIALIZATION -------------- //
	int* R = (int*)malloc((N) * sizeof(int));
	int* pp = (int*)malloc((N) * sizeof(int));
	int* ee = (int*)malloc((N) * sizeof(int));
	int* ww = (int*)malloc((N) * sizeof(int));
	int* nn = (int*)malloc((N) * sizeof(int));
	int* ss = (int*)malloc((N) * sizeof(int));
	int* ii = (int*)malloc((N) * sizeof(int));
	int* jj = (int*)malloc((N) * sizeof(int));
	double* rho = (double*)malloc((N) * sizeof(double));
	double* cp = (double*)malloc((N) * sizeof(double));
	double* ap = (double*)malloc((N) * sizeof(double));
	double* aw = (double*)malloc((N) * sizeof(double));
	double* ae = (double*)malloc((N) * sizeof(double));
	double* an = (double*)malloc((N) * sizeof(double));
	double* as = (double*)malloc((N) * sizeof(double));
	double* Ti = (double*)malloc((N) * sizeof(double));
	double* fi = (double*)malloc((N) * sizeof(double));
	double* sp = (double*)malloc((N) * sizeof(double));
	double* su = (double*)malloc((N) * sizeof(double));
	double* kx = (double*)malloc((N) * sizeof(double));
	double* ky = (double*)malloc((N) * sizeof(double));
	double* f = (double*)malloc((N) * sizeof(double));
	double* b = (double*)malloc((N) * sizeof(double));
	double* T = (double*)malloc((N) * sizeof(double));
	double* pp2 = (double*)malloc((N) * sizeof(double));
	double* ee2 = (double*)malloc((N) * sizeof(double));
	double* ww2 = (double*)malloc((N) * sizeof(double));
	double* nn2 = (double*)malloc((N) * sizeof(double));
	double* ss2 = (double*)malloc((N) * sizeof(double));

	for (int i = 0; i < N; i++) 
	{
		R[i] = 0;
		pp[i] = 0;
		ee[i] = 0;
		ww[i] = 0;
		nn[i] = 0;
		ss[i] = 0;
		pp2[i] = 0.0;
		ee2[i] = 0.0;
		ww2[i] = 0.0;
		nn2[i] = 0.0;
		ss2[i] = 0.0;
		rho[i] = 0.0;
		cp[i] = 0.0;
		ap[i] = 0.0;
		aw[i] = 0.0;
		ae[i] = 0.0;
		an[i] = 0.0;
		as[i] = 0.0;
		Ti[i] = 0.0;
		fi[i] = 0.0;
		sp[i] = 0.0;
		su[i] = 0.0;
		kx[i] = 0.0;
		ky[i] = 0.0;
		f[i] = 0.0;
		b[i] = 0.0;
		T[i] = Tinitial;
	}	

	// -------------------------- INITIALIZATION ----------------------- //
	int i = 0;						// Counter
	int it = 0;						// Iterations counter
	int active = 0;					// Active liquid cooling
	double time = 0.0;				// Current time step [s] 
	double resf = 1.0;				// Initializing residue
	double resmax = 1.0E-4;			// Maximum property residue 
	int previous_time = 0;

	// Simulation timer
	clock_t start, end;

	// Start simulation!
	start = clock();		// Start timer!
	
	// Map mesh
	int o = path(Nx, Ny, pp, ee, ww, nn, ss);
	map(type, o, Nx, Ny, dx, dy, D, t, l, t_fin, t_pcm, kx_bat, ky_bat, k_pcm, k_alu,
		k_cpcm,	kx_gra, ky_gra, rho_bat, rho_pcm, rho_alu, rho_cpcm, rho_gra, cp_bat,
		cp_pcm, cp_alu, cp_cpcm, cp_gra, pp, R, kx, ky, rho, cp);

	for (int i = 0; i < N; i++)
	{
		ww2[i] = 1.0 * ww[i];
		ee2[i] = 1.0 * ee[i];
		nn2[i] = 1.0 * nn[i];
		ss2[i] = 1.0 * ss[i];
		pp2[i] = 1.0 * pp[i];
	}
	
	// Store simulation info
	plot_properties(o, Nx, Ny, dx, dy, pp, R, kx, ky, rho, cp, results);
	plot(o, Nx, Ny, dx, dy, pp, T, time, results, "Temp");
	plot(o, Nx, Ny, dx, dy, pp, f, time, results, "Frac");
	//log(time, o, pp, R, T, f, results);

	while (time < TotalTime) {

		printf("\n// -------------------- Time Step = %5.3fs -------------------- // \n", time);

		// Copy previous time step data and check liquid cooling condition
		for (int i = 0; i < N; i++)
		{
			Ti[i] = T[i];
			fi[i] = f[i];

			if (fi[i] > 0 && active == 0 && type == 3)
			{
				active = 1;
				printf("Liquid cooling activated\n", time);
				break;
			}
		}

		it = 0;
		resf = 1.0;

		while (resf > resmax) {

			assembly(o, pp, type, Nx, Ny, dx, dy, ww, ee, nn, ss, kx, ky, ap, aw, ae, an, as, su, sp, b, Tw, Te, Tn, Ts, qw, qe, qn, qs, T, Ti, rho, cp, L_pcm, L_cpcm, f, fi, R, Q, dt, h_cp, T_cp, h_air, T_air, w, active);

			plot_coef(o, Nx, Ny, dx, dy, pp, aw, ae, an, as, ap, b, time, results);

			if (solver == 1)
			{
				SORt(o, N, pp, ww, ee, nn, ss, ap, aw, ae, an, as, b, T, Ti, Nx, Ny, dx, dy, results);
			}
			else
			{
				SIP(o, N, pp, nn, ss, ee, ww, ap, ae, aw, an, as, b, T, Nx, Ny, dx, dy, results);
			}

			resf = 0.0;

			//SORf(o, pp, ww, ee, nn, ss, ap, aw, ae, an, as, b, T, Ti, R, f, fi, rho, dt, L_pcm, L_cpcm, Tmelt, &resf);
			//it++;
			//printf("\nIteracao = %i\t Residuo = %5.1E\n", it, resf);
		}

		time = time + dt;

		// Store time step info
		//log(time, o, pp, R, T, f, results);

		// Print time step temperature and liquid fraction
		if ((int(time) % 10) == 0 && int(time) != 0 && int(time) != previous_time) 
		{
			//previous_time = int(time);
			//plot(o, Nx, Ny, dx, dy, pp, T, time, results, "Temp");
			//plot(o, Nx, Ny, dx, dy, pp, f, time, results, "Frac");
		}
	}

	end = clock(); // Stop timer!

	// -------------------------------------------------------------------//
	double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
	cout << "\nTime taken by program is : " << time_taken << " s " << endl;
	// -------------------------------------------------------------------//

	average_temperature(o, pp, T);

	return 0;

}
