 #include "Headers/Include.h"

int main() 
{

	int type, pouch;
	int N, Nx, Ny;
	int lines_batt_mod, cols_batt_mod;
	int tr_cell;
	bool tr_active;
	double tr_q_dot;
	double Lx, Ly;
	double TotalTime, dt;
	double dx, dy;
	double D, h, t, l, w, t_fin, t_pcm;
	double Tinitial, Q;
	double Tn, Ts, Te, Tw;
	double qn, qs, qe, qw;
	double rho_bat, kx_bat, ky_bat, cp_bat;
	double rho_pcm, k_pcm, cp_pcm, L_pcm, Tmelt;
	double rho_cpcm, k_cpcm, cp_cpcm, L_cpcm;
	double rho_alu, k_alu, cp_alu;
	double rho_gra, kx_gra, ky_gra, cp_gra;
	double rho_cop, k_cop, cp_cop;
	double porosity;
	double h_cp, T_cp, h_air, T_air;
	double m_dot, cp_liq;
	double tr_time, tr_duration;
	double SOCinit, cell_capacity;
	vector<double> times;
	vector<double> q_dots;
	vector<double> current;
	
	read_data("Input.yaml", type, pouch, TotalTime, dt, dx, dy, D, h, t, l, w, t_fin, t_pcm, Tinitial, Q, Tn, Ts, Te, Tw, qn, qs, qe, qw, rho_bat, kx_bat, ky_bat, cp_bat, rho_pcm, k_pcm, cp_pcm, L_pcm, Tmelt, rho_alu, k_alu, cp_alu, rho_gra, kx_gra, ky_gra, cp_gra, rho_cop, k_cop, cp_cop, porosity, h_cp, T_cp, h_air, T_air, rho_cpcm, k_cpcm, cp_cpcm, L_cpcm, lines_batt_mod, cols_batt_mod, m_dot, cp_liq, times, q_dots, tr_active, tr_cell, tr_q_dot, tr_time, tr_duration, SOCinit, current, cell_capacity);
	
	// --------------------------- INIT SIMULATION DATA ------------------------- //
	set_mesh_problem(type, N, Nx, Ny, Lx, Ly, dx, dy, t, l, t_fin, t_pcm, lines_batt_mod, cols_batt_mod);
	string results = create_results_folder();

	// -------------- MEMORY ALLOCATION AND INITIALIZATION -------------- //
	int* R = (int*)malloc((N) * sizeof(int));
	int* pp = (int*)malloc((N) * sizeof(int));
	int* ee = (int*)malloc((N) * sizeof(int));
	int* ww = (int*)malloc((N) * sizeof(int));
	int* nn = (int*)malloc((N) * sizeof(int));
	int* ss = (int*)malloc((N) * sizeof(int));
	int* batt_pos = (int*)malloc((N) * sizeof(int));
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
	double* L = (double*)malloc((N) * sizeof(double));
	double* T_fluid = (double*)malloc((N) * sizeof(double));;
	
	for (int i = 0; i < N; i++) 
	{
		R[i] = 0;
		pp[i] = 0;
		ee[i] = 0;
		ww[i] = 0;
		nn[i] = 0;
		ss[i] = 0;
		batt_pos[i] = 0;
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
		L[i] = 0.0;
		T[i] = Tinitial;
		T_fluid[i] = T_cp;
	}
	
	// -------------------------- INITIALIZATION ----------------------- //
	double SOCi = 0.0;
	double SOC = SOCinit;
	double time = 0.0;				// Current time step [s] 
	double q_dot = 0.0;
	double i_t = 0.0;				// Cell current [A]
	double rref = 0.0;				// Cell internal resistance [miliOhm]
	int it = 0;						// Iterations counter
	int previous_time = 0;
	int print_every = 10;			// Print every X seconds
	bool print_now = true;
	bool cooling_active = false;	// Liquid cooling active
	bool tr_started = false;		// Thermal runaway event has started

	// Simulation timer
	clock_t start, end;

	// Map mesh
	int o = path(Nx, Ny, pp, ee, ww, nn, ss);
	map_mesh(type, o, Nx, Ny, dx, dy, D, t, l, t_fin, t_pcm, kx_bat, ky_bat, k_pcm, k_alu, k_cpcm, kx_gra, ky_gra, rho_bat, rho_pcm, rho_alu, rho_cpcm, rho_gra, cp_bat, cp_pcm, cp_alu, cp_cpcm, cp_gra, L_pcm, L_cpcm, pp, R, kx, ky, rho, cp, L, lines_batt_mod, cols_batt_mod, batt_pos);

	// Store initial mesh info
	plot_mesh(o, Nx, Ny, dx, dy, pp, ww, ee, nn, ss, R, kx, ky, rho, cp, L, batt_pos, results);
	plot_sim(o, Nx, Ny, dx, dy, pp, T, f, R, time, results);
	plot_coef(o, Nx, Ny, dx, dy, pp, aw, ae, an, as, ap, b, time, results);
	plot_log(time, o, pp, R, T, f, q_dot, SOC, i_t, rref, results);

	// Resíduos (SOR, liquid fraction and conductivity)
	double resmax = 1.0E-4;
	double resf = 1.0;

	// Start simulation!
	start = clock();		// Start timer!

	while (time < TotalTime)
	{

		if (print_now)
		{
			printf("\n// -------------------- Time Step = %5.3fs -------------------- // \n", time);
		}

		// ---------- Copy previous time step data and check liquid cooling condition ---------- //
		for (int i = 0; i < N; i++)
		{	
			Ti[i] = T[i];
			fi[i] = f[i];
			/*
			if (fi[i] > 0 && cooling_active == false && (type == 3 || type == 4))
			{
				cooling_active = true;
				printf("Liquid cooling activated at %5.3fs\n", time);
			}
			*/
		}
		SOCi = SOC;
		// ------------------------------------------------------------------------------------- //

		if ((type == 4 || type == 5) && cooling_active == true)
		{
			fluid_temp(o, pp, N, Ny, T, batt_pos, T_fluid, T_cp, h_cp, m_dot, cp_liq, cols_batt_mod);
		}

		// Init residue and number of iterations
		resf = 1.0;
		it = 0;

		while (resf > resmax)
		{
			// Get Heat generation for module simulation
			if (type == 4 || type == 5)
			{
				get_q_dot(time, times, q_dots, current, SOC, SOCi, cell_capacity, q_dot, i_t, rref, dt);
				
				if (time >= tr_time && time <= tr_time + tr_duration && tr_active == true && tr_started == false)
				{
					tr_started = true;
					printf("Thermal runaway started at %5.3fs\n", time);
				}
				else if((time < tr_time || time > (tr_time + tr_duration)) && tr_started == true)
				{
					tr_started = false;
					printf("Thermal runaway stoped at %5.3fs\n", time);
				}
			}
			
			// Coefficients matrix setup
			assembly(o, pp, type, Nx, Ny, dx, dy, ww, ee, nn, ss, kx, ky, ap, aw, ae, an, as, su, sp, b, Tw, Te, Tn, Ts, qw, qe, qn, qs, T, Ti, rho, cp, L, f, fi, R, Q, dt, h_cp, T_cp, h_air, T_air, w, cooling_active, T_fluid, q_dot, batt_pos, tr_started, tr_cell, tr_q_dot, tr_time, tr_duration);

			// Linear system solver
			SORt(o, pp, ww, ee, nn, ss, ap, aw, ae, an, as, b, T, Ti);
			
			// Compute liquid fraction
			SORf(o, pp, ww, ee, nn, ss, ap, aw, ae, an, as, b, T, Ti, R, f, fi, rho, L, dt, Tmelt, &resf);

			it++;
			if (print_now)
			{
				printf("Iteration = %i\t resf = %5.5f\t SOC = %5.5f\n", it, resf, SOC);
			}
		}

		time = time + dt;

		// Print time step info
		plot_log(time, o, pp, R, T, f, q_dot, SOC, i_t, rref, results);
		
		if ((previous_time != int(time)) && (int(time) % print_every == 0))
		{
			plot_sim(o, Nx, Ny, dx, dy, pp, T, f, R, time, results);
		}

		if ((previous_time != int(time)) && (int(time) % print_every == 0))
		{
			print_now = true;
			previous_time = int(time);
		}
		else
		{
			previous_time = int(time);
			print_now = false;
		}
	}

	end = clock(); // Stop timer!

	// -------------------------------------------------------------------//
	double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
	cout << "\nTime taken by program is : " << time_taken << " s " << endl;
	// -------------------------------------------------------------------//

	free(rho);
	free(pp);
	free(ee);
	free(ww);
	free(nn);
	free(ss);
	free(cp);
	free(ap);
	free(aw);
	free(ae);
	free(an);
	free(as);
	free(Ti);
	free(fi);
	free(sp);
	free(su);
	free(kx);
	free(ky);
	free(R);
	free(b);
	free(T);
	free(f);

	return 0;
}