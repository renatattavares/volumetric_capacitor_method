#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>

using namespace std;
#pragma warning(disable:4996)

// ----------------------- FUNCTIONS DECLARATIONS -------------------------//
void map(int Nx, int Ny, double D, double dx, double dy, double k_bat, double k_pcm, double rho_bat, double rho_pcm, double cp_bat, double cp_pcm, int* ww, int* ee, int* nn, int* ss, int* R, double* k, double* rho, double* cp);
void assembly(int N, int Nx, int Ny, double dx, double dy, int* ww, int* ee, int* nn, int* ss, double* k, double* ap, double* aw, double* ae, double* an, double* as, double* s, double* sp, double* b, double Tw, double Te, double Tn, double Ts, double qw, double qe, double qn, double qs, double* T, double* Ti, double* rho, double* cp, double L_pcm, double* f, double* fi, int* R, double Q, double dt);
void interface_k(int p, int w, int e, int n, int s, double* k, double& kp, double& kw, double& ke, double& kn, double& ks);
void ft(int N, double Tmelt, double* T, double* f, double* fi, int* R);
void kt(int N, double* T, double k_bat, double k_pcm, double* k, double* ki, int* R);
void plot_f(int N, int Nx, int Ny, double dx, double dy, char ffn[20], double* f, int it);
void SORt(int N, int* ww, int* ee, int* nn, int* ss, double* ap, double* aw, double* ae, double* an, double* as, double* b, double* T, double* Ti);
void plot_temperature(int N, int Nx, int Ny, double dx, double dy, char ffn[20], double* T, int it);
void plot_properties(int N, int Nx, int Ny, double dx, double dy, char ffn[20], int* R, double* Ti, double* ki, double* rho, double* cp);
void average_temperature(int N, double* T);

int main() {

	// ---------------------------- DIMENSIONS ---------------------------//
	const int Nx = 200;  // Volumes in x direction
	const int Ny = 200;  // Volumes in y direction 
	double Lx = 0.05;    // Mesh size in x direcition [m] 
	double Ly = 0.05;    // Mesh size in y direcition [m] 
	double dx = Lx / Nx; // CV dimension in x direction [m]
	double dy = Ly / Ny; // CV dimension in y direction [m]
	const int N = (int)Nx * Ny; // Total number of CVs

	// Battery Cell
	double D = 0.018; // Battery cell diameter [m]
	double h = 0.065; // Battery cell height [m]

	// ----------------- BOUNDARY AND INITIAL CONDITIONS ----------------//
	// Initial Condition
	double Tinitial = 20.0; // Temperature in time = 0s [°C]
	double Q = 340074.6;

	// Temperatures for boundary conditions
	double Tn = 0.0; // North
	double Ts = 0.0; // South 
	double Te = 0.0; // East
	double Tw = 0.0; // West

	// Heat flux for boundary conditions
	double qn = 0.0; // North
	double qs = 0.0; // South 
	double qe = 0.0; // East
	double qw = 0.0; // West

	// ------------------------ PHYSICAL PROPERTIES ----------------------//
	// Option 1: RT 18 HC [T_melt = 18°C, Latent heat of fusion = 250 kJ/kg, Cp = 2 kJ/kg. K, k = 0.2 W/kg.K, rho = 880 <-> 770 kg/m³ (solid/liquid)]
	// Option 2: RT 21 HC [T_melt = 21°C, Latent heat of fusion = 250 kJ/kg, Cp = 2 kJ/kg. K, k = 0.2 W/kg.K, rho = 880 <-> 770 kg/m³ (solid/liquid)]

	double rho_bat = 2939;
	double rho_pcm = 825;   // Average value
	double k_bat = 3.4;
	double k_pcm = 0.2;
	double cp_bat = 1280;
	double cp_pcm = 2000;
	double L_pcm = 190000; // PCM latent heat of fusion [J/kg]
	double L_bat = 0.0;    // Battery cell latent heat of fusion [J/kg] -> It doesn't melt!
	double Tmelt = 22;     // PCM fusion temperature

	// ----------------------- MEMORY ALLOCATION ------------------------//
	double* kres = (double*)malloc((N) * sizeof(double));
	double* fres = (double*)malloc((N) * sizeof(double));
	double* rho = (double*)malloc((N) * sizeof(double));
	double* cp = (double*)malloc((N) * sizeof(double));
	double* ap = (double*)malloc((N) * sizeof(double));
	double* aw = (double*)malloc((N) * sizeof(double));
	double* ae = (double*)malloc((N) * sizeof(double));
	double* an = (double*)malloc((N) * sizeof(double));
	double* as = (double*)malloc((N) * sizeof(double));
	double* Ti = (double*)malloc((N) * sizeof(double));
	double* fi = (double*)malloc((N) * sizeof(double));
	double* ki = (double*)malloc((N) * sizeof(double));
	double* sp = (double*)malloc((N) * sizeof(double));
	double* k = (double*)malloc((N) * sizeof(double));
	double* f = (double*)malloc((N) * sizeof(double));
	double* b = (double*)malloc((N) * sizeof(double));
	double* s = (double*)malloc((N) * sizeof(double));
	double* T = (double*)malloc((N) * sizeof(double));

	int* R = (int*)malloc((N) * sizeof(int));
	int* ee = (int*)malloc((N) * sizeof(int));
	int* ww = (int*)malloc((N) * sizeof(int));
	int* nn = (int*)malloc((N) * sizeof(int));
	int* ss = (int*)malloc((N) * sizeof(int));

	for (int i = 0; i < N; i++) {

		// Floats
		kres[i] = 0.0;
		fres[i] = 0.0;
		ap[i] = 0.0;
		ww[i] = 0;
		ee[i] = 0;
		nn[i] = 0;
		ss[i] = 0;
		ap[i] = 0.0;
		aw[i] = 0.0;
		ae[i] = 0.0;
		an[i] = 0.0;
		as[i] = 0.0;
		Ti[i] = Tinitial;
		fi[i] = 0.0;
		ki[i] = 0.0;
		f[i] = 0.0;
		k[i] = 0.0;
		T[i] = 0.0;
		b[i] = 0.0;
		s[i] = 0.0;
		T[i] = 0.0;
		R[i] = 0;

	}

	// Miscellaneous
	char ffn[40];
	clock_t start, end;

	// -------------------------- SIMULATION CODE ------------------------//
	int i = 0;
	int it = 0;
	double C = 0.0;			// SOR property term
	double Y = 0.0;			// SOR property term
	double dt = 1;			// Time step [s] 
	double time = dt;		// Current time step [s] 
	double res = 1.0;		// Initializing residue
	double kappa = 1000;    // Over relaxation factor
	double resmax = 1E-3;	// Maximum property residue 
	double TotalTime = 600; // Total simulation time [s]
	
	start = clock(); // Start timer!

	// Initialization
	map(Nx, Ny, D, dx, dy, k_bat, k_pcm, rho_bat, rho_pcm, cp_bat, cp_pcm, ww, ee, nn, ss, R, k, rho, cp);
	plot_properties(N, Nx, Ny, dx, dy, ffn, R, Ti, k, rho, cp);

	/*for (int i = 0; i < N; i++) {

		T[i] = Ti[i];
		f[i] = fi[i];

	}

	//plot_temperature(N, Nx, Ny, dx, dy, ffn, T, -1);
	//plot_f(N, Nx, Ny, dx, dy, ffn, f, -1);
	*/

	while (time <= TotalTime) {

		printf("// -------------------- Time Step = %5.1fs -------------------- // \n", time);

		while (res > resmax) {

			res = 0.0;

			for (i = 0; i < N; i++) {

				if (R[i] == 1) {

					f[i] = 0.0;
					res = 0.0;
					//printf("Iteracao = %i\t Region = %i \t Residuo = %5.1E\n", it, R[i], res);
					continue;

				}

				fi[i] = f[i];
				Y = aw[i] * T[ww[i]] + ae[i] * T[ee[i]] + an[i] * T[nn[i]] + as[i] * T[ss[i]];
				C = b[i] + (rho[i] * L_pcm * fi[i]) / dt;
				f[i] = (fi[i] + kappa * (dt / (rho[i] * L_pcm)) * (ap[i] * Tmelt - Y - C)) / (1 + kappa);

				res = abs(fi[i] - f[i]);
				//printf("Iteracao = %i\t Region = %i \t Residuo = %5.1E\n", it, R[i], res);

			}

			it++;

		}
		
		assembly(N, Nx, Ny, dx, dy, ww, ee, nn, ss, k, ap, aw, ae, an, as, s, sp, b, Tw, Te, Tn, Ts, qw, qe, qn, qs, T, Ti, rho, cp, L_pcm, f, fi, R, Q, dt);
		SORt(N, ww, ee, nn, ss, ap, aw, ae, an, as, b, T, Ti);
		ft(N, Tmelt, T, f, fi, R);

		time = time + dt;
		res = 1.0;
		it = 0;
		
		if ((int(time) % 30) == 0) { // Print temperature and liquid fraction each 30 s (simulation time)

			plot_temperature(N, Nx, Ny, dx, dy, ffn, T, time / dt);
			plot_f(N, Nx, Ny, dx, dy, ffn, f, time / dt);

		}
				
	}

	end = clock(); // Stop timer!

	// -------------------------------------------------------------------//
	double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
	cout << "\nTime taken by program is : " << time_taken << " s " << endl;
	// -------------------------------------------------------------------//

	average_temperature(N, T);
	return 0;

}

void map(int Nx, int Ny, double D, double dx, double dy, double k_bat, double k_pcm, double rho_bat, double rho_pcm, double cp_bat, double cp_pcm, 
	int* ww, int* ee, int* nn, int* ss, int* R, double* k, double* rho, double* cp) {

	int i, j, o = 0;
	double x_center, y_center, x, y, radius;
	double cell_radius = D / 2;
	double x_cell = (Nx * dx) / 2;
	double y_cell = (Ny * dy) / 2;

	for (i = 0; i < Nx; i++) {
		for (j = 0; j < Ny; j++) {

			// Mesh Navigation
			o = (i * Ny) + j;
			ww[o] = ((i - 1) * Ny) + j; // Id - West neighbor
			ee[o] = ((i + 1) * Ny) + j; // Id - East neighbor
			nn[o] = (i * Ny) + (j + 1); // Id - North neighbor
			ss[o] = (i * Ny) + (j - 1); // Id - South neighbor

			if (i == 0) {
				ww[o] = Nx * Ny;
			}
			if (i == Nx - 1) {
				ee[o] = Nx * Ny;
			}
			if (j == 0) {
				ss[o] = Nx * Ny;
			}
			if (j == Ny - 1) {
				nn[o] = Nx * Ny;
			}

			// Mesh Regions and Properties
			x_center = (dx * 0.5) + i * dx;
			y_center = (dy * 0.5) + j * dx;
			x = x_center - x_cell;
			y = y_center - y_cell;
			radius = sqrt(pow(x, 2) + pow(y, 2));

			if (radius < cell_radius) {

				R[o] = 1; // Battery Region
				k[o] = k_bat;
				rho[o] = rho_bat;
				cp[o] = cp_bat;

			}
			else {

				R[o] = 0; // PCM Region
				k[o] = k_pcm;
				rho[o] = rho_pcm;
				cp[o] = cp_pcm;

			}

		}
	}

}

void assembly(int N, int Nx, int Ny, double dx, double dy, int* ww, int* ee, int* nn, int* ss, double* k, double* ap, double* aw, double* ae, double* an,
	double* as, double* s, double* sp, double* b, double Tw, double Te, double Tn, double Ts, double qw, double qe, double qn, double qs, double* T, double* Ti,
	double* rho, double* cp, double L_pcm, double* f, double* fi, int* R, double Q, double dt) {

	int o = 0;
	int i, x, y;
	double L = 0.0;
	double kw = 0.0;
	double ke = 0.0;
	double kn = 0.0;
	double ks = 0.0;
	double kp = 0.0;

	for (i = 0; i < N; i++) {

		interface_k(i, ww[i], ee[i], nn[i], ss[i], k, kp, kw, ke, kn, ks);
		
		L = 0.0;
		s[i] = 0.0;
		sp[i] = 0.0;
		ae[i] = (2 * ke * kp) / (dx * (ke * dx + kp * dx));
		aw[i] = (2 * kw * kp) / (dx * (kw * dx + kp * dx));
		an[i] = (2 * kn * kp) / (dy * (kn * dy + kp * dy));
		as[i] = (2 * ks * kp) / (dy * (ks * dy + kp * dy));
		

		x = (int)i / Ny;
		y = (int)i - x * Ny;

		// Boundary conditions - West
		if (x == 0) {
			aw[i] = 0;
		}

		// Boundary conditions - East
		if (x == Nx - 1) {
			ae[i] = 0;
			//sp[i] = -(2 * k[i]) / (dx * dx);
			//s[i] = ((2 * k[i]) / (dx * dx)) * Te;
		}

		// Boundary conditions - North
		if (y == Ny - 1) {
			an[i] = 0;
		}

		// Boundary conditions - South
		if (y == 0) {
			as[i] = 0;
		}

		//Boundary conditions - Source term
		if (R[i] == 1) {

			s[i] = Q;
			L = L_pcm;

		}

		//Calculo do ap e do b
		ap[i] = aw[i] + ae[i] + an[i] + as[i] + ((rho[i] * cp[i]) / dt) - sp[i];
		b[i] = s[i] + ((rho[i] * cp[i] * Ti[i]) / dt) - ((rho[i] * cp[i] * L * (f[i] - fi[i])) / dt);

		//printf("CV = %i\t aw = %5.1E\t ae = %5.1E\t an = %5.1E\t as = %5.1E\t ap = %5.1E\t b = %5.1E\n", i, aw[i], ae[i], an[i], as[i], ap[i], b[i]);
		//printf("CV = %i\t kw = %5.1E\t ke = %5.1E\t kn = %5.1E\t ks = %5.1E\t kp = %5.1E\n", i, kw, ke, kn, ks, kp);

	}
}

void SORt(int N, int* ww, int* ee, int* nn, int* ss, double* ap, double* aw, double* ae, double* an, double* as, double* b, double* T, double* Ti) {

	int i, it = 0;
	double R = 0.0;
	double res = 1.0;
	double kappa = 1000;
	double resmax = 1E-6;

	while (res > resmax) {

		res = 0.0;

		for (i = 0; i < N; i++) {

			Ti[i] = T[i];
			R = b[i] + aw[i] * T[ww[i]] + ae[i] * T[ee[i]] + an[i] * T[nn[i]] + as[i] * T[ss[i]];
			T[i] = (T[i] + kappa * (R / ap[i])) / (1 + kappa);

			res = res + pow((Ti[i] - T[i]), 2);

		}

		//printf("Iteracao = %i\t Residuo = %13.5E\n", it, res);
		it++;
	}

}

void interface_k(int p, int w, int e, int n, int s, double* k, double& kp, double& kw, double& ke, double& kn, double& ks) {

	kp = k[p];
	kw = (k[p] + k[w]) / 2;
	ke = (k[p] + k[e]) / 2;
	kn = (k[p] + k[n]) / 2;
	ks = (k[p] + k[s]) / 2;

}

void ft(int N, double Tmelt, double* T, double* f, double* fi, int* R) {

	int i;

	//printf("\n//------------ Evaluating f ------------//\n");

	for (i = 0; i < N; i++) {

		if (R[i] == 0) {
		
			f[i] = 0.5 * T[i] - 10.5;
			
			if (f[i] < 0) {

				f[i] = 0;

			}
			
			if (f[i] >= 1) {

				f[i] = 1;

			}
		}
	}
}

void kt(int N, double* T, double k_bat, double k_pcm, double* k, double* ki, int* R) {

	int i;

	for (i = 0; i < N; i++) {

		ki[i] = k[i];

		if (R[i] == 1) {

			k[i] = k_bat;

		}
		else {

			k[i] = k_pcm;

		}
	}
}

void plot_f(int N, int Nx, int Ny, double dx, double dy, char ffn[20], double* f, int it) {

	double rx{}, ry{};

	cout << "\n_________________________________________________________________________" << endl;
	cout << "\nWriting temperature file..." << endl;

	// Build
	ofstream fout;
	string iteration = to_string(it + 1);
	string filename = "C:/Users/renat/source/repos/volumetric_capacitor_method/VCM2D/LiquidFraction/LiquidFraction_" + iteration + ".vtk";
	fout.open(filename);
	fout << "# vtk DataFile Version 2.0" << endl;
	fout << "Liquid Fraction Mass plot" << endl;
	fout << "ASCII" << endl;
	fout << endl;
	fout << "DATASET RECTILINEAR_GRID" << endl;
	fout << "DIMENSIONS " << Ny << " " << Nx << " 1" << endl;
	fout << endl;
	fout << "Y_COORDINATES " << Nx << " float" << endl;


	for (int i = 0; i < Nx; i++)
	{

		rx = i * dx + 0.5 * dx;
		fout << rx << " ";

	}

	fout << endl;
	fout << endl;
	fout << "X_COORDINATES " << Ny << " float" << endl;

	for (int j = 0; j < Ny; j++)
	{

		ry = j * dy + 0.5 * dy;
		fout << ry << " ";

	}

	fout << endl;
	fout << endl;
	fout << "Z_COORDINATES 1 float" << endl;
	fout << "0" << endl;
	fout << endl;
	fout << "POINT_DATA " << Nx * Ny << endl;

	// Write Region Information
	fout << "SCALARS Liquid_Fraction float" << endl;
	fout << "LOOKUP_TABLE my_table" << endl;

	for (int o = 0; o < N; o++) {

		fout << f[o] << endl;

	}

	// Exit
	fout.close();
	cout << "\nDone." << endl;
	cout << "_________________________________________________________________________" << endl << endl;

}

void plot_temperature(int N, int Nx, int Ny, double dx, double dy, char ffn[20], double* T, int it) {

	double rx{}, ry{};

	cout << "\n_________________________________________________________________________" << endl;
	cout << "\nWriting temperature file..." << endl;

	// Build
	ofstream fout;
	string iteration = to_string(it + 1);
	string filename = "C:/Users/renat/source/repos/volumetric_capacitor_method/VCM2D/Temperature/Temperature_" + iteration + ".vtk";
	fout.open(filename);
	fout << "# vtk DataFile Version 2.0" << endl;
	fout << "Temperature °C plot" << endl;
	fout << "ASCII" << endl;
	fout << endl;
	fout << "DATASET RECTILINEAR_GRID" << endl;
	fout << "DIMENSIONS " << Ny << " " << Nx << " 1" << endl;
	fout << endl;
	fout << "Y_COORDINATES " << Nx << " float" << endl;


	for (int i = 0; i < Nx; i++)
	{

		rx = i * dx + 0.5 * dx;
		fout << rx << " ";

	}

	fout << endl;
	fout << endl;
	fout << "X_COORDINATES " << Ny << " float" << endl;

	for (int j = 0; j < Ny; j++)
	{

		ry = j * dy + 0.5 * dy;
		fout << ry << " ";

	}

	fout << endl;
	fout << endl;
	fout << "Z_COORDINATES 1 float" << endl;
	fout << "0" << endl;
	fout << endl;
	fout << "POINT_DATA " << Nx * Ny << endl;

	// Write Region Information
	fout << "SCALARS Temperature float" << endl;
	fout << "LOOKUP_TABLE my_table" << endl;

	for (int o = 0; o < N; o++) {

		fout << T[o] << endl;

	}

	// Exit
	fout.close();
	cout << "\nDone." << endl;
	cout << "_________________________________________________________________________" << endl << endl;

}

void plot_properties(int N, int Nx, int Ny, double dx, double dy, char ffn[20], int* R, double* Ti, double* k, double* rho, double* cp) {

	double rx{}, ry{};

	cout << "\n_________________________________________________________________________" << endl;
	cout << "\nWriting data files..." << endl;

	// Build
	ofstream fout;
	fout.open("SimulationPlots.vtk");
	fout << "# vtk DataFile Version 2.0" << endl;
	fout << "Physical Properties plot" << endl;
	fout << "ASCII" << endl;
	fout << endl;
	fout << "DATASET RECTILINEAR_GRID" << endl;
	fout << "DIMENSIONS " << Ny << " " << Nx << " 1" << endl;
	fout << endl;
	fout << "Y_COORDINATES " << Nx << " float" << endl;


	for (int i = 0; i < Nx; i++)
	{

		rx = i * dx + 0.5 * dx;
		fout << rx << " ";

	}

	fout << endl;
	fout << endl;
	fout << "X_COORDINATES " << Ny << " float" << endl;

	for (int j = 0; j < Ny; j++)
	{

		ry = j * dy + 0.5 * dy;
		fout << ry << " ";

	}

	fout << endl;
	fout << endl;
	fout << "Z_COORDINATES 1 float" << endl;
	fout << "0" << endl;
	fout << endl;
	fout << "POINT_DATA " << Nx * Ny << endl;

	// Write Region Information
	fout << "SCALARS Region int" << endl;
	fout << "LOOKUP_TABLE my_table" << endl;

	for (int o = 0; o < N; o++) {

		fout << R[o] << endl;

	}

	// Write Conductivity Information
	fout << "SCALARS Thermal_Conductivity float" << endl;
	fout << "LOOKUP_TABLE my_table" << endl;

	for (int o = 0; o < N; o++) {

		fout << k[o] << endl;

	}

	// Write Density Information
	fout << "SCALARS Density float" << endl;
	fout << "LOOKUP_TABLE my_table" << endl;

	for (int o = 0; o < N; o++) {

		fout << rho[o] << endl;

	}

	// Write Specific Heat Information
	fout << "SCALARS Specific_Heat float" << endl;
	fout << "LOOKUP_TABLE my_table" << endl;

	for (int o = 0; o < N; o++) {

		fout << cp[o] << endl;

	}

	// Exit
	fout.close();
	cout << "\nDone." << endl;
	cout << "_________________________________________________________________________" << endl << endl;

}

void average_temperature(int N, double* T) {

	double total = 0.0;

	for (int i = 0; i < N; i++) {

		total = total + T[i];

	}

	printf("Temperature average = %f", total / N);

}