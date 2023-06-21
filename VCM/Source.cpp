#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>

using namespace std;

#pragma warning(disable:4996)

// Processing
//void ft(int N, double* T, double*& f, double*& fi);
//void kt(int N, double* T, double*& k, double*& ki);
void map(int Nx, int Ny, double D, double dx, double dy, int* ww, int* ee, int* nn, int* ss, int* R);
void interface_k(int p, int w, int e, int n, int s, double* k, double& kp, double& kw, double& ke, double& kn, double& ks);
void assembly(int N, int Nx, int Ny, double dx, double dy, int* ww, int* ee, int* nn, int* ss, double* k, double* ap, double* aw, double* ae, double* an, double* as, double* s, double* sp, double* b, double Tw, double Te, double Tn, double Ts, double qw, double qe, double qn, double qs, double* T, double* Ti, double rho, double cp, double L, double* f, double* fi, double dt);
void SOR(int N, int* ww, int* ee, int* nn, int* ss, double* ap, double* aw, double* ae, double* an, double* as, double* b, double* T, double* Ti);


// Post Processing
void plots(int Nx, int Ny, double dx, double dy, char ffn[20], int* R, double* T);
//void output(int Nx, int Ny, double dx, double dy, int* R, char ffn[20]);

int main() {

	// Output file info
	char ffn[40];

	// ---------------------------------- --------------------------------//
	// Cell properties
	double D = 0.018; // Battery cell diameter [m]
	double h = 0.065; // Battery cell height [m]

	// Mesh Properties
	const int Nx = 5; // Volumes in x direction
	const int Ny = 1; // Volumes in y direction 
	double dx = 0.004; // CV dimension in x direction [m]
	double dy = 0.004; // CV dimension in y direction [m]
	const int N = (int)Nx * Ny; // Total number of CVs

	// ---------------------------------- --------------------------------//
	// Allocating memory and cleaning any previously stored value
	double* kres = (double*)malloc((N) * sizeof(double));
	double* fres = (double*)malloc((N) * sizeof(double));
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
	int* ee = (int*)malloc((N) * sizeof(int));     //Coeficiente a Leste
	int* ww = (int*)malloc((N) * sizeof(int));     //Coeficiente a Oeste
	int* nn = (int*)malloc((N) * sizeof(int));     //Coeficiente a Leste
	int* ss = (int*)malloc((N) * sizeof(int));     //Coeficiente a Oeste

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
		Ti[i] = 200;
		fi[i] = 0.0;
		ki[i] = 10.0;
		f[i] = 0.0;
		k[i] = 10.0;
		T[i] = 0.0;
		b[i] = 0.0;
		s[i] = 0.0;
		T[i] = 0.0;
		R[i] = 0;

	}

	// ------------------------ PHYSICAL PROPERTIES ----------------------//
	// Selected PCM: RT 18 HC [T_melt = 18°C, Latent heat of fusion = 250 kJ/kg, Cp = 2 kJ/kgK, k = 2 kJ/kgK, rho = 880 <-> 770 kg/m³ (solid/liquid)]

	double rhocp = 10000000;       // PCM density in solid state [kg/m³]
	//double k_bat = 10.0;   // Battery effective thermal conductivity [W/m.K]
	//double k_pcm = 2.0;    // PCM thermal conductivity [W/m.K]
	double cp = 2;		 // PCM specific heat at constante pressure [J/kg.K]
	double L = 250000;		 // PCM latent heat of fusion [J/kg]

	// ------------------------ BOUNDARY CONDITIONS ----------------------//
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

	// -------------------------- SIMULATION CODE ------------------------//
	//int i;
	bool converged = false; // Property convergence indicator
	double resmax = 1E-3;	// Maximum property residue 
	double dt = 1;			// Time step [s] 
	double time = dt;		// Current time step [s] 
	double TotalTime = 120; // Total simulation time [s]
	int i = 0;
	int l = 0;

	map(Nx, Ny, D, dx, dy, ww, ee, nn, ss, R);
	//assembly(N, Nx, Ny, dx, dy, ww, ee, nn, ss, k, ap, aw, ae, an, as, s, sp, b, Tw, Te, Tn, Ts, qw, qe, qn, qs, T, Ti, rho, cp, L, f, fi, dt);


	while (time <= TotalTime) {

		converged = false;

		while (converged == false) {

			assembly(N, Nx, Ny, dx, dy, ww, ee, nn, ss, k, ap, aw, ae, an, as, s, sp, b, Tw, Te, Tn, Ts, qw, qe, qn, qs, T, Ti, rhocp, cp, L, f, fi, dt);
			SOR(N, ww, ee, nn, ss, ap, aw, ae, an, as, b, T, Ti);

			converged = true;

			for (i = 0; i < N; i++) {

				kres[i] = abs(k[i] - ki[i]);
				fres[i] = abs(f[i] - fi[i]);
				printf("CV = %i\t f = %f\t Residuo = %13.5E\n", i, f[i], fres[i]);

				if ((kres[i] > resmax) && (fres[i] > resmax)) {

					converged = false;

				}
			}


			if (converged == true) {

				printf("\n// ----- Solution for time = %f ----- //\n", time);

				for (i = 0; i < N; i++) {

					printf("T%i = %f\n", i, T[i]);

				}

				time = time + dt;

			}
		}
	}

	//plots(Nx, Ny, dx, dy, ffn, R, T);

	return 0;
}

void map(int Nx, int Ny, double D, double dx, double dy, int* ww, int* ee, int* nn, int* ss, int* R) {

	int i, j, m, o = 0;
	double x_center, y_center, x, y, radius;
	double cell_radius = D / 2;
	double x_cell = (Nx * dx) / 2;
	double y_cell = (Ny * dy) / 2;

	for (i = 0; i < Nx; i++) {
		for (j = 0; j < Ny; j++) {

			// Mesh Regions
			x_center = (dx * 0.5) + i * dx;
			y_center = (dy * 0.5) + j * dx;
			x = x_center - x_cell;
			y = y_center - y_cell;
			radius = sqrt(pow(x, 2) + pow(y, 2));

			if (radius < cell_radius) {

				R[o] = 1;
			}
			else {

				R[o] = 0;
			}

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

		}
	}

}

void assembly(int N, int Nx, int Ny, double dx, double dy, int* ww, int* ee, int* nn, int* ss, double* k, double* ap, double* aw, double* ae, double* an,
	double* as, double* s, double* sp, double* b, double Tw, double Te, double Tn, double Ts, double qw, double qe, double qn, double qs, double* T, double* Ti,
	double rhocp, double cp, double L, double* f, double* fi, double dt) {

	int i, j, x, y;
	int o = 0;
	double kw = 0.0;
	double ke = 0.0;
	double kn = 0.0;
	double ks = 0.0;
	double kp = 0.0;

	for (i = 0; i < N; i++) {

		interface_k(i, ww[i], ee[i], nn[i], ss[i], k, kp, kw, ke, kn, ks);

		ae[i] = (2 * ke * kp) / (dx * (ke * dx + kp * dx));
		aw[i] = (2 * kw * kp) / (dx * (kw * dx + kp * dx));
		an[i] = (2 * kn * kp) / (dy * (kn * dy + kp * dy));
		as[i] = (2 * ks * kp) / (dy * (ks * dy + kp * dy));
		s[i] = 0.0;
		sp[i] = 0.0;

		x = (int)i / Ny;
		y = (int)i - x * Ny;

		// Boundary conditions - West
		if (x == 0) {

			aw[i] = 0;

		}

		// Boundary conditions - East
		if (x == Nx - 1) {

			ae[i] = 0;
			sp[i] = -(2 * k[i]) / (dx * dx);
			s[i] = Te * (2 * k[i]) / (dx * dx);

		}

		// Boundary conditions - North
		if (y == Ny - 1) {

			an[i] = 0;

		}

		// Boundary conditions - South
		if (y == 0) {

			as[i] = 0;

		}

		//Calculo do ap e do b
		ap[i] = aw[i] + ae[i] + ((rhocp) / dt) - sp[i];
		b[i] = s[i] + ((rhocp) / dt) * Ti[i]; //- ((rho * L * (f[pp] - fi[pp])) / dt);

		//printf("CV = %i\t aw = %.2f\t ae = %.2f\t an = %.2f\t as = %.2f\t ap = %.2f\t b = %.2f\n", i, aw[i], ae[i], an[i], as[i], ap[i], b[i]);
		//printf("CV = %i\t kw = %.2f\t ke = %.2f\t kn = %.2f\t ks = %.2f\t kp = %.2f\n", i, kw, ke, kn, ks, k[i]);

	}
}

void SOR(int N, int* ww, int* ee, int* nn, int* ss, double* ap, double* aw, double* ae, double* an, double* as, double* b, double* T, double* Ti) {

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

		printf("Iteracao = %i\t Residuo = %13.5E\n", it, res);
		it++;
	}

}

void plots(int Nx, int Ny, double dx, double dy, char ffn[20], int* R, double* T) {

	int i, j, o;
	double rx{}, ry{};
	cout << "\n_________________________________________________________________________" << endl;
	cout << "\nWriting data files..." << endl;

	// Build
	ofstream fout;
	fout.open("SimulationPlots.vtk");
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
	fout << "SCALARS Region int" << endl;
	fout << "LOOKUP_TABLE my_table" << endl;

	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++) {

			o = (i * Ny) + j;
			fout << R[o] << endl;

		}

	fout << "SCALARS Temperature float" << endl;
	fout << "LOOKUP_TABLE my_table" << endl;

	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++) {

			o = (i * Ny) + j;
			fout << T[o] << endl;

		}

	// Exit
	fout.close();
	cout << "\nDone." << endl;
	cout << "_________________________________________________________________________" << endl << endl;

}

void interface_k(int p, int w, int e, int n, int s, double* k, double& kp, double& kw, double& ke, double& kn, double& ks) {

	kp = k[p];
	kw = (k[p] + k[w]) / 2;
	ke = (k[p] + k[e]) / 2;
	kn = (k[p] + k[n]) / 2;
	ks = (k[p] + k[s]) / 2;

}