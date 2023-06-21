#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include "VCM2D.h"

using namespace std;

#pragma warning(disable:4996)
void plot_regions(int Nx, int Ny, double dx, double dy, char ffn[20], int* R, double* T);
void assembly(int Nx, int Ny, double dx, double dy, int* R, double* ap, double* ae, double* aw, double* an, double* as, double* s,
	double cp, double L, double* T, double* Ti, double* b, double dt);
void output(int Nx, int Ny, double dx, double dy, int* R, char ffn[20]);
void SOR(int Nx, int Ny, double dx, double dy, int* R, double* ap, double* ae, double* aw, double* an, double* as, double* s,
	double k_bat, double k_pcm, double rho_pcm_solid, double rho_pcm_liquid, double cp, double L, double* T, double* Ti, double Tmelt,
	double* b, double dt, int* f, int* fi, double* fres);

int main() {

	// Output file info
	char ffn[40];

	// ---------------------------------- --------------------------------//
	// Cell properties
	double D = 0.018; // Battery cell diameter [m]
	double h = 0.065; // Battery cell height [m]

	// Mesh Properties
	const int Nx = 10; // Volumes in x direction
	const int Ny = 10; // Volumes in y direction 
	double dx = 0.005; // CV dimension in x direction [m]
	double dy = 0.005; // CV dimension in y direction [m]
	const int N = (int)Nx * Ny; // Total number of CVs
	double dt = 2;	   // Time step [s] 
	double time = dt;
	bool fconverged = false;
	double fresmax = 1E-3;
	int i = 0;

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

	for (int i = 0; i < N; i++) {

		// Floats
		kres[i] = 0.0;
		fres[i] = 0.0;
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

	double rho = 1000;       // PCM density in solid state [kg/m³]
	//double k_bat = 10.0;     // Battery effective thermal conductivity [W/m.K]
	//double k_pcm = 2.0;      // PCM thermal conductivity [W/m.K]
	double cp = 2000;			 // PCM specific heat at constante pressure [J/kg.K]
	double L = 250000;			 // PCM latent heat of fusion [J/kg]


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
	bool converged = false; // Property convergence indicator
	double resmax = 1E-3;	// Maximum property residue 
	double dt = 1;			// Time step [s] 
	double time = dt;		// Current time step [s] 
	double TotalTime = 120; // Total simulation time [s]
	int i = 0;

	// mesh_regions(D, dx, dy, Nx, Ny, R);
	plot_regions(Nx, Ny, dx, dy, ffn, R, T);

	while (time <= TotalTime) {

		converged = false;

		while (converged == false) {

			assembly(Nx, Ny, dx, dy, R, ap, ae, aw, an, as, s, cp, L, T, Ti, b, dt);
			SOR(Nx, Ny, dx, dy, R, ap, ae, aw, an, as, s, k_bat, k_pcm, rho, cp, L, T, Ti, b, dt, f, fi, fres);
			kt(N, T, k, ki);
			ft(N, T, f, fi);

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

			return 0;

		}
	}
}

void plot_regions(int Nx, int Ny, double dx, double dy, char ffn[20], int* R, double* T) {

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

void assembly(int Nx, int Ny, double dx, double dy, int* R, double* ap, double* ae, double* aw, double* an, double* as, double* s,
	double cp, double L, double* T, double* Ti, double* b, double dt)
{
	int ow, oe, on, os;
	int i, j, o;
	double rho = 0.0;
	double rhoi = 0.0;
	double kp = 0.0;
	double kw = 0.0;
	double ke = 0.0;
	double kn = 0.0;
	double ks = 0.0;
	int f = 0;
	int fi = 0;


	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++) {

			o = (i * Ny) + j;
			ow = ((i - 1) * Ny) + j; // Id - West neighbor
			oe = ((i + 1) * Ny) + j; // Id - East neighbor
			on = (i * Ny) + (j + 1); // Id - North neighbor
			os = (i * Ny) + (j - 1); // Id - South neighbor

			find_k(dx, dy, R, o, ow, oe, on, os, kp, kw, ke, kn, ks);
			find_rho(o, rho, rhoi, T, Ti);
			find_f(o, f, fi, T, Ti);

			if (R[o] == 1) {

				s[o] = 50000;
			}
			else {

				s[o] = 0;

			}

			ae[o] = (2 * ke * kp) / (dx * (ke * dx + kp * dx));
			aw[o] = (2 * kw * kp) / (dx * (kw * dx + kp * dx));
			an[o] = (2 * kn * kp) / (dy * (kn * dy + kp * dy));
			as[o] = (2 * ks * kp) / (dy * (ks * dy + ks * dy));

			//Condição de contorno - Face norte
			if (j == Ny - 1) {

				an[o] = 0;

			}

			//Condição de contorno - Face sul
			if (j == 0) {

				as[o] = 0;

			}

			//Condição de contorno - Face leste
			if (i == Nx - 1) {

				ae[o] = 0;

			}

			//Condição de contorno - Face oeste
			if (i == 0) {

				aw[o] = 0;
				s[o] = s[o] + 2 * (kw / dy) * dy * 100;

			}

			//Calculo do ap e do b
			ap[o] = aw[o] + ae[o] + an[o] + as[o] + s[o] + ((rho * cp) / dt) + (cp * (rho - rhoi)) / dt;
			b[o] = s[o] + ((rho * cp) / dt) * T[o] - (rho * L * (f - fi)) / dt - (f * L * (rho - rhoi)) / dt;

			//printf("CV = %i\t kn = %5.1E\t ks = %5.1E\t kw = %5.1E\t ke = %5.1E\t kp = %5.1E\t Region = %i\n", o, kn, ks, kw, ke, kp, R[o]);
			//printf("CV = %i\t f = %i\t fi = %i\t rho = %5.1E\t rhoi = %5.1E\t Region = %i\n", o, f, fi, rho, rhoi, R[o]);
			//printf("CV = %i\t an = %5.1E\t as = %5.1E\t aw = %5.1E\t ae = %5.1E\t ap = %5.1E\t b = %5.1E\n", o, an[o], as[o], aw[o], ae[o], ap[o], b[o]);

		}
}

void output(int Nx, int Ny, double dx, double dy, int* R, char ffn[20]) {

	int i, j, o;
	double rx, ry;

	FILE* PFL;

	PFL = fopen(ffn, "w");

	fprintf(PFL, "title = \"tentativa_1\" \n");
	fprintf(PFL, "variables = ");
	fprintf(PFL, "\"x\" \"y\" \"T\"");
	fprintf(PFL, "\nzone i = %i j = %i  T = \"point\" \n", Ny, Nx);

	for (i = 1; i < Nx; i++)
		for (j = 1; j < Ny; j++) {

			o = i * Ny + j;
			rx = i * dx + 0.5 * dx;
			ry = j * dy + 0.5 * dy;

			fprintf(PFL, " %13.5E %13.5E %i", rx, ry, R[o]);


			fprintf(PFL, " \n");

		}
	fprintf(PFL, " \n");
	fclose(PFL);

}

void SOR(int Nx, int Ny, double dx, double dy, int* R, double* ap, double* ae, double* aw, double* an, double* as, double* s,
	double k_bat, double k_pcm, double rho_pcm_solid, double rho_pcm_liquid, double cp, double L, double* T, double* Ti, double Tmelt,
	double* b, double dt, int* f, int* fi, double* fres) {

	int i, j, o, ow, oe, on, os, N = 0;
	double RR = 0.0;
	double res = 1.0;
	double kappa = 1000;
	double resmax = 1E-6;
	int it = 0;

	printf("\nSolver SOR:\n");

	while (res > resmax) {

		res = 0.0;

		for (i = 0; i < Nx; i++)
			for (j = 0; j < Ny; j++) {

				o = (i * Ny) + j;
				ow = ((i - 1) * Ny) + j;
				oe = ((i + 1) * Ny) + j;
				on = (i * Ny) + (j + 1);
				os = (i * Ny) + (j - 1);

				Ti[o] = T[o];
				RR = b[o] + aw[o] * T[ow] + ae[o] * T[oe] + as[o] * T[os] + an[o] * T[on];
				T[o] = (T[o] + kappa * (RR / ap[o])) / (1 + kappa);

				res = res + pow((Ti[o] - T[o]), 2);

			}

		//printf("Iteracao = %i\t Residuo = %13.5E\n", it, res);
		it++;
	}



	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++) {

			o = (i * Ny) + j;

			find_f(o, f[o], fi[o], T, Ti);

			fres[o] = fres[o] + abs(fi[o] - f[o]);

		}
}