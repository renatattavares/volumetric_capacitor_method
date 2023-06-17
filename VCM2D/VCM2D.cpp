#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include "VCM2D.h"

using namespace std;

#pragma warning(disable:4996)

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
	double* fres = (double*)malloc((N) * sizeof(double));
	double* ap = (double*)malloc((N) * sizeof(double));
	double* ae = (double*)malloc((N) * sizeof(double));
	double* aw = (double*)malloc((N) * sizeof(double));
	double* an = (double*)malloc((N) * sizeof(double));
	double* as = (double*)malloc((N) * sizeof(double));
	double* Ti = (double*)malloc((N) * sizeof(double));
	double* b = (double*)malloc((N) * sizeof(double));
	double* s = (double*)malloc((N) * sizeof(double));
	double* T = (double*)malloc((N) * sizeof(double));
	int* fi = (int*)malloc((N) * sizeof(int));
	int* f = (int*)malloc((N) * sizeof(int));
	int* R = (int*)malloc((N) * sizeof(int));

	for (int i = 0; i < N; i++) {

		// Floats
		fres[i] = 0.0;
		ap[i] = 0.0;
		ae[i] = 0.0;
		an[i] = 0.0;
		as[i] = 0.0;
		aw[i] = 0.0;
		Ti[i] = 15;
		T[i] = 15;
		b[i] = 0.0;
		s[i] = 0.0;
		T[i] = 0.0;
		// Integers
		fi[i] = 0;
		f[i] = 0;
		R[i] = 0;

	}
	// ---------------------------------- --------------------------------//
	// Physical properties
	// Selected PCM: RT 18 HC [T_melt = 18°C, Latent heat of fusion = 250 kJ/kg, Cp = 2 kJ/kgK, k = 2 kJ/kgK, rho = 880 <-> 770 kg/m³ (solid/liquid)]

	double Tmelt = 18;			 // PCM melting temperature [°C]
	double rho_pcm_solid = 880;  // PCM density in solid state [kg/m³]
	double rho_pcm_liquid = 770; // PCM density in liquid state [kg/m³]
	double k_bat = 1000.0;       // Battery effective thermal conductivity [W/m.K]
	double k_pcm = 2.0;          // PCM thermal conductivity [W/m.K]
	double cp = 2;				 // PCM specific heat at constante pressure [kJ/kg.K]
	double L = 250;				 // PCM latent heat of fusion [kJ/kg]
	double TotalTime = 1200;		 // Total simulation time [s]
	

	// ---------------------------------- --------------------------------//
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
	// ---------------------------------- --------------------------------//

	// -------------------------- SIMULATION CODE ------------------------//
	mesh_regions(D, dx, dy, Nx, Ny, R);
	plot_regions(Nx, Ny, dx, dy, ffn, R, T);
	
	while (time <= TotalTime) {

		fconverged = false;

		while (fconverged == false) {

			assembly(Nx, Ny, dx, dy, R, ap, ae, aw, an, as, s, cp, L, T, Ti, b, dt);
			SOR(Nx, Ny, dx, dy, R, ap, ae, aw, an, as, s, k_bat, k_pcm, rho_pcm_solid, rho_pcm_liquid, cp, L, T, Ti, Tmelt, b, dt, f, fi, fres);

			fconverged = areAllElementsSmaller(fres, N, fresmax);

		}

		printf("// ----- Solution for Time = %f ----- //\n", time);

		for (i = 0; i < N; i++) {
			
			
			printf("T%i = \t %f\n", i, T[i]);

		}

		time = time + dt;

	}
	
	return 0;

}

void mesh_regions(double D, double dx, double dy, int Nx, int Ny, int* R) {

	// This funtion will define each mesh region based on the distance between the battery cell center and the CV center.

	int i, j, o;
	double x_center, y_center, x, y, radius;
	double cell_radius = D / 2;
	double x_cell = (Nx * dx) / 2;
	double y_cell = (Ny * dy) / 2;

	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++) {

			o = (i * Ny) + j;
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

			//printf("CV = %i\t x = %5.1E\t y = %5.1E\t radius = %5.1E\t Region = %i\n", o, x, y, radius, R[o]);
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

void find_k(double dx, double dy, int* R, int o, int ow, int oe, int on, int os, double& kp, double& kw, double& ke, double& kn, double& ks)
{
	double kWest = 0.0;
	double kEast = 0.0;
	double kNorth = 0.0;
	double kSouth = 0.0;
	double k_bat = 1000.0;       // Battery effective thermal conductivity [W/m.K]
	double k_pcm = 2.0;          // PCM thermal conductivity [W/m.K]

	if (R[o] == 1) {
		kp = k_bat;
	}
	else {
		kp = k_pcm;
	}

	if (R[ow] == 1) {
		kWest = k_bat;
	}
	else {
		kWest = k_pcm;
	}

	if (R[oe] == 1) {
		kEast = k_bat;
	}
	else {
		kEast = k_pcm;
	}

	if (R[on] == 1) {
		kNorth = k_bat;
	}
	else {
		kNorth = k_pcm;
	}

	if (R[os] == 1) {
		kSouth = k_bat;
	}
	else {
		kSouth = k_pcm;
	}

	ke = (kp * kEast * 2 * dx) / (kEast * dx + kp * dx);
	kw = (kp * kWest * 2 * dx) / (kWest * dx + kp * dx);
	kn = (kp * kNorth * 2 * dy) / (kNorth * dy + kp * dy);
	ks = (kp * kSouth * 2 * dy) / (kSouth * dy + kp * dy);

}

void find_rho(int o, double& rho, double& rhoi, double* T, double* Ti)
{

	double Tmelt = 18;			 // PCM melting temperature [°C]
	double rho_pcm_solid = 880;  // PCM density in solid state [kg/m³]
	double rho_pcm_liquid = 770; // PCM density in liquid state [kg/m³]

	if (T[o] > Tmelt) {

		rho = rho_pcm_liquid;

	}
	else {

		rho = rho_pcm_solid;

	}

	if (Ti[o] > Tmelt) {

		rhoi = rho_pcm_liquid;

	}
	else {

		rhoi = rho_pcm_solid;

	}
}

void find_f(int o, int& f, int& fi, double* T, double* Ti)
{

	double Tmelt = 18;

	if (T[o] > Tmelt) {
		f = 1;
	}
	else {
		f = 0;
	}

	if (Ti[o] > Tmelt) {
		fi = 1;
	}
	else {
		fi = 0;
	}

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

bool areAllElementsSmaller(double* ptr, int size, double value) {
	for (int i = 0; i < size; ++i) {
		if (ptr[i] >= value) {
			return false; // If any element is not smaller, return false
		}
	}
	return true; // All elements are smaller
}