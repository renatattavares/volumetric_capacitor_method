#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>

using namespace std;

#pragma warning(disable:4996)

void mesh_region(double D, double dx, double dy, int Nx, int Ny, int* R);
//void output(int Nx, int Ny, double dx, double dy, int* R, char ffn[20]);
void plotRegions(int Nx, int Ny, double dx, double dy, int* R, char ffn[20]);
//void assembly(int Ny, int Nx, double* ap, double* ae, double* aw, double* an, double* as, double* sp,
	//double* su, double k, double dx, double dy, double tk, double qo, double Tn, double Aw, double Ae, double An, double As);
//void SOR(int N, int Nx, int Ny, double* ap, double* ae, double* aw, double* an, double* as, double* su, double* T, double* Ta);

int main() {

	// ---------------------------------- --------------------------------//
	// Cell properties
	double D = 0.018; // Battery cell diameter [m]
	double h = 0.065; // Battery cell height [m]
	
	// Mesh Properties
	const int Nx = 50; // Volumes in x direction
	const int Ny = 50; // Volumes in y direction 
	double dx = 0.001; // CV dimension in x direction [m]
	double dy = 0.001; // CV dimension in y direction [m]
	const int N = (int)Nx * Ny; // Total number of CVs


	// Output file info
	char ffn[40];

	// ---------------------------------- --------------------------------//
	// Allocating memory and cleaning any previously stored value
	double* ap = (double*)malloc((N) * sizeof(double));
	double* ae = (double*)malloc((N) * sizeof(double));
	double* aw = (double*)malloc((N) * sizeof(double));
	double* an = (double*)malloc((N) * sizeof(double));
	double* as = (double*)malloc((N) * sizeof(double));
	double* s =  (double*)malloc((N) * sizeof(double)); 
	double* rx = (double*)malloc((N) * sizeof(double));  
	double* ry = (double*)malloc((N) * sizeof(double));
	double* T = (double*)malloc((N) * sizeof(double));
	int* R = (int*)malloc((N) * sizeof(int));

	for (int i = 0; i < N; i++) {

		ap[i] = 0.0;
		ae[i] = 0.0;
		an[i] = 0.0;
		as[i] = 0.0;
		aw[i] = 0.0;
		s[i] = 0.0;

	}
	// ---------------------------------- --------------------------------//

	// Physical properties
	// Selected PCM: RT 18 HC [T_melt = 18°C, Latent heat of fusion = 250 kJ/kg, Cp = 2 kJ/kgK, k = 2 kJ/kgK, rho = 880 <-> 770 kg/m³ (solid/liquid)]
	double k_bat = 1000.0;       // Battery effective thermal conductivity [W/mK]
	double k_pcm = 2.0;          // PCM thermal conductivity [W/mK]
	double rho_pcm_solid  = 880; // PCM density in solid state [kg/m³]
	double rho_pcm_liquid = 770; // PCM density in liquid state [kg/m³]

	// ---------------------------------- --------------------------------//
	// Temperatures for boundary conditions
	double Tn = 0.0; // North
	double Ts = 0.0; // South 
	double Te = 0.0; // East
	double Tw = 0.0; // WEst
	
	// Heat flux for boundary conditions
	double qn = 0.0; // North
	double qs = 0.0; // South 
	double qe = 0.0; // East
	double qw = 0.0; // West
	// ---------------------------------- --------------------------------//
	
	// -------------------------- SIMULATION CODE ------------------------//
	mesh_region(D, dx, dy, Nx, Ny, R);
	plotRegions(Nx, Ny, dx, dy, R, ffn);
	
	//assembly(Ny, Nx, ap, ae, aw, an, as, sp, su, k, dx, dy, tk, qo, Tn, Aw, Aw, An, As);
	// SOR(N, Nx, Ny, ap, ae, aw, an, as, su, T, Ta);
	// ---------------------------------- --------------------------------//
	return 0;

}

void mesh_region(double D, double dx, double dy, int Nx, int Ny, int *R) {

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
			radius = sqrt(pow(x,2) + pow(y,2));

			if (radius < cell_radius) {

				R[o] = 1;
			}
			else {

				R[o] = 0;
			}

			printf("CV = %i\t x = %5.1E\t y = %5.1E\t radius = %5.1E\t Region = %i\n", o, x, y, radius, R[o]);
		}
		
}

void assembly(int Ny, int Nx, double* ap, double* ae, double* aw, double* an, double* as, double* sp, double* su, double k, double dx, double dy, double tk, double qo, double Tn, double Aw, double Ae, double An, double As)
{
	int i, j, o;

	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++) {

			o = (i * Ny) + j;
			aw[o] = Aw * (k / dx);
			ae[o] = Ae * (k / dx);
			an[o] = An * (k / dy);
			as[o] = As * (k / dy);
			sp[o] = 0.0;
			su[o] = 0.0;

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


			}

			//Calculo do ap
			ap[o] = aw[o] + ae[o] + an[o] + as[o] - sp[o];

			//printf("CV = %i\t i = %i\t j = %i\n", o, i, j);
			//printf("CV = %i\t an = %5.1E\t as = %5.1E\t aw = %5.1E\t ae = %5.1E\t ap = %5.1E\t su = %5.1E\n", o, an[o], as[o], aw[o], ae[o], ap[o], su[o]);

		}
}

void SOR(int N, int Nx, int Ny, double* ap, double* ae, double* aw, double* an, double* as, double* su, double* T, double* Ta) {

	int i, j, o, ow, oe, on, os, it = 0;
	double R;
	double res = 1.0;
	double kappa = 1000.0;
	double resmax = 1E-6;

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

				Ta[o] = T[o];
				R = su[o] + aw[o] * T[ow] + ae[o] * T[oe] + as[o] * T[os] + an[o] * T[on];
				T[o] = (T[o] + kappa * (R / ap[o])) / (1 + kappa);

				res = res + (Ta[o] - T[o]) * (Ta[o] - T[o]);

			}

		printf("Iteracao = %i\t Residuo = %13.5E\n", it, res);
		it++;
	}

	printf("\nSolucao:\n");

	for (i = 0; i < N; i++) {

		printf("T%i = \t %f\n", i, T[i]);

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

void plotRegions(int Nx, int Ny, double dx, double dy, int* R, char ffn[20]) {

	int i, j, o;
	double rx{}, ry{};
	cout << "\n_________________________________________________________________________" << endl;
	cout << "\nWriting data files..." << endl;

	// Build
	ofstream fout;
	fout.open("MeshRegions.vtk");
	fout << "# vtk DataFile Version 2.0" << endl;
	fout << "Temperature °C plot" << endl;
	fout << "ASCII" << endl;
	fout << endl;
	fout << "DATASET RECTILINEAR_GRID" << endl;
	fout << "DIMENSIONS " << Ny << " " << Nx << " 1" << endl;
	fout << endl;
	fout << "Y_COORDINATES " << Nx << " float" << endl;


	for (int i = 0; i < Nx ; i++)
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
	fout << "SCALARS Temperature float" << endl;
	fout << "LOOKUP_TABLE my_table" << endl;

	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++) {

			o = (i * Ny) + j;
			fout << R[o] << endl;		
		
		}

	// Exit
	fout.close();
	cout << "\nDone." << endl;
	cout << "_________________________________________________________________________" << endl << endl;

}