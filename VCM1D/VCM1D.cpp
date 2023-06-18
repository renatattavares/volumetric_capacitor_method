#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>

using namespace std;

#pragma warning(disable:4996)

void assembly(int N, double dx, double* ap, double* ae, double* aw, double* s, double* sp, double* k, double rhocp, double* T, double* Ti, double* b, double dt);
void SOR(int N, double dx, double* ap, double* ae, double* aw, double* s, double* T, double* Ti, double* b, double dt);
void kt(int N, double* T, double*& k, double*& ki);

int main() {

	// Output file info
	//char ffn[40];

	// Mesh Properties
	const int N = 5; // Volumes in x direction
	const double dx = 0.004; // CV dimension in x direction [m]
	double dt = 1;	   // Time step [s] 
	double time = dt;
	bool converged = false;
	double resmax = 1E-6;
	int i = 0;

	// ---------------------------------- --------------------------------//
	// Allocating memory and cleaning any previously stored value
	double* res = (double*)malloc((N) * sizeof(double));
	double* ap = (double*)malloc((N) * sizeof(double));
	double* ae = (double*)malloc((N) * sizeof(double));
	double* aw = (double*)malloc((N) * sizeof(double));
	double* Ti = (double*)malloc((N) * sizeof(double));
	double* sp = (double*)malloc((N) * sizeof(double));
	double* ki = (double*)malloc((N) * sizeof(double));
	double* k = (double*)malloc((N) * sizeof(double));
	double* b = (double*)malloc((N) * sizeof(double));
	double* s = (double*)malloc((N) * sizeof(double));
	double* T = (double*)malloc((N) * sizeof(double));

	for (int i = 0; i < N; i++) {

		// Floats
		res[i] = 0.0;
		ap[i] = 0.0;
		ae[i] = 0.0;
		aw[i] = 0.0;
		Ti[i] = 200.0;
		T[i] = 0.0;
		k[i] = 10.0;
		ki[i] = 10.0;
		b[i] = 0.0;
		s[i] = 0.0;
		T[i] = 0.0;

	}
	// ---------------------------------- --------------------------------//
	// Physical properties

	double rhocp = 10000000;
	//double k = 10.0;
	double TotalTime = 120; // Total simulation time [s]


	// ------------------------ BOUNDARY CONDITIONS ----------------------//
	double Te = 0.0; // East
	double Tw = 0.0; // West

	// -------------------------- SIMULATION CODE ------------------------//

	while (time <= TotalTime) {

		converged = false;

		while (converged == false) {

			assembly(N, dx, ap, ae, aw, s, sp, k, rhocp, T, Ti, b, dt);
			SOR(N, dx, ap, ae, aw, s, T, Ti, b, dt);
			kt(N, T, k, ki);

			converged = true;

			for (i = 0; i < N; i++) {

				res[i] = abs(k[i] - ki[i]);
				printf("CV = %i\t Residuo = %13.5E\n", i, res[i]);

				if (res[i] > resmax) {

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


	/*
	while (time <= TotalTime) {

		converged = false;

		while (converged == false) {

			assembly(N, dx, ap, ae, aw, s, sp, k, rhocp, T, Ti, b, dt);
			SOR(N, dx, ap, ae, aw, s, T, Ti, b, dt);
			kt(N, T, k, ki);

			for (i = 0; i < N; i++) {

				res[i] = abs(k[i] - ki[i]);
				converged = true;

				if (res[i] > resmax) {

					converged = false;

				}
			}

			if(converged == true){

				printf("\n// ----- Solution for time = %f ----- //\n", time);

				for (i = 0; i < N; i++) {

					printf("T%i = %f\n", i, T[i]);

				}

				time = time + dt;

			}		
		}
	}
		*/
	/*for (i = 0; i < N; i++) {
		printf("CV = %i \t k = %f\n", i, k[i]);
	}
	
	kt(N, T, k);

	for (i = 0; i < N; i++) {
		printf("CV = %i \t k = %f\n", i, k[i]);
	}*/
	
	return 0;

}

void assembly(int N, double dx, double* ap, double* ae, double* aw, double* s, double* sp, double* k, double rhocp, double* T, double* Ti, double* b, double dt)
{
	int o;
	double kw = 0.0;
	double ke = 0.0;

	for (o = 0; o < N; o++) {

		if ((o - 1) >= 0)
		{
			kw = (k[o] + k[o - 1]) / 2;
		}
		else
		{
			kw = 0;
		}

		if ((o + 1) < N)
		{
			ke = (k[o] + k[o + 1]) / 2;
		}
		else
		{
			ke = 0;
		}
	
		ae[o] = (ke / (dx * dx));
		aw[o] = (kw / (dx * dx));
		sp[o] = 0;

		//Condição de contorno - Face leste
		if (o == N - 1) {

			ae[o] = 0;
			sp[o] = -(2 * k[o]) / (dx * dx);

		}

		//Condição de contorno - Face oeste
		if (o == 0) {

			aw[o] = 0;

		}

		//Calculo do ap e do b
		ap[o] = aw[o] + ae[o] + ((rhocp) / dt) - sp[o];
		b[o] = s[o] + ((rhocp) / dt) * Ti[o];

		//printf("CV = %i\t kw = %5.1E\t ke = %5.1E\t k = 10.0\n", o, kw, ke);
		//printf("CV = %i\t f = %i\t fi = %i\t rho = %5.1E\t rhoi = %5.1E\t Region = %i\n", o, f, fi, rho, rhoi, R[o]);
		//printf("CV = %i\t aw = %5.1E\t ae = %5.1E\t ap = %5.1E\t b = %5.1E\n", o, aw[o], ae[o], ap[o], b[o]);
	}

}

void SOR(int N, double dx, double* ap, double* ae, double* aw, double* s, double* T, double* Ti, double* b, double dt) {

	int i;
	double R = 0.0;
	double res = 1.0;
	double kappa = 1000;
	double resmax = 1E-6;
	int it = 0;

	//printf("\nSolver SOR:\n");

	while (res > resmax) {

		res = 0.0;

		for (i = 0; i < N; i++) {

			Ti[i] = T[i];
			R = b[i] + aw[i] * T[i - 1] + ae[i] * T[i + 1];
			T[i] = (T[i] + kappa * (R / ap[i])) / (1 + kappa);

			res = res + pow((Ti[i] - T[i]), 2);

		}

		//printf("Iteracao = %i\t Residuo = %13.5E\n", it, res);
		it++;
	}

}

void kt(int N, double* T, double* &k, double* &ki) {

	int i;

	for (i = 0; i < N; i++) {

		ki[i] = k[i];
		k[i] = 10; //+ 0.001 * T[i];

	}
	
}