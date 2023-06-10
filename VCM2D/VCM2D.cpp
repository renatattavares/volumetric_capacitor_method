#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>

using namespace std;

#pragma warning(disable:4996)

void assembly(int Ny, int Nx, double* ap, double* ae, double* aw, double* an, double* as, double* sp,
	double* su, double k, double dx, double dy, double tk, double qo, double Tn, double Aw, double Ae, double An, double As);
void SOR(int N, int Nx, int Ny, double* ap, double* ae, double* aw, double* an, double* as, double* su, double* T, double* Ta);

int main() {

	// Propriedades da malha
	double tk = 0.01; // Espessura da placa [m]
	double Lx = 0.3; // Comprimento na direção x do domínio
	double Ly = 0.4; // Comprimento na direção y do domínio
	const int Nx = 3; // Número de malhas em x
	const int Ny = 4; // Número de malhas em y
	const int N = (int)Nx * Ny; //Número de pontos da malha
	double dx = Lx / Nx; // Comprimento em x do volume de controle
	double dy = Ly / Ny; // Comprimento em y do volume de controle
	double Aw = tk * dy;
	double Ae = tk * dy;
	double An = tk * dx;
	double As = tk * dx;

	printf("Numero de elementos: %i\n", N);

	//Criando vetor de variáveis
	double* ap = (double*)malloc((N) * sizeof(double));
	double* ae = (double*)malloc((N) * sizeof(double));
	double* aw = (double*)malloc((N) * sizeof(double));
	double* an = (double*)malloc((N) * sizeof(double));
	double* as = (double*)malloc((N) * sizeof(double));
	double* su = (double*)malloc((N) * sizeof(double));
	double* sp = (double*)malloc((N) * sizeof(double));
	double* T = (double*)malloc((N) * sizeof(double));
	double* Ta = (double*)malloc((N) * sizeof(double));

	for (int i = 0; i < N; i++) {

		ap[i] = 0.0;
		ae[i] = 0.0;
		an[i] = 0.0;
		as[i] = 0.0;
		aw[i] = 0.0;
		su[i] = 0.0;
		sp[i] = 0.0;
		T[i] = 100.0;
		Ta[i] = 100.0;

	}

	// Propriedades físicas
	double k = 1000.0; // Condutividade térmica [W/mK]

	// Condições de contorno
	// Temperatura prescrita
	double Tn = 100.0; // norte
	double Ts = 0.0; // sul
	double Tl = 0.0; // leste
	double To = 0.0; // oeste

	// Fluxo de calor prescrito
	double qn = 0.0; // norte 
	double qs = 0.0; // sul 
	double ql = 0.0; // leste
	double qo = 500000.0; // oeste

	assembly(Ny, Nx, ap, ae, aw, an, as, sp, su, k, dx, dy, tk, qo, Tn, Aw, Aw, An, As);
	SOR(N, Nx, Ny, ap, ae, aw, an, as, su, T, Ta);

	return 0;

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
				su[o] = su[o] + 2 * (k / dy) * An * Tn; // Temperatura prescrita ?
				sp[o] = sp[o] - 2 * (k / dy) * An; // Temperatura prescrita ?
				//su[o] = su[o] + qn * An; // Fluxo prescrito

			}

			//Condição de contorno - Face sul
			if (j == 0) {

				as[o] = 0;
				//su[o] = su[o] + qs * As; // Fluxo prescrito

			}

			//Condição de contorno - Face leste
			if (i == Nx - 1) {

				ae[o] = 0;
				//su[o] = su[o] + qe * Ae; // Fluxo prescrito
			}

			//Condição de contorno - Face oeste
			if (i == 0) {

				aw[o] = 0;
				su[o] = su[o] + qo * Aw; // Fluxo prescrito

			}

			//Calculo do ap
			ap[o] = aw[o] + ae[o] + an[o] + as[o] - sp[o];

			//printf("CV = %i\t i = %i\t j = %i\n", o, i, j);
			printf("CV = %i\t an = %5.1E\t as = %5.1E\t aw = %5.1E\t ae = %5.1E\t ap = %5.1E\t su = %5.1E\n", o, an[o], as[o], aw[o], ae[o], ap[o], su[o]);

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