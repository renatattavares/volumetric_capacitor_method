#include "..\Headers\Include.h"

void SORt(int o, int* pp, int* ww, int* ee, int* nn, int* ss, double* ap, double* aw, double* ae, double* an, double* as, double* b, double* T, double* Ti) {


	int i, j, it = 0;
	double R, T_old;
	double res = 1.0;
	double kappa = 1000;
	double resmax = 1E-6;

	while (res > resmax)
	{
		res = 0.0;
		#pragma omp parallel for shared (res) private (j, T_old, R)
		for (i = 0; i < o; i++)
		{
			if (pp[i] != 0)
			{
				j = pp[i];
				T_old = T[j];
				R = b[j] + aw[j] * T[ww[j]] + ae[j] * T[ee[j]] + an[j] * T[nn[j]] + as[j] * T[ss[j]];
				T[j] = (T_old + kappa * (R / ap[j])) / (1 + kappa);
				//printf("CV = %i\t w = %5.1E\t e = %5.1E\t n = %5.1E\t s = %5.1E\n", i, aw[i] * T[ww[i]], ae[i] * T[ee[i]], an[i] * T[nn[i]], as[i] * T[ss[i]]);
				res = res + fabs(T[j] - T_old);
			}
		}
		res = res / o;
		//it++;
		//printf("Iteracao = %i\t Residuo = %13.5E\n", it, res);			  
	}
}


void SORf(int o, int* pp, int* ww, int* ee, int* nn, int* ss, double* ap, double* aw, double* ae, double* an, double* as, double* b, double* T, double* Ti, int* R, double* f, double* fi, double* rho, double dt, double L_pcm, double L_cpcm, double Tmelt, double* resf) {

	int i, j;
	int N = 1;
	int can_melt = 1;
	double res = 0.0;
	double total_f = 0.0;
	double kappa = 0.23;			// Under relaxation factor
	double resmax = 1.0E-3;			// Maximum property residue 
	double L, Y, C;
	double f_old;

	for (i = 0; i < o; i++)
	{
		if (pp[i] != 0)
		{
			j = pp[i];

			// Latent Heat of Fusion
			if (R[j] == 0)
			{
				L = L_pcm;
				f_old = f[j];
				Y = aw[j] * T[ww[j]] + ae[j] * T[ee[j]] + an[j] * T[nn[j]] + as[j] * T[ss[j]];
				C = b[j] + ((rho[j] * L) / dt) * f_old;
				f[j] = (f_old + (kappa * (dt / (rho[j] * L)) * (Y + C - ap[j] * Tmelt))) / (1 + kappa);
			}
			else if (R[j] == 5)
			{
				L = L_cpcm;
				f_old = f[j];
				Y = aw[j] * T[ww[j]] + ae[j] * T[ee[j]] + an[j] * T[nn[j]] + as[j] * T[ss[j]];
				C = b[j] + ((rho[j] * L) / dt) * f_old;
				f[j] = (f_old + (kappa * (dt / (rho[j] * L)) * (Y + C - ap[j] * Tmelt))) / (1 + kappa);
			}
			else
			{
				f_old = 0.0;
				f[j] = 0.0;
				continue;
			}

			// Impose physical limit
			if (f[j] > 1.0)
			{
				f[j] = 1.0;
			}
			else if (f[j] < 0.0)
			{
				f[j] = 0.0;
			}

			// Residue
			res = res + fabs(f[j] - f_old);
			if (fabs(f[j] - f_old) > 1.0E-05)
			{
				N++;
			}
			// Behavior of F
			total_f = total_f + f[j];
			can_melt++;
		}
	}
	*resf = res / N;
	total_f = total_f / can_melt;
	printf("Total f = %5.1E\n", total_f);
}


void average_temperature(int o, int* pp, double* T)
{
	int l;
	double total = 0.0;

	for (int i = 0; i < o; i++)
	{
		l = pp[i];
		total = total + T[l];
	}
	printf("Temperature average = %f", total / o);
}