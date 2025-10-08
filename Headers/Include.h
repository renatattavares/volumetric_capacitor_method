// ----------------------------- LIBRARIES ------------------------------ //
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <iomanip>
#include <ctime>
#include <sstream>
#include <direct.h>
#include <omp.h> 
#include <sys/stat.h>

using namespace std;
#pragma warning(disable:4996)

// ----------------------- FUNCTIONS DECLARATIONS ----------------------- //

void map(int type, int o, int Nx, int Ny, double dx, double dy, double D, double t, double l, double t_fin, double t_pcm, double kx_bat, double ky_bat, double k_pcm, double k_alu, double k_cpcm, double kx_gra, double ky_gra, double rho_bat, double rho_pcm, double rho_alu, double rho_cpcm, double rho_gra, double cp_bat, double cp_pcm, double cp_alu, double cp_cpcm, double cp_gra, double L_pcm, double L_cpcm, int* pp, int* R, double* kx, double* ky, double* rho, double* cp, double* L);

void assembly(int o, int* pp, int type, int Nx, int Ny, double dx, double dy, int* ww, int* ee, int* nn, int* ss, double* kx, double* ky, double* ap, double* aw, double* ae, double* an, double* as, double* su, double* sp, double* b, double Tw, double Te, double Tn, double Ts, double qw, double qe, double qn, double qs, double* T, double* Ti, double* rho, double* cp, double* L, double* f, double* fi, int* R, double Q, double dt, double h_cp, double T_cp, double h_air, double T_air, double w, int active);

void SORt(int o, int* pp, int* ww, int* ee, int* nn, int* ss, double* ap, double* aw, double* ae, double* an, double* as, double* b, double* T, double* Ti);

void SORf(int o, int* pp, int* ww, int* ee, int* nn, int* ss, double* ap, double* aw, double* ae, double* an, double* as, double* b, double* T, double* Ti, int* R, double* f, double* fi, double* rho, double* L, double dt, double Tmelt, double* resf);

void SIP(int o, int N, int* pp, int* nn, int* ss, int* ee, int* ww, double* ap, double* ae, double* aw, double* an, double* as, double* b, double* T, double Nx, double Ny, double dx, double dy, double time, string folder);

void average_temperature(int o, int* pp, double* T);

void plot_mesh(int o, int Nx, int Ny, double dx, double dy, int* pp, int* ww, int* ee, int* nn, int* ss, int* R, double* kx, double* ky, double* rho, double* cp, double* L, string dir);

void plot_coef(int o, int Nx, int Ny, double dx, double dy, int* pp, double* aw, double* ae, double* an, double* as, double* ap, double* b, double time, string results_folder);

void plot_sim(int o, int Nx, int Ny, double dx, double dy, int* pp, double* T, double* f, double time, string results_folder);

void set_mesh_problem(int type, int& N, int& Nx, int& Ny, double& Lx, double& Ly,
	double dx, double dy, double t, double l, double t_fin, double t_pcm);

int path(int Nx, int Ny, int* pp, int* ee, int* ww, int* nn, int* ss);

void nonlinear_cond(int o, int* pp, int* R, double* f, double* kx, double* ky, double* resk);

void log(double time, int N, int* R, int* pp, double* T, double* f, string results_folder);

void plot_res(int o, int Nx, int Ny, double dx, double dy, int* pp, double* P, int it, double time, string results_folder);

void log_plot(double time, int o, int* pp, int* R, double* T, double* f);

string create_results_folder();
