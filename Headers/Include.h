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

using namespace std;
#pragma warning(disable:4996)

// ----------------------- FUNCTIONS DECLARATIONS ----------------------- //

void map(int type, int N, int Nx, int Ny, double dx, double dy, double D, double t, double l, double t_fin, double t_pcm, double kx_bat, double ky_bat, double k_pcm, double k_alu, double k_cpcm, double kx_gra, double ky_gra, double rho_bat, double rho_pcm, double rho_alu, double rho_cpcm, double rho_gra, double cp_bat, double cp_pcm, double cp_alu, double cp_cpcm, double cp_gra, int* pp, int* R, double* kx, double* ky, double* rho, double* cp);

void assembly(int o, int* pp, int type, int Nx, int Ny, double dx, double dy, int* ww, int* ee, int* nn, int* ss, double* kx, double* ky, double* ap, double* aw, double* ae, double* an, double* as, double* su, double* sp, double* b, double Tw, double Te, double Tn, double Ts, double qw, double qe, double qn, double qs, double* T, double* Ti, double* rho, double* cp, double L_pcm, double L_cpcm, double* f, double* fi, int* R, double Q, double dt, double h_cp, double T_cp, double h_air, double T_air, double w, int active);

void SORt(int o, int* pp, int* ww, int* ee, int* nn, int* ss, double* ap, double* aw, double* ae, double* an, double* as, double* b, double* T, double* Ti);

void SORf(int o, int* pp, int* ww, int* ee, int* nn, int* ss, double* ap, double* aw, double* ae, double* an, double* as, double* b, double* T, double* Ti, int* R, double* f, double* fi, double* rho, double dt, double L_pcm, double L_cpcm, double Tmelt, double* resf);

void average_temperature(int o, int* pp, double* T);

void plot(int o, int Nx, int Ny, double dx, double dy, int* pp, double* P, double time, string dir,
	string property);

void plot_properties(int o, int Nx, int Ny, double dx, double dy, int* pp, int* R, double* kx,
	double* ky, double* rho, double* cp, string dir);

string create_results_folder(string base_folder);

void set_mesh_problem(int type, int& N, int& Nx, int& Ny, double& Lx, double& Ly,
	double dx, double dy, double t, double l, double t_fin, double t_pcm);

int path(int Nx, int Ny, int* pp, int* ee, int* ww, int* nn, int* ss);

void log(double time, int N, int* R, int* pp, double* T, double* f, string results_folder);