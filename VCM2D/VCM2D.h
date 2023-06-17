#pragma once

void mesh_regions(double D, double dx, double dy, int Nx, int Ny, int* R);

void plot_regions(int Nx, int Ny, double dx, double dy, char ffn[20], int* R, double* T);

void find_k(double dx, double dy, int* R, int o, int ow, int oe, int on, int os, double& kp, double& kw, double& ke, double& kn, double& ks);

void find_rho(int o, double& rho, double& rhoi, double* T, double* Ti);

void find_f(int o, int& f, int& fi, double* T, double* Ti);

void assembly(int Nx, int Ny, double dx, double dy, int* R, double* ap, double* ae, double* aw, double* an, double* as, double* s, double cp, double L, double* T, double* Ti, double* b, double dt);

void output(int Nx, int Ny, double dx, double dy, int* R, char ffn[20]);

void SOR(int Nx, int Ny, double dx, double dy, int* R, double* ap, double* ae, double* aw, double* an, double* as, double* s, double k_bat, double k_pcm, double rho_pcm_solid, double rho_pcm_liquid, double cp, double L, double* T, double* Ti, double Tmelt, double* b, double dt, int* f, int* fi, double* fres);

bool areAllElementsSmaller(double* ptr, int size, double value);
