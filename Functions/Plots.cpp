#include "..\Headers\Include.h"

void plot_properties(int o, int Nx, int Ny, double dx, double dy, int* pp, int* R, 
	double* kx, double* ky, double* rho, double* cp, string dir) {

	double rx{}, ry{};

	cout << "\n_________________________________________________________________________" << endl;
	cout << "\nWriting data files..." << endl;

	// Build
	ofstream fout;
	string filename = dir + "/Properties.vtk";
	fout.open(filename);
	fout << "# vtk DataFile Version 2.0" << endl;
	fout << "Physical Properties plot" << endl;
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

	// Write IDs
	fout << "SCALARS IDs int" << endl;
	fout << "LOOKUP_TABLE my_table" << endl;

	for (int a = 0; a < o; a++)
	{
		fout << pp[a] << endl;
	}

	// Write Region Information
	fout << "SCALARS Region int" << endl;
	fout << "LOOKUP_TABLE my_table" << endl;

	for (int a = 0; a < o; a++) 
	{
		fout << R[pp[a]] << endl;
	}

	// Write Conductivity Information
	fout << "SCALARS Conductivity_X float" << endl;
	fout << "LOOKUP_TABLE my_table" << endl;

	for (int a = 0; a < o; a++) 
	{
		fout << kx[pp[a]] << endl;
	}

	// Write Conductivity Information
	fout << "SCALARS Conductivity_Y float" << endl;
	fout << "LOOKUP_TABLE my_table" << endl;

	for (int a = 0; a < o; a++)
	{
		fout << ky[pp[a]] << endl;
	}

	// Write Density Information
	fout << "SCALARS Density float" << endl;
	fout << "LOOKUP_TABLE my_table" << endl;

	for (int a = 0; a < o; a++)
	{
		fout << rho[pp[a]] << endl;
	}

	// Write Specific Heat Information
	fout << "SCALARS Specific_Heat float" << endl;
	fout << "LOOKUP_TABLE my_table" << endl;

	for (int a = 0; a < o; a++)
	{
		fout << cp[pp[a]] << endl;
	}

	// Exit
	fout.close();
	cout << "\nDone." << endl;
	cout << "_________________________________________________________________________" << endl << endl;

}



void plot(int o, int Nx, int Ny, double dx, double dy, int* pp, double* P, double time, 
	string dir, string property) {

	double rx{}, ry{};
	ofstream fout;
	string iteration = to_string(int(time));
	string filename;
	string title;
	string info;

	// Identify property
	if (property == "Temp") {
		filename = dir + "/Temp/Temp_" + iteration + ".vtk";
		title = "Temperature °C plot";
		info = "SCALARS Temperature float";
	}
	else if (property == "Frac") {
		filename = dir + "/Frac/Frac_" + iteration + ".vtk";
		title = "Liquid Fraction Mass plot";
		info = "SCALARS Liquid_Fraction float";
	}
	else {
		return;
	}

	cout << "\n_________________________________________________________________________" << endl;
	cout << "\nWriting file..." << endl;

	fout.open(filename);
	fout << "# vtk DataFile Version 2.0" << endl;
	fout << title << endl;
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

	// Write Property Information
	fout << info << endl;
	fout << "LOOKUP_TABLE my_table" << endl;

	for (int a = 0; a < o; a++) {

		fout << P[pp[a]] << endl;

	}

	// Exit
	fout.close();
	cout << "\nDone." << endl;
	cout << "_________________________________________________________________________" << endl << endl;

	return;
}


string create_results_folder(string base_folder) {
	// Get the current time
	time_t t = time(nullptr);
	tm* now = localtime(&t);

	// Create a string stream to format the date and time
	ostringstream oss;
	oss << (now->tm_year + 1900) << '-'
		<< setw(2) << setfill('0') << (now->tm_mon + 1) << '-'
		<< setw(2) << setfill('0') << now->tm_mday << '_'
		<< setw(2) << setfill('0') << now->tm_hour << '-'
		<< setw(2) << setfill('0') << now->tm_min;

	// The formatted date-time string
	string dateTimeStr = oss.str();

	// The folder string
	string folder = base_folder + dateTimeStr;
	string tfolder = folder + "/Temp";
	string ffolder = folder + "/Frac";

	if (_mkdir(folder.c_str()) != 0 || _mkdir(tfolder.c_str()) != 0 || _mkdir(ffolder.c_str()) != 0) {
		cout << "Failed creating results folder..." << endl;
		return folder;
	}

	return folder;

}


void log(double time, int o, int* pp, int* R, double* T, double* f, string results_folder) {

	double delta_temp = 0.0;
	double max_temp = 0.0;
	double min_temp = 1000.0;
	double total_f = 0.0;
	int can_melt = 0;
	int l = 0;

	// Log file
	ofstream fout;
	string filename;
	filename = results_folder + "/Log.dat";

	for (int i = 0; i < o; i++)
	{
		l = pp[i];
		if (T[l] > max_temp && R[l] == 1)
		{
			max_temp = T[l];
		}
		if (T[l] < min_temp && R[l] == 1)
		{
			min_temp = T[l];
		}
		if (R[l] == 0 || R[l] == 5)
		{
			total_f = total_f + f[l];
			can_melt++;
		}
	}

	delta_temp = max_temp - min_temp;
	total_f = total_f / can_melt;

	fout.open(filename, std::ios_base::app);
	fout << time << '\t' << max_temp << '\t' << min_temp << '\t' << delta_temp << '\t' << total_f << endl;
	fout.close();

	//cout << "Max temperature:" << max_temp << endl;
	//cout << "Min temperature:" << min_temp << endl;
	//cout << "Delta temperature:" << delta_temp << endl;
	//cout << "Liquid Fraction:" << total_f << endl;

}