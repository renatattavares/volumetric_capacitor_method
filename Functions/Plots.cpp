#include "..\Headers\Include.h"

void plot_properties(int o, int Nx, int Ny, double dx, double dy, int* pp, int* R, double* kx, double* ky, double* rho, double* cp, string dir) {

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

	// Write PPs
	fout << "SCALARS PP int" << endl;
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

void plot_coef(int o, int Nx, int Ny, double dx, double dy, int* pp, double* aw, double* ae, double* an, double* as, double* ap, double* b, double time, string results_folder) {

	double rx{}, ry{};
	string iteration = to_string(int(time));
	string title = "Coefficients �C plot";
	string info = "SCALARS Coefficients float";
	string filename;
	string property = "Coef";

	// Create Results Folder
	string property_folder = results_folder + "\\" + property;

	if (_mkdir(property_folder.c_str()) == 0)
	{
		cout << "Created folder " << property_folder << " succesfully!" << endl;
	}

	// Create filename
	filename = property_folder + "\\" + property + "_" + iteration + ".vtk";

	cout << "\n_________________________________________________________________________" << endl;
	cout << "\nWriting coefficients' files..." << endl;

	// Build
	ofstream fout;
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

	fout << "SCALARS pp int" << endl;
	fout << "LOOKUP_TABLE my_table" << endl;

	for (int a = 0; a < o; a++)
	{
		fout << pp[a] << endl;
	}

	fout << "SCALARS aw float" << endl;
	fout << "LOOKUP_TABLE my_table" << endl;

	for (int a = 0; a < o; a++)
	{
		fout << aw[pp[a]] << endl;
	}

	fout << "SCALARS ae float" << endl;
	fout << "LOOKUP_TABLE my_table" << endl;

	for (int a = 0; a < o; a++)
	{
		fout << ae[pp[a]] << endl;
	}

	fout << "SCALARS an float" << endl;
	fout << "LOOKUP_TABLE my_table" << endl;

	for (int a = 0; a < o; a++)
	{
		fout << an[pp[a]] << endl;
	}

	fout << "SCALARS as float" << endl;
	fout << "LOOKUP_TABLE my_table" << endl;

	for (int a = 0; a < o; a++)
	{
		fout << ap[pp[a]] << endl;
	}

	fout << "SCALARS ap float" << endl;
	fout << "LOOKUP_TABLE my_table" << endl;

	for (int a = 0; a < o; a++)
	{
		fout << b[pp[a]] << endl;
	}

	fout << "SCALARS b float" << endl;
	fout << "LOOKUP_TABLE my_table" << endl;

	for (int a = 0; a < o; a++)
	{
		fout << as[pp[a]] << endl;
	}

	// Exit
	fout.close();
	cout << "\nDone." << endl;
	cout << "_________________________________________________________________________" << endl << endl;

}

void plot(int o, int Nx, int Ny, double dx, double dy, int* pp, double* P, double time, string results_folder, string property) {

	double rx{}, ry{};
	
	string iteration = to_string(int(time));
	string title = "Property �C plot";
	string info = "SCALARS Property float";
	string filename;

	// Create Results Folder
	string property_folder = results_folder + "\\" + property;
		
	if (_mkdir(property_folder.c_str()) == 0) 
	{
		cout << "Created folder " << property_folder << " succesfully!" << endl;
	}

	// Create filename
	filename = property_folder + "\\" + property + "_" + iteration + ".vtk";

	cout << "\n_________________________________________________________________________" << endl;
	cout << "\nWriting file..." << endl;

	// Build VTK
	ofstream fout;
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

	for (int a = 0; a < o; a++)
	{
		fout << P[pp[a]] << endl;
	}

	// Exit
	fout.close();
	cout << "\nDone." << endl;
	cout << "_________________________________________________________________________" << endl << endl;

	return;
}

string create_results_folder() {

	// Get the current time
	time_t t = time(nullptr);
	tm* now = localtime(&t);

	// Create a string stream to format the date and time
	ostringstream oss;
	oss << (now->tm_year + 1900) << '_'
		<< setw(2) << setfill('0') << (now->tm_mon + 1) << '_'
		<< setw(2) << setfill('0') << now->tm_mday << '_'
		<< setw(2) << setfill('0') << now->tm_hour << '_'
		<< setw(2) << setfill('0') << now->tm_min << '_'
		<< setw(2) << setfill('0') << now->tm_sec;

	string dateTimeStr = oss.str();

	// Get current folder
	char cwd[1024];
	getcwd(cwd, 1024);
	string current_folder(cwd);
	
	// Construct results folder
	current_folder = current_folder + "\\Results\\Results_";
	string results_folder = current_folder + dateTimeStr;

	if (_mkdir(results_folder.c_str()) != 0) 
	{
		cout << "Failed creating results folder..." << endl;
		return results_folder;
	}
	return results_folder;
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