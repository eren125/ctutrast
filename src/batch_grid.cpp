#include <local/grid.hpp>
#include <chrono>

int main(int argc, char* argv[]) {
  chrono::high_resolution_clock::time_point t_start = chrono::high_resolution_clock::now();
  string structure_file = argv[1];
  string forcefield_path = argv[2];
  double temperature = stod(argv[3]);
  double cutoff = stod(argv[4]);
  double cutoff_sq = cutoff*cutoff;
  double cutoff_6 = (cutoff_sq)*(cutoff_sq)*(cutoff_sq);
  double inv_cutoff_6 = 1.0/cutoff_6;
  double inv_cutoff_12 = inv_cutoff_6*inv_cutoff_6;
  string element_guest_str = argv[5];
  double approx_spacing = stod(argv[6]);
  double energy_threshold = 40;
  if (argv[7]) {energy_threshold = stod(argv[7]);}
  double access_coeff = 0.85;
  if (argv[8]) {access_coeff = stod(argv[8]);}

  // key constants
  double const R = 8.31446261815324e-3; // kJ/mol/K
  double const N_A = 6.02214076e23;    // part/mol

  // Input
  double molar_mass = 0;
  double boltzmann_energy_lj = 0;
  double sum_exp_energy = 0;
  double volume;
  size_t sample_size;
  gemmi::Grid<double> grid;

  setup_grid(structure_file,forcefield_path,temperature,cutoff,element_guest_str,approx_spacing,energy_threshold, 
  molar_mass, boltzmann_energy_lj, sum_exp_energy, sample_size, volume, grid, access_coeff=access_coeff);

  string structure_name = trim(structure_file);
  double Framework_density = molar_mass/(N_A*volume*1e-30); // g/m3
  double enthalpy_surface = boltzmann_energy_lj/sum_exp_energy - R*temperature;  // kJ/mol
  double henry_surface = sum_exp_energy/(sample_size*R*temperature*Framework_density);    // mol/kg/Pa
  chrono::high_resolution_clock::time_point t_end = chrono::high_resolution_clock::now();
  double elapsed_time_ms = chrono::duration<double, milli>(t_end-t_start).count();
  // Structure name, Enthalpy (kJ/mol), Henry coeff (mol/kg/Pa), Accessible Surface Area (m2/cm3), Time (s)
  cout << structure_name << "," << enthalpy_surface << "," << henry_surface << "," << elapsed_time_ms*0.001 << endl;

}