#include <local/gridcalc.hpp>
#include <chrono>

int main(int argc, char* argv[]) {
  std::chrono::high_resolution_clock::time_point t_start = chrono::high_resolution_clock::now();
  std::string structure_file = argv[1];
  std::string forcefield_path = argv[2];
  double temperature = stod(argv[3]);
  double cutoff = stod(argv[4]);
  double cutoff_sq = cutoff*cutoff;
  double cutoff_6 = (cutoff_sq)*(cutoff_sq)*(cutoff_sq);
  double inv_cutoff_6 = 1.0/cutoff_6;
  double inv_cutoff_12 = inv_cutoff_6*inv_cutoff_6;
  std::string element_guest_str = argv[5];
  double approx_spacing = stod(argv[6]);
  double energy_threshold = 40;
  if (argv[7]) {energy_threshold = stod(argv[7]);}
  double access_coeff = 0.8;
  if (argv[8]) {access_coeff = stod(argv[8]);}

  // Error catch
  if ( temperature < 0 ) {throw std::invalid_argument( "Received negative value for the Temperature" );}
  if ( energy_threshold < 0 ) {throw std::invalid_argument( "Received negative value for the Energy Threshold" );}
  if ( access_coeff > 1 || access_coeff < 0 ) {throw std::invalid_argument( "Accessibility Coefficient out of range (Read the purpose of this coeff)" );}

  // key constants
  double const R = 8.31446261815324e-3; // kJ/mol/K
  double const N_A = 6.02214076e23;    // part/mol

  // Input
  double molar_mass = 0;
  double boltzmann_energy_lj = 0;
  double sum_exp_energy = 0;
  gemmi::Grid<double> grid;

  // read_grid_ccp4(grid, grid_file);

  make_energy_grid_ads(structure_file,forcefield_path,temperature,cutoff,element_guest_str,approx_spacing,energy_threshold,access_coeff, 
  molar_mass, boltzmann_energy_lj, sum_exp_energy, grid, true);

  // save_grid_ccp4(grid, structure_file, approx_spacing, energy_threshold, element_guest_str);

  std::string structure_name = trim(structure_file);
  double Framework_density = molar_mass/(N_A*grid.unit_cell.volume*1e-30); // g/m3
  double enthalpy = boltzmann_energy_lj/sum_exp_energy - R*temperature;  // kJ/mol
  double henry = sum_exp_energy/(grid.data.size()*R*temperature*Framework_density);    // mol/kg/Pa
  auto t_end = std::chrono::high_resolution_clock::now();
  double elapsed_time_ms = std::chrono::duration<double, milli>(t_end-t_start).count();
  // Structure name, Enthalpy (kJ/mol), Henry coeff (mol/kg/Pa), Accessible Surface Area (m2/cm3), Time (s)
  cout << structure_name << "," << enthalpy << "," << henry << "," << elapsed_time_ms*0.001 << endl;
}