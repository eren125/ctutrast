#include <iostream>      // printing
#include <string>

#include <gemmi/grid.hpp>
#include <gemmi/ccp4.hpp>

#define R 8.31446261815324e-3 // kJ/mol/K

double grid_calc_enthalpy(gemmi::Grid<double> grid, double temperature) {
  double boltzmann_energy_lj = 0;
  double sum_exp_energy = 0;
  size_t idx = 0;
  for (int w = 0; w != grid.nw; ++w)
    for (int v = 0; v != grid.nv; ++v)
      for (int u = 0; u != grid.nu; ++u, ++idx) {
        double energy = grid.data[idx];
        if (energy < 100) {
          double exp_energy = exp(-energy/(R*temperature));
          sum_exp_energy += exp_energy;
          boltzmann_energy_lj += exp_energy * energy;
        }
      }
  return boltzmann_energy_lj/sum_exp_energy - R*temperature;  // kJ/mol
}

using namespace std;
int main(int argc, char* argv[]) {
  chrono::high_resolution_clock::time_point t_start = chrono::high_resolution_clock::now();
  string grid_file = argv[1];
  double temperature = stod(argv[2]);
  gemmi::Ccp4<double> map;
  map.read_ccp4_file(grid_file);
  double enthalpy_surface = grid_calc_enthalpy(map.grid, temperature);

  // First identify big channels 
  double energy_threshold = 20; //kJ/mol
  size_t idx = 0;
  for (int w = 0; w != map.grid.nw; ++w)
    for (int v = 0; v != map.grid.nv; ++v)
      for (int u = 0; u != map.grid.nu; ++u, ++idx) {
        double energy = map.grid.data[idx];
        if (energy < energy_threshold) {
          continue;
        }
      }
  // Mask use to remove blocked space

  // Calculate diffusion coefficients

  chrono::high_resolution_clock::time_point t_end = chrono::high_resolution_clock::now();
  double elapsed_time_ms = chrono::duration<double, milli>(t_end-t_start).count();
  cout << grid_file << "," << enthalpy_surface << "," << elapsed_time_ms*0.001 << endl;
}