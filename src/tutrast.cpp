#include <iostream>      // printing
#include <string>

#include <gemmi/grid.hpp>
#include <gemmi/ccp4.hpp>


// Read grid
// Identify tunnel + channel
// Calculate diffusion coefficients
#define R 8.31446261815324e-3 // kJ/mol/K

using namespace std;

int main(int argc, char* argv[]) {
  chrono::high_resolution_clock::time_point t_start = chrono::high_resolution_clock::now();
  string grid_file = argv[1];
  double temperature = stod(argv[2]);

  double boltzmann_energy_lj = 0;
  double sum_exp_energy = 0;

  gemmi::Ccp4<double> map;
  map.read_ccp4_file(grid_file);
  //TODO change it to a fast visit using visited array
  size_t idx = 0;
  for (int w = 0; w != map.grid.nw; ++w)
    for (int v = 0; v != map.grid.nv; ++v)
      for (int u = 0; u != map.grid.nu; ++u, ++idx) {
        double energy = map.grid.data[idx];
        if (energy < 100) {
          double exp_energy = exp(-energy/(R*temperature));
          sum_exp_energy += exp_energy;
          boltzmann_energy_lj += exp_energy * energy;
        }
      }
  double enthalpy_surface = boltzmann_energy_lj/sum_exp_energy - R*temperature;  // kJ/mol
  chrono::high_resolution_clock::time_point t_end = chrono::high_resolution_clock::now();
  double elapsed_time_ms = chrono::duration<double, milli>(t_end-t_start).count();
  cout << grid_file << "," << enthalpy_surface << "," << elapsed_time_ms*0.001 << endl;
}