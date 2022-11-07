#include <iostream>      // printing
#include <string>
#include <chrono>      // timer

#include <algorithm>

#include <gemmi/ccp4.hpp>
#include <gemmi/asumask.hpp>

#define R 8.31446261815324e-3 // kJ/mol/K

using namespace std;

double grid_calc_enthalpy(gemmi::Grid<double> grid, double energy_threshold, double temperature) {
  double boltzmann_energy_lj = 0;
  double sum_exp_energy = 0;
  size_t idx = 0;
  for (int w = 0; w != grid.nw; ++w)
    for (int v = 0; v != grid.nv; ++v)
      for (int u = 0; u != grid.nu; ++u, ++idx) {
        double energy = grid.data[idx];
        if (energy < energy_threshold) {
          double exp_energy = exp(-energy/(R*temperature));
          sum_exp_energy += exp_energy;
          boltzmann_energy_lj += exp_energy * energy;
        }
      }
  return boltzmann_energy_lj/sum_exp_energy - R*temperature;  // kJ/mol
}

double grid_calc_enthalpy_efficient(gemmi::Grid<double> grid, vector<size_t> indexes, double temperature) {
  double boltzmann_energy_lj = 0;
  double sum_exp_energy = 0;
  for (size_t idx: indexes) {
    double energy = grid.data[idx];
    double exp_energy = exp(-energy/(R*temperature));
    sum_exp_energy += exp_energy;
    boltzmann_energy_lj += exp_energy * energy;
    }
  return boltzmann_energy_lj/sum_exp_energy - R*temperature;  // kJ/mol
}

int main(int argc, char* argv[]) {
  chrono::high_resolution_clock::time_point t_start = chrono::high_resolution_clock::now();
  string grid_file = argv[1];
  double temperature = stod(argv[2]);
  double energy_threshold = stod(argv[3]); //kJ/mol

  gemmi::Ccp4<double> map;
  map.read_ccp4_file(grid_file);
  double grid_min = map.hstats.dmin;
  // First identify big channels using an index vector
  vector<size_t> occupied_idx;
  size_t idx = 0;
  for (int w = 0; w != map.grid.nw; ++w)
    for (int v = 0; v != map.grid.nv; ++v)
      for (int u = 0; u != map.grid.nu; ++u, ++idx) {
        double energy = map.grid.data[idx];
        if (energy < energy_threshold) {
          occupied_idx.push_back(idx);
        }
      }

  while (!occupied_idx.empty()) {
    size_t idx = occupied_idx[0];
    gemmi::Fractional fctr = map.grid.point_to_fractional(map.grid.index_to_point(idx));
    double radius = 0.3;
    int du = (int) std::ceil(radius / map.grid.spacing[0]);
    int dv = (int) std::ceil(radius / map.grid.spacing[1]);
    int dw = (int) std::ceil(radius / map.grid.spacing[2]);
    map.grid.use_points_in_box<true>(fctr, du, dv, dw,
        [&](double& ref, const gemmi::Position& delta, int u, int v, int w) {
            map.grid.index_q(u,v,w);
            },true);
    // cout<< point.u << endl;
  }

  // double enthalpy_surface = grid_calc_enthalpy(map.grid, energy_threshold, temperature);
//   double enthalpy_surface = grid_calc_enthalpy_efficient(map.grid, occupied_idx, temperature);
  // Mask use to remove blocked space

  // Calculate diffusion coefficients

  chrono::high_resolution_clock::time_point t_end = chrono::high_resolution_clock::now();
  double elapsed_time_ms = chrono::duration<double, milli>(t_end-t_start).count();
//   cout << grid_file << "," << enthalpy_surface << "," << endl;
  cout << elapsed_time_ms*0.001 << endl;
}