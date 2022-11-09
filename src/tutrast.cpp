#include <iostream>      // printing
#include <string>
#include <chrono>      // timer

#include <algorithm>

#include <gemmi/ccp4.hpp>
#include <gemmi/asumask.hpp>

#include <connected-components-3d/cc3d.hpp>

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

void printDistinct(uint16_t arr[], int n){
    // Pick all elements one by one
    for (size_t i=0; i<n; i++){
        // Check if the picked element is already printed
        size_t j;
        for (j=0; j<i; j++)
           if (arr[i] == arr[j])
               break;
        // If not printed earlier, then print it
        if (i == j)
          cout << arr[i] << " ";
    }
}

using namespace std;
int main(int argc, char* argv[]) {
  chrono::high_resolution_clock::time_point t_start = chrono::high_resolution_clock::now();
  string grid_file = argv[1];
  double temperature = stod(argv[2]);
  double energy_threshold = stod(argv[3]); //kJ/mol
  gemmi::Ccp4<double> map;
  map.read_ccp4_file(grid_file);
  // double grid_min = map.hstats.dmin;

  // Set up a binary 3D array (is channel)
  size_t array_size = map.grid.nu*map.grid.nv*map.grid.nw;
  int* labels = new int[array_size](); 
  size_t idx = 0;
  for (size_t i=0; i<array_size; i++){
    if (map.grid.data[i] < energy_threshold){
      labels[i] = 1;
      // map.grid.data[i]=-40;
    }
  } 
  size_t N = 0;
  uint16_t* cc_labels = cc3d::connected_components3d<int, uint16_t>(
  labels, /*sx=*/map.grid.nu, /*sy=*/map.grid.nv, /*sz=*/map.grid.nw, /*connectivity=*/7, /*N=*/N );
  cout << N << endl;
  // TODO loop over the unique labels and merge them according to pbc (in future change the code to include pbc)
  // for (uint16_t label_1=1; label_1<N-1; label_1++){
  //   for (uint16_t label_2=label_1+1; label_2<N; label_2++){
  //     if ()
  //   }
  // }
  for (size_t i=0; i<array_size; i++) {
    if (cc_labels[i]==0){map.grid.data[i] = 1e3;}
    else {map.grid.data[i] = cc_labels[i];}
  }
  map.write_ccp4_map("grid/KAXQIL_clean_14_0.1_100_l.ccp4");
  // could be done in the code by applying pbc
  // Calculate diffusion coefficients

  chrono::high_resolution_clock::time_point t_end = chrono::high_resolution_clock::now();
  double elapsed_time_ms = chrono::duration<double, milli>(t_end-t_start).count();
//   cout << grid_file << "," << enthalpy_surface << "," << endl;
  cout << elapsed_time_ms*0.001 << endl;
}