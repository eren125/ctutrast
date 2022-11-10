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

string channel_dim(uint16_t* cc_labels, size_t label, size_t nu, size_t nv, size_t nw) {
  // TODO try and replace slice loop by random draws without replacement
  string channels;
  bool present; // checks if a slice contains the label
  // YZ slicing
  for (size_t u = 0; u < nu; u++){
    present = false;
    while (present == false) {
      for (size_t v = 0; v < nv; v++){
        for (size_t w = 0; w < nw; w++){
          size_t loc = size_t(w * nv + v) * nu + u;
          if (cc_labels[loc]==label) {present = true;}
        }
      }
      if (present == false) {break;}
    }
    if (present == false) {break;}
  }
  if (present == true) {channels += 'X';}
  // XZ slicing
  for (size_t v = 0; v < nv; v++){
    present = false;
    while (present == false) {
      for (size_t w = 0; w < nw; w++){
        for (size_t u = 0; u < nu; u++){
          size_t loc = size_t(w * nv + v) * nu + u;
          if (cc_labels[loc]==label) {present = true;}
        }
      }
      if (present == false) {break;}
    }
    if (present == false) {break;}
  }
  if (present == true) {channels += 'Y';}
  // XY slicing
  for (size_t w = 0; w < nw; w++){
    present = false; //if label is in slice
    while (present == false) {
      for (size_t v = 0; v < nv; v++){
        for (size_t u = 0; u < nu; u++){
          size_t loc = size_t(w * nv + v) * nu + u;
          if (cc_labels[loc]==label) {present = true;}
        }
      }
      if (present == false) {break;}
    }
    if (present == false) {break;}
  }
  if (present == true) {channels += 'Z';}
  return channels;
}

using namespace std;
int main(int argc, char* argv[]) {
  chrono::high_resolution_clock::time_point t_start = chrono::high_resolution_clock::now();
  string grid_file = argv[1];
  double temperature = stod(argv[2]);
  double energy_threshold = stod(argv[3]); //kJ/mol
  gemmi::Ccp4<double> map;
  // READ MAP that took about 870 ms to make for AEI
  // (AEI time = 29 ms)
  map.read_ccp4_file(grid_file); 
  double grid_min = map.hstats.dmin;
  // Set up arrays of 0 (if framework) 1 (if void)
  // (AEI time = 11 ms)
  size_t array_size = map.grid.nu*map.grid.nv*map.grid.nw;
  int* channel_labels = new int[array_size](); 
  size_t idx = 0;
  for (size_t i=0; i<array_size; i++){ 
    if (map.grid.data[i] < energy_threshold){
      channel_labels[i] = 1;
    }
  }
  size_t N = 0;
  // Channel labellisation using a 3D connected component algorithm modified to integrate PBC
  // (AEI time = 14ms)
  uint16_t* channel_cc_labels = cc3d::connected_components3d<int, uint16_t, true>(
  channel_labels, /*sx=*/map.grid.nu, /*sy=*/map.grid.nv, /*sz=*/map.grid.nw, /*connectivity=*/26, /*N=*/N );
  cout << N << endl;

  // if channel_dimensions is a null char, it is a pocket  
  // (AEI time = 16ms)
  string channel_dimensions[N+1]; 
  for (uint16_t label=1; label<=N; label++) { 
    channel_dimensions[label] = channel_dim(channel_cc_labels, label, map.grid.nu, map.grid.nv, map.grid.nw);
    cout << channel_dimensions[label] << endl;
  }
  // TODO MAKE 2D array for label indexes to use instead of the grid of labels > then we can replace the previous approach by using the labels instead
  for (size_t i=0; i<array_size; i++){ 

  }

  // remove symmetrical equivalent channels and count them (for the probability calculation)
  // Work with PSI or KAXQIL
  // First approach
  for (uint16_t label_1=1; label_1<N; label_1++) { 
    if (channel_dimensions[label_1] != "") {
      for (uint16_t label_2=label_1+1; label_2<=N; label_2++) {
        if (channel_dimensions[label_2] != "" && 
            channel_dimensions[label_2].length() == channel_dimensions[label_1].length()) {
          // cout << label_1 << " could be equivalent to " << label_2 << endl;
          // function to check if a N randomly selected points of label1 have images in label2
        }
      }
    }
  }

  // Other approach
  // Have an array of sorted index per label for non pocket
  // Make symmetry easily using the grid ops and check if in other label



  // could be done in the code by applying pbc
  // Calculate diffusion coefficients

  chrono::high_resolution_clock::time_point t_end = chrono::high_resolution_clock::now();
  double elapsed_time_ms = chrono::duration<double, milli>(t_end-t_start).count();
//   cout << grid_file << "," << enthalpy_surface << "," << endl;
  cout << elapsed_time_ms*0.001 << endl;
  // // Check the labels
  // for (size_t i=0; i<array_size; i++) {
  //   if (cc_labels[i]==0){map.grid.data[i] = 1e3;}
  //   else {map.grid.data[i] = cc_labels[i];}
  // }
  // map.write_ccp4_map("grid/KAXQIL_clean_14_0.1_100_l.ccp4");
}