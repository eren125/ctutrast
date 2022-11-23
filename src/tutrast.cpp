#include <iostream>      // printing
#include <string>
#include <chrono>      // timer

#include <local/clustering.hpp>

#include <gemmi/ccp4.hpp>
#include <gemmi/asumask.hpp>

#include <connected-components-3d/cc3d.hpp>

#define R 8.31446261815324e-3 // kJ/mol/K

using namespace std;

double grid_calc_enthalpy(gemmi::Grid<double> &grid, double &energy_threshold, double &temperature) {
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


using namespace std;
int main(int argc, char* argv[]) {
  chrono::high_resolution_clock::time_point t_start = chrono::high_resolution_clock::now();
  string grid_file = argv[1];
  double temperature = stod(argv[2]);
  double energy_threshold = stod(argv[3]); //kJ/mol
  gemmi::Ccp4<double> map;
  // READ MAP that took about 870 ms to make for AEI
  // (AEI time = 13 ms)
  map.read_ccp4_file(grid_file); 
  const size_t V = map.grid.nu*map.grid.nv*map.grid.nw;

  // //Breadth first search to get the connected components
  uint8_t* channel_labels = new uint8_t[V]();
  size_t N = 0;
  bfsOfGraph(&channel_labels, vector<bool>(V, false), map.grid, energy_threshold, V, N);

  // Array of the type of channel connectivity in X Y Z but not the dimensionality (BFS to do so)
  // Used to filter out pockets no 
  vector<string> channel_dimensions=channel_dim_array<uint8_t>(channel_labels, N, map.grid.nu, map.grid.nv, map.grid.nw);
  vector<uint8_t> channels;
  for (uint8_t label=0; label!=N; label++) { 
    if (channel_dimensions[label]!="\0") {
      cout << label + 1 << " " << channel_dimensions[label] << endl;
      channels.push_back(label+1);
    }
  }
  cout << channels.size() << " channels out of " << N << " connected clusters" << endl;

  // Vector of channel labels grouped by symmetry
  vector < vector<uint8_t> > channel_unique_labels = sym_unique_labels(map.grid, channel_labels, channels, min(0.0,energy_threshold));
  print_unique_labels(channel_unique_labels);
  // if there are several types of channels: we need to save the weight of each channel and do kMC in each (TODO)

  // Loop over the different energy levels
  double energy_step = R*temperature;
  // double energy_step = 1; // kJ/mol
  for (auto labels: channel_unique_labels){
    auto label = labels[0];
    double min_energy = energy_threshold;
    vector<bool> in_channel = setup_channel_config(channel_labels, label, V, map.grid.data, min_energy);
    double energy_threshold_temp = min_energy + 0.01;
    size_t max_steps = floor((energy_threshold-energy_threshold_temp)/energy_step);
    size_t N_current; size_t N_past=0; 
    // calc weight here later (TODO)
    uint8_t* bassin_labels_current = new uint8_t[V]();
    uint8_t* bassin_labels_past = new uint8_t[V]();
    for (size_t step=0; step<max_steps+1; step++){ 
      energy_threshold_temp += energy_step;
      N_current = 0;
      bool merged = false;
      bfsOfGraph(&bassin_labels_current, in_channel, map.grid, energy_threshold_temp, V, N_current);
      cout << "Step " << (size_t)step << ": Channel " << (size_t)label << " has " << N_current << " components " << energy_threshold_temp << " " ;
      // if N_current changes save it 
      if (N_current == 1){
        vector<string> channel_dimensions_temp=channel_dim_array<uint8_t>(bassin_labels_current, N_current, map.grid.nu, map.grid.nv, map.grid.nw);
        if (!channel_dimensions_temp[0].empty()) {cout << "MERGED " << endl;break;}
      }
      if (step == 0) {
        // Setup the bassin positions & connections & probability
        ;
      }
      else {
        check_merge(bassin_labels_current, bassin_labels_past, N_current, N_past, merged, V);
        if (merged) {cout << "MERGED ";}
        else if (N_current != N_past) {cout << "CHANGED " ;}
      }
      cout << endl;
      N_past = N_current;
      bassin_labels_past = dup<uint8_t>(bassin_labels_current,V); 
    }
    delete [] bassin_labels_current;
    delete [] bassin_labels_past;
  }
  
  // TODO Save the TS and the bassins (displacement+energies)

  //channel dimension > go from one starting point of a channel apply BFS, and if come back to the starting point 
  // (and traversed a PBC save the direction in which it did)
  // Calc the determinent of this matrix to get dimensionality of the channel

  // chrono::high_resolution_clock::time_point t_a = chrono::high_resolution_clock::now();
  // double time = chrono::duration<double, milli>(t_a-t_start).count();
  // cout << time << endl;
  // t_a = chrono::high_resolution_clock::now();
  // time = chrono::duration<double, milli>(t_a-t_start).count();
  // cout << time << endl;

  // TODO we can find the bassin of each channel > min value and check if it is a symmetric image of 
  // another channels
  // These bassins are needed anyway
  // Loop over RT (if number of clusters increase then study in detail to detect TS)
  // We can do a dichotomy algorithm until a given precision on energy is reached (10-1 kJ/mol)
  // Definition TS = least energy point that connects two previous clusters

  delete [] channel_labels;
  // Calculate diffusion coefficients
  chrono::high_resolution_clock::time_point t_end = chrono::high_resolution_clock::now();
  double elapsed_time_ms = chrono::duration<double, milli>(t_end-t_start).count();
  // double enthalpy_surface=grid_calc_enthalpy(map.grid, energy_threshold, temperature);
  // cout << grid_file << "," << enthalpy_surface << "," << endl;
  cout << elapsed_time_ms*0.001 << endl;

  // // Check the labels
  // for (size_t i=0; i<array_size; i++) {
  //   if (cc_labels[i]==0){map.grid.data[i] = 1e3;}
  //   else {map.grid.data[i] = cc_labels[i];}
  // }
  // map.write_ccp4_map("grid/KAXQIL_clean_14_0.1_100_l.ccp4");
}
