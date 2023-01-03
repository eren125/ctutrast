#include <iostream>      // printing
#include <string>
#include <chrono>      // timer

#include <local/gridcalc.hpp>
#include <local/clustering.hpp>


int main(int argc, char* argv[]) {
  std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();
  std::string grid_file = argv[1];
  double temperature = std::stod(argv[2]);
  double energy_threshold = std::stod(argv[3]); //kJ/mol
  double const R = 8.31446261815324e-3; // kJ/mol/K

  gemmi::Grid<double> grid;
  // READ MAP that took about 870 ms to make for AEI
  read_grid_ccp4(grid, grid_file);
  const size_t V = grid.data.size();
  
  // //Breadth first search to get the connected components
  uint8_t* channel_labels = new uint8_t[V]();
  size_t N = 0;
  bfsOfGraph(&channel_labels, vector<bool>(V, false), grid, energy_threshold, V, N);

  // Array of the type of channel connectivity in X Y Z but not the dimensionality (BFS to do so)
  // Used to filter out pockets no 
  vector<std::string> channel_dimensions=channel_dim_array<uint8_t>(channel_labels, N, grid.nu, grid.nv, grid.nw);
  vector<uint8_t> channels;
  for (uint8_t label=0; label!=N; label++) { 
    if (channel_dimensions[label]!="\0") {
      std::cout << label + 1 << " " << channel_dimensions[label] << std::endl;
      channels.push_back(label+1);
    }
  }
  std::cout << channels.size() << " channels out of " << N << " connected clusters" << std::endl;

  // Vector of channel labels grouped by symmetry
  vector < vector<uint8_t> > channel_unique_labels = sym_unique_labels(grid, channel_labels, channels, error, std::min(0.0,energy_threshold));  print_unique_labels(channel_unique_labels);

  // Loop over the different energy levels
  double energy_step = R*temperature;
  // double energy_step = 1; // kJ/mol
  for (auto labels: channel_unique_labels){
    auto label = labels[0];

    double min_energy = energy_threshold; 
    vector<bool> in_channel(V, true); 
    setup_channel_config(in_channel, channel_labels, label, V, grid.data, min_energy);
    
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
      bfsOfGraph(&bassin_labels_current, in_channel, grid, energy_threshold_temp, V, N_current);
      std::cout << "Step " << (size_t)step << ": Channel " << (size_t)label << " has " << N_current << " components " << energy_threshold_temp << " " ;
      // if N_current changes save it 
      if (N_current == 1){
        vector<std::string> channel_dimensions_temp=channel_dim_array<uint8_t>(bassin_labels_current, N_current, grid.nu, grid.nv, grid.nw);
        if (!channel_dimensions_temp[0].empty()) {std::cout << "MERGED " << std::endl;break;}
      }
      if (step == 0) {
        // Setup the bassin positions & connections & probability
        ;
      }
      else {
        check_merge(bassin_labels_current, bassin_labels_past, N_current, N_past, merged, V);
        if (merged) {std::cout << "MERGED ";}
        else if (N_current != N_past) {std::cout << "CHANGED " ;}
      }
      std::cout << std::endl;
      N_past = N_current;
      bassin_labels_past = dup<uint8_t>(bassin_labels_current,V); 
    }
    delete [] bassin_labels_current;
    delete [] bassin_labels_past;
  }
  
  delete [] channel_labels;
  std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
  double elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end-t_start).count();
  std::cout << elapsed_time_ms*0.001 << std::endl;
}