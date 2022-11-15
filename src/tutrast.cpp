#include <iostream>      // printing
#include <string>
#include <chrono>      // timer

#include <algorithm>
#include <vector>

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
  string channels;
  vector<bool> present_X(nu, 0); vector<bool> present_Y(nv, 0); vector<bool> present_Z(nw, 0);
  for (size_t w = 0; w < nw; w++){
    for (size_t v = 0; v < nv; v++){
      for (size_t u = 0; u < nu; u++){
        size_t loc = size_t(w * nv + v) * nu + u;
        if (cc_labels[loc]==label) {present_X[u] = true;present_Y[v] = true;present_Z[w] = true;}
      }
    }
  }
  if (all_of(present_X.begin(),present_X.end(),[](bool v){return v;})) {channels+='X';}
  if (all_of(present_Y.begin(),present_Y.end(),[](bool v){return v;})) {channels+='Y';}
  if (all_of(present_Z.begin(),present_Z.end(),[](bool v){return v;})) {channels+='Z';}
  return channels;
}

vector<string> channel_dim_array(uint16_t* cc_labels, size_t N_label, size_t nu, size_t nv, size_t nw) {
  vector<string> channel_dimensions(N_label,"\0");
  vector< vector<bool> > present_X_2d(N_label, vector<bool>(nu,0)); 
  vector< vector<bool> > present_Y_2d(N_label, vector<bool>(nv,0));
  vector< vector<bool> > present_Z_2d(N_label, vector<bool>(nw,0));
  for (size_t w = 0; w < nw; w++){
    for (size_t v = 0; v < nv; v++){
      for (size_t u = 0; u < nu; u++){
        size_t loc = size_t(w * nv + v) * nu + u;
        size_t label = cc_labels[loc];
        if (label) {
          label -=1;
          present_X_2d[label][u] = true;
          present_Y_2d[label][v] = true;
          present_Z_2d[label][w] = true;
        }
      }
    }
  }
  for (size_t label=0; label!=N_label; label++){
    if (all_of(present_X_2d[label].begin(),present_X_2d[label].end(),[](bool v){return v;})) 
      channel_dimensions[label]+='X';
    if (all_of(present_Y_2d[label].begin(),present_Y_2d[label].end(),[](bool v){return v;})) 
      channel_dimensions[label]+='Y';
    if (all_of(present_Z_2d[label].begin(),present_Z_2d[label].end(),[](bool v){return v;})) 
      channel_dimensions[label]+='Z';
  }
  return channel_dimensions;
}

// string channel_dim_idx(vector<uint16_t> channel_idx_label, size_t nu, size_t nv, size_t nw) {
//   string channels;
//   vector<bool> present_X(nu, 0); vector<bool> present_Y(nv, 0); vector<bool> present_Z(nw, 0);
//   for (uint16_t loc:channel_idx_label){
//     div_t div_temp = div(loc, nu);
//     uint16_t u = div_temp.rem;
//     div_t div_temp_2 = div(div_temp.quot, nv);
//     uint16_t v = div_temp_2.rem;
//     uint16_t w = div_temp_2.quot;
//     present_X[u] = true;present_Y[v] = true;present_Z[w] = true;
//   }
//   if (all_of(present_X.begin(),present_X.end(),[](bool v){return v;})) {channels+='X';}
//   if (all_of(present_Y.begin(),present_Y.end(),[](bool v){return v;})) {channels+='Y';}
//   if (all_of(present_Z.begin(),present_Z.end(),[](bool v){return v;})) {channels+='Z';}
//   return channels;
// }

vector < vector<uint16_t> > sym_unique_labels(gemmi::Grid<double> grid, uint16_t* cc_labels, vector<uint16_t> channels) {
  vector < vector<uint16_t> > unique_labels;

  vector<gemmi::GridOp> grid_ops = grid.get_scaled_ops_except_id();
  vector<uint16_t>::iterator it;
  vector<uint16_t> labels;
  size_t count = 0;
  size_t idx =0;
  for (int w = 0; w != grid.nw; ++w)
    for (int v = 0; v != grid.nv; ++v)
      for (int u = 0; u != grid.nu; ++u, ++idx) {
        uint16_t label = cc_labels[idx];
        double energy = grid.data[idx];
        it = find (channels.begin(), channels.end(), label);
        if (it != channels.end() && energy<0){ // hard coded change the threshold to a meaningful value
          it = find (labels.begin(), labels.end(), label);
          if (it == labels.end()) {
            vector<uint16_t> equiv_label;
            equiv_label.push_back(label);
            labels.push_back(label);count++;
            for (size_t k = 0; k < grid_ops.size(); ++k) {
              array<int,3> t = grid_ops[k].apply(u, v, w);
              size_t mate_idx = grid.index_s(t[0], t[1], t[2]);
              uint16_t mate_label = cc_labels[mate_idx];
              it = find (equiv_label.begin(), equiv_label.end(), mate_label);
              if (it == equiv_label.end()){
                equiv_label.push_back(mate_label);
                labels.push_back(label);count++;
              }
            }
            unique_labels.push_back(equiv_label);
          }
        }
        if (count == channels.size()) {break;}
        else if (count > channels.size()) {cout << "ERROR in unique channel determination" << endl;}
  }
  return unique_labels;
}

void print_unique_labels(vector < vector<uint16_t> > unique_labels){
  cout << "Unique channels" << endl;
  int count = 0;
  for (vector<uint16_t> equiv_labels: unique_labels){
    for (uint16_t label: equiv_labels){
      cout << label << " ";
    }
    count++;
    cout << endl;
  }
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
  for (size_t i=0; i<array_size; i++){ 
    if (map.grid.data[i] < energy_threshold){channel_labels[i] = 1;}
  }
  
  // Channel labellisation using a 3D connected component algorithm modified to integrate PBC
  // (AEI time = 14ms)
  size_t N = 0;
  uint16_t* channel_cc_labels = cc3d::connected_components3d<int,uint16_t,true>(
  channel_labels, /*sx=*/map.grid.nu, /*sy=*/map.grid.nv, /*sz=*/map.grid.nw, /*connectivity=*/26, /*N=*/N );
  delete [] channel_labels;

  // if channel_dimensions is a null char, it is a pocket
  // Array of the dimension for each channel and their direction (X,Y,Z)
  // (AEI time = 3.8ms) TOCHANGE
  vector<string> channel_dimensions=channel_dim_array(channel_cc_labels, N, map.grid.nu, map.grid.nv, map.grid.nw);
  vector<uint16_t> channels;
  for (uint16_t label=0; label!=N; label++) { 
    if (channel_dimensions[label]!="\0") {
      cout << label + 1 << " " << channel_dimensions[label] << endl;
      channels.push_back(label+1);
    }
  }
  cout << channels.size() << " channels out of " << N << " connected clusters" << endl;

  // Map out unique channel label and their multiplicity (only relevant if there are more 
  // than 2 types of channels)
  vector < vector<uint16_t> > unique_labels = sym_unique_labels(map.grid, channel_cc_labels, channels);
  print_unique_labels(unique_labels);
  
  // TODO we can find the bassin of each channel > min value and check if it is a symmetric image of another channels
  // These bassins are needed anyway
  // Loop over RT (if number of clusters increase then study in detail to detect TS)
  // We can do a dichotomy algorithm until a given precision on energy is reached (10-1 kJ/mol)
  // Definition TS = least energy point that connects two previous clusters


  // chrono::high_resolution_clock::time_point t_a = chrono::high_resolution_clock::now();
  // double time = chrono::duration<double, milli>(t_a-t_start).count();
  // cout << time << endl;
  // // Save label index (AEI time = 3.3 ms)
  // vector<uint16_t> channel_idx[N]; // labels are reindexed from 0 to N-1 (instead of 1 to N)
  // for (size_t idx=0; idx<array_size; idx++){
  //   uint16_t label = channel_cc_labels[idx];
  //   if (label!=0){channel_idx[label-1].push_back(idx);}
  // }
  // t_a = chrono::high_resolution_clock::now();
  // time = chrono::duration<double, milli>(t_a-t_start).count();
  // cout << time << endl;

  // // Faster but bug
  // string channel_dimensions[N];
  // for (uint16_t label=0; label!=N; label++) { 
  //   channel_dimensions[label] = channel_dim_idx(channel_idx[label], map.grid.nu, map.grid.nv, map.grid.nw);
  //   cout << channel_dimensions[label] << endl;
  // }

  // remove symmetrical equivalent channels and count them (for the probability calculation)
  // Work with PSI or KAXQIL
  // First approach
  // for (uint16_t label_1=1; label_1<N; label_1++) { 
  //   if (channel_dimensions[label_1] != "") {
  //     for (uint16_t label_2=label_1+1; label_2<=N; label_2++) {
  //       if (channel_dimensions[label_2] != "" && 
  //           channel_dimensions[label_2].length() == channel_dimensions[label_1].length()) {
  //         // cout << label_1 << " could be equivalent to " << label_2 << endl;
  //         // function to check if a N randomly selected points of label1 have images in label2
  //       }
  //     }
  //   }
  // }

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