#include <iostream>      // printing
#include <string>
#include <chrono>      // timer

#include <algorithm>
#include <bits/stdc++.h>
#include <vector>
#include <queue>

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

vector<string> channel_dim_array(uint16_t* cc_labels, size_t &N_label, int &nu, int &nv, int &nw) {
  vector<string> channel_dimensions(N_label,"\0");
  vector< vector<bool> > present_X_2d(N_label, vector<bool>(nu,0)); 
  vector< vector<bool> > present_Y_2d(N_label, vector<bool>(nv,0));
  vector< vector<bool> > present_Z_2d(N_label, vector<bool>(nw,0));
  size_t loc = 0;
  for (size_t w = 0; w < nw; w++){
    for (size_t v = 0; v < nv; v++){
      for (size_t u = 0; u < nu; u++, ++loc){
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

vector < vector<uint16_t> > sym_unique_labels(gemmi::Grid<double> &grid, uint16_t* cc_labels, vector<uint16_t> channels) {
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

void print_unique_labels(vector < vector<uint16_t> > &unique_labels){
  cout << "Unique channels" << endl;
  int count = 0;
  for (vector<uint16_t> equiv_labels: unique_labels){
    for (uint16_t label: equiv_labels){
      std::cout << label << " ";
    }
    count++;
    std::cout << endl;
  }
}

template <typename T = uint32_t >
tuple<T,T,T> index_to_point(T &idx, T &nu, T &nv, T &nw){
    auto d1 = std::div((ptrdiff_t)idx, (ptrdiff_t)nu);
    auto d2 = std::div(d1.quot, (ptrdiff_t)nv);
    T u = (T) d1.rem;
    T v = (T) d2.rem;
    T w = (T) d2.quot;
  return {u, v, w};
}

template <typename T = uint32_t >
vector< vector<T> > bfsOfGraph(vector<bool> &vis, gemmi::Grid<double> &grid, double &energy_threshold, size_t V) {
  vector< vector<T> > bfs_traversal_clusters;
  T nu=(T)grid.nu, nv=(T)grid.nv, nw=(T)grid.nw; 
  for (int idx = 0; idx != V; ++idx) {
    if (!vis[idx]) {
      vis[idx] = true;
      if (grid.data[idx]<energy_threshold) {
        queue<T> q;
        q.push(idx);
        vector<T> bfs_traversal;
        while (!q.empty()) {
          T g_node = q.front();
          q.pop();
          bfs_traversal.push_back(g_node);
          T u_node, v_node, w_node;
          tie(u_node, v_node, w_node) = index_to_point<T>(g_node,nu,nv,nw);
          for(int8_t a=-1; a!=2;a++)
          for(int8_t b=-1; b!=2;b++)
          for(int8_t c=-1; c!=2;c++){
            if (a==0 && b==0 && c==0){continue;}
            T u_temp, v_temp, w_temp;
            if (a==1 && u_node == nu-1) u_temp = 0; 
            else if (a==-1 && u_node == 0) u_temp = nu - 1;
            else u_temp = (T) u_node + a;
            if (b==1 && v_node == nv-1) v_temp = 0; 
            else if (b==-1 && v_node == 0) v_temp = nv - 1;
            else v_temp = (T) v_node + b;
            if (c==1 && w_node == nw-1) w_temp = 0; 
            else if (c==-1 && w_node == 0) w_temp = nw - 1;
            else w_temp = (T) w_node + c;
            T idx_temp = (w_temp * nv + v_temp) * nu + u_temp;
            if (!vis[idx_temp]) {
              vis[idx_temp] = true;
              if (grid.data[idx_temp]<energy_threshold) {
                q.push(idx_temp);
              }
            }
          }
        }
        bfs_traversal_clusters.push_back(bfs_traversal);
      }
    }
  }
  return bfs_traversal_clusters;
}

template <typename T = uint32_t >
string channel_dim_idx(vector<T> &channel_idx_label, T &nu, T &nv, T &nw) {
  string channels;
  vector<bool> present_X(nu, 0); vector<bool> present_Y(nv, 0); vector<bool> present_Z(nw, 0);
  for (T loc:channel_idx_label){
    T u, v, w;
    tie(u, v, w) = index_to_point<T>(loc,nu,nv,nw);
    present_X[u] = true;present_Y[v] = true;present_Z[w] = true;
  }
  if (all_of(present_X.begin(),present_X.end(),[](bool v){return v;})) {channels+='X';}
  if (all_of(present_Y.begin(),present_Y.end(),[](bool v){return v;})) {channels+='Y';}
  if (all_of(present_Z.begin(),present_Z.end(),[](bool v){return v;})) {channels+='Z';}
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
  // (AEI time = 28 ms)
  map.read_ccp4_file(grid_file); 
  const size_t V = map.grid.nu*map.grid.nv*map.grid.nw;
  // Set up arrays of 0 (if framework) 1 (if void)
  // (AEI time = 11 ms)
  //size_t array_size = map.grid.nu*map.grid.nv*map.grid.nw;
  //uint16_t* channel_labels = new uint16_t[array_size](); 
  //double energy_step = R*temperature; // Can be changed
  //for (size_t i=0; i<array_size; i++){ 
  //  double energy = map.grid.data[i];
  //  if (energy<energy_threshold) {
  //    channel_labels[i] = 1;
  //  } 
  //}
  
  // Channel labelisation using a 3D connected component algorithm modified to integrate PBC
  // (AEI time = 13 ms)
  // size_t N = 0;
  // uint16_t* channel_cc_labels = cc3d::connected_components3d<uint16_t,uint16_t,true>(
  // channel_labels, /*sx=*/map.grid.nu, /*sy=*/map.grid.nv, /*sz=*/map.grid.nw, /*connectivity=*/26, /*N=*/N );
  // delete [] channel_labels;
  // cout << N << endl;

  // //Breadth first search to get the connected components
  vector<bool> vis(V, false);
  vector< vector<uint32_t> > channels_idx = bfsOfGraph(vis, map.grid, energy_threshold, V);
  cout << channels_idx.size() << endl;
  size_t N=channels_idx.size();
  vector<string> channel_dimensions(N, "\0");
  for (uint32_t label=0; label!=N; label++) { 
    channel_dimensions[label] = channel_dim_idx(channel_idx[label-1], (uint32_t)map.grid.nu, (uint32_t)map.grid.nv, (uint32_t) map.grid.nw);
    cout << channel_dimensions[label] << endl;
  }

  // Loop over the different energy levels
  // size_t array_size = map.grid.nu*map.grid.nv*map.grid.nw;
  // uint16_t* channel_labels = new uint16_t[array_size](); 
  // double energy_step = R*temperature; // Can be changed
  // uint16_t max_steps = floor((energy_threshold-map.hstats.dmin)/energy_step);
  // for (size_t i=0; i<array_size; i++){ 
  //   double energy = map.grid.data[i];
  //   if (energy<energy_threshold) {
  //     channel_labels[i] = max_steps - floor((energy - map.hstats.dmin)/energy_step) + 1;
  //     if (energy>map.hstats.dmin+max_steps*energy_step){cout << channel_labels[i]<< " ";}
  //   } // initial labels go from 1 to max_steps+1 according to their closeness to  
  // }


  // if channel_dimensions is a null char, it is a pocket
  // Array of the dimension for each channel and their direction (X,Y,Z)
  // (AEI time = 3.6ms) 
  // vector<string> channel_dimensions=channel_dim_array(channel_cc_labels, N, map.grid.nu, map.grid.nv, map.grid.nw);
  // vector<uint16_t> channels;
  // for (uint16_t label=0; label!=N; label++) { 
  //   if (channel_dimensions[label]!="\0") {
  //     cout << label + 1 << " " << channel_dimensions[label] << endl;
  //     channels.push_back(label+1);
  //   }
  // }
  // cout << channels.size() << " channels out of " << N << " connected clusters" << endl;

  //channel dimension > go from one starting point of a channel apply BFS, and if come back to the starting point 
  // (and traversed a PBC save the direction in which it did)
  // Calc the determinent of this matrix to get dimensionality of the channel

  // Map out unique types of channels and their representing labels (a few labels can constitute the same type of channel)
  // (AEI time = 0.5 ms) 
  // vector < vector<uint16_t> > unique_labels = sym_unique_labels(map.grid, channel_cc_labels, channels);
  // print_unique_labels(unique_labels);

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


  // double threshold_temp = map.hstats.dmin + R*temperature;
  // int* channel_labels_temp = new int[array_size](); 
  // for (size_t i=0; i<array_size; i++){ 
  //   double energy = map.grid.data[i];
  //   if (energy < threshold_temp){channel_labels[i] = 1;}
  // }
  // size_t N = 0;
  // uint16_t* channel_cc_labels = cc3d::connected_components3d<int,uint16_t,true>(
  // channel_labels, /*sx=*/map.grid.nu, /*sy=*/map.grid.nv, /*sz=*/map.grid.nw, /*connectivity=*/26, /*N=*/N );
  // delete [] channel_labels_temp;
  // // Faster but bug
  // vector< vector<uint16_t>> channel_idx(N,{});
  // for (uint16_t i=0; i<array_size; i++){ 
  //   channel_idx[channel_cc_labels[i]-1].push_back(i);
  // }


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
