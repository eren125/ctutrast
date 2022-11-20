#include <gemmi/grid.hpp>
#include <algorithm>
#include <vector>
#include <queue>

using namespace std;
template <typename T = uint32_t >
tuple<T,T,T> index_to_point(T &idx, T &nu, T &nv, T &nw){
  auto d1 = std::div((ptrdiff_t)idx, (ptrdiff_t)nu);
  auto d2 = std::div(d1.quot, (ptrdiff_t)nv);
  T u = (T) d1.rem;
  T v = (T) d2.rem;
  T w = (T) d2.quot;
  return {u, v, w};
}

template <typename T = uint32_t, typename OUT = uint8_t>
void bfsOfGraph(OUT *cluster_labels[], vector<bool> vis, gemmi::Grid<double> &grid, const double &energy_threshold, const size_t &V, size_t &N) {
  T nu=(T)grid.nu, nv=(T)grid.nv, nw=(T)grid.nw; 
  T idx = 0;
  for (T w = 0; w != nw; ++w)
  for (T v = 0; v != nv; ++v)
  for (T u = 0; u != nu; ++u, ++idx) {
    if (!vis[idx]) {
      vis[idx] = true;
      if (grid.data[idx]<energy_threshold){
        queue< tuple<T,T,T> > q;
        q.push({u, v, w});
        N = N+1;
        while (!q.empty()) {
          T u_node, v_node, w_node;
          tie(u_node, v_node, w_node) = q.front();
          q.pop();
          T g_node = (T) grid.index_q(u_node, v_node, w_node);
          (*cluster_labels)[g_node] = (OUT) N;
          T temp;
          if (u_node == nu-1) temp = 0; else temp = u_node+1;
          T idx_temp = (T) grid.index_q(temp, v_node, w_node);
          if (!vis[idx_temp]) {
            vis[idx_temp] = true;
            if (grid.data[idx_temp]<energy_threshold) q.push({temp, v_node, w_node});
          }
          if (u_node == 0) temp = nu-1; else temp = u_node-1;
          idx_temp = (T) grid.index_q(temp, v_node, w_node);
          if (!vis[idx_temp]) {
            vis[idx_temp] = true;
            if (grid.data[idx_temp]<energy_threshold) q.push({temp, v_node, w_node});
          }
          if (v_node == nv-1) temp = 0; else temp = v_node+1;
          idx_temp = (T) grid.index_q(u_node, temp, w_node);
          if (!vis[idx_temp]) {
            vis[idx_temp] = true;
            if (grid.data[idx_temp]<energy_threshold) q.push({u_node, temp, w_node});
          }
          if (v_node == 0) temp = nv-1; else temp = v_node-1;
          idx_temp = (T) grid.index_q(u_node, temp, w_node);
          if (!vis[idx_temp]) {
            vis[idx_temp] = true;
            if (grid.data[idx_temp]<energy_threshold) q.push({u_node, temp, w_node});
          }
          if (w_node == nw-1) temp = 0; else temp = w_node+1;
          idx_temp = (T) grid.index_q(u_node, v_node, temp);
          if (!vis[idx_temp]) {
            vis[idx_temp] = true;
            if (grid.data[idx_temp]<energy_threshold) q.push({u_node, v_node, temp});
          }
          if (w_node == 0) temp = nw-1; else temp = w_node-1;
          idx_temp = (T) grid.index_q(u_node, v_node, temp);
          if (!vis[idx_temp]) {
            vis[idx_temp] = true;
            if (grid.data[idx_temp]<energy_threshold) q.push({u_node, v_node, temp});
          }
        }
      }
    }
  }
}

template <typename T = uint8_t >
vector<string> channel_dim_array(T* cc_labels, size_t N_label, int &nu, int &nv, int &nw) {
  vector<string> channel_dimensions(N_label);
  vector< vector<bool> > present_X_2d(N_label, vector<bool>(nu,0)); 
  vector< vector<bool> > present_Y_2d(N_label, vector<bool>(nv,0));
  vector< vector<bool> > present_Z_2d(N_label, vector<bool>(nw,0));
  size_t loc = 0;
  for (size_t w = 0; w < nw; w++){
    for (size_t v = 0; v < nv; v++){
      for (size_t u = 0; u < nu; u++, ++loc){
        T label = cc_labels[loc];
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

template <typename T = uint8_t >
vector < vector<T> > sym_unique_labels(gemmi::Grid<double> &grid, T* cc_labels, vector<T> &channels, double threshold=0) {
  vector < vector<T> > unique_labels;
  vector<gemmi::GridOp> grid_ops = grid.get_scaled_ops_except_id();
  vector<T> labels;
  size_t count = 0;
  size_t idx =0;
  for (int w = 0; w != grid.nw; ++w)
    for (int v = 0; v != grid.nv; ++v)
      for (int u = 0; u != grid.nu; ++u, ++idx) {
        T label = cc_labels[idx];
        double energy = grid.data[idx];
        auto it = find (channels.begin(), channels.end(), label);
        if (it != channels.end() && energy<threshold){ // hard coded change the threshold to a meaningful value
          it = find (labels.begin(), labels.end(), label);
          if (it == labels.end()) {
            vector<T> equiv_label;
            equiv_label.push_back(label);
            labels.push_back(label);count++;
            for (size_t k = 0; k < grid_ops.size(); ++k) {
              array<int,3> t = grid_ops[k].apply(u, v, w);
              size_t mate_idx = grid.index_s(t[0], t[1], t[2]);
              T mate_label = cc_labels[mate_idx];
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

template <typename T = uint8_t >
void print_unique_labels(vector < vector<T> > &unique_labels){
  cout << "Unique channels" << endl;
  int count = 0;
  for (vector<T> equiv_labels: unique_labels){
    for (T label: equiv_labels){
      std::cout << (size_t) label << " ";
    }
    count++;
    std::cout << endl;
  }
}