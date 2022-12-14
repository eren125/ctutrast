#include <local/forcefield.hpp> // define forcefield parameters
#include <local/boxsetting.hpp> // set rectangular shaped box using periodic boundary conditions
#include <gemmi/cif.hpp>        // file -> cif::Document
#include <gemmi/grid.hpp>

#include <chrono>
#define MAX_EXP 1e3 // exp(-1000) = 0

double energy_lj(double &epsilon, double &sigma_6, double &inv_distance_6, double &inv_cutoff_6, double &inv_distance_12, double &inv_cutoff_12, double const R) {
  return 4*R*epsilon*sigma_6*( sigma_6 * (inv_distance_12 - inv_cutoff_12) - inv_distance_6 + inv_cutoff_6 );
}

string trim(string structure_file) {
  string structure_name(structure_file);
  structure_name = structure_name.substr(structure_name.find_last_of("/\\") + 1);
  string::size_type const p(structure_name.find_last_of('.'));
  structure_name = structure_name.substr(0, p);
  return structure_name;
}

using namespace std;
namespace cif = gemmi::cif;

int main(int argc, char* argv[]) {
  // Set up Input Variables
  chrono::high_resolution_clock::time_point t_start = chrono::high_resolution_clock::now();
  string structure_file = argv[1];
  string forcefield_path = argv[2];
  double temperature = stod(argv[3]);
  double cutoff = stod(argv[4]);
  double cutoff_sq = cutoff*cutoff;
  double cutoff_6 = (cutoff_sq)*(cutoff_sq)*(cutoff_sq);
  double inv_cutoff_6 = 1.0/cutoff_6;
  double inv_cutoff_12 = inv_cutoff_6*inv_cutoff_6;
  string element_guest_str = argv[5];
  double approx_spacing = stod(argv[6]);
  double energy_threshold = 40;
  if (argv[7]) {energy_threshold = stod(argv[7]);}
  double access_coeff = 0.8;
  if (argv[8]) {access_coeff = stod(argv[8]);}

  // Error catch
  if ( temperature < 0 ) {throw invalid_argument( "Received negative value for the Temperature" );}
  if ( energy_threshold < 0 ) {throw invalid_argument( "Received negative value for the Energy Threshold" );}
  if ( access_coeff > 1 || access_coeff < 0 ) {throw invalid_argument( "Accessibility Coefficient out of range (Read the purpose of this coeff)" );}

  // key constants
  double const R = 8.31446261815324e-3; // kJ/mol/K
  double const N_A = 6.02214076e23;    // part/mol
  double const MAX_ENERGY = MAX_EXP*R*temperature; // kJ/mol/K

  // Inialize key variables
  double molar_mass = 0;
  double boltzmann_energy_lj = 0;
  double sum_exp_energy = 0;

  // Read Forcefield Infos
  ForceField::Parameters ff_params;
  if (forcefield_path != "DEFAULT") {
    ff_params.read_lj_from_raspa(forcefield_path);
  }
  pair<double,double> epsilon_sigma = ff_params.get_epsilon_sigma(element_guest_str, true);
  double epsilon_guest = epsilon_sigma.first;
  double sigma_guest = epsilon_sigma.second;

  // Cif read
  cif::Document doc = cif::read_file(structure_file);
  cif::Block block = doc.sole_block();
  gemmi::SmallStructure structure = gemmi::make_small_structure_from_block(block);
  int spacegroup_number = 1;
  for (const char* tag : {"_space_group_IT_number",
                          "_symmetry_Int_Tables_number"})
    if (const string* val = block.find_value(tag)) {
      spacegroup_number = (int) cif::as_number(*val);
      break;
    }
  const gemmi::SpaceGroup* sg = gemmi::find_spacegroup_by_number(spacegroup_number);
  structure.cell.set_cell_images_from_spacegroup(sg);
  vector<gemmi::SmallStructure::Site> unique_sites = structure.sites;
  vector<gemmi::SmallStructure::Site> all_sites = structure.get_all_unit_cell_sites();

  // Cell vector
  double a_x = structure.cell.orth.mat[0][0]; double b_x = structure.cell.orth.mat[0][1]; double c_x = structure.cell.orth.mat[0][2];
  double b_y = structure.cell.orth.mat[1][1]; double c_y = structure.cell.orth.mat[1][2];
  double c_z = structure.cell.orth.mat[2][2];
  // Minimal rectangular box that could interact with atoms within the smaller equivalent rectangluar box
  int n_max = int(abs((cutoff + sigma_guest) / a_x)) + 1; 
  int m_max = int(abs((cutoff + sigma_guest) / b_y)) + 1; 
  int l_max = int(abs((cutoff + sigma_guest) / c_z)) + 1; 

  // Grid set-up
  gemmi::Grid<double> grid;
  bool denser = true;
  grid.spacegroup = sg;
  grid.set_unit_cell(structure.cell);
  grid.set_size_from_spacing(approx_spacing, gemmi::GridSizeRounding::Nearest );
  size_t sample_size = grid.data.size();

  // Auxiliary variables
  string element_host_str;
  double dist = 0; double distance_sq = 0;
  double epsilon = 0;
  double sigma = 0; double sigma_sq = 0; double sigma_6 = 0;

  // center position used to reduce the neighbor list
  gemmi::Position center_pos = gemmi::Position(a_x/2,b_y/2,c_z/2);
  double large_cutoff = cutoff + center_pos.length();
  // Creates a list of sites within the cutoff
  vector<array<double,6>> supracell_sites;
  string element_host_str_temp = "X";
  for (auto site: all_sites) {
    element_host_str = site.type_symbol;
    if (element_host_str != element_host_str_temp) {
      epsilon_sigma = ff_params.get_epsilon_sigma(element_host_str, false);
      // Lorentz-Berthelot
      epsilon = sqrt( epsilon_sigma.first * epsilon_guest );
      sigma = 0.5 * ( epsilon_sigma.second+sigma_guest );
    }
    element_host_str_temp = element_host_str;
    sigma_sq = sigma * sigma;
    sigma_6 = sigma_sq * sigma_sq * sigma_sq;
    // Quick calculation in the occupied space
    grid.use_points_around<true>(site.fract, access_coeff*sigma, [&](double& ref, double d2){
      double inv_distance_6 = 1.0 / (d2*d2*d2);
      double inv_distance_12 = inv_distance_6 * inv_distance_6;
      ref += energy_lj(epsilon,sigma_6,inv_distance_6,inv_cutoff_6,inv_distance_12,inv_cutoff_12, R);
    }, false);
    gemmi::Element el(element_host_str.c_str());
    molar_mass += el.weight();
    // neighbor list within rectangular box
    move_rect_box(site.fract,a_x,b_x,c_x,b_y,c_y);
    gemmi::Fractional coord;
    for (int n = -n_max; (n<n_max+1); ++n)
      for (int m = -m_max; (m<m_max+1); ++m) 
        for (int l = -l_max; (l<l_max+1); ++l) {
          // calculate a distance from centre box
          array<double,6> pos_epsilon_sigma;
          coord.x = site.fract.x + n; coord.y = site.fract.y + m; coord.z = site.fract.z + l;
          gemmi::Position pos = gemmi::Position(structure.cell.orthogonalize(coord));
          double delta_x = abs(center_pos.x-pos.x);
          if (delta_x > large_cutoff) {continue;}
          double delta_y = abs(center_pos.y-pos.y);
          if (delta_y > large_cutoff) {continue;}
          double delta_z = abs(center_pos.z-pos.z);
          if (delta_z > large_cutoff) {continue;}
          distance_sq = delta_x*delta_x+delta_y*delta_y+delta_z*delta_z;
          if (distance_sq > large_cutoff*large_cutoff) {continue;}
          pos_epsilon_sigma[0] = pos.x;
          pos_epsilon_sigma[1] = pos.y;
          pos_epsilon_sigma[2] = pos.z;
          pos_epsilon_sigma[3] = epsilon;
          pos_epsilon_sigma[4] = sigma_sq;
          pos_epsilon_sigma[5] = sigma_6;
          supracell_sites.push_back(pos_epsilon_sigma);
        }
  }

  // Symmetry-aware grid construction
  vector<gemmi::GridOp> grid_ops = grid.get_scaled_ops_except_id();
  vector<bool> visited(sample_size, false); 
  size_t idx = 0;
  for (int w = 0; w != grid.nw; ++w)
    for (int v = 0; v != grid.nv; ++v)
      for (int u = 0; u != grid.nu; ++u, ++idx) {
        if (visited[idx]) {continue;}
        else if (grid.data[idx] >= energy_threshold) {continue;}
        gemmi::Fractional V_fract = grid.get_fractional(u,v,w);
        move_rect_box(V_fract,a_x,b_x,c_x,b_y,c_y);
        gemmi::Position pos = gemmi::Position(structure.cell.orthogonalize(V_fract));
        double energy = 0;
        for(array<double,6> pos_epsilon_sigma : supracell_sites) {
          double energy_temp = 0;
          gemmi::Vec3 pos_neigh = gemmi::Vec3(pos_epsilon_sigma[0], pos_epsilon_sigma[1], pos_epsilon_sigma[2]);
          distance_sq = pos.dist_sq(pos_neigh);
          if (distance_sq < cutoff_sq) {
            sigma_sq = pos_epsilon_sigma[4];
            epsilon = pos_epsilon_sigma[3];
            sigma_6 = pos_epsilon_sigma[5];
            double inv_distance_6 = 1.0 / ( distance_sq * distance_sq * distance_sq );
            double inv_distance_12 = inv_distance_6 * inv_distance_6;
            energy += energy_lj(epsilon, sigma_6, inv_distance_6,inv_cutoff_6, inv_distance_12, inv_cutoff_12, R);
            if (energy>MAX_ENERGY) {break;}
          }
        }
        grid.data[idx] = energy;
        // symmetry
        int sym_count = 1;
        for (size_t k = 0; k < grid_ops.size(); ++k) {
          array<int,3> t = grid_ops[k].apply(u, v, w);
          size_t mate_idx = grid.index_s(t[0], t[1], t[2]);
          if (grid.data[mate_idx]!=0.0)
            continue;
          sym_count++;
          grid.data[mate_idx] = energy;
          visited[mate_idx] = true;
        }
        double exp_energy;
        if (energy>MAX_ENERGY) {exp_energy=0;}
        else {exp_energy = exp(-energy/(R*temperature));}
        sum_exp_energy += sym_count * exp_energy;
        boltzmann_energy_lj += sym_count * exp_energy * energy;
      }

  string structure_name = trim(structure_file);
  double Framework_density = molar_mass/(N_A*structure.cell.volume*1e-30); // g/m3
  double enthalpy_surface = boltzmann_energy_lj/sum_exp_energy - R*temperature;  // kJ/mol
  double henry_surface = sum_exp_energy/(sample_size*R*temperature*Framework_density);    // mol/kg/Pa
  chrono::high_resolution_clock::time_point t_end = chrono::high_resolution_clock::now();
  double elapsed_time_ms = chrono::duration<double, milli>(t_end-t_start).count();
  // Structure name, Enthalpy (kJ/mol), Henry coeff (mol/kg/Pa), Accessible Surface Area (m2/cm3), Time (s)
  cout << structure_name << "," << enthalpy_surface << "," << henry_surface << "," << elapsed_time_ms*0.001 << endl;
}