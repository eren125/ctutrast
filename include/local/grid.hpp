#include <local/forcefield.hpp> // define forcefield parameters
#include <local/boxsetting.hpp> // set rectangular shaped box using periodic boundary conditions
#include <gemmi/cif.hpp>        // file -> cif::Document
#include <gemmi/grid.hpp>       // grid construction

#define MAX_EXP 1e3 // exp(-1000) = 1.34*10^434, for argument above we skip calculation

double energy_lj(double &epsilon, double sigma_6, double &inv_distance_6, double &inv_cutoff_6, double &inv_distance_12, double &inv_cutoff_12, double const R) {
  return 4*R*epsilon*sigma_6*( sigma_6 * (inv_distance_12 - inv_cutoff_12) - inv_distance_6 + inv_cutoff_6 );
}

double energy_lj_opt(double &epsilon, double sigma_6, double &inv_distance_6, double &inv_cutoff_6, double &inv_distance_12, double &inv_cutoff_12) {
  return epsilon*sigma_6*( sigma_6 * (inv_distance_12 - inv_cutoff_12) - inv_distance_6 + inv_cutoff_6 );
}

std::string trim(std::string structure_file) {
  string structure_name(structure_file);
  structure_name = structure_name.substr(structure_name.find_last_of("/\\") + 1);
  string::size_type const p(structure_name.find_last_of('.'));
  structure_name = structure_name.substr(0, p);
  return structure_name;
}

gemmi::SmallStructure extract_structure(std::string &structure_file, int &spacegroup_number) {
  gemmi::cif::Document doc = gemmi::cif::read_file(structure_file);
  gemmi::cif::Block block = doc.sole_block();
  for (const char* tag : {"_space_group_IT_number",
                          "_symmetry_Int_Tables_number"})
    if (const std::string* val = block.find_value(tag)) {
      spacegroup_number = (int) gemmi::cif::as_number(*val);
      break;
    }
  return gemmi::make_small_structure_from_block(block);
}

void initialize_grid(gemmi::Grid<double> &grid, vector<gemmi::SmallStructure::Site> &all_sites, std::string &structure_file, double &approx_spacing) {
  // Cif read
  int spacegroup_number = 1;
  gemmi::SmallStructure structure = extract_structure(structure_file, spacegroup_number);
  const gemmi::SpaceGroup* sg = gemmi::find_spacegroup_by_number(spacegroup_number);
  structure.cell.set_cell_images_from_spacegroup(sg);
  all_sites = structure.get_all_unit_cell_sites();
  // Grid set-up
  grid.spacegroup = sg;
  grid.set_unit_cell(structure.cell);
  grid.set_size_from_spacing(approx_spacing, gemmi::GridSizeRounding::Nearest);
  gemmi::check_grid_factors(sg, {{grid.nu, grid.nv, grid.nw}});
}

void set_batch_cells(double &approx_block_size, gemmi::Grid<double> &grid){
  int block_size = (int) approx_block_size/grid.min_spacing();
  int Nu_block = (int) grid.nu/block_size;
  int Ru_block = (int) grid.nu%block_size;
  int Nv_block = (int) grid.nv/block_size;
  int Rv_block = (int) grid.nv%block_size;
  int Nw_block = (int) grid.nw/block_size;
  int Rw_block = (int) grid.nw%block_size;
  std::cout << grid.nu << " " << grid.nv << " " << grid.nw << std::endl;
  std::cout << block_size << std::endl;
  std::cout << Nu_block << " " << Ru_block << std::endl;
  std::cout << Nv_block << " " << Rv_block << std::endl;
  std::cout << Nw_block << " " << Rw_block << std::endl;
}

void energy_grid_opt(std::string &structure_file, std::string &forcefield_path, double &temperature, double &cutoff, std::string &element_guest_str, double &approx_spacing, double &energy_threshold, double &access_coeff, double &molar_mass, double &boltzmann_energy_lj, double &sum_exp_energy, size_t &sample_size, gemmi::Grid<double> &grid, bool calc_block=true) {
  /**
   * Construct the grid while calculating key figures for the Henry constant and Adsorption enthalpy
   * Block by block
   * 
   * @param structure_file [INPUT] PATH to structure file
   * @param forcefield_path [INPUT] PATH to forcefield file
   * @param temperature [INPUT] temperature for Boltzmann averaging
   * @param cutoff [INPUT] cutoff for LJ energy calculation
   * @param element_guest_str [INPUT] Guest element atomic element (only monoatomic now)
   * @param approx_spacing [INPUT] Spacing of the grid
   * @param energy_threshold [INPUT] Energy threshold under which the energy is recalculated and included in the sum
   * @param access_coeff [INPUT] Relative size of the radius that is blocked for each atoms
   * @param molar_mass [OUTPUT] Molar mass of the entire framework (used in density calc)
   * @param boltzmann_energy_lj [OUTPUT] sum of E exp(-E/RT)
   * @param sum_exp_energy [OUTPUT] sum of exp(-E/RT)
   * @param sample_size [OUTPUT] Total number of grid points = Number of samples
   * @param grid [OUTPUT] Interaction Energy Grid in gemmi format
   * 
   */
  double const R = 8.31446261815324e-3; // kJ/mol/K
  double const MAX_ENERGY = MAX_EXP*R*temperature;
  double const MAX_ENERGY_opt = MAX_EXP*temperature/4; 
  double cutoff_sq = cutoff*cutoff;
  double cutoff_6 = cutoff_sq*cutoff_sq*cutoff_sq;
  double inv_cutoff_6 = 1.0/cutoff_6;
  double inv_cutoff_12 = inv_cutoff_6*inv_cutoff_6;
  // Read Forcefield Infos
  LennardJones::Parameters ff_params;
  if (forcefield_path != "DEFAULT_UFF") {ff_params.read_lj_from_raspa(forcefield_path);}
  double epsilon_guest, sigma_guest;
  tie(epsilon_guest,sigma_guest) = ff_params.get_epsilon_sigma(element_guest_str, true);
  // Grid infrastructure
  vector<gemmi::SmallStructure::Site> all_sites;
  initialize_grid(grid, all_sites, structure_file, approx_spacing);
  sample_size = grid.data.size();
  // Cell vector
  double a_x = grid.unit_cell.orth.mat[0][0], b_x = grid.unit_cell.orth.mat[0][1], c_x = grid.unit_cell.orth.mat[0][2];
  double b_y = grid.unit_cell.orth.mat[1][1], c_y = grid.unit_cell.orth.mat[1][2];
  double c_z = grid.unit_cell.orth.mat[2][2];
  // Minimal rectangular box that could interact with atoms within the smaller equivalent rectangluar box
  int n_max = int(std::abs((cutoff + sigma_guest) / a_x)) + 1;
  int m_max = int(std::abs((cutoff + sigma_guest) / b_y)) + 1;
  int l_max = int(std::abs((cutoff + sigma_guest) / c_z)) + 1;
  // center position used to reduce the neighbor list
  gemmi::Position center_pos = gemmi::Position(a_x/2,b_y/2,c_z/2);
  double large_cutoff = cutoff + center_pos.length();
  // Batch_size
  double approx_block_size = 3; // angstrom
  // set_batch_cells(approx_block_size, grid);
  // Creates a list of sites within the cutoff
  std::vector<std::array<double,6>> supracell_sites;
  std::string element_host_str_temp = "X";
  double epsilon = 0, sigma = 0;
  for (auto site: all_sites) {
    std::string element_host_str = site.type_symbol;
    if (element_host_str != element_host_str_temp) {
      double epsilon_host, sigma_host;
      tie(epsilon_host,sigma_host) = ff_params.get_epsilon_sigma(element_host_str, false);
      // Lorentz-Berthelot
      epsilon = sqrt( epsilon_host * epsilon_guest );
      sigma = 0.5 * ( sigma_host + sigma_guest );
    }
    element_host_str_temp = element_host_str;
    double sigma_sq = sigma * sigma;
    double sigma_6 = sigma_sq * sigma_sq * sigma_sq;
    // Remove blocked points
    grid.use_points_around<true>(site.fract, access_coeff*sigma, [&](double& ref, double d2){
      if (calc_block){
        double inv_distance_6 = 1.0 / (d2*d2*d2);
        double inv_distance_12 = inv_distance_6 * inv_distance_6;
        ref += energy_lj(epsilon,sigma_6,inv_distance_6,inv_cutoff_6,inv_distance_12,inv_cutoff_12, R);
      }
      else {ref = MAX_ENERGY;}
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
          std::array<double,6> pos_epsilon_sigma;
          coord.x = site.fract.x + n; coord.y = site.fract.y + m; coord.z = site.fract.z + l;
          gemmi::Position pos = gemmi::Position(grid.unit_cell.orthogonalize(coord));
          double delta_x = abs(center_pos.x-pos.x);
          if (delta_x > large_cutoff) {continue;}
          double delta_y = abs(center_pos.y-pos.y);
          if (delta_y > large_cutoff) {continue;}
          double delta_z = abs(center_pos.z-pos.z);
          if (delta_z > large_cutoff) {continue;}
          double distance_sq = delta_x*delta_x+delta_y*delta_y+delta_z*delta_z;
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
  std::vector<gemmi::GridOp> grid_ops = grid.get_scaled_ops_except_id();
  vector<bool> visited(sample_size, false); 
  size_t idx = 0;
  for (int w = 0; w != grid.nw; ++w)
    for (int v = 0; v != grid.nv; ++v)
      for (int u = 0; u != grid.nu; ++u, ++idx) {
        if (visited[idx]) {continue;}
        else if (grid.data[idx] >= energy_threshold) {continue;}
        gemmi::Fractional V_fract = grid.get_fractional(u,v,w);
        move_rect_box(V_fract,a_x,b_x,c_x,b_y,c_y);
        gemmi::Position pos = gemmi::Position(grid.unit_cell.orthogonalize(V_fract));
        double energy = 0; bool skip = false;
        for(std::array<double,6> pos_epsilon_sigma : supracell_sites) {
          double energy_temp = 0;
          gemmi::Vec3 pos_neigh = gemmi::Vec3(pos_epsilon_sigma[0], pos_epsilon_sigma[1], pos_epsilon_sigma[2]);
          double distance_sq = pos.dist_sq(pos_neigh);
          if (distance_sq < cutoff_sq) {
            double sigma_sq = pos_epsilon_sigma[4];
            epsilon = pos_epsilon_sigma[3];
            double sigma_6 = pos_epsilon_sigma[5];
            double inv_distance_6 = 1.0 / ( distance_sq * distance_sq * distance_sq );
            double inv_distance_12 = inv_distance_6 * inv_distance_6;
            energy += energy_lj_opt(epsilon, sigma_6, inv_distance_6,inv_cutoff_6, inv_distance_12, inv_cutoff_12);
            if (energy>MAX_ENERGY_opt) {skip = true; break;}
          }
        }
        // symmetry
        grid.data[idx] = 4*R*energy;
        visited[idx] = true;
        int sym_count = 1;
        for (size_t k = 0; k < grid_ops.size(); ++k) {
          array<int,3> t = grid_ops[k].apply(u, v, w);
          size_t mate_idx = grid.index_s(t[0], t[1], t[2]);
          if (visited[mate_idx])
            continue;
          sym_count++;
          grid.data[mate_idx] = grid.data[idx];
          visited[mate_idx] = true;
        }
        if (skip==false) {
          double exp_energy = sym_count * exp(-4*energy/(temperature));    
          sum_exp_energy += exp_energy;
          boltzmann_energy_lj += exp_energy * grid.data[idx];
        }
      }
}
