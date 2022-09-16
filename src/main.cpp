#include <local/forcefield.hpp> // define forcefield parameters
#include <local/supracell.hpp> // define supracell size

#include <chrono>

#include <gemmi/cif.hpp>        // file -> cif::Document
#include <gemmi/smcif.hpp>      // cif::Document -> SmallStructure
#include <gemmi/symmetry.hpp>   // Space Group manipulation
#include <gemmi/unitcell.hpp>
#include <gemmi/grid.hpp>
#include <unordered_map>


// Set key constants
#define R 8.31446261815324e-3 // kJ/mol/K
#define sqrt_2 1.414213562373095
#define min_factor 1.122462048309373  // 2^(1/6)
#define N_A 6.02214076e23    // part/mol

using namespace std;
namespace cif = gemmi::cif;

int main(int argc, char* argv[]) {
  // Set up Input Variables
  chrono::high_resolution_clock::time_point t_start = chrono::high_resolution_clock::now();
  auto structure_file = argv[1];
  string forcefield_path = argv[2];
  double temperature = stod(argv[3]);
  double cutoff = stod(argv[4]);
  double cutoff_sq = cutoff*cutoff;
  double cutoff_6 = (cutoff_sq)*(cutoff_sq)*(cutoff_sq);
  double inv_cutoff_6 = 1.0/cutoff_6;
  double inv_cutoff_12 = inv_cutoff_6*inv_cutoff_6;
  string element_guest_str = argv[5];

    // Inialize key variables
  string element_host_str;
  double dist = 0;
  double distance_sq = 0;
  double epsilon = 0;
  double sigma = 0;
  double sigma_sq = 0;
  double sigma_6 = 0;
  double exp_energy = 0;

  // Read Forcefield Infos
  LennardJones::Parameters ff_params;
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
  cout << sg->hm << endl;
  structure.cell.set_cell_images_from_spacegroup(sg);
  vector<gemmi::SmallStructure::Site> unique_sites = structure.sites;
  vector<gemmi::SmallStructure::Site> all_sites = structure.get_all_unit_cell_sites();

  // supercell size
  int n_max=0;int l_max=0; int m_max=0;
  set_supracell(&n_max, &m_max, &l_max, structure, cutoff);

  // Creates a list of sites within the cutoff
  vector<array<double,6>> supracell_sites;
  gemmi::Fractional coord_temp;
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
    gemmi::Fractional coord = site.fract;
    for (int n = -n_max; (n<n_max+1); ++n){
      for (int m = -m_max; (m<m_max+1); ++m) {
        for (int l = -l_max; (l<l_max+1); ++l) {
          // calculate a distance from centre box
          array<double,6> pos_epsilon_sigma;
          coord_temp.x = coord.x + n;
          coord_temp.y = coord.y + m;
          coord_temp.z = coord.z + l;
          gemmi::Position pos = gemmi::Position(structure.cell.orthogonalize(coord_temp));
          pos_epsilon_sigma[0] = pos.x;
          pos_epsilon_sigma[1] = pos.y;
          pos_epsilon_sigma[2] = pos.z;
          pos_epsilon_sigma[3] = epsilon;
          pos_epsilon_sigma[4] = sigma * sigma;
          pos_epsilon_sigma[5] = pos_epsilon_sigma[4] * pos_epsilon_sigma[4] * pos_epsilon_sigma[4];
          supracell_sites.push_back(pos_epsilon_sigma);
        }
      }
    }
  }

  // Grid set-up
  gemmi::Grid<double> grid;
  bool denser = true; double approx_spacing = 0.12;

  grid.spacegroup = sg;
  grid.set_unit_cell(structure.cell);
  grid.set_size_from_spacing(approx_spacing, denser);

  // Symmetry-aware grid construction
  vector<gemmi::GridOp> ops = grid.get_scaled_ops_except_id();
  vector<size_t> mates(ops.size(), 0);
  vector<bool> visited(grid.data.size(), false);
  size_t idx = 0;
  for (int w = 0; w != grid.nw; ++w)
    for (int v = 0; v != grid.nv; ++v)
      for (int u = 0; u != grid.nu; ++u, ++idx) {
        assert(idx == grid.index_q(u, v, w));
        if (visited[idx])
          continue;
        
// HERE

        for (size_t k = 0; k < ops.size(); ++k) {
          std::array<int,3> t = ops[k].apply(u, v, w);
          mates[k] = grid.index_n(t[0], t[1], t[2]);
        }
        double value = grid.data[idx];
        for (size_t k : mates) {
          if (visited[k])
            gemmi::fail("grid size is not compatible with space group");
          // value = func(value, grid.data[k]);
        }
        grid.data[idx] = value;
        visited[idx] = true;
        for (size_t k : mates) {
          grid.data[k] = value;
          visited[k] = true;
        }
      }
  assert(idx == grid.data.size());

  // Save grid in ccp4 binary format
  // visualisation using python

}