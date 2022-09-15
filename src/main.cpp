#include <local/forcefield.hpp> // define forcefield parameters
#include <local/supracell.hpp> // define supracell size

#include <chrono>

#include <gemmi/cif.hpp>        // file -> cif::Document
#include <gemmi/smcif.hpp>      // cif::Document -> SmallStructure
#include <gemmi/symmetry.hpp>   // Space Group manipulation
#include <gemmi/unitcell.hpp>
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

  // Adsorption infos from forcefield dictionary
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

  int n_max=0;int l_max=0; int m_max=0;
  set_supracell(&n_max, &m_max, &l_max, structure, cutoff, sigma_guest);

}
