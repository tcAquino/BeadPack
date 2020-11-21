//
//  main.cpp
//  BeadPack_ReactionMap
//
//  Created by Tomás Aquino on 17/11/2020.
//  Copyright © 2020 Tomás Aquino. All rights reserved.
//

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include "BeadPack/BeadPack.h"
#include "BeadPack/BeadPack_InitialConditions.h"
#include "BeadPack/BeadPack_Input.h"
#include "Field/VectorField_Interpolated.h"
#include "general/Operations.h"
#include "general/Ranges.h"
#include "Geometry/Coordinates.h"
#include "Stochastic/CTRW/Boundary.h"
#include "Stochastic/CTRW/CTRW.h"
#include "Stochastic/CTRW/JumpGenerator.h"
#include "Stochastic/CTRW/Measurer.h"
#include "Stochastic/CTRW/PTRW.h"
#include "Stochastic/CTRW/State.h"
#include "Stochastic/CTRW/StateGetter.h"
#include "Stochastic/CTRW/Transitions_State.h"

template <typename BeadPack>
class Reaction_Decay_Condition_Map_Beadpack_3d
{
public:
  Reaction_Decay_Condition_Map_Beadpack_3d
  (double reaction_rate, double length_discretization,
   std::size_t nr_bins, BeadPack const& bead_pack)
  : reaction_rate{ reaction_rate }
  , length_discretization{ length_discretization }
  , nr_bins_phi{ nr_bins }
  , nr_bins_theta{ std::size_t(nr_bins/2) }
  , bead_pack{ bead_pack }
  {
    map_phi_theta.assign(nr_bins,
      std::vector<double>(std::size_t(nr_bins/2), 0.));
  }
  
  template <typename State>
  void operator()(State& state, double exposure_time)
  {
    auto near = bead_pack.near(state.position, length_discretization);
    if (near.first)
    {
      double old_mass = state.mass;
      state.mass *= std::exp(-reaction_rate*exposure_time);
      
      std::size_t bead = bead_pack.nearest_neighbor(state.position).first;
      auto spherical = geometry::cartesian2spherical(
        operation::minus(state.position, bead_pack.center(bead)));
      std::size_t bin_phi = (spherical[1]+constants::pi)/(2.*constants::pi)*nr_bins_phi;
      if (bin_phi == nr_bins_phi)
        --bin_phi;
      std::size_t bin_theta = spherical[2]/constants::pi*nr_bins_theta;
      if (bin_theta == nr_bins_theta)
        --bin_theta;
      map_phi_theta[bin_phi][bin_theta] += old_mass - state.mass;
    }
  }
  
  void print_map
  (std::string const& filename,
   int precision = 8, std::string delimiter = "\t")
  {
    std::ofstream output{ filename };
    if (!output.is_open())
      throw useful::open_write_error(filename);
    output << std::scientific << std::setprecision(precision);
    
    for (std::size_t ii = 0; ii < nr_bins_phi; ++ii)
    {
      double phi = 2.*constants::pi*(ii+0.5)/nr_bins_phi - constants::pi;
      for (std::size_t jj = 0; jj < nr_bins_theta; ++jj)
      {
        double theta = constants::pi*(jj+0.5)/nr_bins_theta;
        output << phi << delimiter
               << theta << delimiter
               << map_phi_theta[ii][jj] << "\n";
      }
    }
  }
  
private:
  const double reaction_rate;
  const double length_discretization;
  std::size_t nr_bins_phi;
  std::size_t nr_bins_theta;
  BeadPack const& bead_pack;
  std::vector<std::vector<double>> map_phi_theta;
};

int main(int argc, const char * argv[])
{
  const std::size_t dim = 3;

  if (argc != 15 && argc != 16 && argc != 17 && argc != 18)
  {
    throw useful::bad_parameters();
  }

  using BeadPack = beadpack::BeadPack<dim>;
  using Bead = BeadPack::Bead;
  using VelocityField = field::VectorField_LinearInterpolation_UnstructuredGrid<dim>;

  using Boundary_Periodic = boundary::Periodic_WithOutsideInfo;
  using Boundary = boundary::ReflectingBeads_Periodic<BeadPack, Boundary_Periodic>;

  using State = ctrw::State_periodic<std::vector<double>,
    std::vector<int>, double, useful::Empty, std::size_t>;
  using CTRW = ctrw::CTRW<State>;

  using JumpGenerator_Advection = ctrw::JumpGenerator_Velocity_withHint_RK4<VelocityField&, Boundary&>;
  using JumpGenerator_Diffusion = ctrw::JumpGenerator_Diffusion;
  using JumpGenerator = ctrw::JumpGenerator_Add<JumpGenerator_Advection, JumpGenerator_Diffusion>;

  std::size_t arg = 1;
  double domain_side = atof(argv[arg++]);
  double Peclet = atof(argv[arg++]);
  double Damkohler = atof(argv[arg++]);
  double time_step_accuracy_adv = atof(argv[arg++]);
  double time_step_accuracy_diff = atof(argv[arg++]);
  double time_max_diffusion_times = atof(argv[arg++]);
  int measure_type = atoi(argv[arg++]);
  int initial_condition_type = atoi(argv[arg++]);
  double initial_condition_size_domains = atof(argv[arg++]);
  double initial_mass = atof(argv[arg++]);
  std::size_t nr_bins_angle = strtoul(argv[arg++], NULL, 0);
  std::size_t nr_particles = strtoul(argv[arg++], NULL, 0);
  std::size_t run_nr = strtoul(argv[arg++], NULL, 0);
  std::string data_set = argv[arg++];
  std::string filename_input_positions = argc > arg ? argv[arg++] : "";
  std::string input_dir_base = argc > arg ? argv[arg++] : "../input";
  std::string output_dir = argc > arg ? argv[arg++] : "../output";

  std::size_t nr_measures = 50;
  int measure_spacing = 0;

  std::string input_dir = input_dir_base + "/" + data_set;

  std::cout << std::scientific << std::setprecision(2);

  std::vector<std::pair<double, double>> boundaries;
  boundaries.reserve(dim);
  for (std::size_t dd = 0; dd < dim; ++dd)
    boundaries.push_back({ 0., domain_side });
  std::vector<double> domain_dimensions(dim, domain_side);

  std::cout << "Importing beads...\n";
  std::string bead_filename = input_dir + "/" + "spheres.dat";
  BeadPack bead_pack{ beadpack::get_beads<Bead>(dim, bead_filename, 1, domain_side) };
  std::cout << "\tDone!\n";

  std::cout << "Importing contacts...\n";
  std::string contact_filename =  input_dir + "/" + "contacts.dat";
  auto contacts = beadpack::get_contacts<std::vector<double>>(dim, contact_filename, 1, domain_side);
  std::cout << "\tDone!\n";

  std::cout << "Importing grid points and velocities...\n";
  std::string velocity_filename = input_dir + "/" + "velocities.csv";
  auto points_velocities = beadpack::get_points_velocities_velocity_point<
    std::vector<double>,
    std::vector<double>>(dim, velocity_filename, 1);
  std::cout << "\tDone!\n";

  std::cout << "Adding zero velocity grid points at contacts...\n";
  for (auto const& point : contacts)
  {
    points_velocities.first.push_back(point);
    points_velocities.second.emplace_back(dim, 0.);
  }
  std::cout << "\tDone!\n";

  std::cout << "Setting up velocity field...\n";
  VelocityField velocity_field{ points_velocities.first, points_velocities.second };
  std::cout << "\tDone!\n";

  std::cout << "Cleaning up memory...\n";
  decltype(points_velocities.first)().swap(points_velocities.first);
  decltype(points_velocities.second)().swap(points_velocities.second);
  decltype(contacts)().swap(contacts);
  std::cout << "\tDone!\n";

  std::cout << "Importing mean velocity...\n";
  std::vector<double> mean_velocity;
  try
  {
    std::string mean_velocity_filename = input_dir + "/" + "mean_velocity.dat";
    mean_velocity = beadpack::get_mean_velocity<std::vector<double>>(dim, mean_velocity_filename);
  }
  catch (std::runtime_error& err)
  {
    std::cout << "\tFile not available. Computing...\n";
    std::string mean_velocity_filename = output_dir + "/" + "mean_velocity.dat";
    std::ofstream output{ mean_velocity_filename };
    if (!output.is_open())
      throw useful::open_write_error(mean_velocity_filename);
    output << std::setprecision(8)
           << std::scientific;
    std::size_t nr_samples = 1e4;
    mean_velocity = bead_pack.compute_mean_vector(velocity_field, boundaries, nr_samples);
    useful::print(output, mean_velocity);
    output << "\n";
    std::cout << "\t\tDone!\n";
  }
  double mean_velocity_magnitude = operation::abs(mean_velocity);
  std::cout << "\tDone!\n";

  std::cout << "Setting up particles...\n";
  double advection_time = domain_side/mean_velocity_magnitude;
  double diff = domain_side*mean_velocity_magnitude/Peclet;
  double diffusion_time = domain_side*domain_side/(2.*diff);
  double time_step = std::min(time_step_accuracy_adv*advection_time,
    time_step_accuracy_diff*diffusion_time);
  double mass_per_particle = initial_mass/nr_particles;
  double length_discretization = 10.*std::sqrt(2.*diff*time_step);

  double time_min_diffusion_times = time_max_diffusion_times/nr_measures;

  double reaction_rate = Damkohler/diffusion_time*domain_side/length_discretization;

  Boundary_Periodic boundary_periodic{ boundaries };
  Boundary boundary{ bead_pack, boundary_periodic };

  auto state_maker = [&mass_per_particle]()
  { return State{
    std::vector<double>(dim, 0.), std::vector<int>(dim, 0), mass_per_particle }; };

  std::vector<double> domain_midpoint;
  for (auto const& boundary : boundaries)
    domain_midpoint.push_back((boundary.first+boundary.second)/2.);
  std::vector<std::pair<double, double>> initial_box_centered = boundaries;
  for (std::size_t dd = 0; dd < dim; ++dd)
  {
    initial_box_centered[dd].first -= domain_midpoint[dd];
    initial_box_centered[dd].first *= initial_condition_size_domains;
    initial_box_centered[dd].second -= domain_midpoint[dd];
    initial_box_centered[dd].second *= initial_condition_size_domains;
  }

  std::string initial_condition_name =
    beadpack::initial_condition_name(initial_condition_type);

  CTRW ctrw{
    beadpack::make_particles<CTRW::Particle>(
      nr_particles, initial_condition_type,
      domain_midpoint, initial_box_centered,
      velocity_field, mean_velocity,
      bead_pack, boundary_periodic,
      length_discretization,
      filename_input_positions,
      state_maker),
    CTRW::Tag{} };
  std::cout << "\tDone!\n";

  std::cout << "Setting up output...\n";
  double time_min = time_min_diffusion_times*diffusion_time;
  double time_max = time_max_diffusion_times*diffusion_time;
  double dist_min = mean_velocity_magnitude*time_min;
  double dist_max = mean_velocity_magnitude*time_max;
  std::vector<double> measure_times, measure_distances;
  switch (measure_spacing)
  {
    case 0:
      measure_times = range::logspace(time_min, time_max, nr_measures);
      measure_distances = range::logspace(dist_min, dist_max, nr_measures);
      break;
    case 1:
      measure_times = range::linspace(time_min, time_max, nr_measures);
      measure_distances = range::linspace(dist_min, dist_max, nr_measures);
      break;
    case 2:
      break;
    default:
      throw std::runtime_error{ "Undefined measure spacing." };
  }

  std::stringstream stream;
  stream << std::scientific << std::setprecision(2);
  stream << domain_side << "_"
         << initial_condition_size_domains << "_"
         << Peclet << "_"
         << Damkohler << "_"
         << time_step_accuracy_diff << "_"
         << time_step_accuracy_adv << "_"
         << time_max_diffusion_times << "_"
         << nr_bins_angle << "_"
         << nr_particles << "_"
         << run_nr;
  std::string params = stream.str();

  std::string filename_output_map_base = "Data_beadpack_reactive_reactionmap";

  std::string filename_output_map = output_dir + "/" +
    filename_output_map_base + "_" + initial_condition_name +
    "_" + data_set + "_" + params + ".dat";

  std::cout << "\tDone!\n";

  std::cout << "Setting up dynamics...\n";
  auto reaction = Reaction_Decay_Condition_Map_Beadpack_3d{
    reaction_rate,
    length_discretization,
    nr_bins_angle,
    bead_pack
  };
  ctrw::Transitions_PTRW_Transport_Reaction transitions{
    ctrw::Transitions_Position{
      JumpGenerator{
        JumpGenerator_Advection{
          velocity_field,
          time_step,
          nr_particles,
          boundary
        },
        JumpGenerator_Diffusion{
          diff,
          time_step,
          dim
        }
      },
      boundary
    },
    reaction,
    time_step
  };
  ctrw::PTRW ptrw(ctrw, transitions, time_step, 0.);
  std::cout << "\tDone!\n";

  std::cout << "Starting dynamics...\n";
  std::cout << "\tAdvection time = " << advection_time << "\n";
  std::cout << "\tDiffusion time = " << diffusion_time << "\n";
  std::cout << "\tTime step [adv times] = " << time_step/advection_time << "\n";
  std::cout << "\tTime step [diff times] = " << time_step/diffusion_time << "\n";
  std::cout << "\tDiscretization length = " << length_discretization << "\n";
  std::cout << "\tNr of particles = " << ptrw.size() << "\n";

  switch (measure_type)
  {
    case 0:
    {
      for (auto const& time : measure_times)
      {
        ptrw.evolve(time);
        std::cout << "time = " << time
                  << "\ttime_last_measure = " << time_max
                  << "\n";
      }
      reaction.print_map(filename_output_map);
      break;
    }
    default:
      throw std::invalid_argument{ "Undefined measure type" };
  }

  std::cout << "\tDone!\n";

  return 0;
}
