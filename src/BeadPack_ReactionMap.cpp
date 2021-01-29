//
//  BeadPack_ReactionMap.cpp
//
//  Created by Tomás Aquino on 17/11/2020.
//  Copyright © 2020 Tomás Aquino. All rights reserved.
//

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "BeadPack/BeadPack_Models.h"
#include "Field/VectorField_Interpolated.h"
#include "general/Operations.h"
#include "general/Ranges.h"
#include "general/useful.h"
#include "Stochastic/CTRW/CTRW.h"
#include "Stochastic/CTRW/JumpGenerator.h"
#include "Stochastic/CTRW/Measurer.h"
#include "Stochastic/CTRW/PTRW.h"
#include "Stochastic/CTRW/Reaction.h"
#include "Stochastic/CTRW/State.h"
#include "Stochastic/CTRW/StateGetter.h"
#include "Stochastic/CTRW/Transitions_State.h"

int main(int argc, const char * argv[])
{
  using namespace beadpack::model_bcc_cartesian;
  
  if (argc == 1)
  {
    std::cout << "Fluid-phase mass consumed as a function of angle on bead surface\n"
              << "for advective-diffusive particle tracking in 3d beadpacks\n"
              << "with periodic boundary conditions and\n"
              << "decay reaction at constant rate at the bead interfaces.\n"
              << "----------------------------------------------------\n"
              << "Parameters (default value in []):\n"
              << "damkohler : Damkohler number in terms of domain side, reaction rate at interface,\n"
              << "            and diffusion time\n"
              << "peclet : Peclet number in terms of domain side, average velocity,\n"
              << "         and diffusion coefficient\n"
              << "time_step_accuracy_adv : Maximum time step size in units of advection time\n"
              << "time_step_accuracy_diff : Minimum time step size in units of advection time\n"
              << "time_max_diffusion_times : Maximum output time in units of diffusion time\n"
              << "initial_condition_type : 0 - Uniformly random in the void space a plane\n"
              << "                       : 1 - Flux-weighted in the void space on a plane\n"
              << "                       : 2 - Uniformly random in the void space in the periodic domain\n"
              << "                       : 3 - Flux-weighted in the void space in the periodic domain\n"
              << "                       : 4 - Uniformly randomly over twice the discretization distance\n"
              << "                             from the interface in the periodic domain\n"
              << "                       : 5 - Uniformly randomly over twice the discretization distance\n"
              << "                             from the interface of all beads\n"
              << "                       : 6 - Uniformly randomly over the interface in the periodic domain\n"
              << "                       : 7 - Uniformly randomly over the interface of all beads\n"
              << "                       : 8 - Load positions from file\n"
              << "nr_bins_angle : Number of bins to equally discretize the range [0,2pi]"
              << "initial_condition_size_domains : Size of initial condition box or plane in domain sides\n"
              << "initial_mass : Total initial mass\n"
              << "nr_particles : Number of particles to track\n"
              << "run_nr : Nonnegative integer identifier for output files\n"
              << "data_set : Path to input data folder relative to input_dir_base\n"
              << "           (and model name identifier for output files)\n"
              << "filename_input_positions : Filename to read positions from\n"
              << "                           for initial_condition_type = 8 []\n"
              << "input_dir_base : Path to look for input data [../input]\n"
              << "output_dir : Path folder to output to [../output]\n";
    return 0;
  }

  if (argc < 13)
    throw useful::bad_parameters();

  using State = ctrw::State_periodic<std::vector<double>,
    std::vector<int>, double, useful::Empty, std::size_t>;
  using CTRW = ctrw::CTRW<State>;
  using Boundary = Boundaries::Boundary_Reflecting_Periodic;
  using JumpGenerator_Advection
    = ctrw::JumpGenerator_Velocity_withHint_RK4<VelocityField&, Boundary&>;
  using JumpGenerator_Diffusion = ctrw::JumpGenerator_Diffusion;
  using JumpGenerator = ctrw::JumpGenerator_Add<JumpGenerator_Advection, JumpGenerator_Diffusion>;

  std::size_t arg = 1;
  double damkohler = atof(argv[arg++]);
  double peclet = atof(argv[arg++]);
  double time_step_accuracy_adv = atof(argv[arg++]);
  double time_step_accuracy_diff = atof(argv[arg++]);
  double time_max_diffusion_times = atof(argv[arg++]);
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

  std::string input_dir = input_dir_base + "/" + data_set;

  std::cout << std::scientific << std::setprecision(2);

  std::cout << "Making bead pack...\n";
  const BeadPack bead_pack = make_bead_pack(input_dir);
  const Geometry geometry{ bead_pack.radius(0) };
  std::cout << "\tDone!\n";
  
  std::cout << "Setting up boundary conditions...\n";
  Boundaries boundaries{ geometry, bead_pack };
  std::cout << "\tDone!\n";
  
  std::cout << "Setting up velocity field...\n";
  auto velocity_field =
    make_velocity_field(input_dir, output_dir,
                        bead_pack, boundaries.boundary_periodic);
  std::cout << "\tDone!\n";

  std::cout << "Importing mean velocity...\n";
  std::string mean_velocity_filename = input_dir + "/" + "mean_velocity.dat";
  std::vector<double> mean_velocity
    = beadpack::get_mean_velocity(geometry.dim, mean_velocity_filename);
  double magnitude_mean_velocity = operation::abs(mean_velocity);
  std::cout << "\tDone!\n";

  std::cout << "Setting up particles...\n";
  double advection_time = geometry.domain_side/magnitude_mean_velocity;
  double diff = geometry.domain_side*magnitude_mean_velocity/peclet;
  double diffusion_time = geometry.domain_side*geometry.domain_side/(2.*diff);
  double time_step = std::min(time_step_accuracy_adv*advection_time,
    time_step_accuracy_diff*diffusion_time);
  double mass_per_particle = initial_mass/nr_particles;
  double length_discretization = 10.*std::sqrt(2.*diff*time_step);

  double reaction_rate = damkohler/diffusion_time*geometry.domain_side/length_discretization;

  auto state_maker = [&mass_per_particle, &geometry]()
  { return State{
    std::vector<double>(geometry.dim),
    std::vector<int>(geometry.dim), mass_per_particle }; };

  std::vector<double> domain_midpoint;
  for (auto const& boundary : geometry.boundaries)
    domain_midpoint.push_back((boundary.first+boundary.second)/2.);
  auto initial_box_centered = geometry.boundaries;
  for (std::size_t dd = 0; dd < geometry.dim; ++dd)
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
      bead_pack, boundaries.boundary_periodic,
      2.*length_discretization,
      filename_input_positions,
      state_maker),
    CTRW::Tag{} };
  std::cout << "\tDone!\n";

  std::cout << "Setting up output...\n";
  std::size_t check_progress_times = 50;
  double time_max = time_max_diffusion_times*diffusion_time;
  double time_min = time_max/check_progress_times;
  std::vector<double> measure_times = range::logspace(time_min, time_max, check_progress_times);

  std::stringstream stream;
  stream << std::scientific << std::setprecision(2);
  stream << initial_condition_size_domains << "_"
         << peclet << "_"
         << damkohler << "_"
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
  auto reaction = ctrw::Reaction_Decay_Condition_Map_Beadpack_3d{
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
          boundaries.boundary_reflecting_periodic
        },
        JumpGenerator_Diffusion{
          diff,
          time_step,
          geometry.dim
        }
      },
      boundaries.boundary_reflecting_periodic
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

  for (auto const& time : measure_times)
  {
    ptrw.evolve(time);
    std::cout << "time = " << time
              << "\ttime_last_measure = " << time_max
              << "\n";
  }
  reaction.print_map(filename_output_map);

  std::cout << "\tDone!\n";

  return 0;
}
