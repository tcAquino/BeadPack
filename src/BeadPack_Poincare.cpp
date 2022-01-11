//
//  BeadPack_Poincare.cpp
//  BeadPack_Poincare
//
//  Created by Tomás Aquino on 06/01/2022.
//  Copyright © 2022 Tomás Aquino. All rights reserved.
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
    std::cout << "Plane crossings in\n"
              << "advective-diffusive particle tracking in 3d beadpacks\n"
              << "with periodic boundary conditions.\n"
              << "----------------------------------------------------\n"
              << "Parameters (default value in []):\n"
              << "peclet : Peclet number in terms of domain side, average velocity,\n"
              << "         and diffusion coefficient\n"
              << "time_step_accuracy_adv : Maximum time step size in units of advection time\n"
              << "time_step_accuracy_diff : Minimum time step size in units of advection time\n"
              << "nr_measures : Number of measurement planes\n"
              << "initial_condition_type : 0 - Uniformly random in the void space on a plane\n"
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
              << "initial_condition_size_domains : Size of initial condition box or plane in domain sides\n"
              << "                                 (ignored if not applicable)\n"
              << "nr_particles : Number of particles to track\n"
              << "run_nr : Nonnegative integer identifier for output files\n"
              << "data_set : Path to input data folder relative to input_dir_base\n"
              << "           (and model name identifier for output files)\n"
              << "filename_input_positions : Filename to read positions from\n"
              << "                           for initial_condition_type = 8 []\n"
              << "input_dir_base : Path for input [../input]\n"
              << "output_dir : Path for output [../output]\n";
    return 0;
  }
  
  if (argc < 10)
    throw useful::bad_parameters();
  
  using State = ctrw::State_periodic<std::vector<double>,
    std::vector<int>, useful::Empty, useful::Empty, std::size_t>;
  using CTRW = ctrw::CTRW<State>;
  using Boundary = Boundaries::Boundary_Reflecting_Periodic;
  using JumpGenerator_Advection = ctrw::JumpGenerator_Velocity_withHint_RK4<VelocityField&, Boundary&>;
  using JumpGenerator_Diffusion = ctrw::JumpGenerator_Diffusion;
  using JumpGenerator = ctrw::JumpGenerator_Add<JumpGenerator_Advection, JumpGenerator_Diffusion>;
  
  int arg = 1;
  double peclet = atof(argv[arg++]);
  double time_step_accuracy_adv = atof(argv[arg++]);
  double time_step_accuracy_diff = atof(argv[arg++]);
  std::size_t nr_measures = strtoul(argv[arg++], NULL, 0);
  int initial_condition_type = atoi(argv[arg++]);
  double initial_condition_size_domains = atof(argv[arg++]);
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
    make_velocity_field(input_dir, bead_pack, boundaries.boundary_periodic);
  std::cout << "\tDone!\n";
  
  std::cout << "Importing mean velocity...\n";
  std::vector<double> mean_velocity
    = beadpack::get_mean_velocity(input_dir + "/" + "mean_velocity.dat");
  double magnitude_mean_velocity = operation::abs(mean_velocity);
  std::cout << "\tDone!\n";
  
  std::cout << "Setting up particles...\n";
  double advection_time = geometry.domain_side/magnitude_mean_velocity;
  double diff = geometry.domain_side*magnitude_mean_velocity/peclet;
  double diffusion_time = geometry.domain_side*geometry.domain_side/(2.*diff);
  double time_step = std::min(time_step_accuracy_adv*advection_time,
    time_step_accuracy_diff*diffusion_time);
  double length_discretization = 10.*std::sqrt(2.*diff*time_step);
  
  auto state_maker = [&geometry]()
  { return State{ std::vector<double>(geometry.dim),
    std::vector<int>(geometry.dim) }; };
  
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
  
  boundary::Periodic boundary_cubic_cell{ geometry.boundaries };
  CTRW ctrw{
    beadpack::make_particles<CTRW::Particle>(
      nr_particles, initial_condition_type,
      domain_midpoint, initial_box_centered,
      velocity_field, mean_velocity,
      bead_pack,
      [&boundary_cubic_cell, &boundaries]
      (State& state, State const& state_old = {})
      {
        bool b1 = boundary_cubic_cell(state, state_old);
        bool b2 = boundaries.boundary_periodic(state, state_old);
        return b1 || b2;
      },
      length_discretization,
      filename_input_positions,
      state_maker),
    CTRW::Tag{} };
  std::cout << "\tDone!\n";
  
  std::cout << "Setting up output...\n";
  std::stringstream stream;
  stream << std::scientific << std::setprecision(2);
  stream << initial_condition_size_domains << "_"
         << peclet << "_"
         << time_step_accuracy_diff << "_"
         << time_step_accuracy_adv << "_"
         << nr_measures << "_"
         << nr_particles << "_"
         << run_nr;
  std::string params = stream.str();
  
  std::string filename_output_base = "Data_beadpack";
  
  std::string filename_output_positions = output_dir + "/" +
    filename_output_base + "_section_positions_" + initial_condition_name +
    "_" + data_set + "_" + params + ".dat";
  
  // Choose face to monitor crossings
  // closest to perpendicular to mean velocity vector
  std::vector<double> velocity_dot_with_direction(geometry.dim, 0.);
  auto mean_velocity_direction =
    operation::div_scalar(mean_velocity, magnitude_mean_velocity);
  std::size_t crossing_face =
    std::distance(mean_velocity_direction.cbegin(),
                  std::max_element(mean_velocity_direction.cbegin(),
                                   mean_velocity_direction.cend(),
                                   [](double const& first, double const& second)
                                   { return std::abs(first) < std::abs(second); }));
  std::vector<double> crossing_direction(geometry.dim, 0.);
  crossing_direction[crossing_face] = mean_velocity_direction[crossing_face]/
    std::abs(mean_velocity_direction[crossing_face]);
  double crossing_face_position = crossing_direction[crossing_face] > 0
  ? std::abs(geometry.boundaries[crossing_face].second)
  : std::abs(geometry.boundaries[crossing_face].first);
  std::vector<double> crossing_values =
    range::linspace(crossing_face_position,
                    crossing_face_position+(nr_measures-1)*geometry.domain_side,
                    nr_measures);
  std::cout << "\tDone!\n";
  
  std::cout << "Outputting initial info...\n";
  std::string filename_output_initial_positions = output_dir + "/" +
    filename_output_base + "_initial_position_" + initial_condition_name +
    "_" + data_set + "_" + params + ".dat";
  auto output_initial_positions
    = useful::open_write(filename_output_initial_positions);
  output_initial_positions << std::setprecision(8)
                           << std::scientific;
  auto getter_position_real =
    ctrw::Get_position_periodic{ boundaries.boundary_periodic };
  for (auto const& part : ctrw.particles())
  {
    output_initial_positions << part.state_new().tag;
    useful::print(output_initial_positions, getter_position_real(part.state_new()), 1);
    output_initial_positions << "\n";
  }
  output_initial_positions.close();
  
  std::string filename_output_crossing_levels = output_dir + "/" +
    filename_output_base + "_crossing_levels_" + initial_condition_name +
    "_" + data_set + "_" + params + ".dat";
  auto output_crossing_levels
    = useful::open_write(filename_output_crossing_levels);
  output_crossing_levels << std::setprecision(8)
                         << std::scientific;
  output_crossing_levels << crossing_face << "\t"
                         << crossing_direction[crossing_face] << "\t";
  useful::print(output_crossing_levels, crossing_values);
  output_crossing_levels << "\n";
  output_crossing_levels.close();
  std::cout << "\tDone!\n";
  
  std::cout << "Setting up dynamics...\n";
  ctrw::Transitions_Position transitions{
      JumpGenerator{
        JumpGenerator_Advection{
          velocity_field,
          time_step,
          nr_particles,
          boundaries.boundary_reflecting_periodic,
        },
        JumpGenerator_Diffusion{
          diff,
          time_step,
          geometry.dim
        }
      },
      boundaries.boundary_reflecting_periodic
  };
  ctrw::PTRW ptrw(ctrw, transitions, time_step, 0.);
  std::cout << "\tDone!\n";
  
  ctrw::Measurer_Store_FirstCrossing_Tagged<std::vector<double>>
    measurer{ crossing_values, nr_particles };
  auto getter_position_plane =
  [&getter_position_real, &crossing_direction]
  (State const& state)
  {
    return operation::dot(getter_position_real(state),
                          crossing_direction);
  };
  auto getter_position_cubic_cell_plane =
  [&boundary_cubic_cell,&boundaries,&geometry,crossing_face]
  (CTRW::Particle const& particle)
  {
    auto state_copy = particle.state_new();
    boundaries.boundary_periodic.translate(state_copy.position,state_copy.periodicity);
    boundary_cubic_cell(state_copy);
    std::vector<double> position_plane;
    position_plane.reserve(geometry.dim-1);
    for (std::size_t dd = 0; dd < geometry.dim; ++dd)
      if (dd != crossing_face)
        position_plane.push_back(state_copy.position[dd]);
      
    return position_plane;
  };
  
  auto not_finished = [&crossing_values, &measurer]
  (CTRW::Particle part)
  {
    return !(measurer.crossed(crossing_values.size()-1,
                              part.state_new().tag));
  };
  std::cout << "\tDone!\n";
  
  std::cout << "Starting dynamics...\n";
  std::cout << "\tAdvection time = " << advection_time << "\n";
  std::cout << "\tDiffusion time = " << diffusion_time << "\n";
  std::cout << "\tTime step [adv times] = " << time_step/advection_time << "\n";
  std::cout << "\tTime step [diff times] = " << time_step/diffusion_time << "\n";
  std::cout << "\tDiscretization length = " << length_discretization << "\n";
  std::cout << "\tNr of measures = " << nr_measures << "\n";
  
  std::size_t nr_finished = 0;
  std::size_t counter = 0;
  std::size_t notify_every = 1000;
  while(measurer.count(crossing_values.size()-1) < nr_particles)
  {
    ptrw.step(not_finished);
    measurer.update(ptrw,
                    getter_position_cubic_cell_plane,
                    getter_position_plane);
    if (measurer.count(crossing_values.size()-1) > nr_finished)
    {
      nr_finished = measurer.count(crossing_values.size()-1);
      std::cout << nr_finished << " of " << nr_particles << " particles have crossed last level\n";
    }
    
    if (counter % notify_every == 0)
      std::cout
        << "\n"
        << "Last particle at "
        << getter_position_plane(
             std::min_element(ptrw.particles().cbegin(), ptrw.particles().cend(),
                              [&getter_position_plane]
                              (CTRW::Particle const& p1, CTRW::Particle const& p2)
                              { return getter_position_plane(p1.state_new()) < getter_position_plane(p2.state_new()); })->state_new())
        << "\n"
        << "Last level = " << crossing_values.back() << "\n"
        << "Current time [adv times] = " << ptrw.time()/advection_time << "\n"
        << "Current time [diff times] = " << ptrw.time()/diffusion_time << "\n"
        << "\n";
    ++counter;
  }
  
  std::cout << "\tDone!\n";
  
  std::cout << "Outputting...\n";
  measurer.print(filename_output_positions);
  
  std::cout << "\tDone!\n";
  
  return 0;
}
