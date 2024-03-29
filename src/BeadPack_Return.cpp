//
//  BeadPack_Return.cpp
//
//  Created by Tomás Aquino on 13/11/2020.
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
  using namespace beadpack::model_bcc_symmetryplanes;
  
  if (argc == 1)
  {
    std::cout << "Return times and distances to the fluid-solid interface\n"
              << "for advective-diffusive particle tracking in 3d beadpacks\n"
              << "with periodic boundary conditions.\n"
              << "----------------------------------------------------\n"
              << "Parameters (default value in []):\n"
              << "peclet : Peclet number in terms of domain side, average velocity,\n"
              << "         and diffusion coefficient\n"
              << "time_step_accuracy_adv : Maximum time step size in units of advection time\n"
              << "time_step_accuracy_diff : Minimum time step size in units of advection time\n"
              << "nr_measures : Number of measurements\n"
              << "run_nr : Nonnegative integer identifier for output files\n"
              << "data_set : Path to input data folder relative to input_dir_base\n"
              << "           (and model name identifier for output files)\n"
              << "filename_input_positions : Filename to read positions from\n"
              << "                           for initial_condition_type = 8 []\n"
              << "input_dir_base : Path for input [../input]\n"
              << "output_dir : Path for output [../output]\n";
    return 0;
  }
  
  if (argc < 7)
    throw useful::bad_parameters();
  
  using State = ctrw::State_periodic<std::vector<double>,
    std::vector<int>, useful::Empty, useful::Empty, std::size_t>;
  using CTRW = ctrw::CTRW<State>;
  using Boundary = Boundaries::Boundary_Periodic;
  using JumpGenerator_Advection
    = ctrw::JumpGenerator_Velocity_withHint_RK4<VelocityField&, Boundary&>;
  using JumpGenerator_Diffusion = ctrw::JumpGenerator_Diffusion;
  using JumpGenerator
    = ctrw::JumpGenerator_Add<JumpGenerator_Advection, JumpGenerator_Diffusion>;
  
  int arg = 1;
  double peclet = atof(argv[arg++]);
  double time_step_accuracy_adv = atof(argv[arg++]);
  double time_step_accuracy_diff = atof(argv[arg++]);
  std::size_t nr_measures = strtoul(argv[arg++], NULL, 0);
  std::size_t run_nr = strtoul(argv[arg++], NULL, 0);
  std::string data_set = argv[arg++];
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
  
  auto near_wall = [length_discretization, &bead_pack](State const& state)
  { return bead_pack.near(state.position, length_discretization).first; };
  
  auto inside_bead = [&bead_pack](State const& state)
  { return bead_pack.inside(state.position).first; };
  
  auto state_maker = [&geometry]()
  { return State{
    std::vector<double>(geometry.dim),
    std::vector<int>(geometry.dim) }; };
  
  CTRW ctrw{ beadpack::make_particles_random_near_wall_uniform_unit_cell<CTRW::Particle>
    (1, bead_pack, boundaries.boundary_periodic,
     length_discretization, state_maker),
    CTRW::Tag{} };
  std::cout << "\tDone!\n";
  
  std::cout << "Setting up output...\n";
  
  std::stringstream stream;
  stream << std::scientific << std::setprecision(2);
  stream << peclet << "_"
         << time_step_accuracy_diff << "_"
         << time_step_accuracy_adv << "_"
         << nr_measures << "_"
         << run_nr;
  std::string params = stream.str();
  
  std::string filename_output_base = "Data_beadpack_return_wall";
  
  std::string filename_output_time = output_dir + "/" +
    filename_output_base + "_time_" + data_set + "_" + params + ".dat";
  
  std::string filename_output_space = output_dir + "/" +
    filename_output_base + "_space_" + data_set + "_" + params + ".dat";
  
  auto getter_position_longitudinal = ctrw::Get_new_from_particle{
    ctrw::Get_projection{
      ctrw::Get_position_periodic{
        boundaries.boundary_periodic },
    mean_velocity } };
  std::cout << "\tDone!\n";
  
  std::cout << "Setting up dynamics...\n";
  ctrw::Transitions_Position transitions{
      JumpGenerator{
        JumpGenerator_Advection{
          velocity_field,
          time_step,
          1,
          boundaries.boundary_periodic
        },
        JumpGenerator_Diffusion{
          diff,
          time_step,
          geometry.dim
        }
      },
      boundaries.boundary_periodic
  };
  ctrw::PTRW ptrw{ ctrw, transitions, time_step, 0. };
  auto& part = ptrw.particles(0);
  std::cout << "\tDone!\n";
  
  std::cout << "Starting dynamics...\n";
  std::cout << "\tAdvection time = " << advection_time << "\n";
  std::cout << "\tDiffusion time = " << diffusion_time << "\n";
  std::cout << "\tTime step [adv times] = " << time_step/advection_time << "\n";
  std::cout << "\tTime step [diff times] = " << time_step/diffusion_time << "\n";
  std::cout << "\tDiscretization length = " << length_discretization << "\n";
  std::cout << "\tNr of measures = " << nr_measures << "\n";
  
  auto output_time = useful::open_write(filename_output_time);
  output_time << std::setprecision(8)
              << std::scientific;
  
  auto output_space = useful::open_write(filename_output_space);
  output_space << std::setprecision(8)
               << std::scientific;
  
  JumpGenerator_Diffusion jump_generator_diffusion{
    diff,
    time_step,
    geometry.dim
  };
  auto adjust = [&jump_generator_diffusion, &boundaries]
  (State& state)
  {
    auto state_old = state;
    operation::plus_InPlace(state.position, jump_generator_diffusion());
    boundaries.boundary_reflecting_periodic(state, state_old);
  };
  
  std::string delimiter = "";
  for (std::size_t pp = 0; pp < nr_measures; ++pp)
  {
    double start_time = ptrw.time();
    double start_position = getter_position_longitudinal(part);
    std::cout << "Measure " << pp+1 << " of " << nr_measures << "\n";
    while (!inside_bead(part.state_new()))
      ptrw.step();
    output_time << delimiter << ptrw.time()-start_time;
    output_space << delimiter
                 << getter_position_longitudinal(part) - start_position;
    delimiter = "\t";
    
    auto state = part.state_new();
    boundaries.boundary_reflecting_periodic(state, part.state_old());
    while (1)
    {
      // Place the particle at bead surface
      std::size_t bead
        = bead_pack.place_at_closest_surface(state, boundaries.boundary_periodic);
      
      // Check if corresponding position within length_discretization of wall is ok
      auto state_near_wall = state;
      double radius = bead_pack.radius(bead);
      auto radial_vector = operation::minus(state_near_wall.position, bead_pack.center(bead));
      double radius_val = radius + length_discretization;
      operation::times_scalar_InPlace(radius_val/radius, radial_vector);
      operation::plus(bead_pack.center(bead), radial_vector, state_near_wall.position);
      boundaries.boundary_periodic(state);
      
      // If not closer to another bead, accept
      if (!near_wall(state_near_wall))
      {
        ctrw.set(0, state_near_wall);
        break;
      }
      // Otherwise, adjust with diffusive step and try again
      adjust(state);
    }
  }
  output_time << "\n";
  output_space << "n";
  
  output_time.close();
  output_space.close();
  
  std::cout << "\tDone!\n";
  
  std::cout << "Nearest-neighbor interpolation used "
            << velocity_field.nr_interpolation_failures()
            << " times due to numerical issues computing linear interpolation\n";
  
  return 0;
}

