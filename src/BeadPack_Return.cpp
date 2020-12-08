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
#include <utility>
#include <vector>
#include "BeadPack/BeadPack.h"
#include "BeadPack/BeadPack_InitialConditions.h"
#include "BeadPack/BeadPack_Input.h"
#include "Field/VectorField_Interpolated.h"
#include "general/Operations.h"
#include "general/Ranges.h"
#include "Stochastic/CTRW/Boundary.h"
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
  if (argc == 1)
  {
    std::cout << "Return times and distances to the interface for\n"
              << "advective-diffusive particle tracking in 3d beadpacks\n"
              << "with periodic boundary conditions on a cubic domain.\n"
              << "----------------------------------------------------\n"
              << "Parameters (default value in []):\n"
              << "domain_side : Length of domain side or periodic unit cell\n"
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
              << "input_dir_base : Path to look for input data [../input]\n"
              << "output_dir : Path folder to output to [../output]\n";
    return 0;
  }
  
  if (argc != 8 && argc != 9 && argc != 10)
  {
    throw useful::bad_parameters();
  }
  
  const std::size_t dim = 3;
  
  using BeadPack = beadpack::BeadPack<dim>;
  using Bead = BeadPack::Bead;
  using VelocityField = field::VectorField_LinearInterpolation_UnstructuredGrid<dim>;
  
  using Boundary_Periodic = boundary::Periodic_WithOutsideInfo;
  using Boundary = boundary::ReflectingBeads_Periodic<BeadPack, Boundary_Periodic>;
  
  using State = ctrw::State_periodic<std::vector<double>,
    std::vector<int>, useful::Empty, useful::Empty, std::size_t>;
  using CTRW = ctrw::CTRW<State>;
  
  using JumpGenerator_Advection = ctrw::JumpGenerator_Velocity_withHint_RK4<VelocityField&, Boundary&>;
  using JumpGenerator_Diffusion = ctrw::JumpGenerator_Diffusion;
  using JumpGenerator = ctrw::JumpGenerator_Add<JumpGenerator_Advection, JumpGenerator_Diffusion>;
  
  std::size_t arg = 1;
  double domain_side = atof(argv[arg++]);
  double peclet = atof(argv[arg++]);
  double time_step_accuracy_adv = atof(argv[arg++]);
  double time_step_accuracy_diff = atof(argv[arg++]);
  std::size_t nr_measures = strtoul(argv[arg++], NULL, 0);
  std::size_t run_nr = strtoul(argv[arg++], NULL, 0);
  std::string data_set = argv[arg++];
  std::string input_dir_base = argc > arg ? argv[arg++] : "../input/";
  std::string output_dir = argc > arg ? argv[arg++] : "../output/";
  
  std::string input_dir = input_dir_base + data_set + "/";
  std::cout << std::scientific << std::setprecision(2);
    
  std::vector<double> domain_dimensions(dim, domain_side);
  std::vector<std::pair<double, double>> boundaries;
  boundaries.reserve(dim);
  for (std::size_t dd = 0; dd < dim; ++dd)
    boundaries.push_back({ 0., domain_side });
  
  std::cout << "Importing beads...\n";
  std::string bead_filename = input_dir + "spheres.dat";
  
  BeadPack bead_pack{ beadpack::get_beads<Bead>(dim, bead_filename, 1, domain_side) };
  std::cout << "\tDone!\n";
  
  std::cout << "Importing contacts...\n";
  std::string contact_filename = input_dir + "contacts.dat";
  auto contacts = beadpack::get_contacts<std::vector<double>>(dim, contact_filename, 1, domain_side);
  std::cout << "\tDone!\n";
  
  std::cout << "Importing grid points and velocities...\n";
  std::string velocity_filename = input_dir + "velocities.csv";
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
    std::string mean_velocity_filename = input_dir + "mean_velocity.dat";
    mean_velocity = beadpack::get_mean_velocity<std::vector<double>>(dim, mean_velocity_filename);
  }
  catch (std::runtime_error& err)
  {
    std::cout << "\tFile not available. Computing...\n";
    std::size_t nr_samples = 1e4;
    mean_velocity = bead_pack.compute_mean_vector(velocity_field, boundaries, nr_samples);
    std::string mean_velocity_filename = output_dir + "mean_velocity.dat";
    std::ofstream output{ mean_velocity_filename };
    if (!output.is_open())
      throw useful::open_write_error(mean_velocity_filename);
    output << std::setprecision(8)
           << std::scientific;
    useful::print(output, mean_velocity);
    output << "\n";
    std::cout << "\t\tDone!\n";
  }
  double mean_velocity_magnitude = operation::abs(mean_velocity);
  std::cout << "\tDone!\n";
  
  std::cout << "Setting up particles...\n";
  double advection_time = domain_side/mean_velocity_magnitude;
  double diff = domain_side*mean_velocity_magnitude/peclet;
  double diffusion_time = domain_side*domain_side/(2.*diff);
  double time_step = std::min(time_step_accuracy_adv*advection_time,
    time_step_accuracy_diff*diffusion_time);
  double length_discretization = 10.*std::sqrt(2.*diff*time_step);
  
  auto near_wall = [length_discretization, &bead_pack](State const& state)
  { return bead_pack.near(state.position, length_discretization).first; };
  
  Boundary_Periodic boundary_periodic{ boundaries };
  Boundary boundary{ bead_pack, boundary_periodic };
  
  auto state_maker = []()
  { return State{ std::vector<double>(dim, 0.), std::vector<int>(dim, 0) }; };
  
  CTRW ctrw{ beadpack::make_particles_random_near_wall_uniform_unit_cell<CTRW::Particle>
    (1, bead_pack, boundary_periodic, 2.*length_discretization, state_maker), CTRW::Tag{} };
  std::cout << "\tDone!\n";
  
  std::cout << "Setting up output...\n";
  
  std::stringstream stream;
  stream << std::scientific << std::setprecision(2);
  stream << domain_side << "_"
         << peclet << "_"
         << time_step_accuracy_diff << "_"
         << time_step_accuracy_adv << "_"
         << nr_measures << "_"
         << run_nr;
  std::string params = stream.str();
  
  std::string filename_output_base = "Data_beadpack_return_wall";
  
  std::string filename_output_time = output_dir +
    filename_output_base + "_time_" + data_set + "_" + params + ".dat";
  
  std::string filename_output_space = output_dir +
    filename_output_base + "_space_" + data_set + "_" + params + ".dat";
  
  auto getter_position_longitudinal = ctrw::Get_position_periodic_projection{
    domain_dimensions, mean_velocity };
  std::cout << "\tDone!\n";
  
  std::cout << "Setting up dynamics...\n";
  ctrw::Transitions_Position transitions{
      JumpGenerator{
        JumpGenerator_Advection{
          velocity_field,
          time_step,
          1,
          boundary
        },
        JumpGenerator_Diffusion{
          diff,
          time_step,
          dim
        }
      },
      boundary
  };
  ctrw::PTRW ptrw{ ctrw, transitions, time_step, 0. };
  std::cout << "\tDone!\n";
  
  std::cout << "Starting dynamics...\n";
  std::cout << "\tAdvection time = " << advection_time << "\n";
  std::cout << "\tDiffusion time = " << diffusion_time << "\n";
  std::cout << "\tTime step [adv times] = " << time_step/advection_time << "\n";
  std::cout << "\tTime step [diff times] = " << time_step/diffusion_time << "\n";
  std::cout << "\tDiscretization length = " << length_discretization << "\n";
  std::cout << "\tNr of measures = " << nr_measures << "\n";
  
  std::string filename_time{ filename_output_time };
  std::ofstream output_time{ filename_time };
  if (!output_time.is_open())
    throw useful::open_write_error(filename_time);
  output_time << std::setprecision(8)
              << std::scientific;
  
  std::string filename_space{ filename_output_space };
  std::ofstream output_space{ filename_space };
  if (!output_space.is_open())
    throw useful::open_write_error(filename_space);
  output_space << std::setprecision(8)
               << std::scientific;
  
  std::string delimiter = "";
  for (std::size_t pp = 0; pp < nr_measures; ++pp)
  {
    double start_time = ptrw.time();
    double start_position = getter_position_longitudinal(ptrw.particles(0).state_new());
    std::cout << "Measure " << pp+1 << " of " << nr_measures << "\n";
    while (!near_wall(ptrw.particles(0).state_new()))
      ptrw.step();
    output_time << delimiter << ptrw.time()-start_time;
    output_space << delimiter
                 << getter_position_longitudinal(ptrw.particles(0).state_new()) - start_position;
    delimiter = "\t";
    
    while (1)
    {
      // Place the particle at distance 2*length_discretization from nearest interface
      auto state = ptrw.particles(0).state_new();
      std::size_t bead = bead_pack.place_at_closest_surface(state, boundary_periodic);
      double radius = bead_pack.radius(bead);
      auto radial_vector = operation::minus(state.position, bead_pack.center(bead));
      double radius_val = radius + 2.*length_discretization;
      operation::times_scalar_InPlace(radius_val/radius, radial_vector);
      operation::plus(bead_pack.center(bead), radial_vector, state.position);
      boundary_periodic(state);
      
      // If not closer to another bead, accept
      if (!bead_pack.near(state.position, 2.*length_discretization).first)
      {
        ctrw.set(0, state);
        break;
      }
      
      // Othwerise, adjust with a PTRW step and try again
      ptrw.step();
    }
  }
  output_time << "\n";
  output_space << "n";
  
  output_time.close();
  output_space.close();
  
  std::cout << "\tDone!\n";
  
  return 0;
}

