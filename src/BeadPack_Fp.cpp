//
//  Beadpack_Fp.cpp
//
//  Created by Tomás Aquino on 10/11/2020.
//  Copyright © 2020 Tomás Aquino. All rights reserved.
//

#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
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
    std::cout << "First passage times and distances to the fluid-solid interface\n"
              << "for advective-diffusive particle tracking in 3d beadpacks\n"
              << "with periodic boundary conditions on a cubic domain.\n"
              << "----------------------------------------------------\n"
              << "Parameters (default value in []):\n"
              << "peclet : Peclet number in terms of domain side, average velocity,\n"
              << "         and diffusion coefficient\n"
              << "time_step_accuracy_adv : Maximum time step size in units of advection time\n"
              << "time_step_accuracy_diff : Minimum time step size in units of advection time\n"
              << "measure_type : 0 - First passage times and distances (output to file)\n"
              << "               1 - Mean first passage time and distance (output to console)\n"
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
              << "nr_measures : Number of measurements\n"
              << "run_nr : Nonnegative integer identifier for output files\n"
              << "data_set : Path to input data folder relative to input_dir_base\n"
              << "           (and model name identifier for output files)\n"
              << "filename_input_positions : Filename to read positions from for initial_condition_type = 8 []\n"
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
  using JumpGenerator_Advection
    = ctrw::JumpGenerator_Velocity_withHint_RK4<VelocityField&, Boundary&>;
  using JumpGenerator_Diffusion = ctrw::JumpGenerator_Diffusion;
  using JumpGenerator = ctrw::JumpGenerator_Add<JumpGenerator_Advection, JumpGenerator_Diffusion>;
  
  std::size_t arg = 1;
  double peclet = atof(argv[arg++]);
  double time_step_accuracy_adv = atof(argv[arg++]);
  double time_step_accuracy_diff = atof(argv[arg++]);
  int measure_type = atoi(argv[arg++]);
  int initial_condition_type = atoi(argv[arg++]);
  double initial_condition_size_domains = atof(argv[arg++]);
  std::size_t nr_measures = strtoul(argv[arg++], NULL, 0);
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
    = beadpack::get_mean_velocity(geometry.dim, input_dir + "/" + "mean_velocity.dat");
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
  
  auto particles = beadpack::make_particles<CTRW::Particle>(
      nr_measures, initial_condition_type,
      domain_midpoint, initial_box_centered,
      velocity_field, mean_velocity,
      bead_pack, boundaries.boundary_periodic,
      2.*length_discretization,
      filename_input_positions,
      state_maker);
  std::cout << "\tDone!\n";
  
  std::cout << "Setting up output...\n";
  
  std::stringstream stream;
  stream << std::scientific << std::setprecision(2);
  stream << initial_condition_size_domains << "_"
         << peclet << "_"
         << time_step_accuracy_diff << "_"
         << time_step_accuracy_adv << "_"
         << nr_measures << "_"
         << run_nr;
  std::string params = stream.str();
  
  std::string filename_output_base = "Data_beadpack_fp";
  
  std::string filename_output_time = output_dir + "/" +
    filename_output_base + "_time_" + initial_condition_name +
    "_" + data_set + "_" + params + ".dat";
  
  std::string filename_output_space = output_dir + "/" +
    filename_output_base + "_space_" + initial_condition_name +
    "_" + data_set + "_" + params + ".dat";
  
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
          boundaries.boundary_reflecting_periodic
        },
        JumpGenerator_Diffusion{
          diff,
          time_step,
          geometry.dim
        }
      },
      boundaries.boundary_reflecting_periodic
  };
  std::cout << "\tDone!\n";
  
  std::cout << "Starting dynamics...\n";
  std::cout << "\tAdvection time = " << advection_time << "\n";
  std::cout << "\tDiffusion time = " << diffusion_time << "\n";
  std::cout << "\tTime step [adv times] = " << time_step/advection_time << "\n";
  std::cout << "\tTime step [diff times] = " << time_step/diffusion_time << "\n";
  std::cout << "\tDiscretization length = " << length_discretization << "\n";
  std::cout << "\tNr of measures = " << nr_measures << "\n";
  
  switch (measure_type)
  {
    case 0:
    {
      auto output_time = useful::open_write(filename_output_time);
      output_time << std::setprecision(8)
                  << std::scientific;
      
      auto output_space = useful::open_write(filename_output_space);
      output_space << std::setprecision(8)
                   << std::scientific;
      
      std::string delimiter = "";
      for (std::size_t pp = 0; pp < particles.size(); ++pp)
      {
        std::cout << "Measure " << pp+1 << " of " << nr_measures << "\n";
        CTRW ctrw{ { particles[pp] }, CTRW::Tag{} };
        ctrw::PTRW ptrw{ ctrw, transitions, time_step, 0. };
        auto& part = ptrw.particles(0);
        double initial_position = getter_position_longitudinal(part);
        while (!near_wall(part.state_new()))
          ptrw.step();
        output_time << delimiter << ptrw.time();
        output_space << delimiter
                     << getter_position_longitudinal(part)-initial_position;
        delimiter = "\t";
      }
      output_time << "\n";
      output_space << "n";
      
      output_time.close();
      output_space.close();

      break;
    }
    case 1:
    {
      double mean_fpt = 0.;
      double mean_fpd = 0.;
      
      for (std::size_t pp = 0; pp < particles.size(); ++pp)
      {
        CTRW ctrw{ { particles[pp] }, CTRW::Tag{} };
        ctrw::PTRW ptrw{ ctrw, transitions, time_step, 0. };
        auto& part = ptrw.particles(0);
        double initial_position = getter_position_longitudinal(part);
        while (!near_wall(part.state_new()))
          ptrw.step();
        mean_fpt += ptrw.time();
        mean_fpd += getter_position_longitudinal(part)-initial_position;
        
        std::cout << std::setprecision(8) << std::scientific ;
        std::cout << "Mean fpt = " << mean_fpt/nr_measures << "\n";
        std::cout << "Mean fpd = " << mean_fpd/nr_measures << "\n";
        
        if (initial_condition_type == 4)
        {
          std::cout << "Mean fpt * geometry.domain_side/length_discretization = "
                    << mean_fpt/nr_measures
                       *geometry.domain_side/length_discretization << "\n";
          std::cout << "Mean fpd * geometry.domain_side/length_discretization = "
                    << mean_fpd/nr_measures
                       *geometry.domain_side/length_discretization << "\n";
        }
      }
      
      break;
    }
    default:
      throw std::invalid_argument{ "Undefined measure type" };
  }
  
  std::cout << "\tDone!\n";
  
  return 0;
}
