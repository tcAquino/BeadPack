//
//  BeadPack_Conservative.cpp
//
//  Created by Tomás Aquino on 19/11/2020.
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
    std::cout << "Advective-diffusive particle tracking in 3d beadpacks\n"
              << "with periodic boundary conditions.\n"
              << "----------------------------------------------------\n"
              << "Parameters (default value in []):\n"
              << "peclet : Peclet number in terms of domain side, average velocity,\n"
              << "         and diffusion coefficient\n"
              << "time_step_accuracy_adv : Maximum time step size in units of advection time\n"
              << "time_step_accuracy_diff : Minimum time step size in units of advection time\n"
              << "time_min_nondim : Minimum output time (in units specified by measure type)\n"
              << "time_max_nondim : Minimum output time (in units specified by measure type)\n"
              << "nr_measures : Number of measurements\n"
              << "measure_spacing : 0 - Logarithmic spacing between measurements\n"
              << "                : 1 - Linear spacing between measurements\n"
              << "measure_units : 0 - Measure times in units of diffusion time\n"
              << "                1 - Measure times in units of advection time\n"
              << "measure_type :  0 - Output time and all particle positions per file line,\n"
              << "                   at each measurement time\n"
              << "               1 - Output particle positions, one per file line, at final time only\n"
              << "               2 - Output position fluctuations autocorrelation in time\n"
              << "initial_condition_type : 0 - Uniformly random in the void space on a plane\n"
              << "                       : 1 - Flux-weighted in the void space on a plane\n"
              << "                       : 2 - Uniformly random in the void space in bounding box\n"
              << "                       : 3 - Flux-weighted in the void space in bounding box\n"
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
  
  if (argc < 15)
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
  double time_min_nondim = atof(argv[arg++]);
  double time_max_nondim = atof(argv[arg++]);
  std::size_t nr_measures = strtoul(argv[arg++], NULL, 0);
  int measure_spacing = atoi(argv[arg++]);
  int measure_units = atoi(argv[arg++]);
  int measure_type = atoi(argv[arg++]);
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
  
  CTRW ctrw{
    beadpack::make_particles<CTRW::Particle>(
      nr_particles, initial_condition_type,
      domain_midpoint, initial_box_centered,
      velocity_field, mean_velocity,
      bead_pack, boundaries.boundary_periodic,
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
         << time_min_nondim << "_"
         << time_max_nondim << "_"
         << nr_measures << "_"
         << measure_spacing << "_"
         << measure_units << "_"
         << nr_particles << "_"
         << run_nr;
  std::string params = stream.str();
  
  std::string filename_output_base = "Data_beadpack";
  
  std::string filename_output_positions = output_dir + "/" +
    filename_output_base + "_positions_" + initial_condition_name +
    "_" + data_set + "_" + params + ".dat";
  
  double measure_converter = measure_units
  ? advection_time
  : diffusion_time;
  double time_min = time_min_nondim*measure_converter;
  double time_max = time_max_nondim*measure_converter;
  std::vector<double> measure_times;
  switch (measure_spacing)
  {
    case 0:
      measure_times = range::logspace(time_min, time_max, nr_measures);
      break;
    case 1:
      measure_times = range::linspace(time_min, time_max, nr_measures);
      break;
    case 2:
      break;
    default:
      throw std::runtime_error{ "Undefined measure spacing." };
  }
  auto getter_position = ctrw::Get_new_from_particle{
    ctrw::Get_position_periodic{ boundaries.boundary_periodic } };
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
      auto output_positions = useful::open_write(filename_output_positions);
      output_positions << std::setprecision(8)
                       << std::scientific;
      
      for (auto const& time : measure_times)
      {
        ptrw.evolve(time);
        output_positions << time;
        for (auto const& part : ptrw.particles())
          useful::print(output_positions, getter_position(part), 1);
        output_positions << "\n";
        std::cout << "time = " << time
                  << "\ttime_last_measure = " << time_max
                  << "\n";
      }
      output_positions.close();

      break;
    }
    case 1:
    {
      auto output_positions = useful::open_write(filename_output_positions);
      output_positions << std::setprecision(8)
                       << std::scientific;
      
      for (auto const& time : measure_times)
      {
        ptrw.evolve(time);
        std::cout << "time = " << time
                  << "\ttime_last_measure = " << time_max
                  << "\n";
      }
      for (auto const& part : ptrw.particles())
      {
        useful::print(output_positions, getter_position(part));
        output_positions << "\n";
      }
      output_positions.close();
      
      break;
    }
    case 2:
    {
      boundary::Periodic boundary_cubic_cell{ geometry.boundaries };
      auto getter_position_cubic_cell =
      [&boundary_cubic_cell,&boundaries,&geometry]
      (CTRW::Particle const& particle)
      {
        auto state_copy = particle.state_new();
        boundaries.boundary_periodic.translate(state_copy.position,state_copy.periodicity);
        boundary_cubic_cell(state_copy);
          
        return state_copy.position;
      };
      
      std::string filename_output_correlation_time = output_dir + "/" +
        filename_output_base + "_position_fluctuations_autocorrelation_time_" + initial_condition_name +
        "_" + data_set + "_" + params + ".dat";
      auto output_correlation_time = useful::open_write(filename_output_correlation_time);
      output_correlation_time << std::setprecision(8)
                              << std::scientific;
      
      std::cout << std::setprecision(2)
                << std::scientific;
        
      ptrw.evolve(measure_times[0]);
      
      std::cout << "\tTime [adv times] = " << measure_times[0]/advection_time << "\t"
                << "Max time [adv times] = " << measure_times.back()/advection_time << "\n";
      
      std::vector<State::Position> initial_position(nr_particles);
      for (std::size_t pp = 0; pp < nr_particles; ++pp)
        initial_position[pp] = getter_position_cubic_cell(ctrw.particles(pp));
      auto initial_position_mean = State::Position(Geometry::dim);
      for (auto const& part : ctrw.particles())
        operation::plus_InPlace(initial_position_mean, getter_position_cubic_cell(part));
      operation::div_scalar_InPlace(initial_position_mean, nr_particles);
      double autocorrelation = 0.;
      for (std::size_t pp = 0; pp < nr_particles; ++pp)
        autocorrelation +=
          operation::dot(operation::minus(getter_position_cubic_cell(ctrw.particles(pp)),
                                          initial_position_mean),
                         operation::minus(initial_position[pp],
                                          initial_position_mean));
      autocorrelation /= nr_particles;
      output_correlation_time << measure_times[0] << "\t"
                              << autocorrelation << "\n";
      
      for (std::size_t tt = 1; tt < measure_times.size(); ++tt)
      {
        ptrw.evolve(measure_times[tt]);
        
        std::cout << "\tTime [adv times] = " << measure_times[tt]/advection_time << "\t"
                  << "Max time [adv times] = " << measure_times.back()/advection_time << "\n";
        
        auto position_mean = State::Position(Geometry::dim);
        for (auto const& part : ctrw.particles())
          operation::plus_InPlace(position_mean, getter_position_cubic_cell(part));
        operation::div_scalar_InPlace(position_mean, nr_particles);
        double autocorrelation = 0.;
        for (std::size_t pp = 0; pp < nr_particles; ++pp)
          autocorrelation +=
            operation::dot(operation::minus(getter_position_cubic_cell(ctrw.particles(pp)),
                                            position_mean),
                           operation::minus(initial_position[pp],
                                            initial_position_mean));
        autocorrelation /= nr_particles;
        output_correlation_time << measure_times[tt] << "\t"
                                << autocorrelation << "\n";
      }
      break;
    }
    default:
      throw std::invalid_argument{ "Undefined measure type" };
  }
  
  std::cout << "\tDone!\n";
  
  return 0;
}
