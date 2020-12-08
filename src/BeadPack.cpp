//
//  BeadPack.cpp
//
//  Created by Tomás Aquino on 19/11/2020.
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
    std::cout << "Advective-diffusive particle tracking in 3d beadpacks\n"
              << "with periodic boundary conditions on a cubic domain.\n"
              << "----------------------------------------------------\n"
              << "Parameters (default value in []):\n"
              << "domain_side : Length of domain side or periodic unit cell\n"
              << "peclet : Peclet number in terms of domain side, average velocity,\n"
              << "         and diffusion coefficient\n"
              << "time_step_accuracy_adv : Maximum time step size in units of advection time\n"
              << "time_step_accuracy_diff : Minimum time step size in units of advection time\n"
              << "time_min_diffusion_times : Minimum output time in units of diffusion time\n"
              << "time_max_diffusion_times : Maximum output time in units of diffusion time\n"
              << "nr_measures : Number of measurements\n"
              << "measure_spacing : 0 - Logarithmic spacing between measurements\n"
              << "                : 1 - Linear spacing between measurements\n"
              << "measure_type : 0 - Output time and all particle positions per file line,\n"
              << "                   at each measurement time\n"
              << "               1 - Output particle positions, one per file line, at final time only\n"
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
              << "initial_condition_size_domains : Size of initial condition box or plane in domain sides\n"
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
  
  if (argc != 15 && argc != 16 && argc != 17 && argc != 18)
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
  double time_min_diffusion_times = atof(argv[arg++]);
  double time_max_diffusion_times = atof(argv[arg++]);
  std::size_t nr_measures = strtoul(argv[arg++], NULL, 0);
  int measure_spacing = atoi(argv[arg++]);
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
    
  std::vector<double> domain_dimensions(dim, domain_side);
  std::vector<std::pair<double, double>> boundaries;
  boundaries.reserve(dim);
  for (std::size_t dd = 0; dd < dim; ++dd)
    boundaries.push_back({ 0., domain_side });
  
  std::cout << "Importing beads...\n";
  std::string bead_filename = input_dir + "/" + "spheres.dat";
  BeadPack bead_pack{ beadpack::get_beads<Bead>(dim, bead_filename, 1, domain_side) };
  std::cout << "\tDone!\n";
  
  std::cout << "Importing contacts...\n";
  std::string contact_filename = input_dir + "/" + "contacts.dat";
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
    std::size_t nr_samples = 1e4;
    mean_velocity = bead_pack.compute_mean_vector(velocity_field, boundaries, nr_samples);
    std::string mean_velocity_filename = output_dir + "/" + "mean_velocity.dat";
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
  
  Boundary_Periodic boundary_periodic{ boundaries };
  Boundary boundary{ bead_pack, boundary_periodic };
  
  auto state_maker = []()
  { return State{ std::vector<double>(dim, 0.), std::vector<int>(dim, 0) }; };
  
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
  
  std::stringstream stream;
  stream << std::scientific << std::setprecision(2);
  stream << domain_side << "_"
         << initial_condition_size_domains << "_"
         << peclet << "_"
         << time_step_accuracy_diff << "_"
         << time_step_accuracy_adv << "_"
         << time_min_diffusion_times << "_"
         << time_max_diffusion_times << "_"
         << nr_measures << "_"
         << measure_spacing << "_"
         << nr_particles << "_"
         << run_nr;
  std::string params = stream.str();
  
  std::string filename_output_base = "Data_beadpack";
  
  std::string filename_output_positions = output_dir + "/" +
    filename_output_base + "_positions_" + initial_condition_name +
    "_" + data_set + "_" + params + ".dat";
  
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
  
  auto getter_position = ctrw::Get_position_periodic{ domain_dimensions };
  std::cout << "\tDone!\n";
  
  std::cout << "Setting up dynamics...\n";
  ctrw::Transitions_Position transitions{
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
      std::ofstream output_positions{ filename_output_positions };
      if (!output_positions.is_open())
        throw useful::open_write_error(filename_output_positions);
      output_positions << std::setprecision(8)
                       << std::scientific;
      
      for (auto const& time : measure_times)
      {
        ptrw.evolve(time);
        output_positions << time;
        for (auto const& part : ptrw.particles())
          useful::print(output_positions, getter_position(part.state_new()), 1);
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
      std::ofstream output_positions{ filename_output_positions };
      if (!output_positions.is_open())
        throw useful::open_write_error(filename_output_positions);
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
        useful::print(output_positions, getter_position(part.state_new()));
        output_positions << "\n";
      }
      output_positions.close();
      
      break;
    }
    default:
      throw std::invalid_argument{ "Undefined measure type" };
  }
  
  std::cout << "\tDone!\n";
  
  return 0;
}
