//
//  BeadPack_Strips.cpp
//
//  Created by Tomás Aquino on 12/05/2020.
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
#include "Stochastic/CTRW/ParticleCollections.h"
#include "Stochastic/CTRW/PTRW.h"
#include "Stochastic/CTRW/State.h"
#include "Stochastic/CTRW/StateGetter.h"
#include "Stochastic/CTRW/Transitions_State.h"

int main(int argc, const char * argv[])
{
  if (argc == 1)
  {
    std::cout << "Track particle strips for\n"
              << "fully-advective particle tracking in 3d beadpacks\n"
              << "with periodic boundary conditions on a cubic domain.\n"
              << "----------------------------------------------------\n"
              << "Parameters (default value in []):\n"
              << "domain_side : Length of domain side or periodic unit cell\n"
              << "nr_strips : Number of strips to track\n"
              << "max_particles_strip : Maximum number of particles per strip\n"
              << "initial_strip_segment_length_factor : Initial strip length in units of domain sides\n"
              << "time_step_accuracy_adv : Time step size in units of advection time\n"
              << "time_min_advection_times : Minimum output time in units of advection time\n"
              << "time_max_advection_times : Maximum output time in units of advection time\n"
              << "nr_measures : Number of measurements\n"
              << "measure_spacing : 0 - Logarithmic spacing between measurements\n"
              << "                : 1 - Linear spacing between measurements\n"
              << "                : 2 - Output only at last step\n"
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
  
  if (argc != 12 && argc != 13 && argc != 14)
  {
    throw std::runtime_error{ "Inappropriate parameters." };
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
  
  std::size_t arg = 1;
  double domain_side = atof(argv[arg++]);
  std::size_t nr_strips = strtoul(argv[arg++], NULL, 0);
  std::size_t max_particles_strip = strtoul(argv[arg++], NULL, 0);
  double initial_strip_segment_length_factor = atof(argv[arg++]);
  double time_step_accuracy_adv = atof(argv[arg++]);
  double time_min_advection_times = atof(argv[arg++]);
  double time_max_advection_times = atof(argv[arg++]);
  std::size_t nr_measures = strtoul(argv[arg++], NULL, 0);
  int measure_spacing = atoi(argv[arg++]);
  std::size_t run_nr = strtoul(argv[arg++], NULL, 0);
  std::string data_set = argv[arg++];
  std::string input_dir_base = argc > arg ? argv[arg++] : "../input";
  std::string output_dir = argc > arg ? argv[arg++] : "../output";
  
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
  double initial_strip_segment_length = initial_strip_segment_length_factor*domain_side;
  double advection_time = domain_side/mean_velocity_magnitude;
  double time_step = time_step_accuracy_adv*advection_time;
  std::size_t max_particles = max_particles_strip*nr_strips;
  std::size_t particles_strip = 2;
  double max_distance_strip = 1e-2*domain_side;
  Boundary_Periodic boundary_periodic{ boundaries };
  Boundary boundary{ bead_pack, boundary_periodic };
  auto adjust = [&bead_pack, &boundary_periodic](State& state)
  { bead_pack.place_at_closest_surface_if_inside(state, boundary_periodic); };
  auto state_maker = []()
  { return State{ std::vector<double>(dim, 0.), std::vector<int>(dim, 0) }; };
  CTRW ctrw{};
  ctrw::StripHandler strips{
    ctrw,
    ctrw::Get_new_from_particle<ctrw::Get_position_periodic>{ domain_dimensions },
    adjust,
    state_maker };
  beadpack::make_strips_random_uniform_box(nr_strips,
                                           particles_strip, max_particles_strip,
                                           initial_strip_segment_length,
                                           max_distance_strip,
                                           bead_pack, boundary_periodic, boundaries,
                                           ctrw, strips);
  std::cout << "\tDone!\n";
  
  std::cout << "Setting up output...\n";
  double time_min = time_min_advection_times*advection_time;
  double time_max = time_max_advection_times*advection_time;
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
  
  std::stringstream stream;
  stream << std::scientific << std::setprecision(2);
  stream << domain_side << "_"
         << nr_strips << "_"
         << max_particles_strip << "_"
         << initial_strip_segment_length_factor << "_"
         << time_step_accuracy_adv << "_"
         << time_min_advection_times << "_"
         << time_max_advection_times << "_"
         << nr_measures << "_"
         << measure_spacing << "_"
         << run_nr;
  std::string params = stream.str();
  
  std::string filename_output_postions_base = "Data_beadpack_positions";
  std::string filename_output_positions = output_dir + "/"
    + filename_output_postions_base + "_" + data_set + "_" + params + ".dat";
  ctrw::Measurer_Particle measurer_positions{ filename_output_positions };
  
  std::string filename_output_strips_base = "Data_beadpack_strips";
  std::string filename_output_strips = output_dir + "/"
    + filename_output_strips_base + "_" + data_set + "_" + params + ".dat";
  ctrw::Measurer_Collection measurer_strips{ filename_output_strips };
  
  using Get_position = ctrw::Get_new_from_particle<ctrw::Get_position_periodic>;
  std::cout << "\tDone!\n";
  
  std::cout << "Setting up dynamics...\n";
  ctrw::Transitions_Position transitions{
    ctrw::JumpGenerator_Velocity_withHint_RK4{
      velocity_field,
      time_step,
      max_particles,
      boundary
    },
    boundary
  };
  ctrw::PTRW ptrw(ctrw, transitions, time_step, 0.);
  std::cout << "\tDone!\n";
  
  std::cout << "Starting dynamics...\n";
  std::cout << "\tAdvection time = " << advection_time << "\n";
  std::cout << "\ttime step [adv times] = " << time_step/advection_time << "\n";
  std::cout << "\tmax time [adv times] = " << time_max/advection_time << "\n";
  std::cout << "\tnr_particles = " << ptrw.size() << "\n";
  
  if (measure_spacing == 2)
  {
    measurer_positions(ptrw, Get_position{ domain_dimensions }, ptrw.time());
    measurer_strips(strips, ptrw.time());
    while (ptrw.time() < time_max)
    {
      ptrw.step();
      strips.resize();
      measurer_positions(ptrw, Get_position{ domain_dimensions }, ptrw.time());
      measurer_strips(strips, ptrw.time());
      std::cout << "\ttime [adv times] = " << ptrw.time()/advection_time << "\t"
                << "\ttime_max [adv times] = " << time_max/advection_time << "\n";
    }
  }
  else
  {
    for (auto time : measure_times)
    {
      while (ptrw.time() < time)
      {
        ptrw.step();
        strips.resize();
      }
      measurer_positions(ptrw, Get_position{ domain_dimensions }, time);
      measurer_strips(strips, time);
      std::cout << "\ttime [adv times] = " << time/advection_time << "\t"
                << "\ttime_max [adv times]  = " << time_max/advection_time << "\n";
    }
  }
  std::cout << "\tDone!\n";
  
  return 0;
}
