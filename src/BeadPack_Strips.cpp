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
#include <vector>
#include "BeadPack/BeadPack_Models.h"
#include "Field/VectorField_Interpolated.h"
#include "general/Operations.h"
#include "general/Ranges.h"
#include "general/useful.h"
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
  using namespace beadpack::model_bcc_cartesian;
  
  if (argc == 1)
  {
    std::cout << "Track particle strips for\n"
              << "fully-advective particle tracking in 3d beadpacks\n"
              << "with periodic boundary conditions on a cubic domain.\n"
              << "----------------------------------------------------\n"
              << "Parameters (default value in []):\n"
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
  
  if (argc < 11)
    throw std::runtime_error{ "Inappropriate parameters." };
  
  using State = ctrw::State_periodic<std::vector<double>,
    std::vector<int>, useful::Empty, useful::Empty, std::size_t>;
  using CTRW = ctrw::CTRW<State>;
  using Boundary = Boundaries::Boundary_Reflecting_Periodic;
  using JumpGenerator =
    ctrw::JumpGenerator_Velocity_withHint_RK4<VelocityField&, Boundary&>;
  
  std::size_t arg = 1;
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
  double initial_strip_segment_length = initial_strip_segment_length_factor*geometry.domain_side;
  double advection_time = geometry.domain_side/magnitude_mean_velocity;
  double time_step = time_step_accuracy_adv*advection_time;
  std::size_t max_particles = max_particles_strip*nr_strips;
  std::size_t particles_strip = 2;
  double max_distance_strip = 1e-2*geometry.domain_side;
  auto adjust = [&bead_pack, &boundaries](State& state)
  { bead_pack.place_at_closest_surface_if_inside(state, boundaries.boundary_periodic); };
  auto state_maker = [&geometry]()
  { return State{ std::vector<double>(geometry.dim),
    std::vector<int>(geometry.dim) }; };
  CTRW ctrw{};
  ctrw::StripHandler strips{
    ctrw,
    ctrw::Get_new_from_particle{ ctrw::Get_position_periodic{ boundaries.boundary_periodic } },
    adjust,
    state_maker };
  beadpack::make_strips_random_uniform_box(nr_strips,
                                           particles_strip, max_particles_strip,
                                           initial_strip_segment_length,
                                           max_distance_strip,
                                           bead_pack, boundaries.boundary_periodic, geometry.boundaries,
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
  stream << nr_strips << "_"
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
  
  auto getter_position = ctrw::Get_new_from_particle{
    ctrw::Get_position_periodic{ boundaries.boundary_periodic } };
  std::cout << "\tDone!\n";
  
  std::cout << "Setting up dynamics...\n";
  ctrw::Transitions_Position transitions{
    JumpGenerator{
      velocity_field,
      time_step,
      max_particles,
      boundaries.boundary_reflecting_periodic,
    },
    boundaries.boundary_reflecting_periodic
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
    measurer_positions(ptrw, getter_position, ptrw.time());
    measurer_strips(strips, ptrw.time());
    while (ptrw.time() < time_max)
    {
      ptrw.step();
      strips.resize();
      measurer_positions(ptrw, getter_position, ptrw.time());
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
      measurer_positions(ptrw, getter_position, time);
      measurer_strips(strips, time);
      std::cout << "\ttime [adv times] = " << time/advection_time << "\t"
                << "\ttime_max [adv times]  = " << time_max/advection_time << "\n";
    }
  }
  std::cout << "\tDone!\n";
  
  return 0;
}
