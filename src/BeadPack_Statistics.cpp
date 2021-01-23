//
//  BeadPack_Statistics
//
//  Created by Tomás Aquino on 11/01/2021.
//  Copyright © 2021 Tomás Aquino. All rights reserved.
//

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "BeadPack/BeadPack_Models.h"
#include "general/Operations.h"
#include "general/Ranges.h"
#include "general/useful.h"
#include "Stochastic/CTRW/CTRW.h"
#include "Stochastic/CTRW/JumpGenerator.h"
#include "Stochastic/CTRW/PTRW.h"
#include "Stochastic/CTRW/Reaction.h"
#include "Stochastic/CTRW/State.h"
#include "Stochastic/CTRW/StateGetter.h"
#include "Stochastic/CTRW/Transitions_State.h"

int main(int argc, const char * argv[])
{
  using namespace model_beadpack_cartesian_cubic;
  
  if (argc == 1)
  {
    std::cout << "Statistics for 3d beadpacks\n"
              << "with periodic boundary conditions.\n"
              << "----------------------------------------------------\n"
              << "Parameters (default value in []):\n"
              << "nr_samples : Number of statistical samples\n"
              << "measure_type : 0 - Compute Eulerian mean velocity vector\n"
              << "               1 - Compute Eulerian mean velocity magnitude\n"
              << "               2 - Compute porosity\n"
              << "               3 - Compute tortuosity\n"
              << "               4 - Compute uniform Eulerian velocity magnitude samples\n"
              << "               5 - Compute space-Lagrangian velocity magnitude mean\n"
              << "               6 - Compute time-Lagrangian velocity magnitude mean\n"
              << "               7 - Compute space-Lagrangian velocity magnitude autocorrelation\n"
              << "               8 - Compute time-Lagrangian velocity magnitude autocorrelation\n"
              << "               9 - Compute space-Lagrangian velocity magnitude fluctuations autocorrelation\n"
              << "               10 - Compute time-Lagrangian velocity magnitude fluctuations autocorrelation\n"
              << "               11 - Compute space-Lagrangian velocity magnitude series\n"
              << "               12 - Compute time-Lagrangian velocity magnitude series\n"
              << "run_nr : Nonnegative integer identifier for output files\n"
              << "data_set : Path to input data folder relative to input_dir_base\n"
              << "           (and model name identifier for output files)\n"
              << "(if measure_type >= 4) jump_size_domains : Jump size in units of domain side\n"
              << "(if measure_type >= 4) measure_max : Maximum distance along\n"
              << "                                     streamlines in units of domain side\n"
                 "                                     or maximum time in units of\n"
              << "                                     advection time\n"
              << "(if measure_type >= 6) nr_measures : Number of measurements along streamlines\n"
              << "input_dir_base : Path to look for input data [../input]\n"
              << "output_dir : Path folder to output to [../output]\n";
    return 0;
  }
  
  if (argc < 5)
    throw useful::bad_parameters();
  
  using State = ctrw::State_periodic<std::vector<double>,
    std::vector<int>, useful::Empty, double, std::size_t>;
  using CTRW = ctrw::CTRW<State>;
  using Boundary = Boundaries::Boundary_Reflecting_Periodic;
  using JumpGenerator =
    ctrw::JumpGenerator_Velocity_withHint_RK4<VelocityField&, Boundary&>;
  
  std::size_t arg = 1;
  std::size_t nr_samples = strtoul(argv[arg++], NULL, 0);
  int measure_type = atoi(argv[arg++]);
  std::size_t run_nr = strtoul(argv[arg++], NULL, 0);
  std::string data_set = argv[arg++];
  
  double jump_size_domains = 0.;
  double measure_max = 0.;
  if (measure_type >= 5)
  {
    if (argc < 6)
      throw useful::bad_parameters();
    jump_size_domains = atof(argv[arg++]);
    measure_max = atof(argv[arg++]);
  }
  
  std::size_t nr_measures = 0;
  if (measure_type >= 7)
  {
    if (argc < 7)
      throw useful::bad_parameters();
    nr_measures = strtoul(argv[arg++], NULL, 0);
  }
  
  std::string input_dir_base = argc > arg ? argv[arg++] : "../input";
  std::string output_dir = argc > arg ? argv[arg++] : "../output";
  
  std::string input_dir = input_dir_base + "/" + data_set;
  
  Geometry geometry{};
  
  std::cout << "Making bead pack...\n";
  BeadPack bead_pack = make_bead_pack(input_dir, geometry);
  std::cout << "\tDone!\n";
  
  std::cout << "Setting up boundary conditions...\n";
  Boundaries boundaries{ geometry, bead_pack };
  std::cout << "\tDone!\n";
  
  std::cout << "Setting up velocity field...\n";
  auto velocity_field =
    make_velocity_field(input_dir, output_dir,
                        geometry, bead_pack, boundaries.boundary_periodic);
  std::cout << "\tDone!\n";
  
  auto state_maker = [&geometry]()
  { return State{
    std::vector<double>(geometry.dim),
    std::vector<int>(geometry.dim) }; };
  auto getter_position = ctrw::Get_new_from_particle{
    ctrw::Get_position_periodic{ boundaries.boundary_periodic } };
  auto getter_position_old = ctrw::Get_old_from_particle{
    ctrw::Get_position_periodic{ boundaries.boundary_periodic } };
  auto getter_velocity = ctrw::Get_new_from_particle{
    ctrw::Get_position_property{ velocity_field } };
  auto getter_velocity_old = ctrw::Get_old_from_particle{
    ctrw::Get_position_property{ velocity_field } };
  auto getter_velocity_magnitude = [&getter_velocity]
  (CTRW::Particle const& particle)
  { return operation::abs(getter_velocity(particle)); };
  auto getter_velocity_magnitude_old = [&getter_velocity_old]
  (CTRW::Particle const& particle)
  { return operation::abs(getter_velocity_old(particle)); };
  
  switch (measure_type)
  {
    case 0:
    {
      std::cout << "Computing Eulerian mean velocity vector...\n";
      auto particles = beadpack::make_particles_random_uniform_box<CTRW::Particle>(
        nr_samples, bead_pack, boundaries.boundary_periodic, geometry.boundaries, state_maker);
      std::vector<double> mean_velocity(geometry.dim);
      for (auto const& part : particles)
        operation::plus_InPlace(mean_velocity, getter_velocity(part));
      operation::div_scalar_InPlace(mean_velocity, particles.size());
      
      std::stringstream stream_measures;
      stream_measures << nr_samples;
      std::string mean_velocity_filename = output_dir + "/" + "mean_velocity"
        + "_" + data_set + "_" + stream_measures.str() + ".dat";
      std::ofstream output{ mean_velocity_filename };
      if (!output.is_open())
        throw useful::open_write_error(mean_velocity_filename);
      output << std::setprecision(12)
             << std::scientific;
      useful::print(output, mean_velocity);
      output << "\n";
      output.close();
      std::cout << "\tDone!\n";
      break;
    }
    case 1:
    {
      std::cout << "Computing Eulerian mean velocity magnitude...\n";
      auto particles = beadpack::make_particles_random_uniform_box<CTRW::Particle>(
        nr_samples, bead_pack, boundaries.boundary_periodic, geometry.boundaries, state_maker);
      double mean_velocity_magnitude = 0.;
      for (auto const& part : particles)
        mean_velocity_magnitude += getter_velocity_magnitude(part);
      mean_velocity_magnitude /= particles.size();
      std::stringstream stream_measures;
      stream_measures << nr_samples << "_" << run_nr;
      std::string mean_velocity_magnitude_filename = output_dir + "/"
        + "mean_velocity_magnitude" + "_"
        + data_set + "_" + stream_measures.str() + ".dat";
      std::ofstream output{ mean_velocity_magnitude_filename };
      if (!output.is_open())
        throw useful::open_write_error(mean_velocity_magnitude_filename);
      output << std::setprecision(12)
             << std::scientific;
      output << mean_velocity_magnitude << "\n";
      output.close();
      std::cout << "\tDone!\n";
      break;
    }
    case 2:
    {
      std::cout << "Computing porosity...\n";
      double porosity = bead_pack.compute_porosity(geometry.boundaries, nr_samples);
      std::stringstream stream_measures;
      stream_measures << nr_samples << "_" << run_nr;
      std::string filename = output_dir + "/" + "porosity" + "_" + data_set + "_"
        + stream_measures.str() + ".dat";
      std::ofstream output{ filename };
      if (!output.is_open())
        throw useful::open_write_error(filename);
      output << std::setprecision(12)
             << std::scientific;
      output << porosity << "\n";
      output.close();
      std::cout << "\tDone!\n";
      break;
    }
    case 3:
    {
      std::cout << "Computing tortuosity...\n";
      
      std::cout << "\tImporting mean velocity...\n";
      std::string mean_velocity_filename = input_dir + "/" + "mean_velocity.dat";
      std::vector<double> mean_velocity =
        beadpack::get_mean_velocity(geometry.dim, mean_velocity_filename);
      std::ifstream input{ mean_velocity_filename };
      std::string mean_velocity_magnitude_filename = output_dir + "/"
        + "mean_velocity_magnitude.dat";
      if (!input.is_open())
        throw useful::open_read_error(mean_velocity_magnitude_filename);
      double mean_velocity_magnitude;
      input >> mean_velocity_magnitude;
      std::cout << "\t\tDone!\n";
      
      auto getter_velocity_downstream = ctrw::Get_old_from_particle{
        ctrw::Get_projection{
          ctrw::Get_position_property{ velocity_field },
          mean_velocity } };
      auto particles = beadpack::make_particles_random_uniform_box<CTRW::Particle>(
        nr_samples, bead_pack, boundaries.boundary_periodic, geometry.boundaries, state_maker);
      double mean_velocity_downstream = 0.;
      for (auto const& part : particles)
        mean_velocity_downstream += getter_velocity_downstream(part);
      mean_velocity_downstream /= particles.size();
      double tortuosity = mean_velocity_magnitude/mean_velocity_downstream;
      std::stringstream stream_measures;
      stream_measures << nr_samples << "_" << run_nr;
      std::string filename = output_dir + "/" + "tortuosity"
        + "_" + data_set + "_" + stream_measures.str() + ".dat";
      std::ofstream output{ filename };
      if (!output.is_open())
        throw useful::open_write_error(filename);
      output << std::setprecision(12)
             << std::scientific;
      output << tortuosity << "\n";
      output.close();
      std::cout << "\tDone!\n";
      break;
    }
    case 4:
    {
      std::cout << "Computing Eulerian velocity magnitude samples...\n";

      auto particles = beadpack::make_particles_random_uniform_box<ctrw::Particle<State>>(
        nr_samples, bead_pack, boundaries.boundary_periodic, geometry.boundaries, state_maker);
     
      std::stringstream stream_measures;
      stream_measures << nr_samples << "_" << run_nr;
      std::string filename = output_dir + "/"
        + "velocity_magnitude_samples_uniform_unit_cell"
        + "_" + data_set + "_" + stream_measures.str() + ".dat";
      std::ofstream output{ filename };
      if (!output.is_open())
        throw useful::open_write_error(filename);
      output << std::setprecision(12)
             << std::scientific;
      for (auto const& part : particles)
      {
        useful::print(output, getter_velocity_magnitude(part));
        output << "\n";
      }
      output.close();
      std::cout << "\tDone!\n";
      break;
    }
    case 5:
    {
      std::cout << "Computing space-Lagrangian velocity magnitude mean...\n";
      
      std::cout << "\tSetting up trajectories...\n";
      auto particles = beadpack::make_particles_random_uniform_box<CTRW::Particle>(
        nr_samples, bead_pack, boundaries.boundary_periodic, geometry.boundaries, state_maker);
      std::cout << "\t\tDone!\n";
      
      double jump_size = jump_size_domains*geometry.domain_side;
      double distance_max = measure_max*geometry.domain_side;
      
      std::vector<double> velocity_mean(particles.size());
      std::size_t discarded = 0;
      for (std::size_t pp = 0; pp < particles.size(); ++pp)
      {
        std::cout << "\tTrajectory " << pp+1 << " of " << particles.size() << "\n";
        CTRW ctrw{ { particles[pp] }, CTRW::Tag{} };
        ctrw::Transitions_Position_VelocityStep transitions{
            jump_size,
            JumpGenerator{
              velocity_field,
              0.,
              1,
              boundaries.boundary_reflecting_periodic },
            boundaries.boundary_reflecting_periodic
        };
        
        auto const& part = ctrw.particles(0);
        double distance_traveled = 0.;
        while (distance_traveled < distance_max)
        {
          double velocity_magnitude = getter_velocity_magnitude(part);
          if (velocity_magnitude == 0.)
            break;
          ctrw.step(transitions);
          double distance_increment = operation::abs(operation::minus(getter_position(part),
            getter_position_old(part)));
          distance_traveled += distance_increment;
          velocity_mean[pp] += getter_velocity_magnitude_old(part)*distance_increment;
        }
        if (getter_velocity_magnitude(part) == 0.)
          ++discarded;
        if (distance_traveled != 0.)
          velocity_mean[pp] /= distance_traveled;
      }
      std::cout << discarded << " trajectories discarded due to zero velocity\n";
      
      std::stringstream stream;
      stream << std::scientific << std::setprecision(2);
      stream << jump_size_domains << "_"
             << measure_max << "_"
             << nr_samples << "_"
             << run_nr;
      std::string filename = output_dir + "/"
        + "velocity_magnitude_mean_space_lagrangian"
        + "_" + data_set + "_" + stream.str() + ".dat";
      std::ofstream output{ filename };
      if (!output.is_open())
        throw useful::open_write_error(filename);
      output << std::setprecision(12)
             << std::scientific;
      useful::print(output, velocity_mean);
      output << "\n";
      output.close();
      std::cout << "\tDone!\n";
      break;
    }
    case 6:
    {
      std::cout << "Computing time-Lagrangian velocity magnitude mean...\n";
      
      std::cout << "\tImporting mean velocity...\n";
      std::string mean_velocity_filename = input_dir + "/" + "mean_velocity.dat";
      std::vector<double> mean_velocity =
        beadpack::get_mean_velocity(geometry.dim, mean_velocity_filename);
      std::cout << "\t\tDone!\n";
      
      std::cout << "\tSetting up trajectories...\n";
      auto particles = beadpack::make_particles_random_uniform_box<CTRW::Particle>(
        nr_samples, bead_pack, boundaries.boundary_periodic, geometry.boundaries, state_maker);
      std::cout << "\t\tDone!\n";
      
      double jump_size = jump_size_domains*geometry.domain_side;
      double time_max = measure_max*geometry.domain_side/operation::abs(mean_velocity);
      
      std::vector<double> velocity_mean(particles.size());
      std::size_t surviving_trajectories = nr_samples;
      for (std::size_t pp = 0; pp < particles.size(); ++pp)
      {
        std::cout << "\tTrajectory " << pp+1 << " of " << particles.size() << "\n";
        CTRW ctrw{ { particles[pp] }, CTRW::Tag{} };
        ctrw::Transitions_Position_VelocityStep transitions{
            jump_size,
            JumpGenerator{
              velocity_field,
              0.,
              1,
              boundaries.boundary_reflecting_periodic },
            boundaries.boundary_reflecting_periodic
        };
        
        auto const& part = ctrw.particles(0);
        while (part.state_new().time < time_max)
        {
          double velocity_magnitude = getter_velocity_magnitude(part);
          if (velocity_magnitude == 0.)
            break;
          ctrw.step(transitions);
          velocity_mean[pp] += getter_velocity_magnitude_old(part)*transitions.time_step();
        }
        if (getter_velocity_magnitude(part) == 0.)
          --surviving_trajectories;
        if (velocity_mean[pp] != 0.)
          velocity_mean[pp] /= std::min(time_max, part.state_new().time);
      }
      std::cout << nr_samples - surviving_trajectories
                << " trajectories discarded due to zero velocity\n";
      
      std::stringstream stream;
      stream << std::scientific << std::setprecision(2);
      stream << jump_size_domains << "_"
             << measure_max << "_"
             << nr_samples << "_"
             << run_nr;
      std::string filename = output_dir + "/"
        + "velocity_magnitude_mean_time_lagrangian"
        + "_" + data_set + "_" + stream.str() + ".dat";
      std::ofstream output{ filename };
      if (!output.is_open())
        throw useful::open_write_error(filename);
      output << std::setprecision(12)
             << std::scientific;
      useful::print(output, velocity_mean);
      output << "\n";
      output.close();
      std::cout << "\tDone!\n";
      break;
    }
    case 7:
    {
      std::cout << "Computing space-Lagrangian velocity magnitude autocorrelation...\n";
      std::vector<double> distances
        = range::linspace(0., measure_max*geometry.domain_side, nr_measures);
      
      std::cout << "\tSetting up trajectories...\n";
      auto particles = beadpack::make_particles_random_uniform_box<CTRW::Particle>(
        nr_samples, bead_pack, boundaries.boundary_periodic, geometry.boundaries, state_maker);
      std::cout << "\t\tDone!\n";
      
      double jump_size = jump_size_domains*geometry.domain_side;
      
      std::vector<double> velocity_autocorrelation(distances.size());
      std::vector<std::size_t> surviving_trajectories(distances.size(), nr_samples);
      for (std::size_t pp = 0; pp < particles.size(); ++pp)
      {
        std::cout << "\tTrajectory " << pp+1 << " of " << particles.size() << "\n";
        CTRW ctrw{ { particles[pp] }, CTRW::Tag{} };
        ctrw::Transitions_Position_VelocityStep transitions{
            jump_size,
            JumpGenerator{
              velocity_field,
              0.,
              1,
              boundaries.boundary_reflecting_periodic },
            boundaries.boundary_reflecting_periodic
        };
        
        auto const& part = ctrw.particles(0);
        double initial_velocity = getter_velocity_magnitude(part);
        double distance_traveled = 0.;
        for (std::size_t ss = 0; ss < distances.size(); ++ss)
        {
          while (distance_traveled < distances[ss])
          {
            double velocity_magnitude = getter_velocity_magnitude(part);
            if (velocity_magnitude == 0.)
              break;
            ctrw.step(transitions);
            distance_traveled +=
              operation::abs(operation::minus(getter_position(part),
                                              getter_position_old(part)));
          }
          velocity_autocorrelation[ss] +=
            initial_velocity*getter_velocity_magnitude_old(part);
          if (getter_velocity_magnitude(part) == 0.)
          {
            for (std::size_t ss2 = ss+1; ss2 < distances.size(); ++ss2)
              --surviving_trajectories[ss2];
            break;
          }
        }
      }
      operation::div_InPlace(velocity_autocorrelation, surviving_trajectories);
      std::cout << nr_samples - surviving_trajectories.back()
                << " trajectories discarded due to zero velocity\n";
      
      std::stringstream stream;
      stream << std::scientific << std::setprecision(2);
      stream << jump_size_domains << "_"
             << measure_max << "_"
             << nr_measures << "_"
             << nr_samples << "_"
             << run_nr;
      std::string filename = output_dir + "/"
        + "velocity_magnitude_autocorrelation_space_lagrangian"
        + "_" + data_set + "_" + stream.str() + ".dat";
      std::ofstream output{ filename };
      if (!output.is_open())
        throw useful::open_write_error(filename);
      output << std::setprecision(12)
             << std::scientific;
      for (std::size_t ss = 0; ss < distances.size(); ++ss)
      {
        output << distances[ss] << "\t"
               << velocity_autocorrelation[ss] << "\n";
      }
      output.close();
      std::cout << "\tDone!\n";
      break;
    }
    case 8:
    {
      std::cout << "Computing time-Lagrangian velocity magnitude autocorrelation...\n";
      
      std::cout << "\tImporting mean velocity...\n";
      std::string mean_velocity_filename = input_dir + "/" + "mean_velocity.dat";
      std::vector<double> mean_velocity =
        beadpack::get_mean_velocity(geometry.dim, mean_velocity_filename);
      std::cout << "\t\tDone!\n";
      
      std::vector<double> times
        = range::linspace(0., measure_max*geometry.domain_side/operation::abs(mean_velocity),
                          nr_measures);
      
      std::cout << "\tSetting up trajectories...\n";
      auto particles = beadpack::make_particles_random_uniform_box<CTRW::Particle>(
        nr_samples, bead_pack, boundaries.boundary_periodic, geometry.boundaries, state_maker);
      std::cout << "\t\tDone!\n";
      
      double jump_size = jump_size_domains*geometry.domain_side;
      
      std::vector<double> velocity_autocorrelation(times.size());
      std::vector<std::size_t> surviving_trajectories(times.size(), nr_samples);
      for (std::size_t pp = 0; pp < particles.size(); ++pp)
      {
        std::cout << "\tTrajectory " << pp+1 << " of " << particles.size() << "\n";
        CTRW ctrw{ { particles[pp] }, CTRW::Tag{} };
        ctrw::Transitions_Position_VelocityStep transitions{
            jump_size,
            JumpGenerator{
              velocity_field,
              0.,
              1,
              boundaries.boundary_reflecting_periodic },
            boundaries.boundary_reflecting_periodic
        };
        
        auto const& part = ctrw.particles(0);
        double initial_velocity = getter_velocity_magnitude(part);
        for (std::size_t tt = 0; tt < times.size(); ++tt)
        {
          while (part.state_new().time < times[tt])
          {
            double velocity_magnitude = getter_velocity_magnitude(part);
            if (velocity_magnitude == 0.)
              break;
            ctrw.step(transitions);
          }
          velocity_autocorrelation[tt] +=
            initial_velocity*getter_velocity_magnitude_old(part);
          if (getter_velocity_magnitude(part) == 0.)
          {
            for (std::size_t tt2 = tt+1; tt2 < times.size(); ++tt2)
              --surviving_trajectories[tt2];
            break;
          }
        }
      }
      operation::div_InPlace(velocity_autocorrelation, surviving_trajectories);
      std::cout << nr_samples - surviving_trajectories.back()
                << " trajectories discarded due to zero velocity\n";
      
      std::stringstream stream;
      stream << std::scientific << std::setprecision(2);
      stream << jump_size_domains << "_"
             << measure_max << "_"
             << nr_measures << "_"
             << nr_samples << "_"
             << run_nr;
      std::string filename = output_dir + "/"
        + "velocity_magnitude_autocorrelation_time_lagrangian"
        + "_" + data_set + "_" + stream.str() + ".dat";
      std::ofstream output{ filename };
      if (!output.is_open())
        throw useful::open_write_error(filename);
      output << std::setprecision(12)
             << std::scientific;
      for (std::size_t tt = 0; tt < times.size(); ++tt)
      {
        output << times[tt] << "\t"
               << velocity_autocorrelation[tt] << "\n";
      }
      output.close();
      std::cout << "\tDone!\n";
      break;
    }
    case 9:
    {
      std::cout << "Computing space-Lagrangian velocity magnitude fluctuations autocorrelation...\n";
      std::vector<double> distances
        = range::linspace(0., measure_max*geometry.domain_side, nr_measures);
      
      std::cout << "\tSetting up trajectories...\n";
      auto particles = beadpack::make_particles_random_uniform_box<CTRW::Particle>(
        nr_samples, bead_pack, boundaries.boundary_periodic, geometry.boundaries, state_maker);
      std::cout << "\t\tDone!\n";
    
      double jump_size = jump_size_domains*geometry.domain_side;
      
      std::vector<double> velocity_autocorrelation(distances.size());
      std::vector<double> velocity_mean(particles.size());
      std::vector<std::size_t> surviving_trajectories(distances.size(), nr_samples);
      for (std::size_t pp = 0; pp < particles.size(); ++pp)
      {
        std::cout << "\tTrajectory " << pp+1 << " of " << particles.size() << "\n";
        std::vector<double> velocity_magnitude(distances.size());
        CTRW ctrw{ { particles[pp] }, CTRW::Tag{} };
        ctrw::Transitions_Position_VelocityStep transitions{
            jump_size,
            JumpGenerator{
              velocity_field,
              0.,
              1,
              boundaries.boundary_reflecting_periodic },
            boundaries.boundary_reflecting_periodic
        };
        
        auto const& part = ctrw.particles(0);
        double initial_velocity = getter_velocity_magnitude(part);
        double distance_traveled = 0.;
        for (std::size_t ss = 0; ss < distances.size(); ++ss)
        {
          while (distance_traveled < distances[ss])
          {
            double velocity_magnitude = getter_velocity_magnitude(part);
            if (velocity_magnitude == 0.)
              break;
            ctrw.step(transitions);
            double distance_increment = operation::abs(operation::minus(getter_position(part),
              getter_position_old(part)));
            distance_traveled += distance_increment;
            velocity_mean[pp] += velocity_magnitude*distance_increment;
          }
          velocity_magnitude[ss] = getter_velocity_magnitude_old(part);
          if (getter_velocity_magnitude(part) == 0.)
          {
            for (std::size_t ss2 = ss+1; ss2 < distances.size(); ++ss2)
              --surviving_trajectories[ss2];
            break;
          }
        }
        if (distance_traveled != 0.)
          velocity_mean[pp] /= distance_traveled;
        for (std::size_t ss = 0; ss < distances.size(); ++ss)
          velocity_autocorrelation[ss] += (initial_velocity-velocity_mean[pp])
            *(velocity_magnitude[ss]-velocity_mean[pp]);
      }
      operation::div_InPlace(velocity_autocorrelation, surviving_trajectories);
      std::cout << nr_samples - surviving_trajectories.back()
                << " trajectories discarded due to zero velocity\n";
      
      std::stringstream stream;
      stream << std::scientific << std::setprecision(2);
      stream << jump_size_domains << "_"
             << measure_max << "_"
             << nr_measures << "_"
             << nr_samples << "_"
             << run_nr;
      std::string filename = output_dir + "/"
        + "velocity_magnitude_fluctuations_autocorrelation_space_lagrangian"
        + "_" + data_set + "_" + stream.str() + ".dat";
      std::ofstream output{ filename };
      if (!output.is_open())
        throw useful::open_write_error(filename);
      output << std::setprecision(12)
             << std::scientific;
      for (std::size_t ss = 0; ss < distances.size(); ++ss)
      {
        output << distances[ss] << "\t"
               << velocity_autocorrelation[ss] << "\n";
      }
      output.close();
      std::cout << "\tDone!\n";
      break;
    }
    case 10:
    {
      std::cout << "Computing time-Lagrangian velocity magnitude fluctuations autocorrelation...\n";
      
      std::cout << "\tImporting mean velocity...\n";
      std::string mean_velocity_filename = input_dir + "/" + "mean_velocity.dat";
      std::vector<double> mean_velocity =
        beadpack::get_mean_velocity(geometry.dim, mean_velocity_filename);
      std::cout << "\t\tDone!\n";
      
      std::vector<double> times
        = range::linspace(0., measure_max*geometry.domain_side/operation::abs(mean_velocity),
                          nr_measures);
      
      std::cout << "\tSetting up trajectories...\n";
      auto particles = beadpack::make_particles_random_uniform_box<CTRW::Particle>(
        nr_samples, bead_pack, boundaries.boundary_periodic, geometry.boundaries, state_maker);
      std::cout << "\t\tDone!\n";
      
      double jump_size = jump_size_domains*geometry.domain_side;
      
      std::vector<double> velocity_autocorrelation(times.size());
      std::vector<double> velocity_mean(particles.size());
      std::vector<std::size_t> surviving_trajectories(times.size(), nr_samples);
      for (std::size_t pp = 0; pp < particles.size(); ++pp)
      {
        std::cout << "\tTrajectory " << pp+1 << " of " << particles.size() << "\n";
        std::vector<double> velocity_magnitude(times.size());
        CTRW ctrw{ { particles[pp] }, CTRW::Tag{} };
        ctrw::Transitions_Position_VelocityStep transitions{
            jump_size,
            JumpGenerator{
              velocity_field,
              0.,
              1,
              boundaries.boundary_reflecting_periodic },
            boundaries.boundary_reflecting_periodic
        };
        
        auto const& part = ctrw.particles(0);
        double initial_velocity = getter_velocity_magnitude(part);
        for (std::size_t tt = 0; tt < times.size(); ++tt)
        {
          while (part.state_new().time < times[tt])
          {
            double velocity_magnitude = getter_velocity_magnitude(part);
            if (velocity_magnitude == 0.)
              break;
            ctrw.step(transitions);
            velocity_mean[pp] += velocity_magnitude*transitions.time_step();
          }
          velocity_magnitude[tt] = getter_velocity_magnitude_old(part);
          if (getter_velocity_magnitude(part) == 0.)
          {
            for (std::size_t tt2 = tt+1; tt2 < times.size(); ++tt2)
              --surviving_trajectories[tt2];
            break;
          }
        }
        if (velocity_mean[pp] != 0.)
          velocity_mean[pp] /= std::min(times.back(), part.state_new().time);
        for (std::size_t tt = 0; tt < times.size(); ++tt)
          velocity_autocorrelation[tt] += (initial_velocity-velocity_mean[pp])
            *(velocity_magnitude[tt]-velocity_mean[pp]);
      }
      operation::div_InPlace(velocity_autocorrelation, surviving_trajectories);
      std::cout << nr_samples - surviving_trajectories.back()
                << " trajectories discarded due to zero velocity\n";
      
      std::stringstream stream;
      stream << std::scientific << std::setprecision(2);
      stream << jump_size_domains << "_"
             << measure_max << "_"
             << nr_measures << "_"
             << nr_samples << "_"
             << run_nr;
      std::string filename = output_dir + "/"
        + "velocity_magnitude_fluctuations_autocorrelation_time_lagrangian"
        + "_" + data_set + "_" + stream.str() + ".dat";
      std::ofstream output{ filename };
      if (!output.is_open())
        throw useful::open_write_error(filename);
      output << std::setprecision(12)
             << std::scientific;
      for (std::size_t tt = 0; tt < times.size(); ++tt)
      {
        output << times[tt] << "\t"
               << velocity_autocorrelation[tt] << "\n";
      }
      output.close();
      std::cout << "\tDone!\n";
      break;
    }
    case 11:
    {
      std::cout << "Computing space-Lagrangian velocity magnitude series...\n";
      std::vector<double> distances
        = range::linspace(0., measure_max*geometry.domain_side, nr_measures);
      
      std::cout << "\tSetting up trajectories...\n";
      auto particles = beadpack::make_particles_random_uniform_box<CTRW::Particle>(
        nr_samples, bead_pack, boundaries.boundary_periodic, geometry.boundaries, state_maker);
      std::cout << "\t\tDone!\n";
      
      double jump_size = jump_size_domains*geometry.domain_side;

      std::stringstream stream;
      stream << std::scientific << std::setprecision(2);
      stream << jump_size_domains << "_"
             << measure_max << "_"
             << nr_measures << "_"
             << nr_samples << "_"
             << run_nr;
      std::string filename = output_dir + "/"
        + "velocity_magnitude_series_space_lagrangian"
        + "_" + data_set + "_" + stream.str() + ".dat";
      std::ofstream output{ filename };
      if (!output.is_open())
        throw useful::open_write_error(filename);
      output << std::setprecision(12)
             << std::scientific;
      
      std::size_t discarded = 0;
      useful::print(output, distances);
      output << "\n";
      for (std::size_t pp = 0; pp < particles.size(); ++pp)
      {
        std::cout << "\tTrajectory " << pp+1 << " of " << particles.size() << "\n";
        CTRW ctrw{ { particles[pp] }, CTRW::Tag{} };
        ctrw::Transitions_Position_VelocityStep transitions{
            jump_size,
            JumpGenerator{
              velocity_field,
              0.,
              1,
              boundaries.boundary_reflecting_periodic },
            boundaries.boundary_reflecting_periodic
        };
        
        auto const& part = ctrw.particles(0);
        double distance_traveled = 0.;
        std::vector<double> velocity_series(distances.size());
        for (std::size_t ss = 0; ss < distances.size(); ++ss)
        {
          while (distance_traveled < distances[ss])
          {
            double velocity_magnitude = getter_velocity_magnitude(part);
            if (velocity_magnitude == 0.)
              break;
            ctrw.step(transitions);
            distance_traveled +=
              operation::abs(operation::minus(getter_position(part),
                                              getter_position_old(part)));
          }
          velocity_series[ss] = getter_velocity_magnitude_old(part);
          if (getter_velocity_magnitude(part) == 0.)
          {
            ++discarded;
            break;
          }
        }
        useful::print(output, velocity_series);
        output << "\n";
      }
      std::cout << discarded << " trajectories discarded due to zero velocity\n";
      
      output.close();
      std::cout << "\tDone!\n";
      break;
    }
    case 12:
    {
      std::cout << "Computing time-Lagrangian velocity magnitude series...\n";
      
      std::cout << "\tImporting mean velocity...\n";
      std::string mean_velocity_filename = input_dir + "/" + "mean_velocity.dat";
      std::vector<double> mean_velocity =
        beadpack::get_mean_velocity(geometry.dim, mean_velocity_filename);
      std::vector<double> times
        = range::linspace(0., measure_max*geometry.domain_side/operation::abs(mean_velocity),
                          nr_measures);
      std::cout << "\t\tDone!\n";
      
      std::cout << "\tSetting up trajectories...\n";
      auto particles = beadpack::make_particles_random_uniform_box<CTRW::Particle>(
        nr_samples, bead_pack, boundaries.boundary_periodic, geometry.boundaries, state_maker);
      std::cout << "\t\tDone!\n";
      
      double jump_size = jump_size_domains*geometry.domain_side;
      
      std::stringstream stream;
      stream << std::scientific << std::setprecision(2);
      stream << jump_size_domains << "_"
             << measure_max << "_"
             << nr_measures << "_"
             << nr_samples << "_"
             << run_nr;
      std::string filename = output_dir + "/"
        + "velocity_magnitude_series_time_lagrangian"
        + "_" + data_set + "_" + stream.str() + ".dat";
      std::ofstream output{ filename };
      if (!output.is_open())
        throw useful::open_write_error(filename);
      output << std::setprecision(12)
             << std::scientific;
      
      useful::print(output, times);
      output << "\n";
      std::size_t discarded = 0;
      for (std::size_t pp = 0; pp < particles.size(); ++pp)
      {
        std::cout << "\tTrajectory " << pp+1 << " of " << particles.size() << "\n";
        CTRW ctrw{ { particles[pp] }, CTRW::Tag{} };
        ctrw::Transitions_Position_VelocityStep transitions{
            jump_size,
            JumpGenerator{
              velocity_field,
              0.,
              1,
              boundaries.boundary_reflecting_periodic },
            boundaries.boundary_reflecting_periodic
        };
        
        auto const& part = ctrw.particles(0);
        std::vector<double> velocity_series(times.size());
        for (std::size_t tt = 0; tt < times.size(); ++tt)
        {
          while (part.state_new().time < times[tt])
          {
            double velocity_magnitude = getter_velocity_magnitude(part);
            if (velocity_magnitude == 0.)
              break;
            ctrw.step(transitions);
          }
          velocity_series[tt] = getter_velocity_magnitude_old(part);
          if (getter_velocity_magnitude(part) == 0.)
          {
            ++discarded;
            break;
          }
        }
        useful::print(output, velocity_series);
        output << "\n";
      }
      std::cout << discarded << " trajectories discarded due to zero velocity\n";
      output.close();
      std::cout << "\tDone!\n";
      break;
    }
      
    default:
      throw useful::bad_parameters();
  }

  return 0;
}
