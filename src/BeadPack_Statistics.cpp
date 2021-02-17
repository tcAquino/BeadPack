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
  using namespace beadpack::model_bcc_cartesian;
  
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
              << "(if measure_type >= 5) jump_size_domains : Jump size in units of domain side\n"
              << "(if measure_type >= 5) velocity_tolerance : minimum velocity to discard trajectory\n"
              << "                                         in units of mean velocity\n"
              << "(if measure_type >= 5) measure_max : Maximum distance along\n"
              << "                                     streamlines in units of domain side\n"
                 "                                     or maximum time in units of\n"
              << "                                     advection time\n"
              << "(if measure_type >= 7) nr_measures : Number of measurements along streamlines\n"
              << "input_dir_base : Path to look for input data [../input]\n"
              << "output_dir : Path folder to output to [../output]\n";
    return 0;
  }
  
  if (argc < 5)
    throw useful::bad_parameters();
  
  using State = ctrw::State_periodic<std::vector<double>,
    std::vector<int>, useful::Empty, double, std::size_t>;
  using CTRW = ctrw::CTRW<State>;
  
  std::cout << std::setprecision(2)
            << std::scientific;
  
  std::size_t arg = 1;
  std::size_t nr_samples = strtoul(argv[arg++], NULL, 0);
  int measure_type = atoi(argv[arg++]);
  std::size_t run_nr = strtoul(argv[arg++], NULL, 0);
  std::string data_set = argv[arg++];
  
  double jump_size_domains = 0.;
  double measure_max = 0.;
  double velocity_tolerance = 0.;
  if (measure_type >= 5)
  {
    if (argc < 7)
      throw useful::bad_parameters();
    jump_size_domains = atof(argv[arg++]);
    velocity_tolerance = atof(argv[arg++]);
    measure_max = atof(argv[arg++]);
  }
  
  std::size_t nr_measures = 0;
  if (measure_type >= 7)
  {
    if (argc < 8)
      throw useful::bad_parameters();
    nr_measures = strtoul(argv[arg++], NULL, 0);
  }
  
  std::string input_dir_base = argc > arg ? argv[arg++] : "../input";
  std::string output_dir = argc > arg ? argv[arg++] : "../output";
  
  std::string input_dir = input_dir_base + "/" + data_set;
  
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

  // Lightweigth objects for multiple statistics,
  // not always needed but declared here to avoid repetition
  double jump_size = jump_size_domains*geometry.domain_side;
  ctrw::Transitions_Position_VelocityStep transitions{
      jump_size,
      ctrw::JumpGenerator_Velocity_withHint_RK4{
        velocity_field,
        0., 1,
        boundaries.boundary_reflecting_periodic },
      boundaries.boundary_reflecting_periodic
  };

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
  auto getter_velocity_magnitude = [&getter_velocity]
  (CTRW::Particle const& particle)
  { return operation::abs(getter_velocity(particle)); };
  
  std::stringstream stream_samples;
  stream_samples << nr_samples << "_" << run_nr;
  std::string params_samples = stream_samples.str();
  
  std::stringstream stream_lagrangian_no_measures;
  stream_lagrangian_no_measures << std::scientific << std::setprecision(2);
  stream_lagrangian_no_measures << jump_size_domains << "_"
                                << velocity_tolerance << "_"
                                << measure_max << "_"
                                << nr_samples << "_"
                                << run_nr;
  std::string params_lagrangian_no_measures = stream_lagrangian_no_measures.str();
  
  std::stringstream stream_lagrangian_measures;
  stream_lagrangian_measures << std::scientific << std::setprecision(2);
  stream_lagrangian_measures << jump_size_domains << "_"
                             << velocity_tolerance << "_"
                             << measure_max << "_"
                             << nr_measures << "_"
                             << nr_samples << "_"
                             << run_nr;
  std::string params_lagrangian_measures = stream_lagrangian_measures.str();
  
  // Choose statistics to compute
  switch (measure_type)
  {
    case 0:
    {
      std::cout << "Computing Eulerian mean velocity vector...\n";
      
      std::cout << "\tSetting up trajectories...\n";
      auto particles = beadpack::make_particles_random_uniform_box<CTRW::Particle>(
        nr_samples, bead_pack, boundaries.boundary_periodic, geometry.boundaries, state_maker);
      std::cout << "\t\tDone!\n";
      
      std::cout << "\tComputing velocity samples...\n";
      std::vector<double> mean_velocity(geometry.dim);
      std::size_t sample = 0;
      for (auto const& part : particles)
      {
        ++sample;
        if (sample % 10000 == 0)
          std::cout << "\t\tSample " << sample << " of " << nr_samples;
        operation::plus_InPlace(mean_velocity, getter_velocity(part));
      }
      operation::div_scalar_InPlace(mean_velocity, particles.size());
      std::cout << "\t\tDone!\n";
      
      auto output = useful::open_write(output_dir + "/" + "mean_velocity"
        + "_" + data_set + "_" + params_samples + ".dat");
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
      
      std::cout << "\tSetting up trajectories...\n";
      auto particles = beadpack::make_particles_random_uniform_box<CTRW::Particle>(
        nr_samples, bead_pack, boundaries.boundary_periodic, geometry.boundaries, state_maker);
      std::cout << "\t\tDone!\n";
      
      std::cout << "\tComputing velocity magnitude samples...\n";
      double mean_velocity_magnitude = 0.;
      std::size_t sample = 0;
      for (auto const& part : particles)
      {
        ++sample;
        if (sample % 10000 == 0)
          std::cout << "\t\tSample " << sample << " of " << nr_samples << "\n";
        mean_velocity_magnitude += getter_velocity_magnitude(part);
      }
      std::cout << "\t\tDone!\n";
      
      mean_velocity_magnitude /= particles.size();
      auto output = useful::open_write(output_dir + "/"
        + "mean_velocity_magnitude" + "_"
        + data_set + "_" + params_samples + ".dat");
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
      auto output = useful::open_write(output_dir + "/" + "porosity" + "_" + data_set + "_"
        + params_samples + ".dat");
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
      std::vector<double> mean_velocity =
        beadpack::get_mean_velocity(geometry.dim, input_dir + "/" + "mean_velocity.dat");
      std::cout << "\t\tDone!\n";
      
      std::cout << "\tImporting mean velocity magnitude...\n";
      auto input = useful::open_read(input_dir + "/" + "mean_velocity_magnitude.dat");
      double mean_velocity_magnitude;
      input >> mean_velocity_magnitude;
      input.close();
      std::cout << "\t\tDone!\n";
      
      std::cout << "\tSetting up trajectories...\n";
      auto particles = beadpack::make_particles_random_uniform_box<CTRW::Particle>(
        nr_samples, bead_pack, boundaries.boundary_periodic, geometry.boundaries, state_maker);
      std::cout << "\t\tDone!\n";
      
      auto getter_velocity_downstream = ctrw::Get_old_from_particle{
         ctrw::Get_projection{
           ctrw::Get_position_property{ velocity_field },
           mean_velocity } };
      
      std::cout << "\tComputing downstream velocity samples...\n";
      double mean_velocity_downstream = 0.;
      std::size_t sample = 0;
      for (auto const& part : particles)
      {
        ++sample;
        if (sample % 10000 == 0)
          std::cout << "\t\tSample " << sample << " of " << nr_samples << "\n";
        mean_velocity_downstream += getter_velocity_downstream(part);
      }
      mean_velocity_downstream /= particles.size();
      std::cout << "\t\tDone!\n";
      
      double tortuosity = mean_velocity_magnitude/mean_velocity_downstream;
      auto output = useful::open_write(output_dir + "/" + "tortuosity"
        + "_" + data_set + "_" + params_samples + ".dat");
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

      std::cout << "\tSetting up trajectories...\n";
      auto particles = beadpack::make_particles_random_uniform_box<ctrw::Particle<State>>(
        nr_samples, bead_pack, boundaries.boundary_periodic, geometry.boundaries, state_maker);
      std::cout << "\t\tDone!\n";

      auto output = useful::open_write(output_dir + "/"
        + "velocity_magnitude_samples_uniform_unit_cell"
        + "_" + data_set + "_" + params_samples + ".dat");
      output << std::setprecision(12)
             << std::scientific;
      
      std::cout << "\tComputing velocity magnitude samples...\n";
      std::size_t sample = 0;
      for (auto const& part : particles)
      {
        ++sample;
        if (sample % 10000 == 0)
          std::cout << "\t\tSample " << sample << " of " << nr_samples << "\n";
        useful::print(output, getter_velocity_magnitude(part));
        output << "\n";
      }
      std::cout << "\t\tDone!\n";
      
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
      
      std::cout << "\tImporting mean velocity magnitude...\n";
      auto input = useful::open_read(input_dir + "/" + "mean_velocity_magnitude.dat");
      double mean_velocity_magnitude;
      input >> mean_velocity_magnitude;
      input.close();
      std::cout << "\t\tDone!\n";
      
      double distance_max = measure_max*geometry.domain_side;
      
      std::vector<double> velocity_mean(particles.size());
      std::size_t surviving_trajectories = nr_samples;
      double velocity_cutoff = velocity_tolerance*mean_velocity_magnitude;
      for (std::size_t pp = 0; pp < particles.size(); ++pp)
      {
        std::cout << "\tTrajectory " << pp+1 << " of " << particles.size() << "\n";
        CTRW ctrw{ { particles[pp] }, CTRW::Tag{} };
        
        auto const& part = ctrw.particles(0);
        double distance_traveled = 0.;
        double initial_velocity = getter_velocity_magnitude(part);
        double velocity = initial_velocity;
        while (distance_traveled < distance_max)
        {
          if (velocity < velocity_cutoff)
            break;
          ctrw.step(transitions);
          double distance_increment = operation::abs(operation::minus(getter_position(part),
            getter_position_old(part)));
          distance_traveled += distance_increment;
          velocity = distance_increment/transitions.time_step();
          velocity_mean[pp] += velocity*distance_increment;
        }
        if (velocity < velocity_cutoff)
          --surviving_trajectories;
        if (distance_traveled != 0.)
          velocity_mean[pp] /= distance_traveled;
      }
      std::cout << "\t" << nr_samples - surviving_trajectories
                << " trajectories discarded\n";
      
      
      auto output = useful::open_write(output_dir + "/"
        + "velocity_magnitude_mean_space_lagrangian"
        + "_" + data_set + "_" + params_lagrangian_no_measures + ".dat");
      output << std::setprecision(12)
             << std::scientific;
      useful::print(output, velocity_mean);
      output << "\t" << surviving_trajectories << "\n";
      output.close();
      std::cout << "\tDone!\n";
      break;
    }
    case 6:
    {
      std::cout << "Computing time-Lagrangian velocity magnitude mean...\n";
      
      std::cout << "\tImporting mean velocity...\n";
      std::vector<double> mean_velocity =
        beadpack::get_mean_velocity(geometry.dim, input_dir + "/" + "mean_velocity.dat");
      std::cout << "\t\tDone!\n";
      
      std::cout << "\tImporting mean velocity magnitude...\n";
      auto input = useful::open_read(input_dir + "/" + "mean_velocity_magnitude.dat");
      double mean_velocity_magnitude;
      input >> mean_velocity_magnitude;
      input.close();
      std::cout << "\t\tDone!\n";
      
      std::cout << "\tSetting up trajectories...\n";
      auto particles = beadpack::make_particles_random_uniform_box<CTRW::Particle>(
        nr_samples, bead_pack, boundaries.boundary_periodic, geometry.boundaries, state_maker);
      std::cout << "\t\tDone!\n";
      
      double time_max = measure_max*geometry.domain_side/operation::abs(mean_velocity);
      
      std::vector<double> velocity_mean(particles.size());
      std::size_t surviving_trajectories = nr_samples;
      double velocity_cutoff = velocity_tolerance*mean_velocity_magnitude;
      for (std::size_t pp = 0; pp < particles.size(); ++pp)
      {
        std::cout << "\tTrajectory " << pp+1 << " of " << particles.size() << "\n";
        CTRW ctrw{ { particles[pp] }, CTRW::Tag{} };
        
        auto const& part = ctrw.particles(0);
        double velocity = getter_velocity_magnitude(part);
        while (part.state_new().time < time_max)
        {
          if (velocity < velocity_cutoff)
            break;
          ctrw.step(transitions);
          double distance_increment = operation::abs(operation::minus(getter_position(part),
                                                                      getter_position_old(part)));
          velocity = distance_increment/transitions.time_step();
          velocity_mean[pp] += distance_increment;
        }
        if (velocity < velocity_cutoff)
          --surviving_trajectories;
        if (velocity_mean[pp] != 0.)
          velocity_mean[pp] /= std::min(part.state_new().time, time_max);
      }
      std::cout << "\t" << nr_samples - surviving_trajectories
                << " trajectories discarded\n";
      
      auto output = useful::open_write(output_dir + "/"
        + "velocity_magnitude_mean_time_lagrangian"
        + "_" + data_set + "_" + params_lagrangian_no_measures + ".dat");
      output << std::setprecision(12)
             << std::scientific;
      useful::print(output, velocity_mean);
      output << "\t" << surviving_trajectories << "\n";
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
      
      std::cout << "\tImporting mean velocity magnitude...\n";;
      auto input = useful::open_read(input_dir + "/" + "mean_velocity_magnitude.dat");
      double mean_velocity_magnitude;
      input >> mean_velocity_magnitude;
      input.close();
      std::cout << "\t\tDone!\n";
      
      std::vector<double> velocity_autocorrelation(distances.size());
      std::vector<std::size_t> surviving_trajectories(distances.size(), nr_samples);
      double velocity_cutoff = velocity_tolerance*mean_velocity_magnitude;
      for (std::size_t pp = 0; pp < particles.size(); ++pp)
      {
        std::cout << "\tTrajectory " << pp+1 << " of " << particles.size() << "\n";
        CTRW ctrw{ { particles[pp] }, CTRW::Tag{} };
        
        auto const& part = ctrw.particles(0);
        double initial_velocity = getter_velocity_magnitude(part);
        double velocity = initial_velocity;
        double distance_traveled = 0.;
        for (std::size_t ss = 0; ss < distances.size(); ++ss)
        {
          while (distance_traveled < distances[ss])
          {
            if (velocity < velocity_cutoff)
              break;
            ctrw.step(transitions);
            distance_traveled +=
              operation::abs(operation::minus(getter_position(part),
                                              getter_position_old(part)));
          }
          if (velocity < velocity_cutoff)
          {
            for (std::size_t ss2 = ss+1; ss2 < distances.size(); ++ss2)
              --surviving_trajectories[ss2];
            break;
          }
          velocity_autocorrelation[ss] += initial_velocity*velocity;
        }
      }
      operation::div_InPlace(velocity_autocorrelation, surviving_trajectories);
      std::cout << "\t" << nr_samples - surviving_trajectories.back()
                << " trajectories discarded\n";
      
      auto output = useful::open_write(output_dir + "/"
        + "velocity_magnitude_autocorrelation_space_lagrangian"
        + "_" + data_set + "_" + params_lagrangian_measures + ".dat");
      output << std::setprecision(12)
             << std::scientific;
      for (std::size_t ss = 0; ss < distances.size(); ++ss)
        output << distances[ss] << "\t"
               << velocity_autocorrelation[ss] << "\t"
               << surviving_trajectories[ss] << "\n";
      output.close();
      std::cout << "\tDone!\n";
      break;
    }
    case 8:
    {
      std::cout << "Computing time-Lagrangian velocity magnitude autocorrelation...\n";
      
      std::cout << "\tImporting mean velocity...\n";
      std::vector<double> mean_velocity =
        beadpack::get_mean_velocity(geometry.dim, input_dir + "/" + "mean_velocity.dat");
      std::cout << "\t\tDone!\n";
      
      std::cout << "\tImporting mean velocity magnitude...\n";
      auto input = useful::open_read(input_dir + "/" + "mean_velocity_magnitude.dat");
      double mean_velocity_magnitude;
      input >> mean_velocity_magnitude;
      input.close();
      std::cout << "\t\tDone!\n";
      
      std::vector<double> times
        = range::linspace(0., measure_max*geometry.domain_side/operation::abs(mean_velocity),
                          nr_measures);
      
      std::cout << "\tSetting up trajectories...\n";
      auto particles = beadpack::make_particles_random_uniform_box<CTRW::Particle>(
        nr_samples, bead_pack, boundaries.boundary_periodic, geometry.boundaries, state_maker);
      std::cout << "\t\tDone!\n";
      
      std::vector<double> velocity_autocorrelation(times.size());
      std::vector<std::size_t> surviving_trajectories(times.size(), nr_samples);
      double velocity_cutoff = velocity_tolerance*mean_velocity_magnitude;
      for (std::size_t pp = 0; pp < particles.size(); ++pp)
      {
        std::cout << "\tTrajectory " << pp+1 << " of " << particles.size() << "\n";
        CTRW ctrw{ { particles[pp] }, CTRW::Tag{} };
        
        auto const& part = ctrw.particles(0);
        double initial_velocity = getter_velocity_magnitude(part);
        double velocity = initial_velocity;
        for (std::size_t tt = 0; tt < times.size(); ++tt)
        {
          while (part.state_new().time < times[tt])
          {
            if (velocity < velocity_cutoff)
              break;
            ctrw.step(transitions);
            double distance_increment = operation::abs(operation::minus(getter_position(part),
                                                                        getter_position_old(part)));
            velocity = distance_increment/transitions.time_step();
          }
          if (velocity < velocity_cutoff)
          {
            for (std::size_t tt2 = tt+1; tt2 < times.size(); ++tt2)
              --surviving_trajectories[tt2];
            break;
          }
          velocity_autocorrelation[tt] += initial_velocity*velocity;
        }
      }
      operation::div_InPlace(velocity_autocorrelation, surviving_trajectories);
      std::cout << "\t" << nr_samples - surviving_trajectories.back()
                << " trajectories discarded\n";
      
      auto output = useful::open_write(output_dir + "/"
        + "velocity_magnitude_autocorrelation_time_lagrangian"
        + "_" + data_set + "_" + params_lagrangian_measures + ".dat");
      output << std::setprecision(12)
             << std::scientific;
      for (std::size_t tt = 0; tt < times.size(); ++tt)
        output << times[tt] << "\t"
               << velocity_autocorrelation[tt] << "\t"
               << surviving_trajectories[tt] << "\n";
      output.close();
      std::cout << "\tDone!\n";
      break;
    }
    case 9:
    {
      std::cout << "Computing space-Lagrangian velocity magnitude fluctuations autocorrelation...\n";
      std::vector<double> distances
        = range::linspace(0., measure_max*geometry.domain_side, nr_measures);
      
      std::cout << "\tImporting mean velocity magnitude...\n";
      auto input = useful::open_read(input_dir + "/" + "mean_velocity_magnitude.dat");
      double mean_velocity_magnitude;
      input >> mean_velocity_magnitude;
      input.close();
      std::cout << "\t\tDone!\n";

      std::cout << "\tSetting up trajectories...\n";
      auto particles = beadpack::make_particles_random_uniform_box<CTRW::Particle>(
        nr_samples, bead_pack, boundaries.boundary_periodic, geometry.boundaries, state_maker);
      std::cout << "\t\tDone!\n";
      
      std::vector<double> velocity_autocorrelation(distances.size());
      std::vector<double> velocity_mean(particles.size());
      std::vector<std::size_t> surviving_trajectories(distances.size(), nr_samples);
      double velocity_cutoff = velocity_tolerance*mean_velocity_magnitude;
      for (std::size_t pp = 0; pp < particles.size(); ++pp)
      {
        std::cout << "\tTrajectory " << pp+1 << " of " << particles.size() << "\n";
        std::vector<double> velocity_magnitude(distances.size());
        CTRW ctrw{ { particles[pp] }, CTRW::Tag{} };
        
        auto const& part = ctrw.particles(0);
        double initial_velocity = getter_velocity_magnitude(part);
        double velocity = initial_velocity;
        double distance_traveled = 0.;
        std::size_t ss_attained = distances.size();
        for (std::size_t ss = 0; ss < distances.size(); ++ss)
        {
          while (distance_traveled < distances[ss])
          {
            if (velocity < velocity_cutoff)
              break;
            ctrw.step(transitions);
            double distance_increment = operation::abs(operation::minus(getter_position(part),
              getter_position_old(part)));
            velocity = distance_increment/transitions.time_step();
            if (velocity < velocity_cutoff)
              break;
            distance_traveled += distance_increment;
            velocity_mean[pp] += velocity*distance_increment;
          }
          if (velocity < velocity_cutoff)
          {
            ss_attained = ss;
            for (std::size_t ss2 = ss_attained; ss2 < distances.size(); ++ss2)
              --surviving_trajectories[ss2];
            break;
          }
          velocity_magnitude[ss] = velocity;
        }
        if (distance_traveled != 0.)
          velocity_mean[pp] /= distance_traveled;
        for (std::size_t ss = 0; ss < ss_attained; ++ss)
          velocity_autocorrelation[ss] += (initial_velocity-velocity_mean[pp])
            *(velocity_magnitude[ss]-velocity_mean[pp]);
      }
      operation::div_InPlace(velocity_autocorrelation, surviving_trajectories);
      std::cout << "\t" << nr_samples - surviving_trajectories.back()
                << " trajectories discarded\n";

      auto output = useful::open_write(output_dir + "/"
        + "velocity_magnitude_fluctuations_autocorrelation_space_lagrangian"
        + "_" + data_set + "_" + params_lagrangian_measures + ".dat");
      output << std::setprecision(12)
             << std::scientific;
      for (std::size_t ss = 0; ss < distances.size(); ++ss)
        output << distances[ss] << "\t"
               << velocity_autocorrelation[ss] << "\t"
               << surviving_trajectories[ss] << "\n";
      output.close();
      std::cout << "\tDone!\n";
      break;
    }
    case 10:
    {
      std::cout << "Computing time-Lagrangian velocity magnitude fluctuations autocorrelation...\n";
      
      std::cout << "\tImporting mean velocity...\n";
      std::vector<double> mean_velocity =
        beadpack::get_mean_velocity(geometry.dim, input_dir + "/" + "mean_velocity.dat");
      std::cout << "\t\tDone!\n";
      
      std::cout << "\tImporting mean velocity magnitude...\n";
      auto input = useful::open_read(input_dir + "/" + "mean_velocity_magnitude.dat");
      double mean_velocity_magnitude;
      input >> mean_velocity_magnitude;
      input.close();
      std::cout << "\t\tDone!\n";
      
      std::vector<double> times
        = range::linspace(0., measure_max*geometry.domain_side/operation::abs(mean_velocity),
                          nr_measures);
      
      std::cout << "\tSetting up trajectories...\n";
      auto particles = beadpack::make_particles_random_uniform_box<CTRW::Particle>(
        nr_samples, bead_pack, boundaries.boundary_periodic, geometry.boundaries, state_maker);
      std::cout << "\t\tDone!\n";
      
      std::vector<double> velocity_autocorrelation(times.size());
      std::vector<double> velocity_mean(particles.size());
      std::vector<std::size_t> surviving_trajectories(times.size(), nr_samples);
      double velocity_cutoff = velocity_tolerance*mean_velocity_magnitude;
      for (std::size_t pp = 0; pp < particles.size(); ++pp)
      {
        std::cout << "\tTrajectory " << pp+1 << " of " << particles.size() << "\n";
        std::vector<double> velocity_magnitude(times.size());
        CTRW ctrw{ { particles[pp] }, CTRW::Tag{} };
        
        auto const& part = ctrw.particles(0);
        double initial_velocity = getter_velocity_magnitude(part);
        double velocity = initial_velocity;
        std::size_t tt_attained = times.size()-1;
        for (std::size_t tt = 0; tt < times.size(); ++tt)
        {
          while (part.state_new().time < times[tt])
          {
            if (velocity < velocity_cutoff)
              break;
            ctrw.step(transitions);
            double distance_increment = operation::abs(operation::minus(getter_position(part),
                                                                        getter_position_old(part)));
            velocity = distance_increment/transitions.time_step();
            velocity_mean[pp] += distance_increment;
          }
          if (velocity < velocity_cutoff)
          {
            tt_attained = tt;
            for (std::size_t tt2 = tt_attained+1; tt2 < times.size(); ++tt2)
              --surviving_trajectories[tt2];
            break;
          }
          velocity_magnitude[tt] = velocity;
        }
        if (velocity_mean[pp] != 0.)
          velocity_mean[pp] /= std::min(times.back(), part.state_new().time);
        for (std::size_t tt = 0; tt <= tt_attained; ++tt)
          velocity_autocorrelation[tt] += (initial_velocity-velocity_mean[pp])
            *(velocity_magnitude[tt]-velocity_mean[pp]);
      }
      operation::div_InPlace(velocity_autocorrelation, surviving_trajectories);
      std::cout << "\t" << nr_samples - surviving_trajectories.back()
                << " trajectories discarded\n";

      auto output = useful::open_write(output_dir + "/"
        + "velocity_magnitude_fluctuations_autocorrelation_time_lagrangian"
        + "_" + data_set + "_" + params_lagrangian_measures + ".dat");
      output << std::setprecision(12)
             << std::scientific;
      for (std::size_t tt = 0; tt < times.size(); ++tt)
        output << times[tt] << "\t"
               << velocity_autocorrelation[tt] << "\t"
               << surviving_trajectories[tt] << "\n";
      output.close();
      std::cout << "\tDone!\n";
      break;
    }
    case 11:
    {
      std::cout << "Computing space-Lagrangian velocity magnitude series...\n";
      std::vector<double> distances
        = range::linspace(0., measure_max*geometry.domain_side, nr_measures);
      
      std::cout << "\tImporting mean velocity magnitude...\n";
      auto input = useful::open_read(input_dir + "/" + "mean_velocity_magnitude.dat");
      double mean_velocity_magnitude;
      input >> mean_velocity_magnitude;
      input.close();
      std::cout << "\t\tDone!\n";
      
      std::cout << "\tSetting up trajectories...\n";
      auto particles = beadpack::make_particles_random_uniform_box<CTRW::Particle>(
        nr_samples, bead_pack, boundaries.boundary_periodic, geometry.boundaries, state_maker);
      std::cout << "\t\tDone!\n";

      auto output = useful::open_write(output_dir + "/"
        + "velocity_magnitude_series_space_lagrangian"
        + "_" + data_set + "_" + params_lagrangian_measures + ".dat");
      output << std::setprecision(12)
             << std::scientific;
      
      std::size_t surviving_trajectories = nr_samples;
      useful::print(output, distances);
      output << "\n";
      double velocity_cutoff = velocity_tolerance*mean_velocity_magnitude;
      for (std::size_t pp = 0; pp < particles.size(); ++pp)
      {
        std::cout << "\tTrajectory " << pp+1 << " of " << particles.size() << "\n";
        CTRW ctrw{ { particles[pp] }, CTRW::Tag{} };
        
        auto const& part = ctrw.particles(0);
        double distance_traveled = 0.;
        std::vector<double> velocity_series(distances.size());
        double initial_velocity = getter_velocity_magnitude(part);
        double velocity = initial_velocity;
        for (std::size_t ss = 0; ss < distances.size(); ++ss)
        {
          while (distance_traveled < distances[ss])
          {
            if (velocity < velocity_cutoff)
              break;
            ctrw.step(transitions);
            double distance_increment = operation::abs(operation::minus(getter_position(part),
                                                                        getter_position_old(part)));
            distance_traveled += distance_increment;
            velocity = distance_increment/transitions.time_step();
          }
          if (velocity < velocity_cutoff)
          {
            --surviving_trajectories;
            break;
          }
          velocity_series[ss] = velocity;
        }
        useful::print(output, velocity_series);
        output << "\n";
      }
      std::cout << "\t" << nr_samples - surviving_trajectories
                << " trajectories discarded\n";
      output.close();
      std::cout << "\tDone!\n";
      break;
    }
    case 12:
    {
      std::cout << "Computing time-Lagrangian velocity magnitude series...\n";
      
      std::cout << "\tImporting mean velocity...\n";
      std::vector<double> mean_velocity =
        beadpack::get_mean_velocity(geometry.dim, input_dir + "/" + "mean_velocity.dat");
      std::vector<double> times
        = range::linspace(0., measure_max*geometry.domain_side/operation::abs(mean_velocity),
                          nr_measures);
      std::cout << "\t\tDone!\n";
      
      std::cout << "\tImporting mean velocity magnitude...\n";
      auto input = useful::open_read(input_dir + "/" + "mean_velocity_magnitude.dat");
      double mean_velocity_magnitude;
      input >> mean_velocity_magnitude;
      input.close();
      std::cout << "\t\tDone!\n";
      
      std::cout << "\tSetting up trajectories...\n";
      auto particles = beadpack::make_particles_random_uniform_box<CTRW::Particle>(
        nr_samples, bead_pack, boundaries.boundary_periodic, geometry.boundaries, state_maker);
      std::cout << "\t\tDone!\n";
      
      auto output = useful::open_write(output_dir + "/"
        + "velocity_magnitude_series_time_lagrangian"
        + "_" + data_set + "_" + params_lagrangian_measures + ".dat");
      output << std::setprecision(12)
             << std::scientific;
      
      useful::print(output, times);
      output << "\n";
      std::size_t surviving_trajectories = nr_samples;
      double velocity_cutoff = velocity_tolerance*mean_velocity_magnitude;
      for (std::size_t pp = 0; pp < particles.size(); ++pp)
      {
        std::cout << "\tTrajectory " << pp+1 << " of " << particles.size() << "\n";
        CTRW ctrw{ { particles[pp] }, CTRW::Tag{} };
        
        auto const& part = ctrw.particles(0);
        std::vector<double> velocity_series(times.size());
        double initial_velocity = getter_velocity_magnitude(part);
        double velocity = initial_velocity;
        for (std::size_t tt = 0; tt < times.size(); ++tt)
        {
          while (part.state_new().time < times[tt])
          {
            if (velocity < velocity_cutoff)
              break;
            ctrw.step(transitions);
            double distance_increment = operation::abs(operation::minus(getter_position(part),
                                                                        getter_position_old(part)));
            velocity = distance_increment/transitions.time_step();
          }
          if (velocity < velocity_cutoff)
          {
            --surviving_trajectories;
            break;
          }
          velocity_series[tt] = velocity;
        }
        useful::print(output, velocity_series);
        output << "\n";
      }
      std::cout << "\t" << nr_samples - surviving_trajectories
                << " trajectories discarded\n";
      output.close();
      std::cout << "\tDone!\n";
      break;
    }
    default:
      throw useful::bad_parameters();
  }

  return 0;
}
