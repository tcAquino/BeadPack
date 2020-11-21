//
//  BeadPack_InitialConditions.h
//  BeadPack
//
//  Created by Tomás Aquino on 14/05/2020.
//  Copyright © 2020 Tomás Aquino. All rights reserved.
//

#ifndef BeadPack_InitialConditions_h
#define BeadPack_InitialConditions_h

#include <fstream>
#include <cmath>
#include <random>
#include <utility>
#include <vector>
#include "Algebra/Algebra.h"
#include "general/Constants.h"
#include "Stochastic/Random.h"
#include "Stochastic/CTRW/Boundary.h"

namespace beadpack
{
  template
  <typename Particle, typename BeadPack, typename Boundary, typename StateMaker>
  auto make_particles_random_uniform_box
  (std::size_t nr_particles, BeadPack const& bead_pack, Boundary const& boundary,
   std::vector<std::pair<double, double>> const& boundaries, StateMaker state_maker)
  {
    std::vector<Particle> particles;
    particles.reserve(nr_particles);
    
    std::mt19937 rng{ std::random_device{}() };
    std::vector<std::uniform_real_distribution<double>> dist;
    dist.reserve(BeadPack::dim);
    for (std::size_t dd = 0; dd < BeadPack::dim; ++dd)
      dist.emplace_back(boundaries[dd].first, boundaries[dd].second);
    while (particles.size() < nr_particles)
    {
      auto state{ state_maker() };
      for (std::size_t dd = 0; dd < BeadPack::dim; ++dd)
        state.position[dd] = dist[dd](rng);
      boundary(state);
      if (!bead_pack.inside(state.position).first)
        particles.push_back(state);
    }
          
    return particles;
  }
  
  template
  <typename Particle, typename BeadPack, typename FlowField, typename Boundary, typename StateMaker>
  auto make_particles_random_flux_weighted_box
  (std::size_t nr_particles, BeadPack const& bead_pack,
   FlowField const& flow_field, Boundary const& boundary,
   std::vector<std::pair<double, double>> const& boundaries, StateMaker state_maker)
  {
    using State = typename Particle::State;
    std::vector<Particle> particles;
    particles.reserve(nr_particles);
    
    std::mt19937 rng{ std::random_device{}() };
    std::vector<std::uniform_real_distribution<double>> dist;
    dist.reserve(BeadPack::dim);
    for (std::size_t dd = 0; dd < BeadPack::dim; ++dd)
      dist.emplace_back(boundaries[dd].first, boundaries[dd].second);
    std::vector<State> states;
    states.reserve(nr_particles);
    while (states.size() < nr_particles)
    {
      auto state{ state_maker() };
      for (std::size_t dd = 0; dd < BeadPack::dim; ++dd)
        state.position[dd] = dist[dd](rng);
      boundary(state);
      if (!bead_pack.inside(state.position).first)
        states.push_back(state);
    }
    
    std::vector<double> weights;
    weights.reserve(states.size());
    for (auto const& state : states)
      weights.push_back(operation::abs(flow_field(state.position)));
    std::discrete_distribution<std::size_t> fw_dist{
      weights.begin(), weights.end() };
    for (std::size_t ii = 0; ii < nr_particles; ++ii)
      particles.push_back(states[fw_dist(rng)]);
    
//    for (std::size_t ii = 0; ii < particles.size(); ++ii)
//    {
//      useful::print(std::cout, particles[ii].state_new().position);
//      std::cout << "\n";
//    }
          
    return particles;
  }
  
  template <typename BeadPack, typename CTRW, typename Boundary, typename StripHandler>
  void make_strips_random_uniform_box
  (std::size_t nr_strips,
   std::size_t particles_strip, std::size_t max_particles_strip,
   double initial_segment_length, double max_segment_length,
   BeadPack const& bead_pack, Boundary const& boundary,
   std::vector<std::pair<double, double>> const& boundaries,
   CTRW& ctrw,
   StripHandler& strips)
  {
    ctrw.reserve(nr_strips*max_particles_strip);
    strips.reserve(nr_strips);
    
    using State = typename CTRW::State;
    
    std::vector<double> position_first(BeadPack::dim);
    std::vector<double> position_second(BeadPack::dim);
    std::vector<int> periodicity(BeadPack::dim, 0);
    std::mt19937 rng{ std::random_device{}() };
    std::vector<std::uniform_real_distribution<double>> dist;
    dist.reserve(BeadPack::dim);
    
    for (std::size_t dd = 0; dd < BeadPack::dim; ++dd)
      dist.emplace_back(boundaries[dd].first, boundaries[dd].second);
    
    while (ctrw.size() < 2*nr_strips)
    {
      while (1)
      {
        for (std::size_t dd = 0; dd < BeadPack::dim; ++dd)
          position_first[dd] = dist[dd](rng);
        auto state = State{
          position_first, periodicity,
          {}, {}, 2*strips.size() };
        boundary(state);
        if (!bead_pack.inside(position_first).first)
        {
          ctrw.push_back(state);
          break;
        }
      }
      while (1)
      {
        auto displacement = stochastic::isotropic_unit_vector_distribution{
          BeadPack::dim }(rng);
        operation::times_scalar_InPlace(initial_segment_length, displacement);
        operation::plus(position_first, displacement, position_second);
        
        if (boundary::outOfBounds_box(position_second, boundaries))
          continue;
        auto state = State{
          position_second, periodicity,
          {}, {}, 2*strips.size()+1 };
        boundary(state);
        if (!bead_pack.inside(position_second).first)
        {
          ctrw.push_back(state);
          break;
        }
      }
    }
    for (std::size_t pp = 0; pp < nr_strips; pp += 2)
    {
      strips.push_back({ pp, pp + 1 }, max_particles_strip,
        max_segment_length);
    }
    strips.refine(particles_strip);
  }
  
  template
  <typename Particle, typename BeadPack, typename Boundary, typename StateMaker>
  auto make_particles_random_uniform_plane
  (std::size_t nr_particles,
   BeadPack const& bead_pack, Boundary const& boundary,
   std::vector<double> const& plane_normal,
   std::vector<double> const& point_on_plane,
   std::vector<std::pair<double, double>> const& boundaries_on_plane,
   StateMaker state_maker)
  {
    using State = typename Particle::State;
    std::vector<Particle> particles;
    particles.reserve(nr_particles);
    
    auto basis = algebra::gram_schmidt(plane_normal);
    
    std::vector<int> periodicity(BeadPack::dim, 0);
    std::mt19937 rng{ std::random_device{}() };
    std::vector<std::uniform_real_distribution<double>> dist;
    dist.reserve(BeadPack::dim);
    for (std::size_t dd = 0; dd < BeadPack::dim-1; ++dd)
      dist.emplace_back(boundaries_on_plane[dd].first,
                        boundaries_on_plane[dd].second);
    
    while (particles.size() < nr_particles)
    {
      State state{ state_maker() };
      state.position = point_on_plane;
      for (std::size_t dd = 0; dd < BeadPack::dim-1; ++dd)
        operation::plus_InPlace(state.position,
          operation::times_scalar(dist[dd](rng), basis[dd+1]));
      boundary(state);
      
      if (!bead_pack.inside(state.position).first)
        particles.push_back(state);
    }
          
    return particles;
  }
  
  template
  <typename Particle, typename BeadPack, typename FlowField,
  typename Boundary, typename StateMaker>
  auto make_particles_random_fluxweighted_plane
  (std::size_t nr_particles, BeadPack const& bead_pack,
   FlowField const& flow_field, Boundary const& boundary,
   std::vector<double> const& plane_normal,
   std::vector<double> const& point_on_plane,
   std::vector<std::pair<double, double>> const& boundaries_on_plane,
   StateMaker state_maker)
  {
    using State = typename Particle::State;
    std::vector<Particle> particles;
    particles.reserve(nr_particles);
    
    auto basis = algebra::gram_schmidt(plane_normal);
    
    std::vector<int> periodicity(BeadPack::dim, 0);
    std::mt19937 rng{ std::random_device{}() };
    std::vector<std::uniform_real_distribution<double>> dist;
    dist.reserve(BeadPack::dim);
    for (std::size_t dd = 0; dd < BeadPack::dim - 1; ++dd)
      dist.emplace_back(boundaries_on_plane[dd].first,
                        boundaries_on_plane[dd].second);
    
    std::vector<State> states;
    states.reserve(nr_particles);
    while (states.size() < nr_particles)
    {
      State state{ state_maker() };
      state.position = point_on_plane;
      for (std::size_t dd = 0; dd < BeadPack::dim-1; ++dd)
        operation::plus_InPlace(state.position,
          operation::times_scalar(dist[dd](rng), basis[dd+1]));
      boundary(state);
      
      if (!bead_pack.inside(state.position).first)
        states.push_back(state);
    }
    
    std::vector<double> weights;
    for (auto const& state : states)
      weights.push_back(operation::abs(flow_field(state.position)));
    std::discrete_distribution<std::size_t> fw_dist{
      weights.begin(), weights.end() };
    for (std::size_t ii = 0; ii < nr_particles; ++ii)
      particles.push_back(states[fw_dist(rng)]);
          
    return particles;
  }
  
  template
  <typename Particle, typename BeadPack, typename Boundary, typename StateMaker>
  auto make_particles_random_near_wall_uniform_unit_cell
  (std::size_t nr_particles, BeadPack const& bead_pack, Boundary const& boundary,
   double length_near_wall,
   StateMaker state_maker)
  {
    std::vector<Particle> particles;
    particles.reserve(nr_particles);
    
    std::mt19937 rng{ std::random_device{}() };
    std::uniform_real_distribution<double> dist_uniform{ 0., 1. };
    
    std::vector<double> weights;
    weights.reserve(bead_pack.nr_beads());
    for (auto const& bead : bead_pack.beads())
      weights.push_back(bead.radius*bead.radius);
    std::discrete_distribution<std::size_t> dist_beads{ weights.begin(), weights.end() };
    
    while (particles.size() < nr_particles)
    {
      auto state{ state_maker() };
      std::size_t bead = dist_beads(rng);
      double radius = bead_pack.radius(bead)+2.*length_near_wall;
      state.position = bead_pack.bead(bead).center;
      if constexpr (BeadPack::dim == 2)
      {
        double phi = 2.*constants::pi*dist_uniform(rng);
        state.position[0] += radius*std::cos(phi);
        state.position[1] += radius*std::sin(phi);
      }
      if constexpr (BeadPack::dim == 3)
      {
        double phi = 2.*constants::pi*dist_uniform(rng);
        double theta = constants::pi*(2.*dist_uniform(rng)-1.);
        double costheta = std::cos(theta);
        state.position[0] += radius*std::cos(phi)*costheta;
        state.position[1] += radius*std::sin(phi)*costheta;
        state.position[2] += radius*std::sin(theta);
      }
      else
        throw std::invalid_argument( "Dimension not supported" );
      
      if (boundary(state))
        continue;
      if (!bead_pack.near(state.position, 2.*length_near_wall).first)
        particles.push_back(state);
    }
          
    return particles;
  }
  
  template
  <typename Particle, typename BeadPack, typename Boundary, typename StateMaker>
  auto make_particles_random_near_wall_uniform_bead
  (std::size_t nr_particles, BeadPack const& bead_pack, Boundary const& boundary,
   double length_near_wall,
   StateMaker state_maker)
  {
    std::vector<Particle> particles;
    particles.reserve(nr_particles);
    
    std::mt19937 rng{ std::random_device{}() };
    std::uniform_real_distribution<double> dist_uniform{ 0., 1. };
    
    std::vector<double> weights;
    weights.reserve(bead_pack.nr_beads());
    for (auto const& bead : bead_pack.beads())
      weights.push_back(bead.radius*bead.radius);
    std::discrete_distribution<std::size_t> dist_beads{ weights.begin(), weights.end() };
    
    while (particles.size() < nr_particles)
    {
      auto state{ state_maker() };
      std::size_t bead = dist_beads(rng);
      double radius = bead_pack.radius(bead)+2.*length_near_wall;
      state.position = bead_pack.bead(bead).center;
      if constexpr (BeadPack::dim == 2)
      {
        double phi = 2.*constants::pi*dist_uniform(rng);
        state.position[0] += radius*std::cos(phi);
        state.position[1] += radius*std::sin(phi);
      }
      if constexpr (BeadPack::dim == 3)
      {
        double phi = 2.*constants::pi*dist_uniform(rng);
        double theta = constants::pi*(2.*dist_uniform(rng)-1.);
        double sintheta = std::sin(theta);
        state.position[0] += radius*std::cos(phi)*sintheta;
        state.position[1] += radius*std::sin(phi)*sintheta;
        state.position[2] += radius*std::cos(theta);
      }
      else
        throw std::invalid_argument( "Dimension not supported" );
      
      boundary(state);
      if (!bead_pack.near(state.position, 2.*length_near_wall).first)
        particles.push_back(state);
    }
          
    return particles;
  }
  
  template
  <typename Particle, typename BeadPack, typename Boundary, typename StateMaker>
  auto make_particles_random_at_wall_uniform_unit_cell
  (std::size_t nr_particles, BeadPack const& bead_pack, Boundary const& boundary,
   StateMaker state_maker)
  {
    std::vector<Particle> particles;
    particles.reserve(nr_particles);
    
    std::mt19937 rng{ std::random_device{}() };
    std::uniform_real_distribution<double> dist_uniform{ 0., 1. };
    
    std::vector<double> weights;
    weights.reserve(bead_pack.nr_beads());
    for (auto const& bead : bead_pack.beads())
      weights.push_back(bead.radius*bead.radius);
    std::discrete_distribution<std::size_t> dist_beads{ weights.begin(), weights.end() };
    
    while (particles.size() < nr_particles)
    {
      auto state{ state_maker() };
      std::size_t bead = dist_beads(rng);
      double radius = bead_pack.radius(bead);
      state.position = bead_pack.bead(bead).center;
      if constexpr (BeadPack::dim == 2)
      {
        double phi = 2.*constants::pi*dist_uniform(rng);
        state.position[0] += radius*std::cos(phi);
        state.position[1] += radius*std::sin(phi);
      }
      if constexpr (BeadPack::dim == 3)
      {
        double phi = 2.*constants::pi*dist_uniform(rng);
        double theta = constants::pi*(2.*dist_uniform(rng)-1.);
        double sintheta = std::sin(theta);
        state.position[0] += radius*std::cos(phi)*sintheta;
        state.position[1] += radius*std::sin(phi)*sintheta;
        state.position[2] += radius*std::cos(theta);
      }
      else
        throw std::invalid_argument( "Dimension not supported" );
      
      if (boundary(state))
        continue;
      particles.push_back(state);
    }
          
    return particles;
  }
  
  template
  <typename Particle, typename BeadPack, typename Boundary, typename StateMaker>
  auto make_particles_random_at_wall_uniform_bead
  (std::size_t nr_particles, BeadPack const& bead_pack, Boundary const& boundary,
   StateMaker state_maker)
  {
    std::vector<Particle> particles;
    particles.reserve(nr_particles);
    
    std::mt19937 rng{ std::random_device{}() };
    std::uniform_real_distribution<double> dist_uniform{ 0., 1. };
    
    std::vector<double> weights;
    weights.reserve(bead_pack.nr_beads());
    for (auto const& bead : bead_pack.beads())
      weights.push_back(bead.radius*bead.radius);
    std::discrete_distribution<std::size_t> dist_beads{ weights.begin(), weights.end() };
    
    while (particles.size() < nr_particles)
    {
      auto state{ state_maker() };
      std::size_t bead = dist_beads(rng);
      double radius = bead_pack.radius(bead);
      state.position = bead_pack.bead(bead).center;
      if constexpr (BeadPack::dim == 2)
      {
        double phi = 2.*constants::pi*dist_uniform(rng);
        state.position[0] += radius*std::cos(phi);
        state.position[1] += radius*std::sin(phi);
      }
      if constexpr (BeadPack::dim == 3)
      {
        double phi = 2.*constants::pi*dist_uniform(rng);
        double theta = constants::pi*(2.*dist_uniform(rng)-1.);
        double sintheta = std::sin(theta);
        state.position[0] += radius*std::cos(phi)*sintheta;
        state.position[1] += radius*std::sin(phi)*sintheta;
        state.position[2] += radius*std::cos(theta);
      }
      else
        throw std::invalid_argument( "Dimension not supported" );
      
      boundary(state);
      particles.push_back(state);
    }
          
    return particles;
  }
  
  template
  <typename Particle, typename Boundary, typename StateMaker>
  auto make_particles_load_positions
  (std::size_t nr_particles, std::string const& filename,
   Boundary const& boundary, StateMaker state_maker)
  {
    std::vector<Particle> particles;
    particles.reserve(nr_particles);
    
    std::ifstream input{ filename };
    if (!input.is_open())
      throw useful::open_read_error(filename);
    for (std::size_t pp = 0; pp < nr_particles; ++pp)
    {
      auto state{ state_maker() };
      for (std::size_t dd = 0; dd < state.position.size(); ++dd)
        if (!(input >> state.position[dd]))
          throw useful::parse_error_file(filename);
      boundary(state);
      particles.push_back(state);
    }
          
    return particles;
  }
  
  template
  <typename Particle, typename BeadPack, typename FlowField,
  typename Boundary, typename StateMaker>
  auto make_particles
  (std::size_t nr_particles, int initial_condition_type,
   std::vector<double> const& initial_box_midpoint,
   std::vector<std::pair<double, double>> const& initial_box_centered,
   FlowField const &velocity_field,
   std::vector<double> const& mean_velocity,
   BeadPack const& bead_pack, Boundary const& boundary_periodic,
   double length_near_wall,
   std::string const& filename,
   StateMaker state_maker)
  {
    std::vector<std::pair<double, double>> initial_box = initial_box_centered;
    for (std::size_t dd = 0; dd < BeadPack::dim; ++dd)
    {
      initial_box[dd].first += initial_box_midpoint[dd];
      initial_box[dd].second += initial_box_midpoint[dd];
    }
    switch (initial_condition_type)
    {
      case 0:
        return beadpack::make_particles_random_uniform_plane<Particle>
          (nr_particles, bead_pack, boundary_periodic,
           mean_velocity, initial_box_midpoint, initial_box_centered, state_maker);
      case 1:
        return beadpack::make_particles_random_fluxweighted_plane<Particle>
          (nr_particles, bead_pack, velocity_field, boundary_periodic,
           mean_velocity, initial_box_midpoint, initial_box_centered, state_maker);
      case 2:
        return beadpack::make_particles_random_uniform_box<Particle>
          (nr_particles, bead_pack, boundary_periodic,
           initial_box, state_maker);
      case 3:
        return beadpack::make_particles_random_flux_weighted_box<Particle>
          (nr_particles, bead_pack, velocity_field, boundary_periodic,
           initial_box, state_maker);
      case 4:
        return beadpack::make_particles_random_near_wall_uniform_unit_cell<Particle>
          (nr_particles, bead_pack, boundary_periodic, length_near_wall, state_maker);
      case 5:
        return beadpack::make_particles_random_near_wall_uniform_bead<Particle>
          (nr_particles, bead_pack, boundary_periodic, length_near_wall, state_maker);
      case 6:
        return beadpack::make_particles_random_at_wall_uniform_unit_cell<Particle>
          (nr_particles, bead_pack, boundary_periodic, state_maker);
      case 7:
        return beadpack::make_particles_random_at_wall_uniform_bead<Particle>
          (nr_particles, bead_pack, boundary_periodic, state_maker);
      case 8:
        return beadpack::make_particles_load_positions<Particle>
          (nr_particles, filename, boundary_periodic, state_maker);
      default:
        throw std::invalid_argument{ "Undefined initial condition type" };
    }
  }
  
  std::string initial_condition_name(int initial_condition_type)
  {
    switch (initial_condition_type)
    {
      case 0:
        return "uniform_plane";
      case 1:
        return "flux_weighted_plane";
      case 2:
        return "uniform_box";
      case 3:
        return "flux_weighted_box";
      case 4:
        return "near_wall_uniform_unit_cell";
      case 5:
        return "near_wall_uniform_bead";
      case 6:
        return "at_wall_uniform_unit_cell";
      case 7:
        return "at_wall_uniform_bead";
      case 8:
        return "load_positions";
      default:
        throw std::invalid_argument{ "Undefined initial condition type" };
    }
  }
}


#endif /* BeadPack_InitialConditions_h */
