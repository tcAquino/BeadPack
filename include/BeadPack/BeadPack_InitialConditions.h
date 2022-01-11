//
//  BeadPack_InitialConditions.h
//  BeadPack
//
//  Created by Tomás Aquino on 14/05/2020.
//  Copyright © 2020 Tomás Aquino. All rights reserved.
//

#ifndef BeadPack_InitialConditions_h
#define BeadPack_InitialConditions_h

#include <boost/algorithm/string.hpp>
#include <fstream>
#include <cmath>
#include <random>
#include <utility>
#include <vector>
#include "Algebra/Algebra.h"
#include "general/Constants.h"
#include "general/useful.h"
#include "Stochastic/Random.h"
#include "Stochastic/CTRW/Boundary.h"

namespace beadpack
{
  // Functions to initialize particles for CTRW/PTRW in a beadpack
  // StateMaker objects initialize states to a valid state
  // (e.g., positions with appropriate number of components,
  // initial mass, etc)
  
  // Distribute particles uniformly randomly in the void space
  // within specified boundaries along each dimension
  // Each component of the boundaries vector
  // is the lower and upper bound along a dimension
  // Boundary object enforces periodic boundary conditions
  // on beadpack if necessary
  template
  <typename Particle, typename BeadPack,
  typename Boundary, typename StateMaker>
  auto make_particles_random_uniform_box
  (std::size_t nr_particles, BeadPack const& bead_pack,
   Boundary const& boundary,
   std::vector<std::pair<double, double>> const& boundaries,
   StateMaker state_maker)
  {
    std::vector<Particle> particles;
    particles.reserve(nr_particles);
    std::mt19937 rng{ std::random_device{}() };
   
    // To produce uniformly random positions within
    // specified bounding box
    std::vector<std::uniform_real_distribution<double>> dist;
    dist.reserve(BeadPack::dim);
    for (std::size_t dd = 0; dd < BeadPack::dim; ++dd)
      dist.emplace_back(boundaries[dd].first,
                        boundaries[dd].second);
    
    while (particles.size() < nr_particles)
    {
      auto state{ state_maker() };
      for (std::size_t dd = 0; dd < BeadPack::dim; ++dd)
        state.position[dd] = dist[dd](rng);
      boundary(state);
      
      // Discard if inside a bead
      if (!bead_pack.inside(state.position).first)
        particles.push_back(state);
    }
          
    return particles;
  }
  
  // Distribute particles with probability proportional to
  // flow velocity magnitude in the void space
  // within specified boundaries along each dimension
  // Each component of the boundaries vector
  // is the lower and upper bound along a dimension
  // Boundary object enforces periodic boundary conditions
  // on beadpack if necessary
  template
  <typename Particle, typename BeadPack,
  typename FlowField, typename Boundary, typename StateMaker>
  auto make_particles_random_flux_weighted_box
  (std::size_t nr_particles, BeadPack const& bead_pack,
   FlowField const& flow_field, Boundary const& boundary,
   std::vector<std::pair<double, double>> const& boundaries,
   StateMaker state_maker)
  {
    std::vector<Particle> particles;
    particles.reserve(nr_particles);
    std::mt19937 rng{ std::random_device{}() };
    
    // To produce uniformly random positions within
    // specified bounding box
    std::vector<std::uniform_real_distribution<double>> dist;
    dist.reserve(BeadPack::dim);
    for (std::size_t dd = 0; dd < BeadPack::dim; ++dd)
      dist.emplace_back(
        boundaries[dd].first, boundaries[dd].second);
    
    // Find admissible states
    using State = typename Particle::State;
    std::vector<State> states;
    states.reserve(nr_particles);
    while (states.size() < nr_particles)
    {
      auto state{ state_maker() };
      for (std::size_t dd = 0; dd < BeadPack::dim; ++dd)
        state.position[dd] = dist[dd](rng);
      boundary(state);
      
      // Discard if inside a bead
      if (!bead_pack.inside(state.position).first)
        states.push_back(state);
    }
    
    // Choose from admissible states with probability
    // proportional to flow velocity magnitude
    std::vector<double> weights;
    weights.reserve(states.size());
    for (auto const& state : states)
      weights.push_back(operation::abs(
        flow_field(state.position)));
    std::discrete_distribution<std::size_t> fw_dist{
      weights.begin(), weights.end() };
    for (std::size_t ii = 0; ii < nr_particles; ++ii)
      particles.push_back(states[fw_dist(rng)]);
          
    return particles;
  }
  
  // Distribute particles uniformly randomly in the void space
  // on a plane normal to plane_normal and including the point point_on_plane,
  // within specified boundaries along each dimension on the plane
  // Each component of the boundaries_on_plane vector
  // is the lower and upper bound along a dimension on the plane,
  // centered around point_on_plane
  // Boundary object enforces periodic boundary conditions
  // on beadpack if necessary
  template
  <typename Particle, typename BeadPack,
  typename Boundary, typename StateMaker>
  auto make_particles_random_uniform_plane
  (std::size_t nr_particles,
   BeadPack const& bead_pack, Boundary const& boundary,
   std::vector<double> const& plane_normal,
   std::vector<double> const& point_on_plane,
   std::vector<std::pair<double, double>> const& boundaries_on_plane,
   StateMaker state_maker)
  {
    std::vector<Particle> particles;
    particles.reserve(nr_particles);
    std::mt19937 rng{ std::random_device{}() };
    
    // To produce uniformly random positions within
    // specified bounding box
    std::vector<std::uniform_real_distribution<double>> dist;
    dist.reserve(BeadPack::dim);
    for (std::size_t dd = 0; dd < BeadPack::dim-1; ++dd)
      dist.emplace_back(boundaries_on_plane[dd].first,
                        boundaries_on_plane[dd].second);
    
    // Find a set of basis vectors on the plane
    auto basis = algebra::gram_schmidt(plane_normal);
    
    while (particles.size() < nr_particles)
    {
      // Use the basis vectors on the plane to produce
      // uniformly distributed positions within the boundaries
      auto state{ state_maker() };
      state.position = point_on_plane;
      for (std::size_t dd = 0; dd < BeadPack::dim-1; ++dd)
        operation::plus_InPlace(state.position,
          operation::times_scalar(dist[dd](rng), basis[dd+1]));
      boundary(state);
      
      // Discard if inside a bead
      if (!bead_pack.inside(state.position).first)
        particles.push_back(state);
    }
          
    return particles;
  }
  
  // Distribute particles with probability proportional
  // to the flow velocity magnitude in the void space
  // on a plane normal to plane_normal and including the point point_on_plane,
  // within specified boundaries along each dimension on the plane
  // Each component of the boundaries_on_plane vector
  // is the lower and upper bound along a dimension on the plane,
  // centered around point_on_plane
  // Boundary object enforces periodic boundary conditions
  // on beadpack if necessary
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
    std::vector<Particle> particles;
    particles.reserve(nr_particles);
    std::mt19937 rng{ std::random_device{}() };
    
    // To produce uniformly random positions within
    // specified bounding box
    std::vector<std::uniform_real_distribution<double>> dist;
    dist.reserve(BeadPack::dim);
    for (std::size_t dd = 0; dd < BeadPack::dim - 1; ++dd)
      dist.emplace_back(boundaries_on_plane[dd].first,
                        boundaries_on_plane[dd].second);
    
    // Find a set of basis vectors on the plane
    auto basis = algebra::gram_schmidt(plane_normal);
    
    // Find admissible states
    using State = typename Particle::State;
    std::vector<State> states;
    states.reserve(nr_particles);
    while (states.size() < nr_particles)
    {
      auto state{ state_maker() };
      state.position = point_on_plane;
      for (std::size_t dd = 0; dd < BeadPack::dim-1; ++dd)
        operation::plus_InPlace(state.position,
          operation::times_scalar(dist[dd](rng), basis[dd+1]));
      boundary(state);
      
      if (!bead_pack.inside(state.position).first)
        states.push_back(state);
    }
    
    // Choose from admissible states with probability
    // proportional to flow velocity magnitude
    std::vector<double> weights;
    for (auto const& state : states)
      weights.push_back(operation::abs(
        flow_field(state.position)));
    std::discrete_distribution<std::size_t> fw_dist{
      weights.begin(), weights.end() };
    for (std::size_t ii = 0; ii < nr_particles; ++ii)
      particles.push_back(states[fw_dist(rng)]);
          
    return particles;
  }
  
  // Distribute particles uniformly randomly at a
  // distance length_near_wall to the void-bead interface
  // on the periodic unit cell enforced by the Boundary object
  template
  <typename Particle, typename BeadPack,
  typename Boundary, typename StateMaker>
  auto make_particles_random_near_wall_uniform_unit_cell
  (std::size_t nr_particles, BeadPack const& bead_pack,
   Boundary const& boundary, double length_near_wall,
   StateMaker state_maker)
  {
    std::vector<Particle> particles;
    particles.reserve(nr_particles);
    std::mt19937 rng{ std::random_device{}() };
    
    // Probability of being on the surface of a given sphere
    // is proportional to its radius squared (surface area)
    std::vector<double> weights;
    weights.reserve(bead_pack.nr_beads());
    for (auto const& bead : bead_pack.beads())
    {
      double radius = bead.radius+length_near_wall;
      weights.push_back(radius*bead.radius);
    }
    std::discrete_distribution<std::size_t> dist_beads{
      weights.begin(), weights.end() };
    
    while (particles.size() < nr_particles)
    {
      stochastic::isotropic_unit_vector_distribution dist{
        BeadPack::dim };
      std::size_t bead = dist_beads(rng);
      double radius = bead_pack.radius(bead)+length_near_wall;
      auto state{ state_maker() };
      state.position =
        operation::plus(bead_pack.bead(bead).center,
          operation::times_scalar(radius, dist(rng)));
      
      // Discard positions outside unit cell
      if (boundary(state))
        continue;
      // Discard positions too close to another bead
      if (!bead_pack.near(state.position, length_near_wall).first)
        particles.push_back(state);
    }
          
    return particles;
  }
  
  // Distribute particles uniformly randomly at a
  // distance length_near_wall to the void-bead interface
  // over all beads in the beadpack
  // Discard particles outside initial_box
  // Boundary object enforces periodic boundary conditions
  // on beadpack if necessary
  template
  <typename Particle, typename BeadPack,
  typename Boundary, typename StateMaker>
  auto make_particles_random_near_wall_uniform_bead
  (std::size_t nr_particles, BeadPack const& bead_pack,
   std::vector<std::pair<double, double>> const& initial_box,
   Boundary const& boundary, double length_near_wall,
   StateMaker state_maker)
  {
    std::vector<Particle> particles;
    particles.reserve(nr_particles);
    std::mt19937 rng{ std::random_device{}() };
    
    // Probability of being on the surface of a given sphere
    // is proportional to its radius squared (surface area)
    std::vector<double> weights;
    weights.reserve(bead_pack.nr_beads());
    for (auto const& bead : bead_pack.beads())
    {
      double radius = bead.radius+length_near_wall;
      weights.push_back(radius*bead.radius);
    }
    std::discrete_distribution<std::size_t> dist_beads{
      weights.begin(), weights.end() };
    
    while (particles.size() < nr_particles)
    {
      // Produce uniformly distributed positions on a bead
      // surface chosen with appropriate probability
      stochastic::isotropic_unit_vector_distribution dist{
        BeadPack::dim };
      std::size_t bead = dist_beads(rng);
      double radius = bead_pack.radius(bead)+length_near_wall;
      auto state{ state_maker() };
      state.position =
        operation::plus(bead_pack.bead(bead).center,
          operation::times_scalar(radius, dist(rng)));
      for (std::size_t dd = 0; dd < state.position.size(); ++dd)
        if (state.position[dd] < initial_box[dd].first
            || state.position[dd] > initial_box[dd].second)
          continue;
      boundary(state);
      
      // Discard if too close to another bead
      if (!bead_pack.near(state.position,
                          length_near_wall).first)
        particles.push_back(state);
    }
          
    return particles;
  }
  
  // Distribute particles uniformly randomly
  // at the void-bead interface
  // on the periodic unit cell enforced by the Boundary object
  template
  <typename Particle, typename BeadPack,
  typename Boundary, typename StateMaker>
  auto make_particles_random_at_wall_uniform_unit_cell
  (std::size_t nr_particles, BeadPack const& bead_pack,
   Boundary const& boundary, StateMaker state_maker)
  {
    return
      make_particles_random_near_wall_uniform_unit_cell<Particle>(
        nr_particles, bead_pack, boundary, 0., state_maker);
  }
  
  // Distribute particles uniformly randomly at
  // the void-bead interface over all beads in the beadpack
  // Boundary object enforces periodic boundary conditions
  // on beadpack if necessary
  template
  <typename Particle, typename BeadPack,
  typename Boundary, typename StateMaker>
  auto make_particles_random_at_wall_uniform_bead
  (std::size_t nr_particles, BeadPack const& bead_pack,
   std::vector<std::pair<double, double>> const& initial_box,
   Boundary const& boundary, StateMaker state_maker)
  {
    return
      make_particles_random_near_wall_uniform_bead<Particle>(
        nr_particles, bead_pack, initial_box, boundary, 0., state_maker);
  }
  
  // Load particle positions from file
  // Load only the first nr_particles in file
  // File contents are the values of the position components
  // Components of a position must be row-major contiguous
  template
  <typename Particle, typename Boundary, typename StateMaker>
  auto make_particles_load_positions
  (std::size_t nr_particles, std::string const& filename,
   Boundary const& boundary, StateMaker state_maker)
  {
    std::vector<Particle> particles;
    particles.reserve(nr_particles);
    
    auto input = useful::open_read(filename);
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
  
  // Load particle positions from file
  // First dim columns in each line are the values of the position components,
  // any further columns are ignored
  template
  <typename Particle, typename Boundary, typename StateMaker>
  auto make_particles_load_positions
  (std::string const& filename,
   Boundary const& boundary, StateMaker state_maker,
   std::size_t header_lines = 0, double rescale = 1.,
   std::size_t nr_estimate = 0, std::string const& delims = "\t,| ")
  {
    std::vector<Particle> particles;
    particles.reserve(nr_estimate);
    
    auto input = useful::open_read(filename);
    std::string line;
    for (std::size_t ll = 0; ll < header_lines; ++ll)
      getline(input, line);
    
    while (getline(input, line))
    {
      std::vector<std::string> split_line;
      boost::trim_if(line, boost::is_any_of(delims+"\r"));
      boost::algorithm::split(split_line, line, boost::is_any_of(delims),
                              boost::token_compress_on);
      auto state{ state_maker() };
      for (std::size_t dd = 0; dd < state.position.size(); ++dd)
        state.position[dd] = std::stod(split_line[dd]);
      boundary(state);
      particles.push_back(state);
    }
          
    return particles;
  }
  
  // Generic interface to choose an initial condition
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
    std::vector<std::pair<double, double>> initial_box =
      initial_box_centered;
    for (std::size_t dd = 0; dd < BeadPack::dim; ++dd)
    {
      initial_box[dd].first += initial_box_midpoint[dd];
      initial_box[dd].second += initial_box_midpoint[dd];
    }
    switch (initial_condition_type)
    {
      case 0:
        return
          beadpack::make_particles_random_uniform_plane<Particle>(
            nr_particles, bead_pack, boundary_periodic, mean_velocity,
            initial_box_midpoint, initial_box_centered, state_maker);
      case 1:
        return
          beadpack::make_particles_random_fluxweighted_plane<Particle>(
            nr_particles, bead_pack, velocity_field, boundary_periodic,
            mean_velocity, initial_box_midpoint,
            initial_box_centered, state_maker);
      case 2:
        return
          beadpack::make_particles_random_uniform_box<Particle>(
            nr_particles, bead_pack, boundary_periodic,
            initial_box, state_maker);
      case 3:
        return
          beadpack::make_particles_random_flux_weighted_box<Particle>(
            nr_particles, bead_pack, velocity_field,
            boundary_periodic, initial_box, state_maker);
      case 4:
        return
          beadpack::make_particles_random_near_wall_uniform_unit_cell<Particle>(
            nr_particles, bead_pack, boundary_periodic,
            length_near_wall, state_maker);
      case 5:
        return
          beadpack::make_particles_random_near_wall_uniform_bead<Particle>(
            nr_particles, bead_pack, initial_box, boundary_periodic,
            length_near_wall, state_maker);
      case 6:
        return
          beadpack::make_particles_random_at_wall_uniform_unit_cell<Particle>(
            nr_particles, bead_pack, boundary_periodic,
            state_maker);
      case 7:
        return
          beadpack::make_particles_random_at_wall_uniform_bead<Particle>(
            nr_particles, bead_pack, initial_box, boundary_periodic,
            state_maker);
      case 8:
        return
          beadpack::make_particles_load_positions<Particle>(
            nr_particles, filename, boundary_periodic,
            state_maker);
      case 9:
        return
          beadpack::make_particles_load_positions<Particle>(
            filename, boundary_periodic, state_maker);
      default:
        throw std::invalid_argument{ "Undefined initial condition type" };
    }
  }
  
  // Convert integer initial condition type to
  // string initial condition name
  std::string initial_condition_name
  (int initial_condition_type)
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
        throw std::invalid_argument{
          "Undefined initial condition type" };
    }
  }
  
  // Make randomly distributed strips (collections of particles)
  // uniformly randomly distributed on the void space
  // within specified boundaries
  // Parameters:
  //   nr_strips: Number of strips to create
  //   particles_strip: Particles to place in each strip
  //   max_particles_strip: Maximum particles each strip can hold
  //   initial_segment_length: Initial distance between edges of a strip
  //   bead_pack: Beadpack object, strips must be in the void space
  //   boundary: Enforce periodic boundary conditions
  //   boundaries: Lower and upper bounds along each dimension for initial condition
  //   ctrw: Strips are made associated with this CTRW object
  //   strips: Strips are made associated with this strip handler
  template <typename BeadPack, typename CTRW,
  typename Boundary, typename StripHandler>
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
    std::mt19937 rng{ std::random_device{}() };
    
    using State = typename CTRW::State;
    
    // To initialize states, keeps track of which periodic unit
    // cell particles are in
    std::vector<int> periodicity(BeadPack::dim, 0);
    
    // To produce uniformly random positions within
    // specified bounding box
    std::vector<std::uniform_real_distribution<double>> dist;
    dist.reserve(BeadPack::dim);
    for (std::size_t dd = 0; dd < BeadPack::dim; ++dd)
      dist.emplace_back(boundaries[dd].first,
                        boundaries[dd].second);
    
    // Auxiliaries to hold positions of strip edges
    std::vector<double> position_first(BeadPack::dim);
    std::vector<double> position_second(BeadPack::dim);
    
    // Make the two strip edges for each strip
    while (ctrw.size() < 2*nr_strips)
    {
      // Place one end of the strip uniformly randomly
      // in the void space
      while (1)
      {
        for (std::size_t dd = 0; dd < BeadPack::dim; ++dd)
          position_first[dd] = dist[dd](rng);
        auto state = State{
          position_first, periodicity,
          {}, {}, 2*strips.size() };
        boundary(state);
        
        // Discsard if inside bead
        if (!bead_pack.inside(position_first).first)
        {
          ctrw.push_back(state);
          break;
        }
      }
      // Place the other end of the strip
      // at distance initial_segment length in a random direction,
      // within the void space
      while (1)
      {
        auto displacement = stochastic::isotropic_unit_vector_distribution{
          BeadPack::dim }(rng);
        operation::times_scalar_InPlace(initial_segment_length,
                                        displacement);
        operation::plus(position_first, displacement,
                        position_second);
        
        // Discard if outside bounding box
        if (boundary::outOfBounds_box(position_second,
                                      boundaries))
          continue;
        auto state = State{
          position_second, periodicity,
          {}, {}, 2*strips.size()+1 };
        boundary(state);
        // Discard if inside bead
        if (!bead_pack.inside(position_second).first)
        {
          ctrw.push_back(state);
          break;
        }
      }
    }
    
    // Register the 2-particle strips just made with the strip handler
    for (std::size_t pp = 0; pp < nr_strips; pp += 2)
    {
      strips.push_back({ pp, pp + 1 }, max_particles_strip,
        max_segment_length);
    }
    
    // Refine the strips between the edges by adding remaining
    // initial particles to each strip
    strips.refine(particles_strip);
  }
}


#endif /* BeadPack_InitialConditions_h */
