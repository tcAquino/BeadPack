//
//  StateGetter.h
//  BeadPack
//
//  Created by Tomás Aquino on 17/05/2020.
//  Copyright © 2020 Tomás Aquino. All rights reserved.
//

#ifndef StateGetter_h
#define StateGetter_h

#include <utility>
#include "general/Operations.h"

namespace ctrw
{
  // Get info given particle
  // Getter gets info given particle's new state
  template <typename Getter>
  struct Get_new_from_particle
  {
    Getter get;
    
    Get_new_from_particle(Getter&& get = {})
    : get{ std::forward<Getter>(get) }
    {}
    
    template <typename Particle>
    auto operator()(Particle const& particle) const
    {
      return get(particle.state_new());
    }
  };
  template <typename Getter> Get_new_from_particle(Getter&&) -> Get_new_from_particle<Getter>;
  
  // Get info given particle
  // Getter gets info given particle's old state
  template <typename Getter>
  struct Get_old_from_particle
  {
    Getter get;
    
    Get_old_from_particle(Getter&& get = {})
    : get{ std::forward<Getter>(get) }
    {}
    
    template <typename Particle>
    auto operator()(Particle const& particle) const
    {
      return get(particle.state_old());
    }
  };
  template <typename Getter> Get_old_from_particle(Getter&&) -> Get_old_from_particle<Getter>;
  
  // Get info given particle
  // Getter gets info given particle's new state and old state
  template <typename Getter>
  struct Get_from_particle
  {
    Getter get;
    
    Get_from_particle(Getter&& get = {})
    : get{ std::forward<Getter>(get) }
    {}
    template <typename Particle>
    auto operator()(Particle const& particle) const
    {
      return get(particle.state_new(), particle.state_old());
    }
  };
  template <typename Getter> Get_from_particle(Getter&&) -> Get_from_particle<Getter>;
  
  // Get state.position given state
  struct Get_position
  {
    template <typename State>
    auto operator()(State const& state) const
    { return state.position; }
  };
  
  // Get dd component of state.position given state
  template <std::size_t dd>
  struct Get_position_component
  {
    template <typename State>
    auto operator()(State const& state) const
    { return operation::project<dd>(state.position); }
  };
  
  // Get projection of state.position along direction given state
  struct Get_position_projection
  {
    std::vector<double> basis{};
    
    // Construct given direction to project onto
    Get_position_projection(std::vector<double> const& direction)
    : basis{ operation::div_scalar(direction, operation::abs(direction)) }
    {}
    
    template <typename State>
    auto operator()(State const& state) const
    { return operation::dot(state.position, basis); }
  };
  
  // Get dd component squared of state.position given state
  template <std::size_t dd>
  struct Get_position_component_sq
  {
    template <typename State>
    auto operator()(State const& state) const
    {
      double comp = operation::project<dd>(state.position);
      return comp*comp;
    }
  };
  
  // Get dd component of state.position given state
  struct Get_dist_sq
  {
    template <typename State>
    auto operator()(State const& state) const
    { return operation::abs_sq(state.position); }
  };
  
  // Get state.time given state
  struct Get_time
  {
    template <typename State>
    auto operator()(State const& state) const
    { return state.time; }
  };
  
  // Get state.mass given state
  struct Get_mass
  {
    template <typename State>
    auto operator()(State const& state) const
    { return state.mass; }
  };
  
  // Get state.velocity given state
  struct Get_velocity
  {
    template <typename State>
    auto operator()(State const& state) const
    { return state.velocity; }
  };
  
  // Get dd component of velocity given state
  template <std::size_t dd>
  struct Get_velocity_component
  {
    template <typename State>
    auto operator()(State const& state) const
    { return operation::project<dd>(state.velocity); }
  };
  
  // Get function of state.position given state
  template <typename Property>
  struct Get_position_property
  {
    Property property;
    
    // Construct given function of position
    Get_position_property(Property&& property = {})
    : property{ std::forward<Property>(property) }
    {}

    template <typename State>
    auto operator()(State const& state) const
    { return property(state.position); }
  };
  template <typename Property>
  Get_position_property(Property&&) ->
  Get_position_property<Property>;
  
  // Get function of state.velocity given state
  template <typename Property>
  struct Get_velocity_property
  {
    Property property;
    
    // Construct given function of velocity
    Get_velocity_property(Property&& property = {})
    : property{ std::forward<Property>(property) }
    {}

    template <typename State>
    auto operator()(State const& state) const
    { return property(state.velocity); }
  };
  template <typename Property>
  Get_velocity_property(Property&&) ->
  Get_velocity_property<Property>;
  
  // Get absolute position given state
  // from state.position and periodicity information
  // about the unit cell periodic position is in
  // Unit cell faces must be perpendicular to (cartesian) coordinate axes
  struct Get_position_periodic
  {
    std::vector<double> domain_dimensions;
    
    // Construct given size of domain along each dimension
    Get_position_periodic(std::vector<double> domain_dimensions)
    : domain_dimensions{ domain_dimensions }
    {}
    
    template <typename State>
    auto operator()(State const& state) const
    {
      return operation::plus(state.position,
                             operation::times(domain_dimensions,
                                              state.periodicity));
    }
  };
  
  
  // Get dd component of absolute position given state
  // from state.position and periodicity information
  // about the unit cell periodic position is in
  // Unit cell faces must be perpendicular to (cartesian) coordinate axes
  template <std::size_t dd>
  struct Get_position_periodic_component
  {
    std::vector<double> domain_dimensions;
    
    // Construct given size of domain along each dimension
    Get_position_periodic_component(std::vector<double> domain_dimensions)
    : domain_dimensions{ domain_dimensions }
    {}
    
    template <typename State>
    auto operator()(State const& state) const
    {
      return operation::project<dd>(
        operation::plus(state.position,
                        operation::times(domain_dimensions,
                                         state.periodicity)));
    }
  };
  
  // Get projection of absolute position along a direction given state
  // from state.position and periodicity information
  // about the unit cell periodic position is in
  // Unit cell faces must be perpendicular to (cartesian) coordinate axes
  struct Get_position_periodic_projection
  {
    std::vector<double> domain_dimensions;
    std::vector<double> basis;
    
    // Construct given size of domain along each dimension
    // and direction to project onto
    Get_position_periodic_projection
    (std::vector<double> domain_dimensions, std::vector<double> const& direction)
    : domain_dimensions{ domain_dimensions }
    , basis{ operation::div_scalar(direction, operation::abs(direction)) }
    {}
    
    template <typename State>
    auto operator()(State const& state) const
    {
      return operation::dot(
        operation::plus(state.position,
                        operation::times(domain_dimensions,
                                         state.periodicity)),
        basis);
    }
  };
  
  // Get state.tag given state
  struct Get_tag
  {
    template <typename State>
    auto operator()(State const& state) const
    { return state.tag; }
  };
  
  // Get new state given new and old state
  template <typename Getter>
  struct Get_new
  {
    Getter get;
    
    Get_new(Getter&& get = {})
    : get{ std::forward<Getter>(get) }
    {}
    
    template <typename State>
    auto operator()(State const& state_new, State const&) const
    { return get(state_new); }
  };
  template <typename Getter> Get_new(Getter&&) -> Get_new<Getter>;
  
  // Get old state given new and old state
  template <typename Getter>
  struct Get_old
  {
    Getter get;
    
    Get_old(Getter&& get = {})
    : get{ std::forward<Getter>(get) }
    {}
    
    template <typename State>
    auto operator()(State const&, State const& state_old) const
    { return get(state_old); }
  };
  template <typename Getter> Get_old(Getter&&) -> Get_old<Getter>;
  
  // Get linearly interpolated position according to velocity
  // in old state given new state and old state
  template <typename Getter = Get_position,
  typename Getter_velocity = Get_velocity>
  struct Get_position_interp_velocity
  {
    double time;
    Getter get;
    Getter_velocity get_velocity;

    // Construct given current time,
    // getter for position given state and getter for velocity given state
    Get_position_interp_velocity
    (double time, Getter&& get = {}, Getter_velocity&& get_velocity = {})
    : time{ time }
    , get{ std::forward<Getter>(get) }
    , get_velocity{ get_velocity }
    {}

    template <typename State>
    auto operator()
    (State const& state_new, State const& state_old) const
    {
      double delta_t = time-state_old.time;
      auto vel = get_velocity(state_old.velocity);
      return operation::plus(get(state_old.position), vel*delta_t);
    }
  };
  template <typename Getter, typename Getter_velocity>
  Get_position_interp_velocity
  (double, Getter&&, Getter_velocity&&) ->
  Get_position_interp_velocity<Getter, Getter_velocity>;

  // Get linearly interpolated position according to function of state.velocity
  // in old state given new state and old state
  template <typename Getter = Get_time,
  typename VelocityMapper = useful::Forward<double>>
  struct Get_time_interp_velocity
  {
    double position;
    Getter get;
    VelocityMapper velocity_mapper;

    // Construct given current time,
    // getter for position given state and map of state.velocity
    // to actual velocity given state
    Get_time_interp_velocity
    (double position, Getter&& get = {},
     VelocityMapper&& velocity_mapper = {})
    : position{ position }
    , Getter{ std::forward<Getter>(get) }
    , velocity_mapper{ std::forward<VelocityMapper>(velocity_mapper) }
    {}

    template <typename State>
    auto operator()
    (State const& state_new, State const& state_old) const
    {
      double delta_x = position-state_old.position;
      double vel = velocity_mapper(state_old.velocity);
      return state_old.time+delta_x/vel;
    }
  };
  template <typename Getter, typename VelocityMapper>
  Get_time_interp_velocity
  (double, Getter&&, VelocityMapper&&) ->
  Get_time_interp_velocity<Getter, VelocityMapper>;
}

#endif /* StateGetter_h */
