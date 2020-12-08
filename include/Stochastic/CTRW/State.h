//
//  State.h
//  CTRW
//
//  Created by Tomas Aquino on 9/26/17.
//  Copyright © 2017 Tomas Aquino. All rights reserved.
//

#ifndef State_h
#define State_h

#include <valarray>
#include <vector>
#include "general/Operations.h"
#include "general/useful.h"

// States keep track of different quantities
// Some unwanted quantity types can be set to useful::Empty

namespace ctrw
{
  // State for positions including information about periodicity
  // Periodicity information keeps track of
  // the periodic cell the position is in
  // Defines position, periodicity, mass, time, and tag
  template
  <typename Position_t, typename Periodicity_t = std::vector<int>, typename Mass_t = double,
  typename Time_t = double, typename Tag_t = useful::Empty>
  struct State_periodic
  {
    using Position = Position_t;
    using Periodicity = Periodicity_t;
    using Mass = Mass_t;
    using Time = Time_t;
    using Tag = Tag_t;
    
    State_periodic()
    {}
    
    State_periodic
    (Position position, Periodicity periodicity, Mass mass = 1,
     Time_t time = 0, Tag_t tag = {})
    : position{ position }
    , periodicity{ periodicity }
    , mass{ mass }
    , time{ time }
    , tag{ tag }
    {}

    Position position;
    Periodicity periodicity;
    Mass mass;
    Time time;
    Tag tag;
  };
  
  // State defining position, velocity, time, and tag
  template
  <typename Position_t, typename Velocity_t,
  typename Time_t = double, typename Tag_t = useful::Empty>
  struct State_velocity
  {
    using Position = Position_t;
    using Velocity = Velocity_t;
    using Time = Time_t;
    using Tag = Tag_t;
    
    State_velocity()
    {}
    
    State_velocity
    (Position_t position, Velocity_t velocity,
     Time_t time = {}, Tag_t tag = {})
    : position(position)
    , velocity(velocity)
    , time(time)
    , tag(tag)
    {}

    Position position;
    Velocity velocity;
    Time time;
    Tag tag;
  };
  
  // State defining position, velocity, mass, time, and tag
  template
  <typename Position_t, typename Velocity_t, typename Mass_t = double,
  typename Time_t = double, typename Tag_t = useful::Empty>
  struct State_velocity_mass
  {
    using Position = Position_t;
    using Velocity = Velocity_t;
    using Mass = Mass_t;
    using Time = Time_t;
    using Tag = Tag_t;
    
    State_velocity_mass()
    {}
    
    State_velocity_mass
    (Position position, Velocity velocity, Mass mass = {},
     Time time = {}, Tag tag = {})
    : position{ position }
    , velocity{ velocity }
    , mass{ mass }
    , time{ time }
    , tag{ tag }
    {}

    Position position;
    Velocity velocity;
    Mass mass;
    Time time;
    Tag tag;
  };

  // State defining position, time, and tag
  template <typename Position_t, typename Time_t = double, typename Tag_t = useful::Empty>
  struct State_position
  {
    using Position = Position_t;
    using Time = Time_t;
    using Tag = Tag_t;

    State_position()
    {}

    State_position
    (Position position, Time time = {}, Tag tag = {})
    : position(position)
    , time(time)
    {}

    Position position{};
    Time time{};
    Tag tag{};
  };
  
  // State defining position, mass, time, and tag
  template <typename Position_t, typename Mass_t = double, typename Time_t = double, typename Tag_t = useful::Empty>
  struct State_position_mass
  {
    using Position = Position_t;
    using Mass = Mass_t;
    using Time = Time_t;
    using Tag = Tag_t;

    State_position_mass()
    {}

    State_position_mass
    (Position position, Mass mass = 1., Time time = {}, Tag tag = {})
    : position{ position }
    , mass{ mass }
    , time{ time }
    , tag{ tag }
    {}

    Position position;
    Mass mass;
    Time time;
    Tag tag;
  };

  // State defining position
  template <typename Position_t>
  struct State_PTRW_position
  {
    using Position = Position_t;

    State_PTRW_position()
    {}

    State_PTRW_position(Position position)
    : position(position)
    {}

    Position position{};
  };

  // State defining position and type
  template <typename Position_t>
  struct State_PTRW_position_type
  {
    using Position = Position_t;

    State_PTRW_position_type()
    {}

    State_PTRW_position_type(std::size_t type)
    : type(type)
    {}

    State_PTRW_position_type
    (Position position, std::size_t type = 0)
    : position(position)
    , type(type)
    {}

    Position position{};
    std::size_t type{0};
  };

  // State defining node (with position) and time
  template <typename Node_t>
  struct State_Node
  {
    using Node = Node_t;
    
    Node const* node;
    double time{ 0. };

    auto get_position() const
    { return node->position; }
  };

  // State defining node (with position),
  // time, and mass
  template <typename Node_t>
  struct State_Node_Mass
  {
    using Node = Node_t;
    
    Node const* node;
    double time{ 0. };
    double mass{ 1. };

    auto get_position() const
    { return node->position; }
  };

  // State defining node (with position),
  // time, and reaction_time
  template <typename Node_t>
  struct State_Node_ReactionTime
  {
    using Node = Node_t;
    
    Node const* node;
    double time{ 0. };
    double reaction_time{ -1. };

    auto get_position() const
    { return node->position; }
  };

  // State defining node (with position),
  // position, and time
  template <typename Node_t, typename Position_t>
  struct State_Node_Position
  {
    using Node = Node_t;
    using Position = Position_t;
    using Time = double;
    Node const* node;
    Position position;
    double time;

    State_Node_Position(Node const* node, double time = 0.)
    : node(node)
    , position(node->position)
    , time(time)
    {}

    Position get_position() const
    { return position; }
  };

  // State defining position, time, and reaction_time
  template <typename Position_t>
  struct State_Position_ReactionTime
  {
    using Position = Position_t;
    using Time = double;

    State_Position_ReactionTime()
    {}

    State_Position_ReactionTime(Position position, Time time = {})
    : position(position)
    , time(time)
    {}
    
    Position position{};
    Time time{};
    Time reaction_time{ -1. };
  };
  
  // State defining position, orientation, run, time, and tag
  template <typename Orientation_t, typename Tag_t = useful::Empty>
  struct State_RunTumble
  {
    using Position = std::vector<double>;
    using Orientation = Orientation_t;
    using Time = double;
    using Tag = Tag_t;
    
    State_RunTumble()
    {}
    
    State_RunTumble
    (Position position, Orientation orientation, bool run = 0,
     Time time = {}, Tag tag = {})
    : position{ position }
    , orientation{ orientation }
    , run{ run }
    , time{ time }
    , tag{ tag }
    {}
    
    Position position{};
    Orientation orientation{};
    bool run;
    Time time{};
    Tag tag;
  };
  
  // State defining position, orientation, state, and tag
  template <typename Orientation_t = double,
  typename Tag_t = std::size_t,
  typename Position_t = std::vector<double>>
  struct State_RunTumble_PTRW
  {
    using Position = Position_t;
    using Orientation = Orientation_t;
    using Tag = Tag_t;
    
    State_RunTumble_PTRW()
    {}
    
    State_RunTumble_PTRW
    (Position position, Orientation orientation,
     int state, Tag tag = {})
    : position{ position }
    , orientation{ orientation }
    , state{ state }
    , tag{ tag }
    {}
    
    Position position{};
    Orientation orientation{};
    int state;
    Tag tag;
  };
}

#endif /* State_h */
