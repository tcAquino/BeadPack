//
//  ParticleCollections.h
//  BeadPack
//
//  Created by Tomás Aquino on 18/05/2020.
//  Copyright © 2020 Tomás Aquino. All rights reserved.
//

#ifndef ParticleCollections_h
#define ParticleCollections_h

#include <iterator>
#include <list>
#include "general/Operations.h"
#include "general/useful.h"
#include "Stochastic/CTRW/StateGetter.h"

namespace ctrw
{
  // Collection of particles with each particle connected to two others
  // Note: For a closed strip (loop), repeat the 1st particle at the end
  // Warning: Requires particle tags to match position in CTRW particle container.
  // Throughout:
  //   Adjust objects can modify added states on refinement
  //   StateMaker objects return a base particle state
  //   Getter_position objects return a position given a particle
  class Strip
  {
  private:
    using Container = std::list<std::size_t>;  // Type of particle tag container
    using iterator = Container::iterator;      // Type of iterator to particle tag container
    
  public:
    using const_iterator = Container::const_iterator;
    
    // Construct from particle tags, maximum allowed nr of particles,
    // and maximum segment length to trigger refining
    template <typename Container = std::vector<std::size_t>>
    Strip
    (Container const& particle_tags_initial,
     std::size_t max_particles, double max_segment_length)
    : max_particles{ max_particles }
    , max_segment_length{ max_segment_length }
    {
      for (auto const& tag : particle_tags_initial)
        particle_tags.push_back(tag);
    }
    
    // Add particles at the midpoint between connected particles until
    // below segment length,
    // until strip has maximum allowed particles
    template
    <typename CTRW, typename StateMaker,
    typename Getter_position = ctrw::Get_new_from_particle<ctrw::Get_position>,
    typename Adjust = useful::DoNothing>
    void resize
    (CTRW& ctrw, StateMaker state_maker,
     Getter_position get_position = {}, Adjust adjust = {})
    {
      for (auto it = std::next(begin()); it != end();)
      {
        if (full())
          return;
        if (segment_length_sq(it, ctrw, get_position)
            > max_segment_length*max_segment_length)
        {
          insert_midpoint(it, ctrw, state_maker,
                          get_position, adjust);
          std::advance(it, -1);
        }
        else
          ++it;
      }
    }
    
    // Add and connect particles' previous state until
    // below segment length,
    // until strip has maximum allowed particles
    template
    <typename CTRW, typename StateMaker,
    typename Getter_position = ctrw::Get_new_from_particle<ctrw::Get_position>,
    typename Adjust = useful::DoNothing>
    void refine_with_old_state
    (CTRW& ctrw, StateMaker state_maker,
     Getter_position get_position = {}, Adjust adjust = {})
    {
      for (auto it = std::next(begin()); it != end(); ++it)
      {
        if (full())
          return;
        if (segment_length_sq(it, ctrw, get_position)
            > max_segment_length*max_segment_length)
          insert_previous_state(it, ctrw, get_position,
                                adjust, state_maker);
      }
    }
    
    // Add and connect particles' previous state until
    // below segment length,
    // until strip has goal size or maximum allowed particles
    template
    <typename CTRW, typename StateMaker,
       typename Getter_position = ctrw::Get_new_from_particle<ctrw::Get_position>,
       typename Adjust = useful::DoNothing>
    void refine
    (std::size_t goal_size, CTRW& ctrw,
     StateMaker state_maker,
     Getter_position get_position = {}, Adjust adjust = {})
    {
      while (1)
        for (auto it = std::next(begin()); it != end(); ++it)
        {
          if (full() || size() >= goal_size)
            return;
          insert_midpoint(it, ctrw, state_maker,
                          get_position, adjust);
        }
    }
    
    // Get full length of strip as sum of segment lengths
    template <typename CTRW,
    typename Getter_position = ctrw::Get_new_from_particle<ctrw::Get_position>>
    double length
    (CTRW const& ctrw, Getter_position get_position = {}) const
    {
      double len = 0.;
      for (auto it = std::next(cbegin()); it != cend(); ++it)
        len += segment_length(it, ctrw, get_position);
      
      return len;
    }
    
    // Get number of particles in strip
    std::size_t size() const
    { return particle_tags.size(); }
    
    // Get maximum allowed number of particles
    std::size_t maximum_size() const
    { return max_particles; }
    
    // Check if maximum allowed number of particles has been reached
    bool full() const
    { return size() >= maximum_size(); }
    
    // Iterator to beginning of particle tag container
    const_iterator cbegin() const
    { return particle_tags.cbegin(); }
    
    // Iterator to end of particle tag container
    const_iterator cend() const
    { return particle_tags.cend(); }
    
  private:
    std::size_t max_particles;  // Maximum allowed particles in strip
    double max_segment_length;  // Segment lenght above which refinement triggers
    Container particle_tags;    // Container of particle tags in strip, in connected order
    
    // Length of the straight line between two connected particles
    // get_position returns a position from a particle
    template <typename CTRW,
    typename Getter_position = ctrw::Get_new_from_particle<ctrw::Get_position>>
    double segment_length
    (iterator it, CTRW const& ctrw,
     Getter_position get_position = {}) const
    {
      return operation::abs(operation::minus(
        get_position(ctrw.particles(*it)),
        get_position(ctrw.particles(*std::prev(it)))));
    }
    
    // Length squared of the straight line between two connected particles
    template <typename CTRW,
    typename Getter_position = ctrw::Get_new_from_particle<ctrw::Get_position>>
    double segment_length_sq
    (iterator it, CTRW const& ctrw,
     Getter_position get_position = {}) const
    {
      return operation::abs_sq(operation::minus(
        get_position(ctrw.particles(*it)),
        get_position(ctrw.particles(*std::prev(it)))));
    }
    
    // Get midpoint on the straight line between two connected particles
    template <typename CTRW,
    typename Getter_position = ctrw::Get_new_from_particle<ctrw::Get_position>>
    auto segment_midpoint
    (iterator it, CTRW const& ctrw,
     Getter_position get_position = {}) const
    {
      auto midpoint =
        operation::plus(get_position(ctrw.particles(*it)),
                        get_position(ctrw.particles(*std::prev(it))));
      operation::times_scalar_InPlace(0.5, midpoint);
      
      return midpoint;
    }
    
    // Insert and connect particle at midpoint between two particles
    template
    <typename CTRW,
    typename StateMaker,
    typename Getter_position = ctrw::Get_new_from_particle<ctrw::Get_position>,
    typename Adjust>
    void insert_midpoint
    (iterator it, CTRW& ctrw, StateMaker state_maker,
     Getter_position get_position = {}, Adjust adjust = {})
    {
      auto state{ state_maker() };
      state.position = segment_midpoint(it, ctrw, get_position);
      state.tag = ctrw.size();
      adjust(state);
      
      ctrw.push_back(Particle{ state });
      particle_tags.insert(it, state.tag);
    }
    
    // Insert and connect particle with particle's previous state
    template
    <typename CTRW, typename Getter_position,
    typename Adjust, typename StateMaker>
    void insert_previous_state
    (iterator it, CTRW& ctrw, Getter_position get_position,
     Adjust adjust, StateMaker state_maker)
    {
      particle_tags.insert(it, ctrw.size());
      ctrw.push_back(Particle{ ctrw.particles(*it).State_old() });
    }
    
    // Iterator to beginning of particle tag container
    iterator begin()
    { return particle_tags.begin(); }
    
    // Iterator to end of particle tag container
    iterator end()
    { return particle_tags.end(); }
  };
  
  // Handles strips of particles associated with a CTRW object
  template
  <typename CTRW,
  typename Getter_position = ctrw::Get_new_from_particle<Get_position>,
  typename Adjust = useful::DoNothing,
  typename StateMaker = useful::Maker<typename CTRW::State>>
  struct StripHandler
  {
    CTRW& ctrw;                    // CTRW object that handles particle dynamics
    Getter_position get_position;  // Return a position given a particle
    Adjust adjust;                 // Modify added states on refinement
    StateMaker state_maker;        // Return a base particle state
    std::vector<Strip> strips;     // Container of strips to handle
    
    // Construct given the ctrw, position getter, adjuster, and state maker
    StripHandler
    (CTRW& ctrw,
     Getter_position get_position = {},
     Adjust adjust = {}, StateMaker state_maker = {})
    : ctrw{ ctrw }
    , get_position{ get_position }
    , adjust{ adjust }
    , state_maker{ state_maker }
    {}
    
    // Reserve container space for strips
    void reserve(std::size_t nr_strips)
    { strips.reserve(nr_strips); }
    
    // Get number of strips
    std::size_t size() const
    { return strips.size(); }
    
    // Add a strip given particle tags,
    // maximum allowed number of particles,
    // and maximum segment length to trigger refinemenet
    template <typename Container = std::vector<std::size_t>>
    void push_back
    (Container const& particle_tags_initial,
     std::size_t max_particles, double max_segment_length)
    {
      strips.emplace_back(particle_tags_initial,
                          max_particles, max_segment_length);
    }
    
    // Add particles to all strips until goal_size
    // number of particles, or full, or
    // no strips exceed maximum size,
    // adding midpoints between particles
    void refine(std::size_t goal_size)
    {
      for (auto& strip : strips)
        strip.refine(goal_size, ctrw, state_maker,
                     get_position, adjust);
    }
    
    // Add particles to a strip until goal_size
    // number of particles, or full, or
    // no strips exceed maximum size,
    // adding midpoints between particles
    void refine(std::size_t ss, std::size_t goal_size)
    {
      strips[ss].refine(goal_size, ctrw, state_maker,
                        get_position, adjust);
    }
    
    // Add particles to a strip until full or
    // no strips exceed maximum size,
    // adding midpoints between particles
    void resize(std::size_t ss)
    {
      strips[ss].resize(ctrw, state_maker,
                        get_position, adjust);
    }
    
    // Add particles to all strips until full or
    // no strips exceed maximum size,
    // adding midpoints between particles
    void resize()
    {
      for (auto& strip : strips)
        strip.resize(ctrw, state_maker,
                     get_position, adjust);
    }
    
    // Add particles to a strips using previous particle states
    void refine_with_old_state(std::size_t ss)
    {
      strips[ss].refine_with_old_state(ctrw, state_maker,
                                       get_position, adjust);
    }
    
    // Add particles to all strips using previous particle states
    void refine_with_old_state()
    {
      for (auto& strip : strips)
        strip.refine_with_old_state(ctrw, state_maker,
                                    get_position, adjust);
    }
    
    // Get the number of particles in a strip
    std::size_t size(std::size_t ss) const
    { return strips[ss].size(); }
    
    // Get the length of a strip as a sum of segment lengths
    double length(std::size_t ss) const
    { return strips[ss].length(ctrw, get_position); }
    
    // Get the maximum allowed particles in a strip
    std::size_t maximum_size(std::size_t ss) const
    { return strips[ss].maximum_size(); }
    
    // Check if strip has maximum allowed particles
    bool full(std::size_t ss) const
    { return strips[ss].full(); }
    
    // Check if all strips have maximum allowed particles
    bool full() const
    {
      for (auto strip : strips)
        if (!strip.full())
          return 0;
      return 1;
    }
    
    // Iterator to beginning of strip container
    auto cbegin(std::size_t ss) const
    { return strips[ss].cbegin(); }
    
    // Iterator to end of strip container
    auto cend(std::size_t ss) const
    { return strips[ss].cend(); }
  };
}

#endif /* ParticleCollections_h */
