//
//  CTRW.h
//
//  Created by Tomas Aquino on 9/26/17.
//  Copyright © 2017 Tomas Aquino. All rights reserved.
//

#ifndef CTRW_h
#define CTRW_h

#include <functional>
#include <list>
#include <vector>
#include "Particle.h"
#include "general/useful.h"

namespace ctrw
{
  // Handles a set of particles characterized by a state
  // Particles undergo transitions handled by a separate class
  // passed to the evolution methods
  // State_t characterizes each particle's state
  // Note:
  //    Order of particles in container is not preserved by removal methods
  //    Iterators and references can be invalidated when removing or adding
  template <typename State_t>
  class CTRW
  {
  public:
    using State = State_t;                   // Particle state
    using Particle = ctrw::Particle<State>;  // Particle type
    using Container = std::vector<Particle>; // Set of particles
    struct Tag{};                            // To select constructors
                                             // that tag particles

    // Construct empty
    CTRW()
    {}

    // Construct given particles
    CTRW(Container particles)
    : particle_container{ particles }
    {}
    
    // Construct given particles and tag them in ascending order
    // State must define: tag
    CTRW(Container particles, Tag)
    : CTRW{ particles }
    { retag(); }
    
    // Construct given particles and reserve space for given maximum number
    CTRW(Container particles, std::size_t max_nr_particles)
    {
      particle_container.reserve(max_nr_particles);
      particle_container = particles;
    }
    
    // Construct given particles and reserve space for given maximum number
    // Tag them in ascending order
    // State must define: tag
    CTRW(Container particles, std::size_t max_nr_particles, Tag)
    : CTRW(particles, max_nr_particles)
    { retag(); }

    // Handle nr_particles particles
    // Each particle's state is produced by StateMaker
    template <typename StateMaker>
    CTRW(std::size_t nr_particles, StateMaker state_maker)
    { push_back(nr_particles, state_maker); }

    // Add nr_particles particles
    // Each particle's state is produced by StateMaker
    template <typename StateMaker>
    void push_back(std::size_t nr_particles, StateMaker state_maker)
    {
      particle_container.reserve(nr_particles);
      std::generate_n(std::back_inserter(particle_container), nr_particles, state_maker);
    }

    // Add copy of particle
    void push_back(Particle particle)
    { particle_container.push_back(particle); }

    // Add copies of particles in a container
    void push_back(Container particles)
    {
      for (auto const& particle : particles)
        particle_container.push_back(particle);
    }
    
    // Reserve space for nr_particles particles
    void reserve(std::size_t nr_particles)
    { particle_container.reserve(nr_particles); }
    
    // Tag particles in order of container
    void retag()
    {
      std::size_t pp = 0;
      for (auto& particle : particle_container)
      {
        particle.transform([pp](State& state){ state.tag = pp; });
        ++pp;
      }
    }
    
    // Set particle's new and old state to given state
    void set(std::size_t part, State const& state)
    { particle_container[part].set(state); }

    // Remove particles satisfying a criterium
    template<typename Criterium>
    void remove(Criterium criterium)
    {
      for (std::size_t part = size(); part --> 0;)
        remove(part, criterium);
    }

    // Remove particles at given positions in container
    template <typename IntegerType>
    void remove(std::list<IntegerType>& list)
    {
      list.sort(std::greater<IntegerType>{});
      for (auto const& part : list)
        remove(part);
    }

    // Remove all particles
    void clear()
    { particle_container.clear(); }

    // Particles make transitions until their time is >= time_to
    template<typename Transitions_Particle>
    void evolve_time
    (double time_to, Transitions_Particle& transitions_particle)
    {
      evolve
      ([time_to](Particle const& part)
       { return part.state_new().time < time_to; }, transitions_particle);
    }

    // Particles make transitions until their time is >= time_to
    template<typename Transitions_Particle>
    void evolve_space
    (double length_to, Transitions_Particle& transitions_particle)
    {
      evolve
      ([length_to](Particle const& part)
      { return part.state_new().position < length_to; }, transitions_particle);
    }

    // Particles make transitions until a given (particle-based) criterium is met
    template<typename Transitions_Particle, typename Criterium>
    void evolve
    (Criterium criterium, Transitions_Particle& transitions_particle)
    {
      for (auto& part : particle_container)
        while(criterium(part))
          part.transition(transitions_particle);
    }

    // Particles make one transition
    template<typename Transitions_Particle>
    void step(Transitions_Particle& transitions_particle)
    {
      for (auto& part : particle_container)
        part.transition(transitions_particle);
    }
    
    // Particles make one transition if criteritum is satisfied
    template<typename Transitions_Particle, typename Criterium>
    void step(Criterium criterium, Transitions_Particle& transitions_particle)
    {
      for (auto& part : particle_container)
        if (criterium(part))
          part.Transition(transitions_particle);
    }

    // Get reference to particle container
    Container const& particles() const
    { return particle_container; }

    // Get reference to particle
    Particle const& particles(std::size_t part) const
    { return particle_container[part]; }

    // Get number of particles
    std::size_t size() const
    { return particle_container.size(); }

    // Iterator to begining of particle container
    auto cbegin() const
    { return particle_container.cbegin(); }

    // Iterator to end of particle container
    auto cend() const
    { return particle_container.cend(); }
    
    // Iterator to begining of particle container
    auto begin() const
    { return particle_container.cbegin(); }

    // Iterator to end of particle container
    auto end() const
    { return particle_container.cend(); }
    
  private:
    Container particle_container;  // Current particles

    // Remove particle by position in container if criterium is met
    template <typename Criterium>
    void remove(std::size_t part, Criterium criterium)
    { if (criterium(particle_container[part])) remove(part); }
    
    // Remove particle by position in container
    void remove(std::size_t part)
    { useful::swap_erase(particle_container, part); }
  };
}

#endif /* CTRW_h */
