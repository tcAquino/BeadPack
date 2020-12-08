//
//  PTRW.h
//  CTRW
//
//  Created by Tomas Aquino on 1/9/19.
//  Copyright © 2019 Tomas Aquino. All rights reserved.
//

#ifndef PTRW_h
#define PTRW_h

#include <list>

namespace ctrw
{
  // Interface for fixed-time-step steps PTRW using CTRW object framework
	template <typename CTRW, typename Transitions, typename Time_t = double>
	class PTRW
	{
	public:
		using Particle = typename CTRW::Particle;    // Type of particle
    using State = typename CTRW::State;          // Type of particle state
		using Container = typename CTRW::Container;  // Type of particle container
    using Time = Time_t;                         // Type of time variable

    // Construct with given CTRW, transition handler, and initial time
    // time step obtained from transitions
    PTRW(CTRW& ctrw, Transitions& transitions, Time time = 0)
    : ctrw{ ctrw }
    , transitions{ transitions }
    , current_time{ time }
		{ time_step(transitions.timestep()); }
    
    // Construct with given CTRW, transition handler,
    // time step, and initial time
    PTRW(CTRW& ctrw, Transitions& transitions, Time dt, Time time)
    : ctrw{ ctrw }
    , transitions{ transitions }
    , current_time{ time }
    , dt{ dt }
    {}

    // Add copy of particle
		void push_back(Particle particle)
		{ ctrw.push_back(particle); }

    // Remove particles satisfying a criterium
		template <typename Criterium>
		void remove(Criterium criterium)
		{ ctrw.remove(criterium); }

    // Remove particles at given positions in container
		template <typename IntegerType>
		void remove(std::list<IntegerType>& list)
		{ ctrw.remove(list); }
    
    // Remove all particles
    void clear()
    { ctrw.clear(); }

    // Set time step
		void step(Time dt)
		{
			transitions.timestep(dt);
			ctrw.step(transitions);
      current_time += dt;
		}

    // Evolve all particles one time step
    void step()
    {
      ctrw.step(transitions);
      current_time += dt;
    }
    
    // Set time step
    void time_step(Time dt)
    {
      this->dt = dt;
      transitions.time_step(dt);
    }

    // Step each particle directly up to time_max
    // Requires particle time to exist and be updated
    void evolve_by_particle(Time time_max)
    {
      ctrw.evolve_time(time_max, transitions);
      current_time = time_max;
    }
    
    // Step all particles step by step up to time_max
    void evolve(Time time_max)
    {
      while (current_time < time_max)
        step();
      current_time = time_max;
    }

    // Get current time
    double time() const
    { return current_time; }

    // Get reference to particle container
		auto const& particles() const
		{ return ctrw.particles(); }

    // Get reference to particle
		auto const& particles(std::size_t part) const
		{ return ctrw.particles(part); }
    
    // Iterator to beginning of particle container
    auto cbegin() const
    { return ctrw.cbegin(); }

    // Iterator to end of particle container
    auto cend() const
    { return ctrw.cend(); }
    
    // Iterator to beginning of particle container
    auto begin() const
    { return ctrw.begin(); }

    // Iterator to end of particle container
    auto end() const
    { return ctrw.end(); }
    
    // Get time step
    void time_step() const
    {
      return dt;
    }
    
    // Get number of particles
    std::size_t size() const
    { return ctrw.size(); }

	private:
		CTRW& ctrw;                // Wrapped CTRW object
		Transitions& transitions;  // CTRW transition handler with fixed time step

    Time current_time;         // System time
    Time dt;                   // Time step
  };
}

#endif /* PTRW_h */
