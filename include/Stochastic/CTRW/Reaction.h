//
//  Reaction.h
//  rCTRW
//
//  Created by Tomas Aquino on 4/30/18.
//  Copyright © 2018 Tomas Aquino. All rights reserved.
//

#ifndef Reaction_CTRW_h
#define Reaction_CTRW_h

#include <cmath>
#include <random>

namespace ctrw
{
	template <int stoichiometry_type>
	class Reaction_MassAction_SingleSpecies
	{
	public:
		const double stoichiometry;
		const double rate;

		Reaction_MassAction_SingleSpecies(double stoichiometry, double rate)
		: stoichiometry(stoichiometry)
		, rate(rate)
		{}

		template <typename Link, typename State>
		void operator() (Link const& link, State& state, double exposure_time) const
		{
			if (link.trap.params.reactivity)
				state.mass = AnalyticalSol(exposure_time, state.mass);
		}

		double AnalyticalSol(double exposure_time, double mass) const;

	private:
		double coeff = rate/std::tgamma(stoichiometry - 1);
		double exponent = 1. / (1. - stoichiometry);
	};

	class Reaction_Decay_Mass
	{
	public:
		const double rate;

		Reaction_Decay_Mass(double rate = 1.)
		: rate(rate)
		{}

		template <typename Link, typename State>
		void operator() (Link const& link, State& state, double exposure_time) const
		{
			state.mass *= std::exp(-link.trap.params.reactivity*rate*exposure_time);
		}
	};

	class Reaction_Decay_Inst
	{
	public:
		template <typename Link, typename State>
		void operator() (Link const& link, State& state, double = 0.)
		{
			if (link.trap.params.reactivity)
				state.reaction_time = state.time;
		}
	};

	class Reaction_Decay_Prob_Rate
	{
    std::mt19937 rng{ std::random_device{}() };
    std::exponential_distribution<double> exp_dist;

	public:
		Reaction_Decay_Prob_Rate(double rate)
		: exp_dist{ rate }
		{}

		template <typename Link, typename State>
		void operator() (Link const& link, State& state, double exposure_time)
		{
			if (link.trap.params.reactivity)
			{
				double time = exp_dist(rng);
				if (time < exposure_time)
					state.reaction_time = state.time + time;
			}
		}
	};

	class Reaction_Decay_Prob_Het
	{
    std::mt19937 rng{ std::random_device{}() };
    std::exponential_distribution<double> exp_dist{ 1. };

	public:
		template <typename Link, typename State>
		void operator() (Link const& link, State& state, double exposure_time)
		{
			double time = exp_dist(rng)/link.trap.params.reactivity;
			if (time < exposure_time)
				state.reaction_time = state.time + time;
		}
	};

	class Reaction_CTRW_Decay_Inst
	{
    std::mt19937 rng{ std::random_device{}() };
    std::bernoulli_distribution reactivity;

	public:
		Reaction_CTRW_Decay_Inst(double reactive_prob)
		: reactivity{ reactive_prob }
		{}

		template <typename State>
		void operator() (State& state, double = 0.)
		{
			if (reactivity(rng))
				state.reaction_time = state.time;
		}
	};

	class Reaction_CTRW_Decay_Prob_Rate
	{
    std::mt19937 rng{ std::random_device{}() };
    std::bernoulli_distribution reactivity;
		std::exponential_distribution<double> exp_dist;

	public:

		Reaction_CTRW_Decay_Prob_Rate(double reactive_prob, double rate)
		: reactivity{ reactive_prob }
		, exp_dist{ rate }
		{}

		template <typename State>
		void operator() (State& state, double exposure_time)
		{
			if (reactivity(rng))
			{
				double time = exp_dist(rng);
				if (time < exposure_time)
					state.reaction_time = state.time + time;
			}
		}
	};

	template <typename Rate_generator>
	class Reaction_CTRW_Decay_Prob_Het
	{
    std::mt19937 rng{ std::random_device{}() };
    std::exponential_distribution<double> exp_dist{ 1. };
		Rate_generator rate_generator;

	public:
		Reaction_CTRW_Decay_Prob_Het(Rate_generator rate_generator)
		: rate_generator{ rate_generator }
		{}

		template <typename State>
		void operator() (State& state, double exposure_time)
		{
			double time = exp_dist(rng)/rate_generator();
			if (time < exposure_time)
				state.reaction_time = state.time + time;
		}
	};

  template <typename Condition>
  class Reaction_CTRW_Decay_Condition
  {
  public:
    const double reaction_rate;
    
    Reaction_CTRW_Decay_Condition(double reaction_rate, Condition condition)
    : reaction_rate{ reaction_rate }
    , condition{ condition }
    {}
    
    template <typename State>
    void operator()(State& state, double exposure_time)
    {
      if (condition(state))
        state.mass *= std::exp(-reaction_rate*exposure_time);
    }
    
  private:
    Condition condition;
  };
  
	template <>
	double Reaction_MassAction_SingleSpecies<-1>::
	AnalyticalSol(double exposure_time, double mass) const
	{
		double val = mass*std::pow(1. + coeff*exposure_time, exponent);
		return val < 0. ? 0. : val;
	}

	template <>
	double Reaction_MassAction_SingleSpecies<0>::
	AnalyticalSol(double exposure_time, double mass) const
	{
		return mass*std::exp( -rate * exposure_time );
	}

	template <>
	double Reaction_MassAction_SingleSpecies<1>::
	AnalyticalSol(double exposure_time, double mass) const
	{
		return mass*std::pow(1. + coeff*exposure_time, exponent);
	}
  
  // Decay grid concentration
  // Second-order rate based on particle number concentration and grid concentration
  class Reaction_GridConcentrationDecay_ParticleNr_GridConcentration
  {
  public:
    double rate;
    double time_step;
    
    Reaction_GridConcentrationDecay_ParticleNr_GridConcentration
    (double rate, double time_step)
    : rate{ rate }
    , time_step{ time_step }
    {}
    
    template <typename Concentration, typename GridParticles, typename Grid>
    void operator()
    (Concentration& concentration, GridParticles const& particle_grid, Grid const& grid)
    {
      for (auto const& val : particle_grid.particles_by_grid_idx())
      {
        std::size_t grid_idx = val.first;
        double concentration_particles = val.second.size()/grid.cell_volume(grid_idx);
        concentration[grid_idx] *= std::exp(-time_step*rate*concentration_particles);
      }
    }
  };
}

#endif /* Reaction_CTRW_h */
