//
//  Reaction.h
//  CTRW
//
//  Created by Tomas Aquino on 4/30/18.
//  Copyright © 2018 Tomas Aquino. All rights reserved.
//

#ifndef Reaction_CTRW_h
#define Reaction_CTRW_h

#include <cmath>
#include <random>
#include <vector>
#include "Geometry/Coordinates.h"

namespace ctrw
{
  // Single-species mass action reaction stoichiometry A -> \varnothing
  // using analytical solution
  // stoiochiometry_type should be set to -1, 0, 1
  // for stoiohiometry < 1, stoichiometry = 1, and stoichiometry > 1,
  // respectively, to implement appropriate analytical solution
  // at compile time
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
  
  // Analytical solution for stoichiometry A -> \varnothing
  // for stoichiometry < 1, for which mass becomes 0 in finite time
  template <>
  double Reaction_MassAction_SingleSpecies<-1>::
  AnalyticalSol(double exposure_time, double mass) const
  {
    double val = mass*std::pow(1. + coeff*exposure_time, exponent);
    return val < 0. ? 0. : val;
  }

  // Analytical solution for A -> \varnothing
  template <>
  double Reaction_MassAction_SingleSpecies<0>::
  AnalyticalSol(double exposure_time, double mass) const
  {
    return mass*std::exp( -rate * exposure_time );
  }

  // Analytical solution for stoichiometry A -> \varnothing
  // for stoichiometry > 1
  template <>
  double Reaction_MassAction_SingleSpecies<1>::
  AnalyticalSol(double exposure_time, double mass) const
  {
    return mass*std::pow(1. + coeff*exposure_time, exponent);
  }

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

  // Instantaneous decay along links between trips
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

  // Exponential decay along links between traps
  // with quenched reactive or non-reactive traps
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

  // Exponential decay along links between traps
  // with quenched heterogeneous trap reaction rate
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

  // Instantaneous decay with annealed reactive or non-reactive steps
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

  // Exponential decay with annealed reactive or non-reactive steps
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

  // Exponential decay with annealed reaction rate in each step
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

  // Exponential decay when state meets condition
  // Condition must take a state and returns true if condition
  // for reaction is met, and false otherwise
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
  
  // Reaction to track reaction prevalence on bead surfaces
  // in a beadpack
  template <typename BeadPack>
  class Reaction_Decay_Condition_Map_Beadpack_3d
  {
  public:
    // Construct for given reaction rate, discretization length for reaction,
    // number of bins to equally discretize [0, 2pi] for reaction map angles,
    // and beadpack
    Reaction_Decay_Condition_Map_Beadpack_3d
    (double reaction_rate, double length_discretization,
     std::size_t nr_bins, BeadPack const& bead_pack)
    : reaction_rate{ reaction_rate }
    , length_discretization{ length_discretization }
    , nr_bins_phi{ nr_bins }
    , nr_bins_theta{ std::size_t(nr_bins/2) }
    , bead_pack{ bead_pack }
    {
      // Equally-discretized angle bins for azimuthal and elevation angles phi and theta
      // Note: the solid angle associated with each bin is not constant
      map_phi_theta.assign(nr_bins,
        std::vector<double>(std::size_t(nr_bins/2), 0.));
    }
    
    // Consume state.mass during exposure time according to reaction
    // Store amount of reaction in appropriate angle bin
    template <typename State>
    void operator()(State& state, double exposure_time)
    {
      auto near = bead_pack.near(state.position, length_discretization);
      if (near.first)
      {
        double old_mass = state.mass;
        state.mass *= std::exp(-reaction_rate*exposure_time);
        
        std::size_t bead = bead_pack.nearest_neighbor(state.position).first;
        auto spherical = geometry::cartesian2spherical(
          operation::minus(state.position, bead_pack.center(bead)));
        std::size_t bin_phi = (spherical[1]+constants::pi)/(2.*constants::pi)*nr_bins_phi;
        if (bin_phi == nr_bins_phi)
          --bin_phi;
        std::size_t bin_theta = spherical[2]/constants::pi*nr_bins_theta;
        if (bin_theta == nr_bins_theta)
          --bin_theta;
        map_phi_theta[bin_phi][bin_theta] += old_mass - state.mass;
      }
    }
    
    // Print angle values and corresponding amount of mass consumed to file
    // file columns: phi theta consumption
    void print_map
    (std::string const& filename,
     int precision = 8, std::string delimiter = "\t")
    {
      std::ofstream output{ filename };
      if (!output.is_open())
        throw useful::open_write_error(filename);
      output << std::scientific << std::setprecision(precision);
      
      for (std::size_t ii = 0; ii < nr_bins_phi; ++ii)
      {
        double phi = 2.*constants::pi*(ii+0.5)/nr_bins_phi - constants::pi;
        for (std::size_t jj = 0; jj < nr_bins_theta; ++jj)
        {
          double theta = constants::pi*(jj+0.5)/nr_bins_theta;
          output << phi << delimiter
                 << theta << delimiter
                 << map_phi_theta[ii][jj] << "\n";
        }
      }
    }
    
  private:
    const double reaction_rate;          // Reaction rate in discretization reactive region
    const double length_discretization;  // Reaction occurs within this length of bead surfaces
    std::size_t nr_bins_phi;             // Number of bins to equally discretize azimuthal angle
    std::size_t nr_bins_theta;           // Number of bins to equally discretize elevation angle
    BeadPack const& bead_pack;           // Beadpack whose bead surfaces are reactive
    std::vector<std::vector<double>>
      map_phi_theta;                     // Mass consumed at each angle bin
  };
}

#endif /* Reaction_CTRW_h */
