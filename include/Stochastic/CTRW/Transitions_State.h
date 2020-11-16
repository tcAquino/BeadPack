//
//  Transitions_State.h
//  CTRW
//
//  Created by Tomas Aquino on 8/6/18.
//  Copyright © 2018 Tomas Aquino. All rights reserved.
//

#ifndef Transitions_State_h
#define Transitions_State_h

#include <utility>
#include "Boundary.h"
#include "JumpGenerator.h"
#include "TimeGenerator.h"
#include "general/Operations.h"

namespace ctrw
{
	template <typename TimeGenerator, typename JumpGenerator,
  typename Boundary = boundary::DoNothing>
	class Transitions_Time_Position
	{
	public:
		Transitions_Time_Position
    (TimeGenerator&& time_generator, JumpGenerator&& jump_generator,
     Boundary&& boundary = {})
    : time_generator{ std::forward<TimeGenerator>(time_generator) }
    , jump_generator{ std::forward<JumpGenerator>(jump_generator) }
    , boundary{ std::forward<Boundary>(boundary) }
		{}

		template <typename State>
		void operator() (State& state)
		{
      auto state_old = state;
      operation::plus_InPlace(state.position, jump_generator(state));
      boundary(state, state_old);
      state.time += time_generator(state);
		}

	private:
		TimeGenerator time_generator;
		JumpGenerator jump_generator;
    Boundary boundary;
	};
  template <typename TimeGenerator, typename JumpGenerator,
  typename Boundary>
  Transitions_Time_Position
  (TimeGenerator&&, JumpGenerator&&, Boundary&&) ->
  Transitions_Time_Position<TimeGenerator, JumpGenerator, Boundary>;
  
  template <typename JumpGenerator, typename Boundary = boundary::DoNothing>
  class Transitions_Position
  {
  public:
    Transitions_Position
    (JumpGenerator&& jump_generator, Boundary&& boundary = {})
    : jump_generator{ std::forward<JumpGenerator>(jump_generator) }
    , boundary{ std::forward<Boundary>(boundary) }
    {}

    template <typename State>
    void operator()(State& state)
    {
      auto state_old = state;
      operation::plus_InPlace(state.position, jump_generator(state));
      boundary(state, state_old);
    }

  private:
    JumpGenerator jump_generator;
    Boundary boundary;
  };
  template <typename JumpGenerator, typename Boundary>
  Transitions_Position(JumpGenerator&&, Boundary&&) ->
  Transitions_Position<JumpGenerator, Boundary>;

	template <typename FlowField, typename Boundary = boundary::DoNothing>
	class Transitions_PTRW_FlowField_Diff
	{
  private:
    auto make_diff_generators(std::vector<double> const& diff, double dt)
    {
      std::vector<JumpGenerator_Diffusion_1d> diff_generators;
      diff_generators.reserve(diff.size());
      for (auto const& diff_val : diff)
        diff_generators.emplace_back(diff_val, dt);

      return diff_generators;
    }
    
	public:
    Transitions_PTRW_FlowField_Diff
    (std::vector<double> const& diff, double dt,
     FlowField&& flow_field, Boundary&& boundary = {})
    : time_generator{ dt }
    , diff_generators{ make_diff_generators(diff, dt) }
    , flow_field{ std::forward<FlowField>(flow_field) }
    , boundary{ std::forward<Boundary>(boundary) }
		{}
    
		template <typename State>
		void operator()(State& state)
		{
      auto state_old = state;
      std::vector<double> jump(state.position.size());
      for (std::size_t dd = 0; dd < jump.size(); ++dd)
        jump[dd] = diff_generators[dd]();
      jump[0] += flow_field(state.position)*time_generator.time_step();

      for (std::size_t dd = 0; dd < jump.size(); ++dd)
        state.position[dd] += jump[dd];

			boundary(state, state_old);
		}

		void timestep(double time_step)
		{
      time_generator.time_step(time_step);
      for (auto& gen : diff_generators)
        gen.time_step(time_step);
    }

    void diff(std::size_t dd, double diff)
    { return diff_generators[dd].diff(diff); }

    double timestep() const
    { return time_generator.time_step(); }

    double diff(std::size_t dd) const
    { return diff_generators[dd].diff(); }

	private:
    TimeGenerator_Step<double> time_generator;
    std::vector<JumpGenerator_Diffusion_1d> diff_generators;
    FlowField flow_field;
		Boundary boundary;
	};
  template <typename FlowField, typename Boundary>
  Transitions_PTRW_FlowField_Diff
  (std::vector<double> const&, double,
   FlowField&&, Boundary&&) ->
  Transitions_PTRW_FlowField_Diff<FlowField, Boundary>;

  template <typename Boundary = boundary::DoNothing>
  class Transitions_PTRW_Diffusion_1d
  {
  public:
    Transitions_PTRW_Diffusion_1d
    (double diff, double time_step, Boundary&& boundary = {})
    : time_generator{ time_step }
    , jump_generator{ diff, time_step }
    , boundary{ std::forward<Boundary>(boundary) }
    {}

    template <typename State>
    void operator()(State& state)
    {
      auto state_old = state;
      state.position += jump_generator(state);
      boundary(state, state_old);
      state.time += timestep();
    }

    void timestep(double time_step)
    {
      time_generator.time_step(time_step);
      jump_generator.time_step(time_step);
    }

    void diff(double diff)
    { return jump_generator.diff(diff); }

    double timestep() const
    { return time_generator.time_step(); }

    double diff() const
    { return jump_generator.diff(); }

  private:
    TimeGenerator_Step<double> time_generator;
    JumpGenerator_Diffusion_1d jump_generator;
    Boundary boundary;
  };
  template <typename Boundary>
  Transitions_PTRW_Diffusion_1d
  (double, double, Boundary&&) ->
  Transitions_PTRW_Diffusion_1d<Boundary>;
  
  template <typename Transitions_Transport, typename Reaction>
  class Transitions_PTRW_Transport_Reaction
  {
  public:
    Transitions_PTRW_Transport_Reaction
    (Transitions_Transport&& transitions_transport,
     Reaction&& reaction)
    : transitions_transport{
      std::forward<Transitions_Transport>(transitions_transport) }
    , reaction{ std::forward<Reaction>(reaction) }
    , time_step{ transitions_transport.timestep() }
    {}
    
    Transitions_PTRW_Transport_Reaction
    (Transitions_Transport&& transitions_transport,
     Reaction reaction, double time_step)
    : transitions_transport{
      std::forward<Transitions_Transport>(transitions_transport) }
    , reaction{ std::forward<Reaction>(reaction) }
    , time_step{ time_step }
    {}

    template <typename State>
    void operator()(State& state)
    {
      transitions_transport(state);
      reaction(state, time_step);
    }
    
    void timestep(double time_step)
    {
      transitions_transport.timestep(time_step);
      this->time_step = time_step;
    }
    
    double timestep()
    { return time_step; }

  private:
    Transitions_Transport transitions_transport;
    Reaction reaction;
    double time_step;
  };
  template <typename Transitions_Transport, typename Reaction>
  Transitions_PTRW_Transport_Reaction
  (Transitions_Transport&&, Reaction&&) ->
  Transitions_PTRW_Transport_Reaction<Transitions_Transport, Reaction>;
  template <typename Transitions_Transport, typename Reaction>
  Transitions_PTRW_Transport_Reaction
  (Transitions_Transport&&, Reaction&&, double) ->
  Transitions_PTRW_Transport_Reaction<Transitions_Transport, Reaction>;
  
  template <typename Transitions_Transport, typename Reaction>
  class Transitions_CTRW_Transport_Reaction
  {
  public:
    Transitions_CTRW_Transport_Reaction
    (Transitions_Transport&& transitions_transport,
     Reaction&& reaction)
    : transitions_transport{
      std::forward<Transitions_Transport>(transitions_transport) }
    , reaction{ std::forward<Reaction>(reaction) }
    {}
    
    Transitions_CTRW_Transport_Reaction
    (Transitions_Transport&& transitions_transport,
     Reaction reaction, double time_step)
    : transitions_transport{
      std::forward<Transitions_Transport>(transitions_transport) }
    , reaction{ std::forward<Reaction>(reaction) }
    {}

    template <typename State>
    void operator()(State& state)
    {
      auto time_old = state.time;
      transitions_transport(state);
      reaction(state, state.time-time_old);
    }

  private:
    Transitions_Transport transitions_transport;
    Reaction reaction;
  };
  template <typename Transitions_Transport, typename Reaction>
  Transitions_CTRW_Transport_Reaction
  (Transitions_Transport&&, Reaction&&) ->
  Transitions_CTRW_Transport_Reaction<Transitions_Transport, Reaction>;
  template <typename Transitions_Transport, typename Reaction>
  Transitions_CTRW_Transport_Reaction
  (Transitions_Transport&&, Reaction&&, double) ->
  Transitions_CTRW_Transport_Reaction<Transitions_Transport, Reaction>;

	template <typename TimeGenerator,
  typename JumpGenerator, typename Reaction>
	class Transitions_Reaction_Position
	{
	public:
		Transitions_Reaction_Position
		(TimeGenerator&& time_generator,
     JumpGenerator&& jump_generator,
     Reaction reaction)
    : time_generator{
      std::forward<TimeGenerator>(time_generator) }
    , jump_generator{
      std::forward<JumpGenerator>(jump_generator) }
    , reaction{ reaction }
		{}

		template <typename State>
		void operator()(State& state)
		{
			double exposure_time = time_generator();
			reaction(state, exposure_time);
			state.position += jump_generator();
			state.time += exposure_time;
		}

	private:
		TimeGenerator time_generator;
		JumpGenerator jump_generator;
		Reaction reaction;
	};
  template <typename TimeGenerator,
  typename JumpGenerator, typename Reaction>
  Transitions_Reaction_Position
  (TimeGenerator&&, JumpGenerator&&, Reaction&&) ->
  Transitions_Reaction_Position<TimeGenerator, JumpGenerator, Reaction>;
  
  template
  <typename TimeGenerator_run, typename TimeGenerator_tumble,
  typename JumpGenerator, typename OrientationGenerator>
  class Transitions_RunTumble
  {
  public:
    Transitions_RunTumble
    (TimeGenerator_run&& time_generator_run,
     TimeGenerator_tumble&& time_generator_tumble,
     JumpGenerator&& jump_generator,
     OrientationGenerator&& orientation_generator)
    : time_generator_run{
      std::forward<TimeGenerator_run>(time_generator_run) }
    , time_generator_tumble{
      std::forward<TimeGenerator_tumble>(time_generator_tumble) }
    , jump_generator{
      std::forward<JumpGenerator>(jump_generator) }
    , orientation_generator{
      std::forward<OrientationGenerator>(orientation_generator) }
    {}

    template <typename State>
    void operator()(State& state)
    {
      if (state.run)
      {
        state.position += jump_generator(state);
        state.time += time_generator_run(state);
      }
      else
      {
        state.orientation += orientation_generator(state);
        state.time += time_generator_tumble(state);
      }
      state.run = !state.run;
    }

  private:
    TimeGenerator_run time_generator_run;
    TimeGenerator_tumble time_generator_tumble;
    JumpGenerator jump_generator;
    OrientationGenerator orientation_generator;
  };
  template
  <typename TimeGenerator_run, typename TimeGenerator_tumble,
  typename JumpGenerator, typename OrientationGenerator>
  Transitions_RunTumble
  (TimeGenerator_run&&, TimeGenerator_tumble&&,
   JumpGenerator&&, OrientationGenerator&&) ->
  Transitions_RunTumble<TimeGenerator_run, TimeGenerator_tumble,
  JumpGenerator, OrientationGenerator>;
  
  template
  <typename Boundary, typename StateSwitcher,
  typename JumpGenerator, typename OrientationGenerator,
  typename OrientationGenerator_Wall = OrientationGenerator>
  class Transitions_RunAndTumble_PTRW
  {
  public:
    Transitions_RunAndTumble_PTRW
    (JumpGenerator&& jump_generator,
     OrientationGenerator&& orientation_generator,
     OrientationGenerator_Wall&& orientation_generator_wall,
     StateSwitcher&& state_switcher, Boundary&& boundary)
    : boundary{ std::forward<Boundary>(boundary) }
    , state_switcher{
      std::forward<StateSwitcher>(state_switcher) }
    , jump_generator{
      std::forward<JumpGenerator>(jump_generator) }
    , orientation_generator{
      std::forward<OrientationGenerator>(orientation_generator) }
    , orientation_generator_wall{
      std::forward<OrientationGenerator_Wall>(orientation_generator_wall) }
    {}
    
    Transitions_RunAndTumble_PTRW
    (JumpGenerator&& jump_generator,
     OrientationGenerator&& orientation_generator,
     StateSwitcher&& state_switcher, Boundary&& boundary)
    : Transitions_RunAndTumble_PTRW{
      std::forward<JumpGenerator>(jump_generator),
      std::forward<OrientationGenerator>(orientation_generator),
      std::forward<OrientationGenerator>(orientation_generator),
      std::forward<StateSwitcher>(state_switcher),
      std::forward<Boundary>(boundary) }
    {}
    
    template <typename State>
    void operator()(State& state)
    {
      auto const& state_old = state;
      switch (state.state)
      {
        case 0:
          operation::plus_InPlace(state.position,
                                  jump_generator(state));
          if(!boundary(state, state_old))
            state.state = state_switcher.run(state);
          break;
        case 1:
          state.state = state_switcher.tumble(state);
          if (state.state == 0)
            operation::plus_InPlace(state.orientation,
                                    orientation_generator(state));
          break;
        case 2:
          state.state = state_switcher.wall_tumble(state);
          if (state.state == 0)
            operation::plus_InPlace(state.orientation,
                                    orientation_generator_wall(state));
          break;
        default:
          break;
      }
    }
    
  private:
    Boundary boundary;
    StateSwitcher state_switcher;
    JumpGenerator jump_generator;
    OrientationGenerator orientation_generator;
    OrientationGenerator_Wall orientation_generator_wall;
  };
  template
  <typename Boundary, typename StateSwitcher,
  typename JumpGenerator, typename OrientationGenerator,
  typename OrientationGenerator_Wall>
  Transitions_RunAndTumble_PTRW
  (JumpGenerator&&, OrientationGenerator&&, OrientationGenerator_Wall&&,
  StateSwitcher&&, Boundary&&) ->
  Transitions_RunAndTumble_PTRW
  <Boundary, StateSwitcher, JumpGenerator,
  OrientationGenerator, OrientationGenerator_Wall>;
  template
  <typename Boundary, typename StateSwitcher,
  typename JumpGenerator, typename OrientationGenerator>
  Transitions_RunAndTumble_PTRW
  (JumpGenerator&&, OrientationGenerator&&, StateSwitcher&&, Boundary&&) ->
  Transitions_RunAndTumble_PTRW
  <Boundary, StateSwitcher, JumpGenerator,
  OrientationGenerator, OrientationGenerator>;
}

#endif /* Transitions_State_h */
