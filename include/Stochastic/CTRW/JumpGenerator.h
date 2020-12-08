//
//  JumpGenerator.h
//  CTRW
//
//  Created by Tomas Aquino on 1/20/18.
//  Copyright © 2018 Tomas Aquino. All rights reserved.
//

#ifndef JumpGenerator_h
#define JumpGenerator_h

#include <algorithm>
#include <cmath>
#include <random>
#include <type_traits>
#include <utility>
#include "Stochastic/CTRW/Boundary.h"
#include "general/Constants.h"
#include "general/Operations.h"
#include "general/useful.h"

// A jump generator must implement the following basic functionality:
// class JumpGenerator
// {
//   // Return the jump increment
//   // with type of state element to be incremented
//   template <typename State>
//   auto operator() (State const& state)
//   {
//     // Return the jump increment
//   }
//}

namespace ctrw
{
  // Jump a fixed step
  template <typename Step = double>
  class JumpGenerator_Step
  {
  public:
    // Construct given the fixed step
    JumpGenerator_Step(Step step)
    : step{ step }
    {}
    
    // Return the jump increment
    template <typename State = useful::Empty>
    auto operator() (State const& state = {})
    { return step; }
    
    const Step step;  // The fixed step to jump
  };
  
  // Jump according to the sum of two jump generators
  template <typename JumpGenerator_1, typename JumpGenerator_2>
  class JumpGenerator_Add
  {
  public:
    // Construct given the two jump generators
    JumpGenerator_Add
    (JumpGenerator_1 jump_generator_1,
     JumpGenerator_2 jump_generator_2)
    : jump_generator_1{ jump_generator_1 }
    , jump_generator_2{ jump_generator_2 }
    {}
    
    // Return the jump increment
    template <typename State = useful::Empty>
    auto operator() (State const& state = {})
    {
      return operation::plus(jump_generator_1(state),
                             jump_generator_2(state));
    }
    
  private:
    JumpGenerator_1 jump_generator_1;
    JumpGenerator_2 jump_generator_2;
  };
  
  // Jump according to a velocity field and time step (forward Euler)
  // State must define: position
  template <typename VelocityField>
  class JumpGenerator_Velocity
  {
  public:
    // Construct given the velocity field as a function of position
    // and the time step
    JumpGenerator_Velocity
    (VelocityField&& velocity, double dt)
    : velocity{ std::forward<VelocityField>(velocity) }
    , dt{ dt }
    {}

    // Set the time step
    void time_step(double dt)
    { this->dt = dt; }

    // Get the time step
    double time_step() const
    { return dt; }

    // Return the jump increment
    template <typename State>
    auto operator() (State const& state)
    {
      return operation::times_scalar(dt, velocity(state.position));
    }
    
  private:
    VelocityField velocity;   // Velocity as a function of position
    double dt;                // Time step
  };
  template <typename VelocityField>
  JumpGenerator_Velocity(VelocityField&&, double)
  -> JumpGenerator_Velocity<VelocityField>;
  
  // Jump according to a velocity field and time step (forward Euler)
  // with the previous position as a hint to locate
  // the current position for velocity interpolation
  // State must define: position
  //                    tag (to identify particles for hints)
  // Note: Particle number and tags must remain fixed
  template <typename VelocityField>
  class JumpGenerator_Velocity_withHint
  {
  public:
    // Construct given the velocity field as a function of position,
    // the time step, and the nr of particles, for which hints will be stored
    JumpGenerator_Velocity_withHint
    (VelocityField&& velocity, double dt, std::size_t nr_particles)
    : velocity{ std::forward<VelocityField>(velocity) }
    , dt{ dt }
    , hint(nr_particles)
    {}

    // Set the time step
    void time_step(double dt)
    { this->dt = dt; }

    // Get the time step
    double time_step() const
    { return dt; }

    // Return the jump increment
    template <typename State>
    auto operator() (State const& state)
    {
      hint[state.tag] = velocity.locate(state.position, hint[state.tag]);
      return operation::times_scalar(dt,
        velocity(state.position, hint[state.tag]));
    }
    
  private:
    VelocityField velocity;  // Velocity as a function of position
    double dt;               // Time step
    
    using Hint =
      typename std::remove_reference<VelocityField>::
      type::Hint;            // Hint type
    std::vector<Hint> hint;  // Hints for each particle
  };
  template <typename VelocityField>
  JumpGenerator_Velocity_withHint
  (VelocityField&&, double, std::size_t)
  -> JumpGenerator_Velocity_withHint<VelocityField>;
  
  // Jump according to a velocity field and time step, using RK4 scheme
  // Boundary conditions are enforced in predictor-corrector steps
  // State must define: position
  template <typename VelocityField, typename Boundary = boundary::DoNothing>
  class JumpGenerator_Velocity_RK4
  {
  public:
    // Construct given the velocity field as a function of position,
    // the time step, and the boundary to enforce boundary conditions
    JumpGenerator_Velocity_RK4
    (VelocityField&& velocity, double dt, Boundary&& boundary = {})
    : velocity{ std::forward<VelocityField>(velocity) }
    , dt{ dt }
    , boundary{ std::forward<Boundary>(boundary) }
    {}

    // Set the time step
    void time_step(double dt)
    { this->dt = dt; }

    // Get the time step
    double time_step() const
    { return dt; }

    // Return the jump increment
    template <typename State>
    auto operator() (State const& state)
    {
      auto state_intermediate = state;
      
      auto k1 = velocity(state.position);
      operation::linearOp(1., state.position,
                          dt/2., k1, state_intermediate.position);
      boundary(state_intermediate, state);
      auto k2 = velocity(state_intermediate.position);
      operation::linearOp(1., state.position,
                          dt/2., k2, state_intermediate.position);
      boundary(state_intermediate, state);
      auto k3 = velocity(state_intermediate.position);
      operation::linearOp(1., state.position,
                          dt, k3, state_intermediate.position);
      boundary(state_intermediate, state);
      auto k4 = velocity(state_intermediate.position);
      
      auto jump = operation::sum(k1, k4);
      operation::linearOp_InPlace(1., jump, 2., k2);
      operation::linearOp_InPlace(1., jump, 2., k3);
      operation::times_scalar_InPlace(dt/6., jump);
      
      return jump;
    }
    
  private:
    VelocityField velocity;
    double dt;
    Boundary boundary;
  };
  template <typename VelocityField, typename Boundary>
  JumpGenerator_Velocity_RK4
  (VelocityField&&, double, Boundary&&)
  -> JumpGenerator_Velocity_RK4<VelocityField, Boundary>;
  
  // Jump according to a velocity field and time step, using RK4 scheme,
  // with the previous position as a hint to locate
  // the current position for velocity interpolation
  // State must define: position
  //                    tag (to identify particles for hints)
  // Note: Particle number and tags must remain fixed
  template
  <typename VelocityField, typename Boundary = boundary::DoNothing>
  class JumpGenerator_Velocity_withHint_RK4
  {
  public:
    // Construct given the velocity field as a function of position,
    // the time step, the nr of particles, for which hints will be stored,
    // and the boundary to enforce boundary conditions
    JumpGenerator_Velocity_withHint_RK4
    (VelocityField&& velocity, double dt,
     std::size_t nr_particles, Boundary&& boundary = {})
    : velocity{ std::forward<VelocityField>(velocity) }
    , dt{ dt }
    , hint(nr_particles)
    , boundary{ std::forward<Boundary>(boundary) }
    {}

    // Set the time step
    void time_step(double dt)
    { this->dt = dt; }

    // Get the time step
    double time_step() const
    { return dt; }
    
    // Return the jump increment
    template <typename State>
    auto operator() (State const& state)
    {
      auto state_intermediate = state;
      hint[state.tag] = velocity.locate(state.position, hint[state.tag]);
      
      auto k1 = velocity(state.position, hint[state.tag]);
      operation::linearOp(dt/2., k1, state.position,
                          state_intermediate.position);
      boundary(state_intermediate, state);
      auto k2 = velocity(state_intermediate.position, hint[state.tag]);
      operation::linearOp(dt/2., k2, state.position,
                          state_intermediate.position);
      boundary(state_intermediate, state);
      auto k3 = velocity(state_intermediate.position, hint[state.tag]);
      operation::linearOp(dt, k3, state.position,
                          state_intermediate.position);
      boundary(state_intermediate, state);
      auto k4 = velocity(state_intermediate.position, hint[state.tag]);
      
      auto jump = operation::plus(k1, k4);
      operation::linearOp_InPlace(1., jump, 2., k2);
      operation::linearOp_InPlace(1., jump, 2., k3);
      operation::times_scalar_InPlace(dt/6., jump);
      
      return jump;
    }
    
  private:
    VelocityField velocity;  // Velocity as a function of position
    double dt;               // Time step
    
    using Hint =
      typename std::remove_reference<VelocityField>::type::Hint; // Hint type
    std::vector<Hint> hint;
    
    Boundary boundary;
  };
  template <typename VelocityField, typename Boundary>
  JumpGenerator_Velocity_withHint_RK4
  (VelocityField&&, double, std::size_t, Boundary&&) ->
  JumpGenerator_Velocity_withHint_RK4<VelocityField, Boundary>;
  
  // One dimensional diffusion jumps
  // using Gaussian distribution
  class JumpGenerator_Diffusion_1d
  {
    double dt;          // Time step
    double diff_coeff;  // Diffusion coefficient
    
    double diff_aux{ std::sqrt(2.*diff_coeff*dt) };  // Typicall jump size
    std::mt19937 rng{ std::random_device{}() };      // Random number generator
    std::normal_distribution<double> normal_dist{
      0., 1. };                                      // Unit normal distribution

  public:

    // Construct given diffusion coefficient and time step
    JumpGenerator_Diffusion_1d(double diff_coeff, double dt)
    : dt{ dt }
    , diff_coeff{ diff_coeff }
    {}

    // Set time step
    void time_step(double dt)
    {
      this->dt = dt;
      diff_aux = std::sqrt(2.*diff_coeff*dt);
    }

    // Set diffusion coefficient
    void diff(double diff_coeff)
    {
      this->diff_coeff = diff_coeff;
      diff_aux = std::sqrt(2.*diff_coeff*dt);
    }

    // Get time step
    double time_step() const
    { return dt; }
    
    // Get diffusion coefficient
    double diff() const
    { return diff_coeff; }

    // Return the jump increment
    template <typename State = useful::Empty>
    double operator() (State const& = {})
    { return diff_aux*normal_dist(rng); }
  };
  
  // Diffusion jumps along each dimension
  // using Gaussian distributions
  class JumpGenerator_Diffusion
  {
  public:
    // Construct given diffusion coefficients along each dimension
    // and time step
    JumpGenerator_Diffusion
    (std::vector<double> const& diff, double dt)
    {
      jump_generator.reserve(diff.size());
      for (auto val : diff)
        jump_generator.emplace_back(val, dt);
    }
    
    // Construct given same diffusion coefficient along all dimensions,
    // time step, and spatial dimension
    JumpGenerator_Diffusion
    (double diff, double dt, std::size_t dim)
    {
      jump_generator.reserve(dim);
      for (std::size_t dd = 0; dd < dim; ++dd)
        jump_generator.emplace_back(diff, dt);
    }
    
    // Set time step
    void time_step(double dt)
    {
      for (std::size_t dd = 0; dd < dim(); ++dd)
        jump_generator[dd].time_step(dt);
    }

    // Set same diffusion coefficient along all dimensions
    void diff(double diff_coeff)
    {
      for (std::size_t dd = 0; dd < dim(); ++dd)
        jump_generator[dd].diff(diff_coeff);
    }
    
    // Set diffusion coefficient along each dimension
    void diff(std::vector<double> const& diff_coeff)
    {
      for (std::size_t dd = 0; dd < dim(); ++dd)
        jump_generator[dd].diff(diff_coeff[dd]);
    }
    
    // Get time step
    double time_step()
    { return jump_generator[0].time_step(); }
    
    // Get diffusion coefficient along dimension dd
    double diff(std::size_t idx)
    { return jump_generator[0].diff(); }
    
    // Get spatial dimension
    std::size_t dim()
    { return jump_generator.size(); }
    
    // Return the jump increment
    template <typename State = useful::Empty>
    auto operator() (State const& state = {})
    {
      decltype(State::position) jump(dim(), 0.);
      for (std::size_t dd = 0; dd < dim(); ++dd)
        jump[dd] = jump_generator[dd](state);
      
      return jump;
    }
    
  private:
    std::vector<JumpGenerator_Diffusion_1d>
      jump_generator;  // Diffusion jump generators along each dimension
  };
    
  // One-dimensional random walk jumps with fixed step sizes
  // to left or right
  template <typename Distance_type = double>
  class JumpGenerator_RW_1d
  {
    const Distance_type jump_size;               // Fixed jump size
    const double prob_right;                     // Probability to jump right
    std::mt19937 rng{ std::random_device{}() };  // Random number generator
    std::bernoulli_distribution bernoulli_dist{
      prob_right };                              // Bernoulli distribution

  public:
    // Construct given fixed jump size and probability to jump right
    JumpGenerator_RW_1d
    (Distance_type jump_size = 1, double prob_right = 0.5)
    : jump_size(jump_size)
    , prob_right(prob_right)
    {}

    // Return the jump increment
    template <typename State = useful::Empty>
    double operator() (State const& = {})
    { return jump_size*(2*bernoulli_dist(rng)-1); }
  };
  
  // One-dimensional random walk jumps with uniformly-distributed step sizes
  class JumpGenerator_Uniform_1d
  {
  public:
    const double max_jump_size;  // Jumps uniformly distributed between
                                 // this maximum size and its negative
    
    // Construct given maximum jump size
    JumpGenerator_Uniform_1d(double max_jump_size)
    : max_jump_size{ max_jump_size }
    {}
    
    // Return the jump increment
    template <typename State = useful::Empty>
    double operator() (State const& = {})
    { return max_jump_size*(2.*dist(rng) - 1.); }
    
  private:
    std::mt19937 rng{ std::random_device{}() };     // Random number generator
    std::uniform_real_distribution<double> dist{};  // Unit uniform distribution
  };
  
  // Random walk jumps along current orientation
  // with fixed step size
  // State must implement: orientation (scalar)
  class JumpGenerator_JumpAngle_2d
  {
  public:
    const double jump_size;  // Fixed jump size
    
    // Construct given maximum jump size
    JumpGenerator_JumpAngle_2d(double jump_size)
    : jump_size{ jump_size }
    {}
    
    // Return the jump increment
    template <typename State>
    std::vector<double> operator()(State const& state)
    {
      return { jump_size*std::cos(state.orientation),
        jump_size*std::sin(state.orientation) };
    }
  };
  
  // Gaussian-distributed angle jumps
  // around local gradient
  // towards preferred concentration value
  // State must define: orientation (scalar)
  template
  <typename Concentration, typename Gradient>
  class OrientationGenerator_Gradient_1d
  {
  public:
    double variance;                 // Variance of jumps around local gradient
    double preferred_concentration;  // Concentration to jump towards
    
    // Construct given concentration as a function of state,
    // gradient as a function of state,
    // variance of jumps around local gradient,
    // and preferred concentration to seek
    OrientationGenerator_Gradient_1d
    (Concentration const& concentration,
     Gradient const& gradient,
     double variance, double preferred_concentration)
    : variance{ variance }
    , preferred_concentration{ preferred_concentration }
    , concentration{ concentration }
    , gradient{ gradient }
    {}
    
    // Avoid temporary concentration object
    OrientationGenerator_Gradient_1d
    (Concentration&& concentration,
     Gradient const& gradient,
     double variance, double preferred_concentration) = delete;
    
    // Avoid temporary gradient object
    OrientationGenerator_Gradient_1d
    (Concentration const& concentration,
     Gradient&& gradient,
     double variance, double preferred_concentration) = delete;
    
    // Avoid temporary concentration and gradient objects
    OrientationGenerator_Gradient_1d
    (Concentration&& concentration,
     Gradient&& gradient,
     double variance, double preferred_concentration) = delete;
    
    // Return the orientation increment
    template <typename State>
    double operator() (State const& state)
    {
      double concentration_val = concentration(state);
      std::vector<double> gradient_val = gradient(state);
      
      // If local gradient is zero jump in a random direction
      if (operation::abs(gradient_val) == 0.)
        return constants::pi*(2.*dist(rng) - 1.);
      
      // Orient mean jump direction along gradient if below preferred concentration,
      // or along opposite direction of gradient if above
      double reference_angle =
        concentration_val < preferred_concentration
        ? std::atan2(gradient_val[1], gradient_val[0])
        : std::atan2(-gradient_val[1], -gradient_val[0]);
      
      // Compute angle jump and bound it to ]-pi,pi]
      double angle = variance*dist(rng) + reference_angle;
      angle = std::abs(angle) > constants::pi
      ? constants::pi
      : angle;
      
      return angle - state.orientation;
    }
    
  private:
    std::mt19937 rng{ std::random_device{}() };       // Random number generator
    std::normal_distribution<double> dist{ 0., 1. };  // Unit normal distribution
    Concentration concentration;                      // Concentration as a function of state
    Gradient gradient;                                // Gradient as a function of state
  };
    
  // Flip orientation angle
  class OrientationGenerator_Flip
  {
  public:
    // Return the orientation increment
    template <typename State>
    double operator() (State const& state)
    { return constants::pi; }
  };
}


#endif /* JumpGenerator_h */
