//
//  JumpGenerator.h
//  CTRW
//
//  Created by Tomas Aquino on 1/20/18.
//  Copyright © 2018 Tomas Aquino. All rights reserved.
//

// Jump step generators for CTRW

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

namespace ctrw
{
  template <typename Step = double>
  class JumpGenerator_Step
  {
  public:
    JumpGenerator_Step(Step step)
    : step{ step }
    {}
    
    template <typename State>
    auto operator() (State const& state)
    { return step; }
    
    const Step step;
  };
  
  template <typename JumpGenerator_1, typename JumpGenerator_2>
  class JumpGenerator_Add
  {
  public:
    JumpGenerator_Add
    (JumpGenerator_1 jump_generator_1,
     JumpGenerator_2 jump_generator_2)
    : jump_generator_1{ jump_generator_1 }
    , jump_generator_2{ jump_generator_2 }
    {}
    
    template <typename State>
    auto operator() (State const& state)
    {
      return operation::plus(jump_generator_1(state), jump_generator_2(state));
    }
    
  private:
    JumpGenerator_1 jump_generator_1;
    JumpGenerator_2 jump_generator_2;
  };
  
  template <typename VelocityField>
  class JumpGenerator_Velocity
  {
  public:
    JumpGenerator_Velocity
    (VelocityField&& velocity, double dt)
    : velocity{ std::forward<VelocityField>(velocity) }
    , dt{ dt }
    {}

    void time_step(double dt)
    { this->dt = dt; }

    double time_step() const
    { return dt; }

    template <typename State>
    auto operator() (State const& state)
    {
      return operation::times_scalar(dt, velocity(state.position));
    }
    
  private:
    VelocityField velocity;
    double dt;
  };
  template <typename VelocityField>
  JumpGenerator_Velocity(VelocityField&&, double) -> JumpGenerator_Velocity<VelocityField>;
  
  template <typename VelocityField>
  class JumpGenerator_Velocity_withHint
  {
  public:
    JumpGenerator_Velocity_withHint
    (VelocityField&& velocity, double dt, std::size_t nr_particles)
    : velocity{ std::forward<VelocityField>(velocity) }
    , dt{ dt }
    , hint(nr_particles)
    {}

    void time_step(double dt)
    { this->dt = dt; }

    double time_step() const
    { return dt; }

    template <typename State>
    auto operator() (State const& state)
    {
      hint[state.tag] = velocity.locate(state.position, hint[state.tag]);
      return operation::times_scalar(dt, velocity(state.position, hint[state.tag]));
    }
    
  private:
    VelocityField velocity;
    double dt;
    
    using Hint = typename VelocityField::Hint;
    std::vector<Hint> hint;
  };
  template <typename VelocityField>
  JumpGenerator_Velocity_withHint
  (VelocityField&&, double, std::size_t) -> JumpGenerator_Velocity_withHint<VelocityField>;
  
  template <typename VelocityField, typename Boundary = boundary::DoNothing>
  class JumpGenerator_Velocity_RK4
  {
  public:
    JumpGenerator_Velocity_RK4
    (VelocityField&& velocity, double dt, Boundary&& boundary = {})
    : velocity{ std::forward<VelocityField>(velocity) }
    , dt{ dt }
    , boundary{ std::forward<Boundary>(boundary) }
    {}

    void time_step(double dt)
    { this->dt = dt; }

    double time_step() const
    { return dt; }

    template <typename State>
    auto operator() (State const& state)
    {
      auto state_intermediate = state;
      
      auto k1 = velocity(state.position);
      operation::linearOp(1., state.position, dt/2., k1, state_intermediate.position);
      boundary(state_intermediate, state);
      auto k2 = velocity(state_intermediate.position);
      operation::linearOp(1., state.position, dt/2., k2, state_intermediate.position);
      boundary(state_intermediate, state);
      auto k3 = velocity(state_intermediate.position);
      operation::linearOp(1., state.position, dt, k3, state_intermediate.position);
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
  (VelocityField&&, double, Boundary&&) -> JumpGenerator_Velocity_RK4<VelocityField, Boundary>;
  
  template <typename VelocityField, typename Boundary = boundary::DoNothing>
  class JumpGenerator_Velocity_withHint_RK4
  {
  public:
    JumpGenerator_Velocity_withHint_RK4
    (VelocityField&& velocity, double dt, std::size_t nr_particles, Boundary&& boundary = {})
    : velocity{ std::forward<VelocityField>(velocity) }
    , dt{ dt }
    , hint(nr_particles)
    , boundary{ std::forward<Boundary>(boundary) }
    {}

    void time_step(double dt)
    { this->dt = dt; }

    double time_step() const
    { return dt; }
    
    template <typename State>
    auto operator() (State const& state)
    {
      auto state_intermediate = state;
      hint[state.tag] = velocity.locate(state.position, hint[state.tag]);
      
      auto k1 = velocity(state.position, hint[state.tag]);
      operation::linearOp(1., state.position, dt/2., k1, state_intermediate.position);
      boundary(state_intermediate, state);
      auto k2 = velocity(state_intermediate.position, hint[state.tag]);
      operation::linearOp(1., state.position, dt/2., k2, state_intermediate.position);
      boundary(state_intermediate, state);
      auto k3 = velocity(state_intermediate.position, hint[state.tag]);
      operation::linearOp(1., state.position, dt, k3, state_intermediate.position);
      boundary(state_intermediate, state);
      auto k4 = velocity(state_intermediate.position, hint[state.tag]);
      
      auto jump = operation::plus(k1, k4);
      operation::linearOp_InPlace(1., jump, 2., k2);
      operation::linearOp_InPlace(1., jump, 2., k3);
      operation::times_scalar_InPlace(dt/6., jump);
      
      return jump;
    }
    
  private:
    VelocityField velocity;
    double dt;
    
    using Hint = typename std::remove_reference<VelocityField>::type::Hint;
    std::vector<Hint> hint;
    
    Boundary boundary;
  };
  template <typename VelocityField, typename Boundary>
  JumpGenerator_Velocity_withHint_RK4
  (VelocityField&&, double, std::size_t, Boundary&&) ->
  JumpGenerator_Velocity_withHint_RK4<VelocityField, Boundary>;
  
  class JumpGenerator_Diffusion_1d
  {
    double dt;
    double diff_coeff;
    double diff_aux{ std::sqrt(2.*diff_coeff*dt) };
    std::mt19937 rng{ std::random_device{}() };
    std::normal_distribution<double> normal_dist{ 0., 1. };

  public:

    JumpGenerator_Diffusion_1d(double diff_coeff, double dt)
    : dt{ dt }
    , diff_coeff{ diff_coeff }
    {}

    void time_step(double dt)
    {
      this->dt = dt;
      diff_aux = std::sqrt(2.*diff_coeff*dt);
    }

    void diff(double diff_coeff)
    {
      this->diff_coeff = diff_coeff;
      diff_aux = std::sqrt(2.*diff_coeff*dt);
    }

    double diff() const
    { return diff_coeff; }

    double time_step() const
    { return dt; }

    template <typename State = useful::Empty>
    double operator() (State const& = {})
    { return diff_aux*normal_dist(rng); }
  };
  
  class JumpGenerator_Diffusion
  {
  public:
    JumpGenerator_Diffusion(std::vector<double> const& diff, double dt)
    {
      jump_generator.reserve(diff.size());
      for (auto val : diff)
        jump_generator.emplace_back(val, dt);
    }
    
    JumpGenerator_Diffusion(double diff, double dt, std::size_t dim)
    {
      jump_generator.reserve(dim);
      for (std::size_t dd = 0; dd < dim; ++dd)
        jump_generator.emplace_back(diff, dt);
    }
    
    void time_step(double dt)
    {
      for (std::size_t dd = 0; dd < dim(); ++dd)
        jump_generator[dd].time_step(dt);
    }

    void diff(double diff_coeff)
    {
      for (std::size_t dd = 0; dd < dim(); ++dd)
        jump_generator[dd].diff(diff_coeff);
    }
    
    void diff(std::vector<double> const& diff_coeff)
    {
      for (std::size_t dd = 0; dd < dim(); ++dd)
        jump_generator[dd].diff(diff_coeff[dd]);
    }
    
    std::size_t dim()
    { return jump_generator.size(); }
    
    template <typename State>
    auto operator() (State const& state)
    {
      decltype(State::position) jump(dim(), 0.);
      for (std::size_t dd = 0; dd < dim(); ++dd)
        jump[dd] = jump_generator[dd](state);
      
      return jump;
    }
    
  private:
    std::vector<JumpGenerator_Diffusion_1d> jump_generator;
  };
    
  template <typename Distance_type = double>
  class JumpGenerator_RW_1d
  {
    const Distance_type jump_size;
    const double prob_right;
    std::mt19937 rng{ std::random_device{}() };
    std::bernoulli_distribution bernoulli_dist{ prob_right };

  public:
    JumpGenerator_RW_1d(Distance_type jump_size = 1, double prob_right = 0.5)
    : jump_size(jump_size)
    , prob_right(prob_right)
    {}

    template <typename State = useful::Empty>
    double operator() (State const& = {})
    { return jump_size*(2*bernoulli_dist(rng)-1); }
  };
    
  class JumpGenerator_Uniform_1d
  {
  public:
    const double max_jump_size;
    
    JumpGenerator_Uniform_1d(double max_jump_size)
    : max_jump_size{ max_jump_size }
    {}
    
    
    template <typename State = useful::Empty>
    double operator() (State const& = {})
    { return max_jump_size*(2.*dist(rng) - 1.); }
    
  private:
    std::mt19937 rng{ std::random_device{}() };
    std::uniform_real_distribution<double> dist{};
  };
  
  class JumpGenerator_JumpAngle_2d
  {
  public:
    const double jump_size;
    
    JumpGenerator_JumpAngle_2d(double jump_size)
    : jump_size{ jump_size }
    {}
    
    template <typename State>
    std::vector<double> const& operator()(State const& state)
    {
      jump[0] = jump_size*std::cos(state.orientation);
      jump[1] = jump_size*std::sin(state.orientation);
      
      return jump;
    }
    
  private:
    std::vector<double> jump{ 0., 0. };
  };
  
  template <typename Reference_angle>
  class OrientationGenerator_AngleClasses_Gaussian_1d
  {
  public:
    std::vector<double> angle_classes;
    std::vector<double> variances;
    
    OrientationGenerator_AngleClasses_Gaussian_1d
    (std::vector<double> angle_classes, std::vector<double> variances,
     Reference_angle const& reference_angle)
    : angle_classes{ angle_classes }
    , variances{ variances }
    , reference_angle{ reference_angle }
    {}
    
    template <typename State>
    double operator() (State const& state)
    {
      std::size_t angle_idx = angle_class(state.angle - reference_angle(state));
      double angle_increment = variances[angle_idx]*dist(rng);
      
      return std::abs(angle_increment) > constants::pi
      ? constants::pi
      : angle_increment;
    }
    
    std::size_t angle_class(double angle)
    {
      return std::lower_bound(angle_classes.begin(), angle_classes.end(), angle) -
        angle_classes.begin();
    }
    
  private:
    std::mt19937 rng{ std::random_device{}() };
    std::normal_distribution<double> dist{ 0., 1. };
    Reference_angle const& reference_angle;
  };
    
  template <typename Concentration_particle, typename Gradient_particle>
  class OrientationGenerator_Gradient_1d
  {
  public:
    double variance;
    double preferred_concentration;
    
    OrientationGenerator_Gradient_1d
    (Concentration_particle const& concentration, Gradient_particle const& gradient,
     double variance, double preferred_concentration)
    : variance{ variance }
    , preferred_concentration{ preferred_concentration }
    , concentration{ concentration }
    , gradient{ gradient }
    {}
    
    OrientationGenerator_Gradient_1d
    (Concentration_particle&& concentration, Gradient_particle const& gradient,
     double variance, double preferred_concentration) = delete;
    
    OrientationGenerator_Gradient_1d
    (Concentration_particle const& concentration, Gradient_particle&& gradient,
     double variance, double preferred_concentration) = delete;
    
    OrientationGenerator_Gradient_1d
    (Concentration_particle&& concentration, Gradient_particle&& gradient,
     double variance, double preferred_concentration) = delete;
    
    template <typename State>
    double operator() (State const& state)
    {
      double concentration_val = concentration(state);
      std::vector<double> gradient_val = gradient(state);
      if (operation::abs(gradient_val) == 0.)
        return constants::pi*(2.*dist(rng) - 1.);
      
      double reference_angle = concentration_val < preferred_concentration
      ? std::atan2(gradient_val[1], gradient_val[0])
      : std::atan2(-gradient_val[1], -gradient_val[0]);
      
      double angle = variance*dist(rng) + reference_angle;
      angle = std::abs(angle) > constants::pi
      ? constants::pi
      : angle;
      
      return angle - state.orientation;
    }
    
  private:
    std::mt19937 rng{ std::random_device{}() };
    std::normal_distribution<double> dist{ 0., 1. };
    Concentration_particle const& concentration;
    Gradient_particle const& gradient;
  };
    
  class OrientationGenerator_Flip
  {
  public:
    template <typename State>
    double operator() (State const& state)
    { return constants::pi; }
  };
}


#endif /* JumpGenerator_h */
