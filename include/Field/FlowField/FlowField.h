//
//  FlowField.h
//  FlowField
//
//  Created by Tomas Aquino on 10/23/19.
//  Copyright © 2019 Tomas Aquino. All rights reserved.
//

#ifndef FlowField_h
#define FlowField_h

#include <cmath>
#include "general/Constants.h"
#include "general/useful.h"

namespace flowfield
{
  double poiseuille(double radial_position, double coeff, double radius)
  {
    return coeff*(radius-radial_position)*(radius+radial_position);
  }

  template <std::size_t dim>
  class FlowField_Poiseuille
  {
  public:
    const struct Params
    {
      double average_velocity;
      double half_width;
      double cross_section;
      double vmax;
      double vmin;
    } params;

    FlowField_Poiseuille(double average_velocity, double radius)
    : params{ make_params(average_velocity, radius, useful::Selector<std::size_t, dim>{}) }
    , coeff{ params.vmax/(radius*radius) }
    {}

    template <typename Position_type>
    double radial_position(Position_type const& position) const
    {
      double radial_pos = 0.;
      for (std::size_t dd = 1; dd < position.size(); ++dd)
        radial_pos += position[dd]*position[dd];
      return std::sqrt(radial_pos);
    }

    template <typename Position_type>
    double operator()(Position_type const& position) const
    {
      return poiseuille(radial_position(position), coeff, params.half_width);
    }

  private:
    Params make_params(double average_velocity, double radius, useful::Selector<std::size_t, 3>) const
    {
      Params params;
      params.average_velocity = average_velocity;
      params.half_width = radius;
      params.cross_section = constants::pi*radius*radius;
      params.vmax = 2.*average_velocity;
      params.vmin = 0.;

      return params;
    }

    Params make_params(double average_velocity, double radius, useful::Selector<std::size_t, 2>) const
    {
      Params params = make_params(average_velocity, radius, useful::Selector<std::size_t, 3>{});
      params.vmax *= 3./4.;

      return params;
    }

    const double coeff;
  };

  class FlowField_Poiseuille_2d
  {
  public:
    const struct Params
    {
      double average_velocity;
      double half_width;
      double cross_section;
      double vmax;
      double vmin;
    } params;

    FlowField_Poiseuille_2d(double average_velocity, double half_width)
    : params{ average_velocity, half_width, 2.*half_width, 3./2.*average_velocity, 0. }
    , coeff{ params.vmax/(params.half_width*params.half_width) }
    {}

    template <typename Position_type>
    double operator()(Position_type const& position) const
    {
      return poiseuille(std::abs(position[1]), coeff, params.half_width);
    }

  private:
    const double coeff;
  };

  class FlowField_M_2d
  {
  public:
    const struct Params
    {
      double cross_section;
      double half_width;
      double average_velocity;
      double average_deg;
      double vmax;
      double vmin;
      double critical_velocity;
      double shear;
    } params;

    FlowField_M_2d(double average_velocity, double half_width, double average_deg)
    : params{ make_params(average_velocity, half_width, average_deg) }
    , length_inner{ (average_deg-2.)*half_width/4. }
    {}

    template <typename Position_type>
    double operator()(Position_type const& position) const
    {
      double pos = std::abs(position[1]);
      if (pos > length_inner)
        return params.shear*(params.half_width - pos);
      else
        return params.critical_velocity + params.shear*pos;
    }

  private:
    Params make_params(double average_velocity, double half_width, double average_deg) const
    {
      Params params;
      params.cross_section = 2.*half_width;
      params.half_width = half_width;
      params.average_velocity = average_velocity;
      params.average_deg = average_deg;

      double aux6 = 6. - average_deg;
      double aux4 = 4. - average_deg;
      params.vmax = 4.*aux6/(aux6*aux6 - 2.*aux4*aux4);
      params.vmin = 0.;
      params.critical_velocity = 2*params.vmax*aux4/aux6;
      params.shear = 4.*params.vmax/(aux6*half_width);

      return params;
    }

    const double length_inner;
  };

  class FlowField_Couette_2d
  {
  public:
    const struct Params
    {
      double cross_section;
      double half_width;
      double average_velocity;
      double vmax;
      double vmin;
      double shear;
    } params;

    FlowField_Couette_2d(double average_velocity, double half_width)
    : params{ 2.*half_width, half_width, average_velocity,
      2.*average_velocity, 0., average_velocity/half_width }
    {}

    template <typename Position_type>
    double operator()(Position_type const& position) const
    { return params.shear*(params.half_width+position[1]); }
  };

  class FlowField_Triangular_2d
  {
  public:
    const struct Params
    {
      double cross_section;
      double half_width;
      double average_velocity;
      double vmax;
      double vmin;
      double shear;
    } params;

    FlowField_Triangular_2d(double average_velocity, double half_width, double average_deg)
    : params{ 2.*half_width, half_width, average_velocity,
      2.*average_velocity, 0., 2.*average_velocity/half_width }
    {}

    template <typename Position_type>
    double operator()(Position_type const& position) const
    { return params.shear*(params.half_width - std::abs(position[1])); }
  };
  
  template <typename FlowField, typename VelocityClasses>
  struct Velocity_idx_from_position
  {
    FlowField const& flow_field;
    VelocityClasses const& velocity_classes;
    
    template <typename Position>
    std::size_t operator()(Position const& position) const
    {
      return velocity_classes.idx(flow_field(position));
    }
  };
}

#endif /* FlowField_h */
