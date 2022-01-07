//
//  Shape.h
//  RunAndTumble
//
//  Created by Tomás Aquino on 07/04/2020.
//  Copyright © 2020 Tomás Aquino. All rights reserved.
//

#ifndef Shape_h
#define Shape_h

#include <cmath>
#include <vector>
#include "general/Operations.h"

namespace geometry
{
  // Rectangle in any dimension
  template <typename Position_t = std::vector<double>>
  struct Parallelepiped
  {
    using Position = Position_t;  // Type of coordinate container
    
    Position corner;              // Lower left corner
    Position dimensions;          // Sizes of rectangle along each dimension
    
    // Construct undefined rectangle
    Parallelepiped()
    {}
    
    // Construct rectangle given lowe left corner and dimensions
    Parallelepiped(Position corner, Position dimensions)
    : corner{ corner }
    , dimensions{ dimensions }
    {}
      
    // Spatial dimension of rectangle
    // Note: Zero if dimensions have not been assigned
    std::size_t dim() const
    { return dimensions.size(); }
    
    // Return true if position is inside rectangle, false otherwise
    template <typename Position>
    bool inside(Position const& position) const
    {
      auto position_helper = operation::minus(position, get_center());
      auto half_dim = get_half_dimensions();
      for (std::size_t dd = 0; dd < dim(); ++dd)
        if (std::abs(position_helper[dd]) > half_dim[dd])
          return false;
      return true;
    }
    
    // Half the size along each dimension
    Position get_half_dimensions() const
    { return operation::div_scalar(dimensions, 2.); }
      
    // Position of rectangle center
    Position get_center() const
    { return operation::plus(corner, get_half_dimensions()); };
  };
  
  // Sphere in any dimension
  template <typename Position_t = std::vector<double>>
  struct Sphere
  {
    using Position = Position_t;  // Type of center position
    Position center;              // Sphere center position
    double radius;                // Sphere radius
    
    // Construct undefined sphere
    Sphere()
    {}

    // Construct with given center and radius
    Sphere(Position center, double radius)
    : center{ center }
    , radius{ radius }
    {}

    // Spatial dimension of sphere
    // Note: Zero if center has not been assigned
    std::size_t dim() const
    { return center.size(); }
    
    // Position of rectangle center
    Position get_center() const
    { return center; };

    // Return true if position is inside sphere, false otherwise
    template <typename Position>
    bool inside(Position const& position) const
    {
      auto position_helper = operation::minus(position, center);
      if (operation::abs_sq(position_helper) > radius*radius)
        return false;
      return true;
    }
  };
  
  // Return true if any ctrw particle position
  // is outside domain, false otherwise
  template <typename CTRW, typename Domain>
  bool out_of_bounds(CTRW const& ctrw, Domain const& domain)
  {
    for (auto const& part : ctrw.particles())
    if (domain.out_of_bounds(part.state_new().position))
      return 1;
    return 0;
  }
  
  // Domain with shape Shape,
  // with rectangular and/or spherical inclusions in any spatial dimension
  template <typename Shape>
  struct Domain
  {
    using DomainShape = Shape;  // Bounding domain shape type
    
    // Make empty domain
    Domain()
    {}
    
    // Make empty domain with given shape box
    Domain(DomainShape box)
    : box{ box }
    {}
    
    // Make domain with given shape box
    // and given rectangular and spherical inclusions
    Domain
    (DomainShape box,
     std::vector<geometry::Parallelepiped<>> parallelepipeds,
     std::vector<geometry::Sphere<>> sphere)
    : box{ box }
    {}
    
    Shape box;  // Bounding domain shape object
    std::vector<geometry::Parallelepiped<>>
      parallelepipeds; // Rectangular inclusions
    std::vector<geometry::Sphere<>> spheres; // Spherical inclusions
    
    std::vector<double> dimensions() const
    {
      if constexpr (std::is_same<Shape,geometry::Sphere<>>::value)
        return std::vector<double>(box.dim(), 2.*box.radius);
      else
        return box.dimensions;
    }
    
    // Return true if position is outside domain, false otherwise
    template <typename Position>
    bool out_of_bounds(Position const& position) const
    {
      if (!box.inside(position))
        return 1;
      for (auto const& shape : parallelepipeds)
        if (shape.inside(position))
          return 1;
      for (auto const& shape : spheres)
        if (shape.inside(position))
          return 1;
      
      return 0;
    }
  };
  
  // Reflect off sphere exterior (in any dimension)
  // with old position outside and new position inside
  // Put reflected position in reflected and contact point with sphere in contacg_point
  // Note: it is safe to pass the same container for new and reflected quantities
  template <typename Position, typename Sphere>
  void reflectOffSphere_outsideToInside
  (Position const& position_new, Position const& position_old,
   Sphere const& sphere, Position& reflected, Position& contact_point)
  {
    // Vector from outside position to inside position
    auto insideminusoutside = operation::minus(position_new, position_old);
    
    // Vector from outside position to sphere center
    auto outsideminuscenter = operation::minus(position_old, sphere.center);
    
    // Find contact point with sphere
    double coeff_a = operation::abs_sq(insideminusoutside);
    double coeff_b = 2.*operation::dot(insideminusoutside, outsideminuscenter);
    double coeff_c = operation::dot(outsideminuscenter, outsideminuscenter)
      - sphere.radius*sphere.radius;
    double fraction_to_intersection =
      -(coeff_b + std::sqrt(coeff_b*coeff_b - 4.*coeff_a*coeff_c))/(2.*coeff_a);
    operation::linearOp(fraction_to_intersection, insideminusoutside,
                        position_old, contact_point);
    
    // Vector from contact point to inside position
    auto leftover = operation::times_scalar(1.-fraction_to_intersection,
                                            insideminusoutside);
    
    // Compute unit normal to sphere surface at contact
    auto normal = operation::minus(contact_point, sphere.center);
    operation::div_scalar_InPlace(normal, sphere.radius);
    
    // Compute reflected position
    operation::linearOp(-2.*operation::dot(leftover, normal), normal,
                        position_new, reflected);
  }
  
  // Reflect off sphere exterior (in any dimension)
  // with old position outside and new position inside
  // Put reflected position in reflected
  // Note: it is safe to pass the same container for new and reflected quantities
  template <typename Position, typename Sphere>
  void reflectOffSphere_outsideToInside
  (Position const& position_new, Position const& position_old,
   Sphere const& sphere, Position& reflected)
  {
    auto contact_point(position_new.size());
    reflectOffSphere_outsideToInside(
      position_new, position_old, sphere, reflected, contact_point);
  }
  
  
  // Reflect off sphere exterior (in any dimension)
  // based on position and velocity inside state
  // with old position outside and new position inside
  // Trajectory is approximated as a straight line
  // Put reflected position and velocity in position_reflected and velocity_reflected
  // and contact point with sphere in contact_point
  // Note: it is safe to pass the same container for new and reflected quantities
  template <typename Position, typename Sphere>
  void reflectOffSphere_velocity_outsideToInside
  (Position const& position_new, Position const& position_old,
   Position const& velocity_old,
   Sphere const& sphere,
   Position& position_reflected, Position& velocity_reflected,
   Position& contact_point, Position& contact_velocity, double& time_to_contact,
   double time_step, Position const& acceleration,
   double restitution_coeff_norm = 1.,
   double restitution_coeff_tang = 1., double radius = 0.)
  {
    double effective_radius = sphere.radius + radius;
    
    // Vector from outside position to inside position
    auto insideminusoutside = operation::minus(position_new, position_old);
    
    // Vector from outside position to sphere center
    auto outsideminuscenter = operation::minus(position_old, sphere.center);
    
    // Find contact point with sphere
    double coeff_a = operation::abs_sq(insideminusoutside);
    double coeff_b = 2.*operation::dot(insideminusoutside, outsideminuscenter);
    double coeff_c = operation::dot(outsideminuscenter, outsideminuscenter)
      - effective_radius*effective_radius;
    double fraction_to_intersection =
      -(coeff_b + std::sqrt(coeff_b*coeff_b - 4.*coeff_a*coeff_c))/(2.*coeff_a);
    operation::linearOp(fraction_to_intersection, insideminusoutside,
                        position_old, contact_point);
    
    // Vector from outside position to contact point
    auto to_contact = operation::times_scalar(fraction_to_intersection,
                                              insideminusoutside);
    
    // Compute unit normal to sphere surface at contact
    auto normal = operation::minus(contact_point, sphere.center);
    operation::div_scalar_InPlace(normal, effective_radius);

    time_to_contact = std::sqrt(operation::abs_sq(to_contact)/
                                operation::abs_sq(velocity_old));
    operation::linearOp(time_to_contact, acceleration,
                        contact_velocity,
                        contact_velocity);
    operation::linearOp(-2.*operation::dot(contact_velocity, normal), normal,
                        velocity_old,
                        contact_velocity);
    
    auto normal_vel = operation::times_scalar(operation::dot(contact_velocity, normal), normal);
    auto tangential_vel = operation::minus(contact_velocity, normal_vel);
    operation::times_scalar_InPlace(restitution_coeff_norm, normal_vel);
    operation::times_scalar_InPlace(restitution_coeff_tang, tangential_vel);
    operation::plus(normal_vel, tangential_vel, contact_velocity);
    
    // Compute reflected position and velocity
    double time_leftover = time_step - time_to_contact;
    operation::linearOp(time_leftover, contact_velocity,
                        contact_point,
                        position_reflected);
    operation::linearOp(time_leftover, acceleration,
                        contact_velocity,
                        velocity_reflected);
  }
  
  // Reflect off sphere exterior (in any dimension)
  // based on position and velocity inside state
  // with old position outside and new position inside
  // Trajectory is approximated as a straight line
  // Put reflected position and velocity in position_reflected and velocity_reflected
  // Note: it is safe to pass the same container for new and reflected quantities
  template <typename Position, typename Sphere, typename State,
  typename Acceleration >
  void reflectOffSphere_velocity_outsideToInside
  (Position const& position_new, Position const& position_old,
  Position const& velocity_old,
  Sphere const& sphere,
  Position& position_reflected, Position& velocity_reflected,
  double time_step, Position acceleration,
  double restitution_coeff = 1., double radius = 0.)
  {
    auto contact_point(position_new.size());
    auto contact_velocity(position_new.size());
    double time_to_contact;
    reflectOffSphere_velocity_outsideToInside(
      position_new, position_old, velocity_old,
      sphere,
      position_reflected, velocity_reflected,
      contact_point, contact_velocity, time_to_contact,
      time_step, acceleration, restitution_coeff, radius);
  }
}

#endif /* Shape_h */
