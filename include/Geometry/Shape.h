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
    using Position = Position_t;  // Type of center position
    Position corner;              // Lower left corner
    Position dimensions;          // Sizes of rectangle along each dimension
    
    Position half_dimensions{
      operation::div_scalar(dimensions, 2.) };     // Half the size along each dimension
    Position center{
      operation::plus(corner, half_dimensions) };  // Position of rectangle center
    
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
      auto position_helper = operation::minus(position, center);
      for (std::size_t dd = 0; dd < dim(); ++dd)
        if (std::abs(position_helper[dd]) > half_dimensions[dd])
          return false;
      return true;
    }
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
  // Put reflected position in reflected and contact point with sphere on
  // Note: It is safe to reuse position_new and position_old
  //       for reflected and contact_point
  template <typename Position, typename Sphere>
  void ReflectOffSphere_OutsideToInside
  (Position const& position_new, Position const& position_old,
   Sphere const& sphere, Position& reflected, Position& contact_point)
  {
    // Vector from outside position to outside position
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
  // Note: It is safe to reuse position_new and position_old for reflected
  template <typename Position, typename Sphere>
  void ReflectOffSphere_OutsideToInside
  (Position const& position_new, Position const& position_old,
   Sphere const& sphere, Position& reflected)
  {
    auto contact_point(position_new.size());
    ReflectOffSphere_OutsideToInside(
      position_new, position_old, sphere, reflected, contact_point);
  }
}

#endif /* Shape_h */
