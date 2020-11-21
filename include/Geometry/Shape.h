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
#include "Geometry/Shape.h"

namespace geometry
{
  template <typename Container_t = std::vector<double>>
  struct Parallelepiped
  {
    using Position_type = Container_t;
    Position_type corner;
    Position_type dimensions;
    
    Position_type half_dimensions{ operation::div_scalar(dimensions, 2.) };
    Position_type center{ operation::plus(corner, half_dimensions) };
    
    Parallelepiped()
    {}
    
    Parallelepiped(Position_type corner, Position_type dimensions)
    : corner{ corner }
    , dimensions{ dimensions }
    {}
      
    std::size_t dim() const
    { return dimensions.size(); }
    
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
  
  template <typename Container_t = std::vector<double>>
  struct Sphere
  {
    using Position_type = Container_t;
    Position_type center;
    double radius;
    
    Sphere()
    {}

    Sphere(Position_type center, double radius)
    : center{ center }
    , radius{ radius }
    {}

    std::size_t dim() const
    { return center.size(); }

    template <typename Position>
    bool inside(Position const& position) const
    {
      auto position_helper = operation::minus(position, center);
      if (operation::abs_sq(position_helper) > radius*radius)
        return false;
      return true;
    }
  };
  
  template <typename CTRW, typename Domain>
  bool out_of_bounds(CTRW const& ctrw, Domain const& domain)
  {
    for (auto const& part : ctrw.particles())
    if (domain.out_of_bounds(part.state_new().position))
      return 1;
    return 0;
  }
  
  template <typename Shape>
  struct Domain
  {
    using DomainShape = Shape;
    
    Domain()
    {}
    
    Domain(DomainShape box)
    : box{ box }
    {}
    
    Domain
    (DomainShape box, std::vector<geometry::Parallelepiped<>> parallelepipeds,
     std::vector<geometry::Sphere<>> sphere)
    : box{ box }
    {}
    
    Shape box;
    std::vector<geometry::Parallelepiped<>> parallelepipeds;
    std::vector<geometry::Sphere<>> spheres;
    
    std::vector<double> dimensions() const
    {
      if constexpr (std::is_same<Shape,geometry::Sphere<>>::value)
        return std::vector<double>(box.dim(), 2.*box.radius);
      else
        return box.dimensions;
    }
    
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
  
  template <typename Position, typename Sphere>
  void ReflectOffSphere_OutsideToInside
  (Position const& position_new, Position const& position_old, Sphere const& sphere,
   Position& reflected, Position& contact_point)
  {
    auto insideminusoutside = operation::minus(position_new, position_old);
    auto outsideminuscenter = operation::minus(position_old, sphere.center);
    
    double coeff_a = operation::abs_sq(insideminusoutside);
    double coeff_b = 2.*operation::dot(insideminusoutside, outsideminuscenter);
    double coeff_c = operation::dot(outsideminuscenter, outsideminuscenter) - sphere.radius*sphere.radius;
    double fraction_to_intersection = -(coeff_b +
                                        std::sqrt(coeff_b*coeff_b - 4.*coeff_a*coeff_c))/(2.*coeff_a);
    operation::linearOp(fraction_to_intersection, insideminusoutside, position_old, contact_point);
    
    auto leftover = operation::times_scalar(1.-fraction_to_intersection, insideminusoutside);
    auto normal = operation::minus(contact_point, sphere.center);
    operation::div_scalar_InPlace(normal, sphere.radius);
    
    operation::linearOp(-2.*operation::dot(leftover, normal), normal, position_new, reflected);
  }
  
  template <typename Position, typename Sphere>
  void ReflectOffSphere_OutsideToInside
  (Position const& position_new, Position const& position_old, Sphere const& sphere,
   Position& reflected)
  {
    auto contact_point(position_new.size());
    ReflectOffSphere_OutsideToInside(position_new, position_old, sphere, reflected, contact_point);
  }
}

#endif /* Shape_h */
