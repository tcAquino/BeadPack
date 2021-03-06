//
//  SymmetryPlanes.h
//
//  Created by Tomás Aquino on 13/05/2020.
//  Copyright © 2020 Tomás Aquino. All rights reserved.
//

#ifndef SymmetryPlanes_h
#define SymmetryPlanes_h

#include <cmath>
#include "general/Operations.h"

namespace geometry
{
  // A SymmetryPlanes object must define:
  // dim : the spatial dimension
  // length : the length of the unit cell along the normal vector associated with each plane
  // normal : a vector of unit normal vectors associated with each plane
  // translation : a vector of translation vectors associated with each plane
  
  // Symmetry plane for bcc packing of unit radius
  struct SymmetryPlanes_Bcc
  {
    // Construct given bead radius
    SymmetryPlanes_Bcc()
    : normal{
        operation::times_scalar(1./std::sqrt(2.),
                                std::vector<double>{ 0., 1., 1. }),
        operation::times_scalar(1./std::sqrt(2.),
                                std::vector<double>{ 1., 0., 1. }),
        operation::times_scalar(1./std::sqrt(2.),
                                std::vector<double>{ 1., 1., 0. })
    }
    , translation{
        operation::times_scalar(2./std::sqrt(3.),
                                std::vector<double>{ -1., 1., 1. }),
        operation::times_scalar(2./std::sqrt(3.),
                                std::vector<double>{ 1., -1., 1. }),
        operation::times_scalar(2./std::sqrt(3.),
                                std::vector<double>{ 1., 1., -1. })
    }
    {}
    
    static constexpr std::size_t dim{ 3 };
    const std::vector<double> length{
      std::vector<double>(dim, 4./std::sqrt(6.)) };
    const std::vector<std::vector<double>> normal;
    const std::vector<std::vector<double>> translation;
  };
  
  // Symmetry plane for sc packing of unit radius
  struct SymmetryPlanes_Sc
  {
    // Construct given bead radius
    SymmetryPlanes_Sc()
    : normal{
        { 1., 0., 0. },
        { 0., 1., 0. },
        { 0., 0., 1. }
    }
    , translation{
        std::vector<double>{ 2., 0., 0. },
        std::vector<double>{ 0., 2., 0. },
        std::vector<double>{ 0., 0., 2. }
    }
    {}
    
    static constexpr std::size_t dim{ 3 };
    const std::vector<double> length{
      std::vector<double>(dim, 1.) };
    const std::vector<std::vector<double>> normal;
    const std::vector<std::vector<double>> translation;
  };
  
  template <typename Position, typename SymmetryPlanes>
  auto project
  (Position const& position, SymmetryPlanes const& symmetry_planes,
   std::size_t dd, double scale)
  {
    return operation::dot(position, symmetry_planes.normal[dd])
      /(scale*symmetry_planes.length[dd]);
  }
  
  template <typename Position, typename SymmetryPlanes>
  auto project
  (Position const& position, SymmetryPlanes const& symmetry_planes,
   double scale)
  {
    Position projections;
    projections.reserve(symmetry_planes.dim);
    for (std::size_t dd = 0; dd < symmetry_planes.dim; ++dd)
      projections.push_back(project(position, symmetry_planes,
                                    dd, scale));
    return projections;
  }
  
  template <typename Position, typename SymmetryPlanes, typename Origin>
  auto project
  (Position const& position, SymmetryPlanes const& symmetry_planes,
   std::size_t dd, double scale, Origin const& origin)
  {
    return project(operation::minus(position, origin),
                   symmetry_planes, dd, scale);
  }
  
  template <typename Position, typename SymmetryPlanes, typename Origin>
  auto project
  (Position const& position, SymmetryPlanes const& symmetry_planes,
   double scale, Origin const& origin)
  {
    return project(operation::minus(position, origin),
                   symmetry_planes, scale);
  }
  
  template <typename Position, typename SymmetryPlanes,
  typename Projections = std::vector<double>>
  auto translate
  (Position& position, SymmetryPlanes const& symmetry_planes,
   Projections const& projections,
   double scale)
  {
    for (std::size_t dd = 0; dd < symmetry_planes.dim; ++dd)
      operation::plus_InPlace(position,
        operation::times_scalar(scale*projections[dd],
                                symmetry_planes.translation[dd]));
    
    return projections;
  }
  
  template <typename Position, typename SymmetryPlanes,
  typename Projections = std::vector<double>>
  auto translate_back
  (Position& position, SymmetryPlanes const& symmetry_planes,
   Projections const& projections,
   double scale)
  {
    for (std::size_t dd = 0; dd < symmetry_planes.dim; ++dd)
      operation::minus_InPlace(position,
        operation::times_scalar(scale*projections[dd],
                                symmetry_planes.translation[dd]));
    
    return projections;
  }
  
  template <typename Position, typename SymmetryPlanes>
  auto place_in_unit_cell
  (Position& position, SymmetryPlanes const& symmetry_planes,
   double scale)
  {
    std::vector<int> projections;
    projections.reserve(symmetry_planes.dim);
    for (std::size_t dd = 0; dd < symmetry_planes.dim; ++dd)
      projections.push_back(std::floor(project(position, symmetry_planes,
                                               dd, scale)));
    translate_back(position, symmetry_planes, projections, scale);

    return projections;
  }
  
  template <typename Position, typename SymmetryPlanes, typename Origin>
  auto place_in_unit_cell
  (Position& position, SymmetryPlanes const& symmetry_planes,
   double scale, Origin const& origin)
  {
    std::vector<int> projections;
    projections.reserve(symmetry_planes.dim);
    for (std::size_t dd = 0; dd < symmetry_planes.dim; ++dd)
      projections.push_back(std::floor(project(position, symmetry_planes,
                                               dd, scale, origin)));
    translate_back(position, symmetry_planes, projections, scale);

    return projections;
  }
  
  // Add periodic velocity point images outside domain
  // if the projection of a point is within a given fraction of
  // the boundary along each symmetry plane, add an image along that plane
  template
  <typename Points = std::vector<std::vector<double>>,
   typename Velocities = std::vector<std::vector<double>>,
   typename Boundary_Periodic>
  void add_periodic_images
  (Points& points, Velocities& velocities,
   double fraction, Boundary_Periodic const& boundary_periodic)
  {
    std::size_t nr_points = points.size();
    if (nr_points == 0)
      return;
    
    std::size_t dim = points.back().size();
    std::vector<std::vector<int>>
      projections_forward(dim, std::vector<int>(dim));
    std::vector<std::vector<int>>
      projections_backward(dim, std::vector<int>(dim));
    for (std::size_t dd = 0; dd < dim; ++dd)
    {
      projections_forward[dd][dd] = 1.;
      projections_backward[dd][dd] = -1.;
    }
    
    for (std::size_t pp = 0; pp < nr_points; ++pp)
    {
      auto projections =
        geometry::project(points[pp],
                          boundary_periodic.symmetry_planes,
                          boundary_periodic.scale,
                          boundary_periodic.origin);
      for (std::size_t dd = 0; dd < dim; ++dd)
      {
        // If close to lower bound, make an image forward
        if (projections[dd] < fraction)
        {
          points.push_back(points[pp]);
          velocities.push_back(velocities[pp]);
          boundary_periodic.translate(points.back(),
                                      projections_forward[dd]);
        }
        // If close to upper bound, make an image backward
        if (projections[dd] > 1.-fraction)
        {
          points.push_back(points[pp]);
          velocities.push_back(velocities[pp]);
          boundary_periodic.translate(points.back(),
                                      projections_backward[dd]);
        }
      }
    }
  }
}


#endif /* SymmetryPlanes_h */
