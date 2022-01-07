//
//  BeadPack.h
//  BeadPack
//
//  Created by Tomás Aquino on 12/05/2020.
//  Copyright © 2020 Tomás Aquino. All rights reserved.
//

#ifndef BeadPack_h
#define BeadPack_h

#include <algorithm>
#include <nanoflann.hpp>
#include <random>
#include <vector>
#include "general/Constants.h"
#include "general/Operations.h"
#include "Geometry/Coordinates.h"
#include "Geometry/Shape.h"

namespace beadpack
{
  // Handles a beadpack (set of spherical beads) in any spatial dimension
  template <std::size_t dimension>
  class BeadPack
  {
  private:
    using KDTreeAdaptor = nanoflann::KDTreeSingleIndexAdaptor<
      nanoflann::L2_Simple_Adaptor<double, BeadPack>,
      BeadPack, dimension, std::size_t>;           // Type for KDTree searching
    using KDTreeParams =
      nanoflann::KDTreeSingleIndexAdaptorParams;   // Parameter type for KDTree algorithm
    using SearchParams = nanoflann::SearchParams;  // Parameter type for KDTree searching
    const SearchParams kdtree_search_params;       // Parameter object for KDTree searching
    
  public:
    using Bead = geometry::Sphere<>;                // Spherical bead type (any spatial dimension)
    static constexpr std::size_t dim{ dimension };  // Spatial dimension
    
    // Construct object given bead centers, all of same given radius
    // Optional KDTree search and algorithm parameters may be passed
    BeadPack
    (std::vector<std::vector<double>> const& centers, double radius,
     SearchParams kdtree_search_params = { 32, 0., false },
     KDTreeParams params = { 10 })
    : max_radius{ radius }
    , min_radius{ radius }
    , kdtree{ dimension, *this, params }
    {
      bead_container.reserve(centers.size());
      for (auto const& center : centers)
        bead_container.push_back({ center, radius });
        
      kdtree.buildIndex();
    }
    
    // Construct object given bead objects
    // Optional KDTree search and algorithm parameters may be passed
    BeadPack
    (std::vector<Bead> beads,
     SearchParams kdtree_search_params = { 32, 0., false },
     std::size_t kdtree_leaf_max_size = 10)
    : bead_container{ beads }
    , kdtree{ dimension, *this, KDTreeParams{ kdtree_leaf_max_size } }
    {
      if (!beads.empty())
      {
        auto radius_comp = [](Bead const& s1, Bead const& s2)
        { return s1.radius < s2.radius; };
        max_radius = std::max_element(bead_container.cbegin(),
          bead_container.cend(), radius_comp)->radius;
        min_radius = std::min_element(bead_container.cbegin(),
          bead_container.cend(), radius_comp)->radius;
      }
        
      kdtree.buildIndex();
    }
    
    // Get number of beads
    std::size_t nr_beads() const
    { return bead_container.size(); }
    
    // Get reference to bead
    Bead const& bead(std::size_t idx) const
    { return bead_container[idx]; }
    
    // Get reference to bead container
    auto const& beads() const
    { return bead_container; }
    
    // Get radius of bead
    double radius(std::size_t idx) const
    { return bead_container[idx].radius; }
    
    // Get center of bead
    auto center(std::size_t idx) const
    { return bead_container[idx].center; }
    
    // Check whether a position is inside a bead
    // Returns a pair, first element is 1 if inside a bead and 0 otherwise,
    // second element is the first bead it is inside of or 0 if not inside any
    template <typename Position>
    std::pair<bool, std::size_t> inside(Position const& position) const
    {
      std::vector<std::pair<size_t,double>> near_center;
      radiusSearch(position, max_radius*max_radius, near_center);
      
      for (auto const& point : near_center)
      {
        double radius_val = radius(point.first);
        if (point.second < radius_val*radius_val)
          return { 1, point.first };
      }
      
      return { 0, 0 };
    }
    
    // Check whether a position is near (within distance) of a bead
    // Returns a pair, first element is 1 if near any bead and 0 otherwise,
    // second element is the first bead it is near of or 0 if not near any
    template <typename Position>
    std::pair<bool, std::size_t> near
    (Position const& position, double distance) const
    {
      std::vector<std::pair<size_t,double>> near_center;
      double max_distance = max_radius + distance;
      radiusSearch(position, max_distance*max_distance, near_center);
      
      for (auto const& point : near_center)
      {
        double distance_from_center = radius(point.first) + distance;
        if (point.second < distance_from_center*distance_from_center)
          return { 1, point.first };
      }
      
      return { 0, 0 };
    }
    
    // Place position at the closest point on the given bead's surface
    template <typename Position>
    void place_at_surface
    (Position& position, std::size_t bead) const
    {
      auto radial_vector = operation::minus(position, center(bead));
      double distance_to_center = operation::abs(radial_vector);
      if (distance_to_center != 0.)
      {
        double distance_fraction = radius(bead)/distance_to_center;
        operation::times_scalar_InPlace(distance_fraction, radial_vector);
        operation::plus(center(bead), radial_vector, position);
      }
      // If at bead center, place at surface along first dimension
      else
      {
        std::vector<double> displacement(position.size(), 0.);
        displacement[0] = radius(bead);
        operation::plus(center(bead), displacement, position);
      }
    }
    
    // Place position at given distance from the closest point on the given bead's surface
    // Warning: Does not check if inside another bead
    template <typename Position>
    void place_near_surface
    (Position& position, std::size_t bead, double distance) const
    {
      auto radial_vector = operation::minus(position, center(bead));
      double distance_to_center = operation::abs(radial_vector);
      if (distance_to_center != 0.)
      {
        double distance_fraction = (radius(bead)+distance)/distance_to_center;
        operation::times_scalar_InPlace(distance_fraction, radial_vector);
        operation::plus(center(bead), radial_vector, position);
      }
      // If at bead center, place at surface along first dimension
      else
      {
        std::vector<double> displacement(position.size(), 0.);
        displacement[0] = radius(bead)+distance;
        operation::plus(center(bead), displacement, position);
      }
    }
    
    // Place state.position at the closest point on the given bead's surface
    // Applies boundary condition (e.g. periodic)
    template <typename State, typename Boundary>
    void place_at_surface
    (State& state, Boundary const& boundary, std::size_t bead) const
    {
      place_at_surface(state.position);
      boundary(state);
    }
    
    // Check if position is inside a bead
    // If it is, place at the closest surface of that bead
    // Return a pair with
    // first element true if position was inside bead, false otherwise
    // and second element bead index if inside and zero othwerwise
    template <typename Position>
    std::pair<std::size_t, std::size_t> place_at_closest_surface_if_inside
    (Position& position) const
    {
      auto inside_bead = inside(position);
      
      if (inside_bead.first)
        place_at_surface(position, inside_bead.second);
      
      return { inside_bead.first, inside_bead.second };
    }
    
    // Check if state.position is inside a bead
    // If it is, place at the closest surface of that bead
    // Applies boundary condition (e.g. periodic)
    // Return a pair with
    // first element true if position was inside bead, false otherwise
    // and second element bead index if inside and zero othwerwise
    template <typename State, typename Boundary>
    std::pair<std::size_t, std::size_t> place_at_closest_surface_if_inside
    (State& state, Boundary const& boundary) const
    {
      auto inside = place_at_closest_surface_if_inside(state.position);
      if (inside.first)
        boundary(state);
      return inside;
    }
    
    // Place position at the closest bead surface
    // Return the bead index
    // Applies boundary condition (e.g. periodic)
    template <typename Position>
    std::size_t place_at_closest_surface(Position& position) const
    {
      auto inside_bead = inside(position);
      if (inside_bead.first)
      {
        place_at_surface(position, inside_bead.second);
        return inside_bead.second;
      }
      else
      {
        std::size_t nearest = nearest_surface(position).first;
        place_at_surface(position, nearest);
        return nearest;
      }
    }
    
    // Place state.position at the closest bead surface
    // Return the bead index
    // Applies boundary condition (e.g. periodic)
    template <typename State, typename Boundary>
    std::size_t place_at_closest_surface
    (State& state, Boundary const& boundary) const
    {
      std::size_t bead = place_at_closest_surface(state.position);
      boundary(state);
      
      return bead;
    }
    
    // Compute mean value of a vector field outside of beads
    // in region with boundaries along each dimension
    template
    <typename Field, typename Boundary, typename Container = std::vector<double>>
    Container compute_mean_vector
    (Field const& field,
     std::vector<std::pair<double,double>> const& boundaries,
     std::size_t nr_samples) const
    {
      Container mean(dim);
      compute_mean(field, boundaries, nr_samples, mean);
      return mean;
    }
    
    // Compute mean value of a scalar field outside of beads
    // in region with boundaries along each dimension
    template <typename Field, typename Value = double>
    Value compute_mean
    (Field const& field,
     std::vector<std::pair<double,double>> const& boundaries,
     std::size_t nr_samples) const
    {
      Value mean{};
      compute_mean(field, boundaries, nr_samples, mean);
      return mean;
    }
    
    // Compute mean value of a field outside of beads
    // in region with boundaries along each dimension
    // Place result in mean
    template <typename Field, typename Container>
    void compute_mean
    (Field const& field,
     std::vector<std::pair<double,double>> const& boundaries,
     std::size_t nr_samples, Container& mean) const
    {
      std::mt19937 rng{ std::random_device{}() };
      std::vector<std::uniform_real_distribution<double>> dist;
      dist.reserve(dim);
      for (std::size_t dd = 0; dd < dim; ++dd)
        dist.emplace_back(boundaries[dd].first,
                          boundaries[dd].second);
      
      std::vector<double> position(dim);
      std::size_t samples = 0;
      while (samples < nr_samples)
      {
        for (std::size_t dd = 0; dd < dim; ++dd)
          position[dd] = dist[dd](rng);
        if (!inside(position).first)
        {
          operation::plus_InPlace(mean, field(position));
          ++samples;
        }
      }
      
      operation::div_scalar_InPlace(mean, double(nr_samples));
    }
    
    // Compute porosity in region with boundaries along each dimension
    double compute_porosity
    (std::vector<std::pair<double,double>> const& boundaries,
     std::size_t nr_samples) const
    {
      std::mt19937 rng{ std::random_device{}() };
      std::vector<std::uniform_real_distribution<double>> dist;
      dist.reserve(dim);
      for (std::size_t dd = 0; dd < dim; ++dd)
        dist.emplace_back(boundaries[dd].first,
                          boundaries[dd].second);
      
      std::vector<double> position(dim);
      std::size_t samples = 0;
      std::size_t void_samples = 0;
      while (samples < nr_samples)
      {
        for (std::size_t dd = 0; dd < dim; ++dd)
          position[dd] = dist[dd](rng);
        if (!inside(position).first)
          ++void_samples;
        ++samples;
      }
      
      return void_samples/double(samples);
    }
    
    // Find bead with center nearest to position
    // Return a pair, first element is the bead index,
    // second element is the distance squared
    template <typename Position>
    std::pair<std::size_t, double> nearest_neighbor
    (Position const& position) const
    {
      std::pair<std::size_t, double> idx_dist_sq;
      kdtree.knnSearch(&(position[0]), 1,
                       &idx_dist_sq.first, &idx_dist_sq.second);
      
      return idx_dist_sq;
    }
    
    // Find beads with center with center within given squared distance
    // Results placed a vector of pairs, first elements are the bead index,
    // second elements are the distance squared
    // Returns the number of beads within the distance
    template <typename Position>
    std::size_t radiusSearch
    (Position const& position, double distance_sq,
     std::vector<std::pair<std::size_t, double>> &indices_dists_sq) const
    {
      return kdtree.radiusSearch(&(position[0]), distance_sq,
                                 indices_dists_sq, kdtree_search_params);
    }
    
    // Find bead with surface nearest to position
    // Return a pair, first element is the bead index,
    // second element is the distance squared to its surface
    template <typename Position>
    std::pair<std::size_t, double> nearest_surface
    (Position const& position) const
    {
      // Find distance to surface of bead with closest center
      auto closest_center = nearest_neighbor(position);
      double dist_closest_center = std::sqrt(closest_center.second);
      double dist_surface = dist_closest_center
        - radius(closest_center.first);
      
      // If all radii are the same we're done
      if (max_radius == min_radius)
        return { closest_center.first, dist_surface*dist_surface };
      
      // Find all beads with centers closer than nearest neighbor plus max radius
      double max_dist = dist_surface + max_radius;
      std::vector<std::pair<std::size_t, double>> neighbors;
      radiusSearch(&(position[0]), max_dist*max_dist, neighbors);
      
      // Find beads with closest surface and distance squared to its center
      auto closest_surface =
        std::min_element(neighbors.cbegin(), neighbors.cend(),
          [this](std::pair<std::size_t, double> neighbor1,
                 std::pair<std::size_t, double> neighbor2)
          { return std::sqrt(neighbor1.second)-radius(neighbor1.first)
            < std::sqrt(neighbor2.second)-radius(neighbor2.first); });
      std::size_t closest_bead = closest_surface->first;
      
      // Get distance to closest surface
      double dist_surface_closest = std::sqrt(closest_surface->second)
        - radius(closest_bead);
      
      return { closest_surface->first,
        dist_surface_closest*dist_surface_closest };
    }
    
    // Method required by kdtree class, see nanoflann documentation
    std::size_t kdtree_get_point_count() const
    { return nr_beads(); }

    // Method required by kdtree class, see nanoflann documentation
    double kdtree_get_pt(const size_t idx, std::size_t dd) const
    { return bead_container[idx].center[dd]; }

    // Method required by kdtree class, see nanoflann documentation
    template <class BBOX>
    bool kdtree_get_bbox(BBOX &bb) const
    { return false; }
    
  private:
    std::vector<Bead> bead_container;  // Container of beads
    double max_radius;                 // Maximum bead radius
    double min_radius;                 // Minimum bead radius
    KDTreeAdaptor kdtree;              // KDTree for bead searching
  };
}

#endif /* BeadPack_h */
