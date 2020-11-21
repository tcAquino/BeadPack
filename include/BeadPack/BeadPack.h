//
//  BeadPack.h
//  BeadPack
//
//  Created by Tomás Aquino on 12/05/2020.
//  Copyright © 2020 Tomás Aquino. All rights reserved.
//

#ifndef BeadPack_h
#define BeadPack_h

#include "general/Operations.h"
#include "Geometry/Shape.h"
#include <nanoflann.hpp>
#include <random>
#include <vector>

namespace beadpack
{
  template <std::size_t dimension>
  class BeadPack
  {
  private:
    using KDTreeAdaptor = nanoflann::KDTreeSingleIndexAdaptor<
      nanoflann::L2_Simple_Adaptor<double, BeadPack>,
      BeadPack, dimension, std::size_t>;
    using KDTreeParams = nanoflann::KDTreeSingleIndexAdaptorParams;
    using SearchParams = nanoflann::SearchParams;
    const SearchParams kdtree_search_params;
    
  public:
    using Bead = geometry::Sphere<>;
    static constexpr std::size_t dim{ dimension };
    
    BeadPack
    (std::vector<std::vector<double>> const& centers, double radius,
     SearchParams kdtree_search_params = { 32, 0., false },
     std::size_t kdtree_leaf_max_size = 10)
    : max_radius{ radius }
    , min_radius{ radius }
    , kdtree{ dimension, *this, KDTreeParams{ kdtree_leaf_max_size } }
    {
      bead_container.reserve(centers.size());
      for (auto const& center : centers)
        bead_container.push_back({ center, radius });
        
      kdtree.buildIndex();
    }
    
    BeadPack
    (std::vector<Bead> beads,
     SearchParams kdtree_search_params = { 32, 0., false },
     std::size_t kdtree_leaf_max_size = 10)
    : bead_container{ beads }
    , kdtree{ dimension, *this, KDTreeParams{ kdtree_leaf_max_size } }
    {
      auto radius_comp = [](Bead const& s1, Bead const& s2)
      { return s1.radius < s2.radius; };
      max_radius = std::max_element(bead_container.cbegin(),
                                    bead_container.cend(), radius_comp)->radius;
      min_radius = std::min_element(bead_container.cbegin(),
                                    bead_container.cend(), radius_comp)->radius;
        
      kdtree.buildIndex();
    }
    
    std::size_t nr_beads() const
    { return bead_container.size(); }
    
    Bead const& bead(std::size_t idx) const
    { return bead_container[idx]; }
    
    auto const& beads() const
    { return bead_container; }
    
    double radius(std::size_t idx) const
    { return bead_container[idx].radius; }
    
    auto center(std::size_t idx) const
    { return bead_container[idx].center; }
    
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
    
    template <typename Position>
    bool place_at_closest_surface(Position& position) const
    {
      auto inside_bead = inside(position);
      
      if (inside_bead.first)
      {
        auto radial_vector =
          operation::minus(position, center(inside_bead.second));
        double distance_to_center = operation::abs(radial_vector);
        double radius_val = radius(inside_bead.second);
        if (distance_to_center != 0.)
        {
          operation::times_scalar_InPlace(radius_val/distance_to_center, radial_vector);
          operation::plus(bead_container[inside_bead.second].center, radial_vector, position);
        }
        else
        {
          std::vector<double> displacement(position.size(), 0.);
          displacement[0] = radius_val;
          operation::plus_InPlace(position, displacement);
        }
      }
      
      return inside_bead.first;
    }
    
    template <typename State, typename Boundary>
    bool place_at_closest_surface
    (State& state, Boundary const& boundary) const
    {
      bool inside_bead = 1;
      bool outofbounds = 1;
      
      int counter = 0;
      
      while (inside_bead || outofbounds)
      {
        outofbounds = boundary(state);
        inside_bead = place_at_closest_surface(state.position);
        counter += outofbounds || inside_bead;
      }
      
      return counter;
    }
    
    template <typename Field, typename Container = std::vector<double>>
    Container compute_mean_vector
    (Field const& field,
     std::vector<std::pair<double,double>> const& boundaries,
     std::size_t nr_samples) const
    {
      Container mean(dim);
      compute_mean(field, boundaries, nr_samples, mean);
      return mean;
    }
    
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
    
    template <typename Field, typename Container>
    Container compute_mean
    (Field const& field,
     std::vector<std::pair<double,double>> const& boundaries,
     std::size_t nr_samples, Container& mean) const
    {
      std::mt19937 rng{ std::random_device{}() };
      std::vector<std::uniform_real_distribution<double>> dist;
      dist.reserve(dim);
      for (std::size_t dd = 0; dd < dim; ++dd)
        dist.emplace_back(boundaries[dd].first, boundaries[dd].second);
      
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
      
      return mean;
    }
    
    // Interface methods for kdtree searching
    
    template <typename Position>
    std::pair<std::size_t, double> nearest_neighbor
    (Position const& position) const
    {
      std::pair<std::size_t, double> idx_dist_sq;
      kdtree.knnSearch(&(position[0]), 1,
                       &idx_dist_sq.first, &idx_dist_sq.second);
      
      return idx_dist_sq;
    }
    
    template <typename Position>
    std::size_t radiusSearch
    (Position const& position, double radius_sq,
     std::vector<std::pair<std::size_t, double>> &indices_dists_sq) const
    {
      return kdtree.radiusSearch(&(position[0]), radius_sq,
                                 indices_dists_sq, kdtree_search_params);
    }
    
    // Methods required by kdtree class
    
    size_t kdtree_get_point_count() const
    { return nr_beads(); }

    double kdtree_get_pt(const size_t idx, std::size_t dd) const
    { return bead_container[idx].center[dd]; }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX &bb) const
    { return false; }
    
  private:
    std::vector<Bead> bead_container;
    double max_radius;
    double min_radius;
    KDTreeAdaptor kdtree;
  };
}

#endif /* BeadPack_h */
