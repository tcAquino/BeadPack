//
//  Grid.h
//  Grid
//
//  Created by Tomas Aquino on 1/11/19.
//  Copyright © 2019 Tomas Aquino. All rights reserved.
//

#ifndef Grid_h
#define Grid_h

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <nanoflann.hpp>
#include "general/MultiArray.h"
#include "general/Operations.h"
#include "general/useful.h"
#include "Geometry/Shape.h"

namespace grid
{
  template
  <typename value_t,
  typename Position_t = std::vector<double>,
  typename Index_t = std::vector<std::size_t>>
  class RegularGrid : public useful::MultiArray<value_t>
  {
  public:
    using value_type = value_t;
    using index_type = Index_t;
    using position_type = Position_t;

    const position_type lengths;
    const position_type corner;

    RegularGrid
    (index_type const& nr_grid_points, position_type lengths,
     position_type corner, value_type value = {})
    : useful::MultiArray<value_t>(nr_grid_points, value)
    , lengths{ lengths }
    , corner{ corner }
    , cell_centers{ compute_cell_centers() }
    {}

    std::size_t dim() const
    { return this->rank(); }

    std::size_t cell(position_type const& position) const
    {
      operation::minus(position, corner, position_helper);
      operation::div(position_helper, lengths, index);
      return this->computeIndex(index);
    }
    
    template <typename Index = useful::Empty>
    double cell_volume(Index idx = {}) const
    { return volume; }

    position_type cell_center(std::size_t cell) const
    { return cell_centers[cell]; }

    double cell_center(std::size_t cell, std::size_t dd) const
    { return cell_centers[cell][dd]; }

    std::size_t neighbor_cell_up
    (index_type const& indexes, std::size_t dd,
     std::size_t offset = 1) const
    {
      index = indexes;
      index[dd] = index[dd] + offset >= this->size(dd)
      ? offset - (this->size(dd) - index[dd])
      : index[dd] + offset;
      
      return this->computeIndex(indexes);
    }

    std::size_t neighbor_cell_down
    (index_type const& indexes,
     std::size_t dd, std::size_t offset = 1) const
    {
      index = indexes;
      index[dd] = index[dd] - offset > this->size(dd)
      ? this->size(dd) - (offset - index[dd])
      : index[dd] - offset;
      
      return this->computeIndex(indexes);
    }

    std::size_t neighbor_cell_up
    (position_type const& position, std::size_t dim,
     std::size_t offset = 1) const
    {
      operation::minus(position, corner, position_helper);
      operation::div(position_helper, lengths, index);
      return neighbor_cell_up(index, dim, offset);
    }

    std::size_t neighbor_cell_down
    (position_type const& position, std::size_t dim,
     std::size_t offset = 1) const
    {
      operation::minus(position, corner, position_helper);
      operation::div(position_helper, lengths, index);
      return neighbor_cell_down(index, dim, offset);
    }
    
    std::size_t neighbor_cell_up
    (std::size_t idx, std::size_t dim,
     std::size_t offset = 1) const
    {
      this->computeIndexes(idx, index);
      return neighbor_cell_up(index, dim, offset);
    }

    std::size_t neighbor_cell_down
    (std::size_t idx,
     std::size_t dim, std::size_t offset = 1) const
    {
      this->computeIndexes(idx, index);
      return neighbor_cell_down(index, dim, offset);
    }

  private:
    mutable index_type index{ index_type(corner.size()) };
    mutable position_type position_helper{
      position_type(corner.size()) };
    const std::vector<position_type> cell_centers;
    double volume{ operation::prod(lengths) };

    auto compute_cell_centers()
    {
      auto centers = decltype(cell_centers){};
      std::vector<double> corner_plus_half_cell{
        operation::plus(corner,
          operation::div_scalar(lengths, 2.)) };
      for (std::size_t cell = 0; cell < this->size(); ++cell)
      {
        this->computeIndexes(cell, index);
        centers.push_back(operation::plus(operation::times(
          lengths, index), corner_plus_half_cell));
      }
      return centers;
    }
  };
  
  template <typename Grid, std::size_t dim>
  class KDTree_Grid
  {
  private:
    Grid const& grid;
    
  public:
    using KDTreeAdaptor = nanoflann::KDTreeSingleIndexAdaptor<
      nanoflann::L2_Simple_Adaptor<double, KDTree_Grid>,
      KDTree_Grid, dim, std::size_t>;
    using KDTreeParams = nanoflann::KDTreeSingleIndexAdaptorParams;
    using SearchParams = nanoflann::SearchParams;
    const SearchParams kdtree_search_params;
    KDTreeAdaptor kdtree;
    
    using Pair = std::pair<std::size_t, double>;
    
    KDTree_Grid
    (Grid const& grid,
     SearchParams kdtree_search_params = { 32, 0., false },
     std::size_t kdtree_leaf_max_size = 10)
    : grid{ grid }
    , kdtree_search_params{ kdtree_search_params }
    , kdtree{ dim, *this, KDTreeParams{ kdtree_leaf_max_size } }
    { kdtree.buildIndex(); }
    
    std::size_t size() const
    { return grid.size(); }
    
    template <typename Position>
    Pair nearest_neighbor
    (Position const& position) const
    {
      std::pair<std::size_t, double> idx_dist_sq;
      kdtree.knnSearch(&(position[0]), 1, &idx_dist_sq.first,
        &idx_dist_sq.second);
      
      return idx_dist_sq;
    }
    
    template <typename Position>
    std::size_t radiusSearch
    (Position const& position, double radius_sq,
     std::vector<Pair> &indices_dists_sq) const
    {
      return kdtree.radiusSearch(&(position[0]), radius_sq,
        indices_dists_sq, kdtree_search_params);
    }
    
    template <typename Shape>
    void inside
    (Shape const& shape, std::vector<std::size_t>& indices) const
    {
      inside(shape, indices, useful::Selector_t<Shape>{});
    }
    
    template <typename Shape>
    void outside
    (Shape const& shape, std::vector<std::size_t>& indices) const
    {
      return outside(shape, indices, useful::Selector_t<Shape>{});
    }
    
    template <typename Shape>
    void inside
    (Shape const& shape, std::vector<std::size_t>& indices,
     useful::Selector_t<geometry::Parallelepiped<>>) const
    {
      std::vector<Pair> near_center;
      double radius_sq = operation::abs_sq(shape.half_dimensions);
      radiusSearch(shape.center, radius_sq, near_center);
      
      for (auto const& point : near_center)
        if (shape.inside(grid.cell_center(point.first)))
          indices.push_back(point.first);
    }
     
    template <typename Shape>
    void outside
    (Shape const& shape, std::vector<std::size_t>& indices,
     useful::Selector_t<geometry::Parallelepiped<>>) const
    {
      for (std::size_t idx = 0; idx < grid.size(); ++idx)
        if (!shape.inside(grid.cell_center(idx)))
          indices.push_back(idx);
    }
    
    template <typename Shape>
    void inside
    (Shape const& shape, std::vector<std::size_t>& indices,
     useful::Selector_t<geometry::Sphere<>>) const
    {
      std::vector<std::pair<size_t,double>> near_center;
      double radius_sq = shape.radius*shape.radius;
      radiusSearch(shape.center,
                   radius_sq, near_center);
      std::sort(near_center.begin(), near_center.end());
      
      for (auto const& point : near_center)
        indices.push_back(point.first);
    }
    
    template <typename Shape>
    void outside
    (Shape const& shape, std::vector<std::size_t>& indices,
     useful::Selector_t<geometry::Sphere<>>) const
    {
      std::vector<std::pair<size_t,double>> near_center;
      double radius_sq = shape.radius*shape.radius;
      radiusSearch(shape.center, radius_sq, near_center);
      std::sort(near_center.begin(), near_center.end());

      for (std::size_t idx = 0; idx < grid.size(); ++idx)
        if (!useful::contains(near_center, idx,
            [](std::pair<size_t,double> const& elem, std::size_t val)
            { return elem.first < val; },
            [](std::pair<size_t,double> const& elem, std::size_t val)
            { return elem.first == val; }))
          indices.push_back(idx);
    }
    
    // Methods for kdtree class
    
    size_t kdtree_get_point_count() const
    { return grid.size(); }

    double kdtree_get_pt(const size_t idx, std::size_t dd) const
    { return grid.cell_center(idx, dd); }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX &bb) const
    { return false; }
  };
  
  template<typename Grid, std::size_t dim>
  class KDTree_Grid_Mask
  {
  private:
    Grid const& grid;
    std::vector<std::size_t> const& mask;
    
  public:
    using KDTreeAdaptor = nanoflann::KDTreeSingleIndexAdaptor<
      nanoflann::L2_Simple_Adaptor<double, KDTree_Grid_Mask>,
      KDTree_Grid_Mask, dim, std::size_t>;
    using KDTreeParams = nanoflann::KDTreeSingleIndexAdaptorParams;
    using SearchParams = nanoflann::SearchParams;
    const SearchParams kdtree_search_params{};
    KDTreeAdaptor kdtree{ dim, *this, KDTreeParams{} };
    
    KDTree_Grid_Mask
    (Grid const& grid, std::vector<std::size_t> const& mask,
     SearchParams kdtree_search_params = { 0, 0., false })
    : grid{ grid }
    , mask{ mask }
    { kdtree.buildIndex(); }
    
    std::size_t size() const
    { return mask.size(); }
    
    auto const& points() const
    { return mask; }
    
    std::size_t points(std::size_t idx) const
    { return mask[idx]; }
    
    std::pair<std::size_t, double> nearest_neighbor
    (std::vector<double> const& position) const
    {
      std::pair<std::size_t, double> idx_dist_sq;
      nanoflann::KNNResultSet<double> resultSet(1);
      resultSet.init(&idx_dist_sq.first, &idx_dist_sq.second);
      kdtree.findNeighbors(resultSet, &position[0],
        kdtree_search_params);
      
      return idx_dist_sq;
    }
    
    std::size_t radiusSearch
    (std::vector<double> const& position, double radius_sq,
     std::vector<std::pair<std::size_t, double>> &indices_dists_sq) const
    {
      return kdtree.radiusSearch(&(position[0]),
        radius_sq, indices_dists_sq, kdtree_search_params);
    }
    
    // Methods for kdtree class
    
    size_t kdtree_get_point_count() const
    { return mask.size(); }

    double kdtree_get_pt(const size_t idx, std::size_t dd) const
    { return grid.cell_center(mask[idx], dd); }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX &bb) const
    { return false; }
  };
  
  // Identify boundary positions along a dimension
  // Values according to neighbors to "left" and "right"
  // 0: void - void
  // 1: solid - void
  // 2: void - solid
  // 3: solid - solid
  template <typename Grid>
  int void_type
  (Grid const& grid,
   std::vector<std::size_t> const& solid_idx,
   std::size_t idx,
   std::size_t dd)
  {
    std::size_t idx_left = grid.neighbor_cell_down(idx, dd);
    std::size_t idx_right = grid.neighbor_cell_up(idx, dd);
    
    bool left_is_solid = useful::contains(solid_idx, idx_left);
    bool right_is_solid = useful::contains(solid_idx, idx_right);
    
    if (!left_is_solid && !right_is_solid)
      return 0;
    if (left_is_solid && !right_is_solid)
      return 1;
    if (!left_is_solid && right_is_solid)
      return 2;
    return 3;
  }

  template <typename Grid>
  auto void_type
  (Grid const& grid,
   std::vector<std::size_t> const& solid_idx,
   std::size_t idx)
  {
    std::vector<int> type;
    type.reserve(grid.dim());
    for (std::size_t dd = 0; dd < grid.dim(); ++dd)
      type.push_back(void_type(grid, solid_idx, idx, dd));
      
    return type;
  }
  
  template <typename Grid>
  auto void_type
  (Grid const& grid,
   std::vector<std::size_t> const& void_idx,
   std::vector<std::size_t> const& solid_idx)
  {
    std::vector<std::vector<int>> type(grid.size());
    
    for (auto idx : void_idx)
      type[idx] = void_type(grid, solid_idx, idx);
    
    return type;
  }
  
  template <typename Grid>
  void print_centers
  (Grid const& grid, std::vector<std::size_t> const& points,
   std::string const& filename, int precision = 8,
   std::string delimiter = "\t")
  {
    std::ofstream output{ filename };
    if (!output.is_open())
      throw useful::open_write_error(filename);
    output << std::setprecision(precision);
    output << std::scientific;
    for (auto idx : points)
    {
      std::string delim = "";
      for (std::size_t dd = 0; dd < grid.dim(); ++dd)
      {
        output << delim << grid.cell_center(idx, dd);
        delim = delimiter;
      }
      output << "\n";
    }
    output.close();
  }
  
  class Grid_void_solid
  {
  public:
    template <typename Domain, typename KDTree>
    Grid_void_solid(Domain const& domain, KDTree const& kdtree_grid)
    {
      for (auto const& shape : domain.parallelepipeds)
        kdtree_grid.inside(shape, solid_idx);
      for (auto const& shape : domain.spheres)
        kdtree_grid.inside(shape, solid_idx);
      kdtree_grid.outside(domain.box, solid_idx);
      std::sort(solid_idx.begin(), solid_idx.end());
      
      for (std::size_t idx = 0; idx < kdtree_grid.size(); ++idx)
        if (!useful::contains(solid_idx, idx))
          void_idx.push_back(idx);
    }
    
    std::size_t nr_voids() const
    { return void_idx.size(); }
    
    std::size_t nr_solids() const
    { return solid_idx.size(); }
    
    std::size_t voids(std::size_t idx) const
    { return void_idx[idx]; }
    
    std::size_t solids(std::size_t idx) const
    { return solid_idx[idx]; }
    
    auto const& voids() const
    { return void_idx; }
    
    auto const& solids() const
    { return solid_idx; }
    
  private:
    std::vector<std::size_t> void_idx;
    std::vector<std::size_t> solid_idx;
  };
}

#endif /* Grid_h */

