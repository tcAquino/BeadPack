//
//  VectorField_Interpolated.h
//  Field
//
//  Created by Tomás Aquino on 13/05/2020.
//  Copyright © 2020 Tomás Aquino. All rights reserved.
//

#ifndef VectorField_Interpolated_h
#define VectorField_Interpolated_h

#include <iterator>
#include <type_traits>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/interpolation_functions.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/natural_neighbor_coordinates_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <nanoflann.hpp>

namespace field
{
  // KDTree class for nearest-neighbor searching
  template <typename Data, std::size_t dim>
  class KDTree
  {
  public:
    using KDTreeAdaptor = nanoflann::KDTreeSingleIndexAdaptor<
      nanoflann::L2_Simple_Adaptor<double, KDTree>,
      KDTree, dim, std::size_t>;
    using KDTreeParams = nanoflann::KDTreeSingleIndexAdaptorParams;
    using SearchParams = nanoflann::SearchParams;
    const SearchParams kdtree_search_params;
    KDTreeAdaptor kdtree;
    
  private:
    using Pair = std::pair<std::size_t, double>;
    using Point = std::vector<double>;
    using Value = std::vector<double>;
    std::vector<Point> points;
    std::vector<Value> values;
    
  public:
    // Build and set up kdtree
    // Copy of data is made
    // because CGAL interpolators use a map,
    // which is very slow for this purpose
    KDTree
    (Data const& data,
     SearchParams kdtree_search_params = { 32, 0., false },
     std::size_t kdtree_leaf_max_size = 10)
    : kdtree_search_params{ kdtree_search_params }
    , kdtree{ dim, *this, KDTreeParams{ kdtree_leaf_max_size } }
    {
      points.reserve(data.size());
      values.reserve(data.size());
      for (auto const& val : data)
      {
        points.push_back({});
        points.back().reserve(dim);
        values.push_back({});
        values.back().reserve(dim);
        for (std::size_t dd = 0; dd < dim; ++dd)
        {
          points.back().push_back(val.first[int(dd)]);
          values.back().push_back(val.second[int(dd)]);
        }
      }
      kdtree.buildIndex();
    }
    
    std::size_t size() const
    { return points.size(); }
    
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
  
    template <typename Position>
    auto nearest_neighbor_value
    (Position const& position) const
    {
      return values[nearest_neighbor(position).first];
    }
    
    // Methods for kdtree class
    
    auto kdtree_get_point_count() const
    { return size(); }

    auto kdtree_get_pt(std::size_t idx, std::size_t dd) const
    { return points[idx][dd]; }
                          
    auto
    kdtree_distance(const double *p1, const std::size_t idx_p2, std::size_t size) const
    {
      double dist_sq = 0.;
      for (std::size_t dd = 0; dd < dim; ++dd)
        dist_sq += (p1[dd]-points[idx_p2])*(p1[dd]-points[idx_p2]);

      return dist_sq;
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX &bb) const
    { return false; }
  };
  
  // Template class for linear vector interpolation on unstructured grid
  // to be specialized for different dimensions
  // If the exact parameter is 1, the interpolator uses
  // exact (arbitrary-precision) constructions. This is slower, but necessary
  // for some meshes. If it is 0 (default), inexact (double-precision) constructions are
  // used, which is faster but may result in NaN interpolation values
  // for some meshes
  // If check_interpolation is 1 (meant to be used with inexact construction),
  // the interpolator checks if the construction is valid. If not,
  // it uses the nearest point for velocity instead of interpolation
  template
  <std::size_t dim, bool exact = 0, bool check_interpolation = 0,
  typename Value_type = std::vector<double>>
  class VectorField_LinearInterpolation_UnstructuredGrid;
  
  // Linear vector interpolation on unstructured grid in 3d
  // Uses CGAL library for
  // Delaunay triangulation and Sibson natural neighbor coordinates
  // Value_type is the type of the interpolated field values
  template <bool exact, bool check_interpolation, typename Value_type>
  class VectorField_LinearInterpolation_UnstructuredGrid
  <3, exact, check_interpolation, Value_type>
  {
  private:
    // Types for interpolation, see CGAL library documentation
    using Kernel = typename std::conditional<exact,
      CGAL::Exact_predicates_exact_constructions_kernel,
      CGAL::Exact_predicates_inexact_constructions_kernel>::type;
    using Triangulation =
      CGAL::Delaunay_triangulation_3<Kernel, CGAL::Fast_location>;
    using Cell_handle = typename Triangulation::Cell_handle;
    using Coord_type = typename Kernel::FT;
    using Vertex_handle = typename Triangulation::Vertex_handle;
    
  public:
    static constexpr std::size_t dim = 3;  // Spatial dimension
    
    using Hint = Cell_handle;                  // Hint type to speed up locating a point
    using Point = typename Kernel::Point_3;    // CGAL spatial point type
    using Vector = typename Kernel::Vector_3;  // CGAL Field vector type
    using value_type = Value_type;             // Type of the interpolated field values
    
  private:
    // More types for interpolation, see CGAL library documentation
    using Coord_map = std::map<Point, Vector, typename Kernel::Less_xyz_3>;
    using Value_access = CGAL::Data_access<Coord_map>;
    using Vertex_coordinate_vector =
      std::vector<std::pair<Vertex_handle, Coord_type>>;
    using Point_coordinate_vector =
      std::vector<std::pair<Point, Coord_type>>;
    
    // Tree to find nearest neighbors when check_interpolation is on
    using Tree = typename std::conditional<check_interpolation,
      KDTree<Coord_map, dim>, useful::Empty>::type;
    
  public:
    
    // Construct object given grid points and associated vector field values
    // from vectors of points and vector values in random-access containers
    template <typename Container_point, typename Container_vector>
    VectorField_LinearInterpolation_UnstructuredGrid
    (Container_point const& grid_points,
     Container_vector const& function_values)
    : function_values{ make_coord_map(grid_points, function_values) }
    , tree{ this->function_values }
    {}
    
    // Vector field interpolation value at position given in internal Point type
    // Optional Hint can be given to speed up point triangulation
    auto operator()(Point const& point, Hint const& hint = {}) const
    {
      Vertex_coordinate_vector coords_vertex;
      Coord_type norm;
      CGAL::sibson_natural_neighbor_coordinates_3(
        triangulation, point, std::back_inserter(coords_vertex), norm, hint);
      
      if constexpr (check_interpolation)
      {
        if (norm == 0.)
        {
          count_interpolation_failures++;
          if constexpr(std::is_same<value_type, Vector>::value == true)
            return tree.nearest_neighbor_value(point);
          else
          {
            auto val = tree.nearest_neighbor_value(point);
            return value_type{
              CGAL::to_double(val[0]), CGAL::to_double(val[1]), CGAL::to_double(val[2]) };
          }
        }
      }
      
      Point_coordinate_vector coords_point;
      coords_point.reserve(coords_vertex.size());
      for (auto const& val : coords_vertex)
        coords_point.push_back(std::make_pair(val.first->point(), val.second));
      
      if constexpr(std::is_same<value_type, Vector>::value == true)
        return CGAL::to_double(CGAL::linear_interpolation(
          coords_point.begin(), coords_point.end(), norm,
          Value_access(function_values)));
      else
      {
        auto val = CGAL::linear_interpolation(
          coords_point.begin(), coords_point.end(), norm,
          Value_access(function_values));
        return value_type{
          CGAL::to_double(val[0]), CGAL::to_double(val[1]), CGAL::to_double(val[2]) };
      }
    }
    
    // Vector field interpolation value at position given in random-access container
    // Optional Hint can be given to speed up point triangulation
    template <typename Position>
    auto operator()
    (Position const& position, Hint const& hint = {}) const
    {
      return this->operator()(make_point(position), hint);
    }
    
    // Triangulate a position given in internal Point type
    auto locate
    (Point const& point, Hint const& hint = {}) const
    {
      return triangulation.locate(point, hint);
    }
    
    // Triangulate a position given in a different container
    template <typename Position>
    auto locate
    (Position const& position, Hint const& hint = {}) const
    {
      return locate(make_point(position), hint);
    }
    
    // Convert a random-access container to internal Point type
    template <typename Position>
    Point make_point(Position const& position) const
    {
      return { position[0], position[1], position[2] };
    }
    
    // Convert a random-access container to internal Vector type
    template <typename Value>
    Vector make_vector(Value const& value) const
    {
      return { value[0], value[1], value[2] };
    }
    
    // Return how many times nearest-neighbor interpolation had to be used
    std::size_t nr_interpolation_failures()
    {
      return count_interpolation_failures;
    }
    
  private:
    Triangulation triangulation;  // Triangulation object for interpolation
    Coord_map function_values;    // Vector field values at points
    Tree tree;                    // Tree to find nearest neighbors if needed
    mutable std::size_t count_interpolation_failures{ 0 };
    
    template <typename Container_point, typename Container_vector>
    Coord_map make_coord_map
    (std::vector<Point> const& grid_points,
     std::vector<Vector> const& function_values_data)
    {
      Coord_map function_values;
      for (std::size_t idx = 0; idx < grid_points.size(); ++idx)
      function_values.insert(
        std::make_pair(grid_points[idx],
                       function_values_data[idx]));
      
      triangulation.insert(grid_points.begin(), grid_points.end());
      
      return function_values;
    }
    
    template <typename Container_point, typename Container_vector>
    Coord_map make_coord_map
    (Container_point const& grid_points_data,
    Container_vector const& function_values_data)
    {
      Coord_map function_values;
      std::vector<Point> grid_points;
      grid_points.reserve(grid_points_data.size());
      for (auto const& point : grid_points_data)
        grid_points.push_back(make_point(point));
      
      for (std::size_t idx = 0; idx < grid_points.size(); ++idx)
        function_values.insert(
          std::make_pair(grid_points[idx],
                         make_vector(function_values_data[idx])));
      
      triangulation.insert(grid_points.begin(), grid_points.end());
      
      return function_values;
    }
    
    template <typename Container>
    Coord_map make_coord_map
    (std::vector<Point> const& grid_points,
    Container const& function_values_data)
    {
      Coord_map function_values;
      for (std::size_t idx = 0; idx < grid_points.size(); ++idx)
      function_values.insert(
        std::make_pair(grid_points[idx],
                       make_vector(function_values_data[idx])));
      
      triangulation.insert(grid_points.begin(), grid_points.end());
      
      return function_values;
    }
  };
  
  // Linear vector interpolation on unstructured grid in 2d
  // Uses CGAL library for
  // Delaunay triangulation and natural neighbor coordinates
  // Type of the interpolated field values
  template <bool exact, bool check_interpolation, typename Value_type>
  class VectorField_LinearInterpolation_UnstructuredGrid
  <2, exact, check_interpolation, Value_type>
  {
  private:
    // Types for interpolation, see CGAL library documentation
    using Kernel = typename std::conditional<exact,
      CGAL::Exact_predicates_exact_constructions_kernel,
      CGAL::Exact_predicates_inexact_constructions_kernel>::type;
    using Triangulation = CGAL::Delaunay_triangulation_2<Kernel>;
    using Face_handle = typename Triangulation::Face_handle;
    using Coord_type = typename Kernel::FT;
    
  public:
    static constexpr std::size_t dim = 2;  // Spatial dimension
      
    using Hint = Face_handle;                  // Hint type to speed up locating a point
    using Point = typename Kernel::Point_2;    // CGAL spatial point type
    using Vector = typename Kernel::Vector_2;  // CGAL Field vector type
    using value_type = Value_type;             // Type of the interpolated field values
    
  private:
    // More types for interpolation, see CGAL library documentation
    using Coord_map = std::map<Point, Vector, typename Kernel::Less_xy_2>;
    using Value_access = CGAL::Data_access<Coord_map>;
    using Point_coordinate_vector =
      std::vector<std::pair<Point, Coord_type>>;
    
    // Tree to find nearest neighbors when check_interpolation is on
    using Tree = typename std::conditional<check_interpolation,
      KDTree<Coord_map, dim>, useful::Empty>::type;
    
  public:
    // Construct object given grid points and associated vector field values
    // from vectors of points and vector values in random-access containers
    template <typename Container_point, typename Container_vector>
    VectorField_LinearInterpolation_UnstructuredGrid
    (Container_point const& grid_points,
     Container_vector const& function_values)
    : function_values{ make_coord_map(grid_points, function_values) }
    , tree{ this->function_values }
    {}
    
    // Vector field interpolation value at position given in internal Point type
    // Optional Hint can be given to speed up point triangulation
    auto operator()
    (Point const& point, Hint const& hint = {}) const
    {
      Point_coordinate_vector coords_point;
      Coord_type norm = CGAL::natural_neighbor_coordinates_2(
        triangulation, point, std::back_inserter(coords_point), hint).second;
      
      if constexpr (check_interpolation)
      {
        if (norm == 0.)
        {
          count_interpolation_failures++;
          if constexpr(std::is_same<value_type, Vector>::value == true)
            return tree.nearest_neighbor_value(point);
          else
          {
            auto val = tree.nearest_neighbor_value(point);
            return value_type{
              CGAL::to_double(val[0]), CGAL::to_double(val[1]) };
          }
        }
      }
      
      if constexpr(std::is_same<value_type, Vector>::value == true)
        return CGAL::to_double(CGAL::linear_interpolation(
          coords_point.begin(), coords_point.end(), norm,
          Value_access(function_values)));
      else
      {
        auto val = CGAL::linear_interpolation(
          coords_point.begin(), coords_point.end(), norm,
          Value_access(function_values));
        return value_type{ CGAL::to_double(val[0]), CGAL::to_double(val[1]) };
      }
    }
    
    // Vector field interpolation value at position using a different container
    // Optional Hint can be given to speed up point triangulation
    template <typename Position>
    auto operator()
    (Position const& position, Hint const& hint = {}) const
    {
      return this->operator()(make_point(position), hint);
    }
    
    // Triangulate a position given in internal Point structure
    auto locate(Point const& point, Hint const& hint = {}) const
    {
      return triangulation.locate(point, hint);
    }
    
    // Triangulate a position given in a random-access container
    template <typename Position>
    auto locate(Position const& position, Hint const& hint = {}) const
    {
      return locate(make_point(position), hint);
    }
    
    // Convert a random-access container to internal Point type
    template <typename Position>
    Point make_point(Position const& position) const
    {
      return { position[0], position[1] };
    }
    
    // Convert a random-access container to internal Vector type
    template <typename Value>
    Vector make_vector(Value const& value) const
    {
      return { value[0], value[1] };
    }
    
    // Return how many times nearest-neighbor interpolation had to be used
    std::size_t nr_interpolation_failures()
    {
      return count_interpolation_failures;
    }
    
  private:
    Triangulation triangulation;  // Triangulation object for interpolation
    Coord_map function_values;    // Vector field values at points
    Tree tree;                    // Tree to find nearest neighbors if needed
    mutable std::size_t count_interpolation_failures{ 0 };
    
    template <typename Container_point, typename Container_vector>
    Coord_map make_coord_map
    (std::vector<Point> const& grid_points,
     std::vector<Vector> const& function_values_data)
    {
      Coord_map function_values;
      for (std::size_t idx = 0; idx < grid_points.size(); ++idx)
      function_values.insert(
        std::make_pair(grid_points[idx],
                       function_values_data[idx]));
      
      triangulation.insert(grid_points.begin(), grid_points.end());
      
      return function_values;
    }
    
    template <typename Container_point, typename Container_vector>
    Coord_map make_coord_map
    (Container_point const& grid_points_data,
    Container_vector const& function_values_data)
    {
      Coord_map function_values;
      std::vector<Point> grid_points;
      grid_points.reserve(grid_points_data.size());
      for (auto const& point : grid_points_data)
        grid_points.push_back(make_point(point));
      
      for (std::size_t idx = 0; idx < grid_points.size(); ++idx)
        function_values.insert(
          std::make_pair(grid_points[idx],
                         make_vector(function_values_data[idx])));
      
      triangulation.insert(grid_points.begin(), grid_points.end());
      
      return function_values;
    }
    
    template <typename Container>
    Coord_map make_coord_map
    (std::vector<Point> const& grid_points,
    Container const& function_values_data)
    {
      Coord_map function_values;
      for (std::size_t idx = 0; idx < grid_points.size(); ++idx)
      function_values.insert(
        std::make_pair(grid_points[idx],
                       make_vector(function_values_data[idx])));
      
      triangulation.insert(grid_points.begin(), grid_points.end());
      
      return function_values;
    }
  };
}


#endif /* VectorField_Interpolated_h */
