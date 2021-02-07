//
//  VectorField_Interpolated.h
//  Field
//
//  Created by Tomás Aquino on 13/05/2020.
//  Copyright © 2020 Tomás Aquino. All rights reserved.
//

#ifndef VectorField_Interpolated_h
#define VectorField_Interpolated_h

#include <type_traits>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/interpolation_functions.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/natural_neighbor_coordinates_3.h>

namespace field
{
  // Template class for linear vector interpolation on unstructured grid
  // to be specialized for different dimensions
  // If the exact parameter is 1, the interpolator uses
  // exact (arbitrary-precision) constructions. This is slower, but necessary
  // for some meshes. If it is 0 (default), inexact (double-precision) constructions are
  // used, which is faster but may result in NaN interpolation values
  // for some meshes
  template
  <std::size_t dim, bool exact = 0, typename Value_type = std::vector<double>>
  class VectorField_LinearInterpolation_UnstructuredGrid;
  
  // Linear vector interpolation on unstructured grid in 3d
  // Uses CGAL library for
  // Delaunay triangulation and Sibson natural neighbor coordinates
  // Value_type is the type of the interpolated field values
  template <bool exact, typename Value_type>
  class VectorField_LinearInterpolation_UnstructuredGrid<3, exact, Value_type>
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
    
  public:
    // Construct object given grid points and associated vector field values
    // from vectors of internal Point and Vector types
    VectorField_LinearInterpolation_UnstructuredGrid
    (std::vector<Point> const& grid_points,
     std::vector<Vector> const& function_values_data)
    {
      for (std::size_t idx = 0; idx < grid_points.size(); ++idx)
        function_values.insert(
          std::make_pair(grid_points[idx],
                         function_values_data[idx]));
      
      triangulation.insert(grid_points.begin(), grid_points.end());
    }
    
    // Construct object given grid points and associated vector field values
    // from vectors of internal Point type
    // and vector values in a random-access container
    template <typename Container>
    VectorField_LinearInterpolation_UnstructuredGrid
    (std::vector<Point> const& grid_points,
     Container const& function_values_data)
    {
      for (std::size_t idx = 0; idx < grid_points.size(); ++idx)
        function_values.insert(
          std::make_pair(grid_points[idx],
                         make_vector(function_values_data[idx])));
      
      triangulation.insert(grid_points.begin(), grid_points.end());
    }
    
    // Construct object given grid points and associated vector field values
    // from vectors of points and vector values in a random-access container
    template <typename Container>
    VectorField_LinearInterpolation_UnstructuredGrid
    (Container const& grid_points_data,
     Container const& function_values_data)
    {
      std::vector<Point> grid_points;
      grid_points.reserve(grid_points_data.size());
      for (auto const& point : grid_points_data)
        grid_points.push_back(make_point(point));
      
      for (std::size_t idx = 0; idx < grid_points.size(); ++idx)
        function_values.insert(
          std::make_pair(grid_points[idx],
                         make_vector(function_values_data[idx])));
      
      triangulation.insert(grid_points.begin(), grid_points.end());
    }
    
    // Vector field interpolation value at position given in internal Point type
    // Optional Hint can be given to speed up point triangulation
    auto operator()(Point const& point, Hint const& hint = {}) const
    {
      Vertex_coordinate_vector coords_vertex;
      Coord_type norm;
      CGAL::sibson_natural_neighbor_coordinates_3(
        triangulation, point, std::back_inserter(coords_vertex), norm, hint);
      
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
    
  private:
    Coord_map function_values;    // Vector field values at points
    Triangulation triangulation;  // Triangulation object for interpolation
  };
  
  // Linear vector interpolation on unstructured grid in 2d
  // Uses CGAL library for
  // Delaunay triangulation and natural neighbor coordinates
  // Type of the interpolated field values
  template <bool exact, typename Value_type>
  class VectorField_LinearInterpolation_UnstructuredGrid<2, exact, Value_type>
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
    
  public:
    // Construct object given grid points and associated vector field values
    // from vectors of internal Point and Vector types
    VectorField_LinearInterpolation_UnstructuredGrid
    (std::vector<Point> const& grid_points,
     std::vector<Vector> const& function_values_data)
    {
      for (std::size_t idx = 0; idx < grid_points.size(); ++idx)
        function_values.insert(
          std::make_pair(grid_points[idx], function_values_data[idx]));
      
      triangulation.insert(grid_points.begin(), grid_points.end());
    }
    
    // Construct object given grid points and associated vector field values
    // from vectors of internal Point type
    // and vector values in a random-access container
    template <typename Container>
    VectorField_LinearInterpolation_UnstructuredGrid
    (std::vector<Point> const& grid_points,
     Container const& function_values_data)
    {
      for (std::size_t idx = 0; idx < grid_points.size(); ++idx)
        function_values.insert(
          std::make_pair(grid_points[idx],
                         make_vector(function_values_data[idx])));
      
      triangulation.insert(grid_points.begin(), grid_points.end());
    }
    
    // Construct object given grid points and associated vector field values
    // from vectors of points and vector values in a random-access container
    template <typename Container>
    VectorField_LinearInterpolation_UnstructuredGrid
    (Container const& grid_points_data, Container const& function_values_data)
    {
      std::vector<Point> grid_points;
      grid_points.reserve(grid_points_data.size());
      for (auto const& point : grid_points_data)
        grid_points.push_back(make_point(point));
      
      for (std::size_t idx = 0; idx < grid_points.size(); ++idx)
        function_values.insert(
          std::make_pair(grid_points[idx],
                         make_vector(function_values_data[idx])));
      
      triangulation.insert(grid_points.begin(), grid_points.end());
    }
    
    // Vector field interpolation value at position given in internal Point type
    // Optional Hint can be given to speed up point triangulation
    auto operator()
    (Point const& point, Hint const& hint = {}) const
    {
      Point_coordinate_vector coords_point;
      Coord_type norm = CGAL::natural_neighbor_coordinates_2(
        triangulation, point, std::back_inserter(coords_point), hint).second;
      
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
    
  private:
    Coord_map function_values;    // Vector field values at points
    Triangulation triangulation;  // Triangulation object for interpolation
  };
}


#endif /* VectorField_Interpolated_h */
