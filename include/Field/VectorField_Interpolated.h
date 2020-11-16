//
//  VectorField_Interpolated.h
//  BeadPack
//
//  Created by Tomás Aquino on 13/05/2020.
//  Copyright © 2020 Tomás Aquino. All rights reserved.
//

#ifndef VectorField_Interpolated_h
#define VectorField_Interpolated_h

#include <type_traits>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/natural_neighbor_coordinates_3.h>
#include <CGAL/interpolation_functions.h>

namespace field
{
  template <std::size_t dim, typename Value_type = std::vector<double>>
  class VectorField_LinearInterpolation_UnstructuredGrid;
  
  template <typename Value_type>
  class VectorField_LinearInterpolation_UnstructuredGrid<3, Value_type>
  {
  private:
    using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Triangulation = CGAL::Delaunay_triangulation_3<Kernel, CGAL::Fast_location>;
    using Cell_handle = Triangulation::Cell_handle;
    using Coord_type = Kernel::FT;
    using Vertex_handle = Triangulation::Vertex_handle;
    
  public:
    const std::size_t dim = 3;
    
    using Hint = Cell_handle;
    using Point = Kernel::Point_3;
    using Vector = Kernel::Vector_3;
    using value_type = Value_type;
    
  private:
    using Coord_map = std::map<Point, Vector, Kernel::Less_xyz_3>;
    using Value_access = CGAL::Data_access<Coord_map>;
    using Vertex_coordinate_vector = std::vector<std::pair<Vertex_handle, Coord_type>>;
    using Point_coordinate_vector = std::vector<std::pair<Point, Coord_type>>;
    
  public:
    template <typename Container>
    VectorField_LinearInterpolation_UnstructuredGrid
    (Container const& grid_points_data, Container const& function_values_data)
    {
      std::vector<Point> grid_points;
      grid_points.reserve(grid_points_data.size());
      for (auto const& point : grid_points_data)
        grid_points.push_back(make_point(point));
      
      for (std::size_t idx = 0; idx < grid_points.size(); ++idx)
        function_values.insert(std::make_pair(grid_points[idx],
                                              make_vector(function_values_data[idx])));
      
      triangulation.insert(grid_points.begin(), grid_points.end());
    }
    
    template <typename Container>
    VectorField_LinearInterpolation_UnstructuredGrid
    (std::vector<Point> const& grid_points, Container const& function_values_data)
    {
      for (std::size_t idx = 0; idx < grid_points.size(); ++idx)
        function_values.insert(std::make_pair(grid_points[idx],
                                              make_vector(function_values_data[idx])));
      
      triangulation.insert(grid_points.begin(), grid_points.end());
    }
    
    VectorField_LinearInterpolation_UnstructuredGrid
    (std::vector<Point> const& grid_points, std::vector<Vector> const& function_values_data)
    {
      for (std::size_t idx = 0; idx < grid_points.size(); ++idx)
        function_values.insert(std::make_pair(grid_points[idx], function_values_data[idx]));
      
      triangulation.insert(grid_points.begin(), grid_points.end());
    }
    
    auto operator()(Point const& point, Hint const& hint = {}) const
    {
      Vertex_coordinate_vector coords_vertex;
      Coord_type norm;
      CGAL::sibson_natural_neighbor_coordinates_3(triangulation,
                                                  point, std::back_inserter(coords_vertex),
                                                  norm, hint);
      
      Point_coordinate_vector coords_point;
      coords_point.reserve(coords_vertex.size());
      for (auto const& val : coords_vertex)
        coords_point.push_back(std::make_pair(val.first->point(), val.second));
      
      if constexpr(std::is_same<value_type, Vector>::value == true)
        return CGAL::linear_interpolation(coords_point.begin(), coords_point.end(), norm,
                                          Value_access(function_values));
      else
      {
        auto val = CGAL::linear_interpolation(coords_point.begin(), coords_point.end(), norm,
                                              Value_access(function_values));
        return value_type{ val[0], val[1], val[2] };
      }
    }
    
    template <typename Position>
    auto operator()(Position const& position, Hint const& hint = {}) const
    {
      return this->operator()(make_point(position), hint);
    }
    
    auto locate(Point const& point, Hint const& hint = {}) const
    {
      return triangulation.locate(point, hint);
    }
    
    template <typename Position>
    auto locate(Position const& position, Hint const& hint = {}) const
    {
      return locate(make_point(position), hint);
    }
    
    template <typename Position>
    Point make_point(Position const& position) const
    {
      return { position[0], position[1], position[2] };
    }
    
    template <typename Value>
    Vector make_vector(Value const& value) const
    {
      return { value[0], value[1], value[2] };
    }
    
  private:
    Coord_map function_values;
    Triangulation triangulation;
  };
  
  template <typename Value_type>
  class VectorField_LinearInterpolation_UnstructuredGrid<2, Value_type>
  {
  private:
    using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Triangulation = CGAL::Delaunay_triangulation_2<Kernel>;
    using Face_handle = Triangulation::Face_handle;
    using Coord_type = Kernel::FT;
    
  public:
    const std::size_t dim = 3;
    
    using Hint = Face_handle;
    using Point = Kernel::Point_2;
    using Vector = Kernel::Vector_2;
    using value_type = Value_type;
    
  private:
    using Coord_map = std::map<Point, Vector, Kernel::Less_xy_2>;
    using Value_access = CGAL::Data_access<Coord_map>;
    using Point_coordinate_vector = std::vector<std::pair<Point, Coord_type>>;
    
  public:
    template <typename Container>
    VectorField_LinearInterpolation_UnstructuredGrid
    (Container const& grid_points_data, Container const& function_values_data)
    {
      std::vector<Point> grid_points;
      grid_points.reserve(grid_points_data.size());
      for (auto const& point : grid_points_data)
        grid_points.push_back(make_point(point));
      
      for (std::size_t idx = 0; idx < grid_points.size(); ++idx)
        function_values.insert(std::make_pair(grid_points[idx],
                                              make_vector(function_values_data[idx])));
      
      triangulation.insert(grid_points.begin(), grid_points.end());
    }
    
    template <typename Container>
    VectorField_LinearInterpolation_UnstructuredGrid
    (std::vector<Point> const& grid_points, Container const& function_values_data)
    {
      for (std::size_t idx = 0; idx < grid_points.size(); ++idx)
        function_values.insert(std::make_pair(grid_points[idx],
                                              make_vector(function_values_data[idx])));
      
      triangulation.insert(grid_points.begin(), grid_points.end());
    }
    
    VectorField_LinearInterpolation_UnstructuredGrid
    (std::vector<Point> const& grid_points, std::vector<Vector> const& function_values_data)
    {
      for (std::size_t idx = 0; idx < grid_points.size(); ++idx)
        function_values.insert(std::make_pair(grid_points[idx], function_values_data[idx]));
      
      triangulation.insert(grid_points.begin(), grid_points.end());
    }
    
    auto operator()(Point const& point, Hint const& hint = {}) const
    {
      Point_coordinate_vector coords_point;
      Coord_type norm = CGAL::natural_neighbor_coordinates_2(triangulation,
                                                             point,
                                                             std::back_inserter(coords_point),
                                                             hint).second;
      
      if constexpr(std::is_same<value_type, Vector>::value == true)
        return CGAL::linear_interpolation(coords_point.begin(), coords_point.end(), norm,
                                          Value_access(function_values));
      else
      {
        auto val = CGAL::linear_interpolation(coords_point.begin(), coords_point.end(), norm,
                                              Value_access(function_values));
        return value_type{ val[0], val[1] };
      }
    }
    
    template <typename Position>
    auto operator()(Position const& position, Hint const& hint = {}) const
    {
      return this->operator()(make_point(position), hint);
    }
    
    auto locate(Point const& point, Hint const& hint = {}) const
    {
      return triangulation.locate(point, hint);
    }
    
    template <typename Position>
    auto locate(Position const& position, Hint const& hint = {}) const
    {
      return locate(make_point(position), hint);
    }
    
    template <typename Position>
    Point make_point(Position const& position) const
    {
      return { position[0], position[1] };
    }
    
    template <typename Value>
    Vector make_vector(Value const& value) const
    {
      return { value[0], value[1] };
    }
    
  private:
    Coord_map function_values;
    Triangulation triangulation;
  };
}


#endif /* VectorField_Interpolated_h */
