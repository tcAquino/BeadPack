//
//  BeadPack_Input.h
//  BeadPack
//
//  Created by Tomás Aquino on 14/05/2020.
//  Copyright © 2020 Tomás Aquino. All rights reserved.
//

#ifndef BeadPack_Input_h
#define BeadPack_Input_h

#include <boost/algorithm/string.hpp>
#include <fstream>
#include <string>
#include "general/useful.h"

namespace beadpack
{
  // Get bead positions and radii from file
  // File columns are bead center position components, radius
  // (any further columns ignored)
  // Parameters:
  //   dim : Spatial dimension (number of center components)
  //   filename : Name of file to be read
  //   header_lines : Number of lines to skip at top of file
  //   rescale : Multiply centers and radii by this factor
  //   nr_estimate : Estimate the number of beads for efficiency
  //   delims : possible delimiter characters separating data
  template <typename Bead>
  std::vector<Bead> get_beads_centers_radius
  (std::size_t dim, std::string const& filename,
   std::size_t header_lines = 0, double rescale = 1.,
   std::size_t nr_estimate = 0, std::string const& delims = "\t,| ")
  {
    std::vector<Bead> beads;
    beads.reserve(nr_estimate);
    
    auto file = useful::open_read(filename);
    std::string line;
    for (std::size_t ll = 0; ll < header_lines; ++ll)
      getline(file, line);
    
    while (getline(file, line))
    {
      std::vector<std::string> split_line;
      boost::trim_if(line, boost::is_any_of(delims+"\r"));
      boost::algorithm::split(split_line, line, boost::is_any_of(delims),
                              boost::token_compress_on);
      
      beads.push_back({ typename Bead::Position(dim),
        rescale*std::stod(split_line[dim]) });
      for (std::size_t dd = 0; dd < dim; ++dd)
        beads.back().center[dd] = rescale*std::stod(split_line[dd]);
    }
    file.close();
    
    return beads;
  }
  
  // Get bead positions and radii from file
  // File columns are bead center position components, diameter
  // (any further columns ignored)
  // Parameters:
  //   dim : Spatial dimension (number of center components)
  //   filename : Name of file to be read
  //   header_lines : Number of lines to skip at top of file
  //   rescale : Multiply centers and radii by this factor
  //   nr_estimate : Estimate the number of beads for efficiency
  //   delims : possible delimiter characters separating data
  template <typename Bead>
  std::vector<Bead> get_beads_centers_diameter
  (std::size_t dim, std::string const& filename,
   std::size_t header_lines = 0, double rescale = 1.,
   std::size_t nr_estimate = 0, std::string const& delims = "\t,| ")
  {
    std::vector<Bead> beads;
    beads.reserve(nr_estimate);
    
    auto file = useful::open_read(filename);
    std::string line;
    for (std::size_t ll = 0; ll < header_lines; ++ll)
      getline(file, line);
    
    while (getline(file, line))
    {
      std::vector<std::string> split_line;
      boost::trim_if(line, boost::is_any_of(delims+"\r"));
      boost::algorithm::split(split_line, line, boost::is_any_of(delims),
                              boost::token_compress_on);
      
      beads.push_back({ typename Bead::Position(dim),
        rescale*std::stod(split_line[dim])/2. });
      for (std::size_t dd = 0; dd < dim; ++dd)
        beads.back().center[dd] = rescale*std::stod(split_line[dd]);
    }
    file.close();
    
    return beads;
  }
  
  // Get positions of contact points betweem beads from file
  // File columns are contact point position components (any further columns ignored)
  // Parameters:
  //   dim : Spatial dimension (number of position components)
  //   filename : Name of file to be read
  //   header_lines : Number of lines to skip at top of file
  //   rescale : Multiply positions by this factor
  //   nr_estimate : Estimate the number of contacts for efficiency
  //   delims : possible delimiter characters separating data
  template <typename Point = std::vector<double>>
  std::vector<Point> get_contacts
  (std::size_t dim, std::string const& filename,
   std::size_t header_lines = 0, double rescale = 1.,
   std::size_t nr_estimate = 0, std::string const& delims = "\t,| ")
  {
    std::vector<Point> contacts;
    contacts.reserve(nr_estimate);
    
    auto file = useful::open_read(filename);
    std::string line;
    for (std::size_t ll = 0; ll < header_lines; ++ll)
      getline(file, line);
    
    while (getline(file, line))
    {
      std::vector<std::string> split_line;
      boost::trim_if(line, boost::is_any_of(delims+"\r"));
      boost::algorithm::split(split_line, line, boost::is_any_of(delims),
                              boost::token_compress_on);
      
      contacts.push_back(Point(dim));
      for (std::size_t dd = 0; dd < dim; ++dd)
        contacts.back()[dd] = rescale*std::stod(split_line[dd]);
    }
    file.close();
    
    return contacts;
  }
  
  // Get positions of contact points betweem beads from file
  // with cylinders centered at contacts
  // File columns are the positions components of the center of each cylinder base
  // (any further columns ignored)
  // Parameters:
  //   dim : Spatial dimension (number of position components)
  //   filename : Name of file to be read
  //   header_lines : Number of lines to skip at top of file
  //   rescale : Multiply positions by this factor
  //   nr_estimate : Estimate the number of contacts for efficiency
  //   delims : possible delimiter characters separating data
  template <typename Point = std::vector<double>>
  std::vector<Point> get_contacts_from_cylinders
  (std::size_t dim, std::string const& filename,
   std::size_t header_lines = 0, double rescale = 1.,
   std::size_t nr_estimate = 0, std::string const& delims = "\t,| ")
  {
    std::vector<Point> contacts;
    contacts.reserve(nr_estimate);
    
    auto file = useful::open_read(filename);
    std::string line;
    for (std::size_t ll = 0; ll < header_lines; ++ll)
      getline(file, line);
    
    while (getline(file, line))
    {
      std::vector<std::string> split_line;
      boost::trim_if(line, boost::is_any_of(delims+"\r"));
      boost::algorithm::split(split_line, line, boost::is_any_of(delims),
                              boost::token_compress_on);
      
      contacts.push_back(Point(dim));
      for (std::size_t dd = 0; dd < dim; ++dd)
        contacts.back()[dd] = rescale*
          (std::stod(split_line[dim+dd])+std::stod(split_line[dd]))/2.;
    }
    file.close();
    
    return contacts;
  }
  
  // Get velocity vectors (to interpolate) and corresponding positions from file
  // File columns are velocity components, point position components (any further columns ignored)
  // Parameters:
  //   dim : Spatial dimension (number of position components)
  //   filename : Name of file to be read
  //   header_lines : Number of lines to skip at top of file
  //   rescale_point : Multiply point positions by this factor
  //   rescale_velocities : Multiply velocity vectors by this factor
  //   nr_estimate : Estimate the number of points for efficiency
  template <typename Point = std::vector<double>,
  typename Vector = std::vector<double>>
  std::pair<std::vector<Point>, std::vector<Vector>>
  get_points_velocities_velocity_point
  (std::size_t dim, std::string const& filename,
   std::size_t header_lines = 0,
   double rescale_points = 1., double rescale_velocities = 1.,
   std::size_t nr_estimate = 0, std::string const& delims = "\t,| ")
  {
    std::pair<std::vector<Point>, std::vector<Vector>> points_velocities;
    points_velocities.first.reserve(nr_estimate);
    points_velocities.second.reserve(nr_estimate);
    
    auto file = useful::open_read(filename);
    std::string line;
    for (std::size_t ll = 0; ll < header_lines; ++ll)
      getline(file, line);
    
    while (getline(file, line))
    {
      std::vector<std::string> split_line;
      boost::trim_if(line, boost::is_any_of(delims+"\r"));
      boost::algorithm::split(split_line, line, boost::is_any_of(delims),
                              boost::token_compress_on);
      
      points_velocities.first.push_back(Point(dim));
      points_velocities.second.push_back(Vector(dim));
      for (std::size_t dd = 0; dd < dim; ++dd)
        points_velocities.first.back()[dd] =
          rescale_points*std::stod(split_line[dim+dd]);
      for (std::size_t dd = 0; dd < dim; ++dd)
        points_velocities.second.back()[dd] =
          rescale_velocities*std::stod(split_line[dd]);
    }
    file.close();
    
    return points_velocities;
  }
  
  // Get velocity vectors (to interpolate) and corresponding positions from file
  // File columns are point position components, velocity components (any further columns ignored)
  // Parameters:
  //   dim : Spatial dimension (number of position components)
  //   filename : Name of file to be read
  //   header_lines : Number of lines to skip at top of file
  //   rescale_point : Multiply point positions by this factor
  //   rescale_velocities : Multiply velocity vectors by this factor
  //   nr_estimate : Estimate the number of points for efficiency
  template <typename Point = std::vector<double>,
  typename Vector = std::vector<double>>
  std::pair<std::vector<Point>, std::vector<Vector>>
  get_points_velocities_point_velocity
  (std::size_t dim, std::string const& filename,
   std::size_t header_lines = 0,
   double rescale_points = 1., double rescale_velocities = 1.,
   std::size_t nr_estimate = 0, std::string const& delims = "\t,| ")
  {
    std::pair<std::vector<Point>, std::vector<Vector>> points_velocities;
    points_velocities.first.reserve(nr_estimate);
    points_velocities.second.reserve(nr_estimate);
    
    auto file = useful::open_read(filename);
    std::string line;
    for (std::size_t ll = 0; ll < header_lines; ++ll)
      getline(file, line);
    
    while (getline(file, line))
    {
      std::vector<std::string> split_line;
      boost::trim_if(line, boost::is_any_of(delims+"\r"));
      boost::algorithm::split(split_line, line, boost::is_any_of(delims),
                              boost::token_compress_on);
      
      points_velocities.first.push_back(Point(dim));
      points_velocities.second.push_back(Vector(dim));
      for (std::size_t dd = 0; dd < dim; ++dd)
        points_velocities.first.back()[dd] =
          rescale_points*std::stod(split_line[dd]);
      for (std::size_t dd = 0; dd < dim; ++dd)
        points_velocities.second.back()[dd] =
          rescale_velocities*std::stod(split_line[dim+dd]);
    }
    file.close();
    
    return points_velocities;
  }
  
  // Get velocity vectors (to interpolate) and corresponding positions from file
  // File columns are:
  // $v_1$ $v_2$ ... $\del_{x_1} v_1$ $\del_{x_2} v_1$ ... $\del_{x_1} v_2$ $\del_{x_2} v_2$ ... $x_1$ $x_2$ ...
  // (any further columns ignored)
  // Parameters:
  //   dim : Spatial dimension (number of position components)
  //   filename : Name of file to be read
  //   header_lines : Number of lines to skip at top of file
  //   rescale_point : Multiply point positions by this factor
  //   rescale_velocities : Multiply velocity vectors by this factor
  //   nr_estimate : Estimate the number of points for efficiency
  template
  <typename Point = std::vector<double>,
  typename Vector = std::vector<double>>
  std::pair<std::vector<Point>, std::vector<Vector>>
  get_points_velocities_velocity_gradient_point
  (std::size_t dim, std::string const& filename,
   std::size_t header_lines = 0,
   double rescale_points = 1., double rescale_velocities = 1.,
   std::size_t nr_estimate = 0, std::string const& delims = "\t,| ")
  {
    std::pair<std::vector<Point>, std::vector<Vector>> points_velocities;
    points_velocities.first.reserve(nr_estimate);
    points_velocities.second.reserve(nr_estimate);
    
    auto file = useful::open_read(filename);
    std::string line;
    for (std::size_t ll = 0; ll < header_lines; ++ll)
      getline(file, line);
    
    while (getline(file, line))
    {
      std::vector<std::string> split_line;
      boost::trim_if(line, boost::is_any_of(delims+"\r"));
      boost::algorithm::split(split_line, line, boost::is_any_of(delims),
                              boost::token_compress_on);
      
      points_velocities.first.push_back(Point(dim));
      points_velocities.second.push_back(Vector(dim));
      for (std::size_t dd = 0; dd < dim; ++dd)
        points_velocities.first.back()[dd] =
          rescale_points*std::stod(split_line[dim*(dim+1)+dd]);
      for (std::size_t dd = 0; dd < dim; ++dd)
        points_velocities.second.back()[dd] = rescale_velocities*std::stod(split_line[dd]);
    }
    file.close();
    
    return points_velocities;
  }
  
  // Get velocity vectors (to interpolate) and corresponding positions from file
  // File columns are:
  // $x_1$ $x_2$ ...$v_1$ $v_2$ ... $\del_{x_1} v_1$ $\del_{x_2} v_1$ ... $\del_{x_1} v_2$ $\del_{x_2} v_2$ ...
  // (any further columns ignored)
  // Parameters:
  //   dim : Spatial dimension (number of position components)
  //   filename : Name of file to be read
  //   header_lines : Number of lines to skip at top of file
  //   rescale_point : Multiply point positions by this factor
  //   rescale_velocities : Multiply velocity vectors by this factor
  //   nr_estimate : Estimate the number of points for efficiency
  template <typename Point = std::vector<double>, typename Vector = std::vector<double>>
  std::pair<std::vector<Point>, std::vector<Vector>>
  get_points_velocities_point_velocity_gradient
  (std::size_t dim, std::string const& filename,
   std::size_t header_lines = 0,
   double rescale_points = 1., double rescale_velocities = 1.,
   std::size_t nr_estimate = 0, std::string const& delims = "\t,| ")
  {
    std::pair<std::vector<Point>, std::vector<Vector>> points_velocities;
    points_velocities.first.reserve(nr_estimate);
    points_velocities.second.reserve(nr_estimate);
    
    auto file = useful::open_read(filename);
    std::string line;
    for (std::size_t ll = 0; ll < header_lines; ++ll)
      getline(file, line);
    
    while (getline(file, line))
    {
      std::vector<std::string> split_line;
      boost::trim_if(line, boost::is_any_of(delims+"\r"));
      boost::algorithm::split(split_line, line, boost::is_any_of(delims),
                              boost::token_compress_on);
      
      points_velocities.first.push_back(Point(dim));
      points_velocities.second.push_back(Vector(dim));
      for (std::size_t dd = 0; dd < dim; ++dd)
        points_velocities.first.back()[dd] =
          rescale_points*std::stod(split_line[dd]);
      for (std::size_t dd = 0; dd < dim; ++dd)
        points_velocities.second.back()[dd] =
          rescale_velocities*std::stod(split_line[dim+dd]);
    }
    file.close();
    
    return points_velocities;
  }
  
  // Get velocity vectors (to interpolate) and corresponding positions from file
  // File columns are:
  // $x_1$ $x_2$ ... $v_1$ $\del_{x_1} v_1$ $\del_{x_2} v_1$ ... $v_2$ ...
  // (any further columns ignored)
  // Parameters:
  //   dim : Spatial dimension (number of position components)
  //   filename : Name of file to be read
  //   header_lines : Number of lines to skip at top of file
  //   rescale_point : Multiply point positions by this factor
  //   rescale_velocities : Multiply velocity vectors by this factor
  //   nr_estimate : Estimate the number of points for efficiency
  template <typename Point = std::vector<double>,
  typename Vector = std::vector<double>>
  std::pair<std::vector<Point>, std::vector<Vector>>
  get_points_velocities_point_velocity_gradient_2
  (std::size_t dim, std::string const& filename,
   std::size_t header_lines = 0,
   double rescale_points = 1., double rescale_velocities = 1.,
   std::size_t nr_estimate = 0, std::string const& delims = "\t,| ")
  {
    std::pair<std::vector<Point>, std::vector<Vector>> points_velocities;
    points_velocities.first.reserve(nr_estimate);
    points_velocities.second.reserve(nr_estimate);
    
    auto file = useful::open_read(filename);
    std::string line;
    for (std::size_t ll = 0; ll < header_lines; ++ll)
      getline(file, line);
    
    while (getline(file, line))
    {
      std::vector<std::string> split_line;
      boost::trim_if(line, boost::is_any_of(delims+"\r"));
      boost::algorithm::split(split_line, line, boost::is_any_of(delims),
                              boost::token_compress_on);
      
      points_velocities.first.push_back(Point(dim));
      points_velocities.second.push_back(Vector(dim));
      for (std::size_t dd = 0; dd < dim; ++dd)
        points_velocities.first.back()[dd] =
          rescale_points*std::stod(split_line[dd]);
      for (std::size_t dd = 0; dd < dim; ++dd)
        points_velocities.second.back()[dd] =
          rescale_velocities*std::stod(split_line[dim+dd*(dim+1)]);
    }
    file.close();
    
    return points_velocities;
  }
  
  // Get velocity vectors (to interpolate) and corresponding positions from file
  // File columns are:
  // $v_1$ $\del_{x_1} v_1$ $\del_{x_2} v_1$ ... $v_2$ ... $x_1$ $x_2$ ...
  // (any further columns ignored)
  // Parameters:
  //   dim : Spatial dimension (number of position components)
  //   filename : Name of file to be read
  //   header_lines : Number of lines to skip at top of file
  //   rescale_point : Multiply point positions by this factor
  //   rescale_velocities : Multiply velocity vectors by this factor
  //   nr_estimate : Estimate the number of points for efficiency
  template
  <typename Point = std::vector<double>,
  typename Vector = std::vector<double>>
  std::pair<std::vector<Point>, std::vector<Vector>>
  get_points_velocities_velocity_gradient_point_2
  (std::size_t dim, std::string const& filename,
   std::size_t header_lines = 0,
   double rescale_points = 1., double rescale_velocities = 1.,
   std::size_t nr_estimate = 0, std::string const& delims = "\t,| ")
  {
    std::pair<std::vector<Point>, std::vector<Vector>> points_velocities;
    points_velocities.first.reserve(nr_estimate);
    points_velocities.second.reserve(nr_estimate);
    
    auto file = useful::open_read(filename);
    std::string line;
    for (std::size_t ll = 0; ll < header_lines; ++ll)
      getline(file, line);
    
    while (getline(file, line))
    {
      std::vector<std::string> split_line;
      boost::trim_if(line, boost::is_any_of(delims+"\r"));
      boost::algorithm::split(split_line, line, boost::is_any_of(delims),
                              boost::token_compress_on);
      
      points_velocities.first.push_back(Point(dim));
      points_velocities.second.push_back(Vector(dim));
      for (std::size_t dd = 0; dd < dim; ++dd)
        points_velocities.first.back()[dd] =
          rescale_points*std::stod(split_line[dim*(dim+1)+dd]);
      for (std::size_t dd = 0; dd < dim; ++dd)
        points_velocities.second.back()[dd] = rescale_velocities*std::stod(split_line[dd*(dim+1)]);
    }
    file.close();
    
    return points_velocities;
  }
  
  // Get mean velocity vector from file
  // First dim values in file (row major) are  mean velocity components
  // (any further values ignored)
  // Parameters:
  //   dim : Spatial dimension (number of position components)
  //   filename : Name of file to be read
  //   header_lines : Number of lines to skip at top of file
  //   nr_estimate : Estimate the number of points for efficiency
  template <typename Vector = std::vector<double>>
  Vector get_mean_velocity
  (std::size_t dim, std::string const& filename,
   std::size_t header_lines = 0, double rescale = 1.,
   std::string const& delims = "\t,| ")
  {
    Vector mean_velocity;
    
    auto file = useful::open_read(filename);
    std::string line;
    for (std::size_t ll = 0; ll < header_lines; ++ll)
      getline(file, line);
    
    while (getline(file, line))
    {
      std::vector<std::string> split_line;
      boost::trim_if(line, boost::is_any_of(delims+"\r"));
      boost::algorithm::split(split_line, line, boost::is_any_of(delims),
                              boost::token_compress_on);
      
      for (std::size_t ii = 0; ii < split_line.size(); ++ii)
        mean_velocity.push_back(std::stod(split_line[ii]));
    }
    file.close();
    
    return mean_velocity;
  }
}


#endif /* BeadPack_Input_h */
