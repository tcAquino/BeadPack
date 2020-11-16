//
//  BeadPack_Input.h
//  BeadPack
//
//  Created by Tomás Aquino on 14/05/2020.
//  Copyright © 2020 Tomás Aquino. All rights reserved.
//

#ifndef BeadPack_Input_h
#define BeadPack_Input_h

#include <fstream>
#include <string>
#include <boost/algorithm/string.hpp>
#include "../general/useful.h"

namespace beadpack
{
  template <typename Bead>
  std::vector<Bead> get_beads
  (std::size_t dim, std::string const& filename,
   std::size_t header_lines = 0, double rescale = 1., std::size_t nr_estimate = 0, std::string const& delims = "\t,| ")
  {
    std::vector<Bead> beads;
    beads.reserve(nr_estimate);
    
    std::ifstream file(filename);
    if (!file.is_open())
      throw useful::open_read_error(filename);
    std::string line;
    for (std::size_t ll = 0; ll < header_lines; ++ll)
      getline(file, line);
    
    while (getline(file, line))
    {
      std::vector<std::string> split_line;
      boost::algorithm::split(split_line, line, boost::is_any_of(delims));
      
      beads.push_back({ typename Bead::Position_type(dim), rescale*std::stod(split_line[dim]) });
      for (std::size_t dd = 0; dd < dim; ++dd)
        beads.back().center[dd] = rescale*std::stod(split_line[dd]);
    }
    file.close();
    
    return beads;
  }
  
  template <typename Point>
  std::vector<Point> get_contacts
  (std::size_t dim, std::string const& filename,
   int header_lines = 0, double rescale = 1., std::size_t nr_estimate = 0, std::string const& delims = "\t,| ")
  {
    std::vector<Point> contacts;
    contacts.reserve(nr_estimate);
    
    std::ifstream file(filename);
    if (!file.is_open())
      throw useful::open_read_error(filename);
    std::string line;
    for (std::size_t ll = 0; ll < header_lines; ++ll)
      getline(file, line);
    
    while (getline(file, line))
    {
      std::vector<std::string> split_line;
      boost::algorithm::split(split_line, line, boost::is_any_of(delims));
      
      contacts.push_back(Point(dim));
      for (std::size_t dd = 0; dd < dim; ++dd)
        contacts.back()[dd] = rescale*std::stod(split_line[dd]);
    }
    file.close();
    
    return contacts;
  }
  
  template <typename Point, typename Vector>
  std::pair<std::vector<Point>, std::vector<Vector>> get_points_velocities_velocity_point
  (std::size_t dim, std::string const& filename,
   int header_lines = 0, double rescale_points = 1., double rescale_velocities = 1.,
   std::size_t nr_estimate = 0, std::string const& delims = "\t,| ")
  {
    std::pair<std::vector<Point>, std::vector<Vector>> points_velocities;
    points_velocities.first.reserve(nr_estimate);
    points_velocities.second.reserve(nr_estimate);
    
    std::ifstream file(filename);
    if (!file.is_open())
      throw useful::open_read_error(filename);
    std::string line;
    for (std::size_t ll = 0; ll < header_lines; ++ll)
      getline(file, line);
    
    while (getline(file, line))
    {
      std::vector<std::string> split_line;
      boost::algorithm::split(split_line, line, boost::is_any_of(delims));
      
      points_velocities.first.push_back(Point(dim));
      points_velocities.second.push_back(Vector(dim));
      for (std::size_t dd = 0; dd < dim; ++dd)
        points_velocities.first.back()[dd] = rescale_points*std::stod(split_line[dd+dim]);
      for (std::size_t dd = 0; dd < dim; ++dd)
        points_velocities.second.back()[dd] = rescale_velocities*std::stod(split_line[dd]);
    }
    file.close();
    
    return points_velocities;
  }
  
  template <typename Vector>
  Vector get_mean_velocity
  (std::size_t dim, std::string const& filename,
   int header_lines = 0, double rescale = 1., std::string const& delims = "\t,| ")
  {
    Vector mean_velocity;
    
    std::ifstream file(filename);
    if (!file.is_open())
      throw useful::open_read_error(filename);
    std::string line;
    for (std::size_t ll = 0; ll < header_lines; ++ll)
      getline(file, line);
    
    while (getline(file, line))
    {
      std::vector<std::string> split_line;
      boost::algorithm::split(split_line, line, boost::is_any_of(delims));
      
      for (std::size_t ii = 0; ii < split_line.size(); ++ii)
        mean_velocity.push_back(std::stod(split_line[ii]));
    }
    file.close();
    
    return mean_velocity;
  }
}


#endif /* BeadPack_Input_h */
