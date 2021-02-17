//
//  BeadPack_Models.h
//
//  Created by Tomás Aquino on 20/01/2021.
//  Copyright © 2021 Tomás Aquino. All rights reserved.
//

#ifndef BeadPack_Models_h
#define BeadPack_Models_h

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "BeadPack/BeadPack.h"
#include "BeadPack/BeadPack_InitialConditions.h"
#include "BeadPack/BeadPack_Input.h"
#include "Field/VectorField_Interpolated.h"
#include "general/Operations.h"
#include "general/Ranges.h"
#include "general/useful.h"
#include "Geometry/SymmetryPlanes.h"
#include "Stochastic/CTRW/Boundary.h"

namespace beadpack
{
  struct Tag{
    struct Cylinder{};
    struct Diameter{};
    struct None{};
    struct Point{};
    struct Radius{};
    struct Sphere{};
  };
  
  template <bool centered>
  struct Geometry_Bcc
  {
    static constexpr std::size_t dim{ 3 };
    double radius;
    double domain_side{ 4./std::sqrt(3.)*radius };
    std::vector<double> domain_dimensions;
    std::vector<std::pair<double, double>> boundaries;
    
    Geometry_Bcc(double radius = 1.)
    : radius{ radius }
    , domain_dimensions(dim, domain_side)
    , boundaries{ centered
      ? std::vector<std::pair<double, double>>(dim,
        { -domain_side/2., domain_side/2. })
      : std::vector<std::pair<double, double>>(dim,
        { 0., domain_side }) }
    {}
  };
  
  template <typename BeadPack, typename BuildFrom>
  BeadPack make_bead_pack
  (std::size_t dim, std::string const& filename,
   std::size_t header_lines = 0, double rescale = 1.,
   std::size_t nr_estimate = 0, std::string const& delims = "\t,| ")
  {
    if constexpr (std::is_same_v<BuildFrom, Tag::Radius>)
      return beadpack::get_beads_centers_radius<typename BeadPack::Bead>(
        dim, filename, header_lines, rescale, nr_estimate, delims);
    if constexpr (std::is_same_v<BuildFrom, Tag::Diameter>)
      return beadpack::get_beads_centers_diameter<typename BeadPack::Bead>(
        dim, filename, header_lines, rescale, nr_estimate, delims);
    throw useful::bad_parameters();
  }
    
  template
  <typename VelocityField,
  bool add_periodic_images = 1,
  bool velocity_file_points_first = 1,
  bool velocity_file_has_gradients = 0,
  bool velocity_file_all_velocities_before_gradients = 1,
  bool add_zero_centers = 1,
  typename Contacts = Tag::Point,
  typename BeadPack,
  typename Boundary_Periodic = boundary::DoNothing>
  VelocityField make_velocity_field
  (BeadPack const& bead_pack,
   std::string const& velocity_filename,
   std::string const& contact_filename = {},
   Boundary_Periodic const& boundary_periodic = {},
   double fraction_domain_image_periodic = 1e-3,
   std::size_t header_lines_velocity = 0,
   std::size_t header_lines_contact = 0,
   std::size_t nr_estimate_velocity = 0,
   std::size_t nr_estimate_contact = 0,
   double rescale_point = 1.,
   double rescale_velocity = 1.,
   double rescale_contact = 1.,
   std::string const& delims = "\t,| ")
  {
    auto points_velocities =
      velocity_file_has_gradients
      ? velocity_file_points_first
        ? velocity_file_all_velocities_before_gradients
          ? get_points_velocities_point_velocity_gradient(
              bead_pack.dim,
              velocity_filename, header_lines_velocity,
              rescale_point, rescale_velocity,
              nr_estimate_velocity, delims)
          : get_points_velocities_point_velocity_gradient_2(
              bead_pack.dim,
              velocity_filename, header_lines_velocity,
              rescale_point, rescale_velocity,
              nr_estimate_velocity, delims)
        : velocity_file_all_velocities_before_gradients
          ? get_points_velocities_velocity_gradient_point(
              bead_pack.dim,
              velocity_filename, header_lines_velocity,
              rescale_point, rescale_velocity,
              nr_estimate_velocity, delims)
          : get_points_velocities_velocity_gradient_point_2(
              bead_pack.dim,
              velocity_filename, header_lines_velocity,
              rescale_point, rescale_velocity,
              nr_estimate_velocity, delims)
      : velocity_file_points_first
        ? get_points_velocities_point_velocity(
            bead_pack.dim,
            velocity_filename, header_lines_velocity,
            rescale_point, rescale_velocity,
            nr_estimate_velocity, delims)
        : get_points_velocities_velocity_point(
            bead_pack.dim,
            velocity_filename, header_lines_velocity,
            rescale_point, rescale_velocity,
            nr_estimate_velocity, delims);
   
    if constexpr (std::is_same_v<Contacts, Tag::None>)
    {}
    else if constexpr (std::is_same_v<Contacts, Tag::Cylinder>)
    {
      auto contacts = beadpack::get_contacts_from_cylinders(
        bead_pack.dim, contact_filename, header_lines_contact,
        rescale_contact, nr_estimate_contact, delims);
      for (auto const& point : contacts)
      {
        points_velocities.first.push_back(point);
        points_velocities.second.emplace_back(bead_pack.dim, 0.);
      }
    }
    else
    {
      auto contacts = beadpack::get_contacts(
        bead_pack.dim, contact_filename, header_lines_contact,
        rescale_contact, nr_estimate_contact, delims);
      for (auto const& point : contacts)
      {
        points_velocities.first.push_back(point);
        points_velocities.second.emplace_back(bead_pack.dim, 0.);
      }
    }
    
    if constexpr (add_zero_centers)
      for (auto const& bead : bead_pack.beads())
      {
        points_velocities.first.push_back(bead.center);
        points_velocities.second.emplace_back(bead_pack.dim, 0.);
      }
    
    if constexpr (add_periodic_images)
      geometry::add_periodic_images(points_velocities.first,
                                    points_velocities.second,
                                    fraction_domain_image_periodic,
                                    boundary_periodic);
    
    return { points_velocities.first, points_velocities.second };
  }

  namespace model_bcc_cartesian
  {
    using Geometry = Geometry_Bcc<0>;
    using BeadPack = beadpack::BeadPack<Geometry::dim>;
    using VelocityField =
      field::VectorField_LinearInterpolation_UnstructuredGrid<Geometry::dim>;
    
    struct Boundaries
    {
      using Boundary_Periodic = boundary::Periodic_WithOutsideInfo;
      using Boundary_Reflecting_Periodic
        = boundary::ReflectingBeads_Periodic<BeadPack, Boundary_Periodic>;
      
      Boundaries(Geometry const& geometry, BeadPack const& bead_pack)
      : boundary_periodic{ geometry.boundaries }
      , boundary_reflecting_periodic{ bead_pack, boundary_periodic }
      {}
      
      Boundary_Periodic boundary_periodic;
      Boundary_Reflecting_Periodic boundary_reflecting_periodic;
    };
    
    auto make_bead_pack
    (std::string const& input_dir)
    {
      return beadpack::make_bead_pack<BeadPack, Tag::Radius>(
        Geometry::dim, input_dir + "/" + "beads.dat", 1);
    }
    
    auto make_velocity_field
    (std::string const& input_dir, BeadPack const& bead_pack)
    {
      return beadpack::make_velocity_field
        <VelocityField, 0, 0, 0, 0, 1, Tag::Point>(
        bead_pack,
        input_dir + "/" + "velocities.dat",
        input_dir + "/" + "contacts.dat",
        {}, 0, 1, 1);
    }
    
    auto make_velocity_field
    (std::string const& input_dir, BeadPack const& bead_pack,
     Boundaries::Boundary_Periodic const& boundary_periodic)
    {
      return make_velocity_field(input_dir, bead_pack);
    }
  }

  namespace model_bcc_symmetryplanes
  {
    using Geometry = Geometry_Bcc<1>;
    using BeadPack = beadpack::BeadPack<Geometry::dim>;
    using VelocityField =
      field::VectorField_LinearInterpolation_UnstructuredGrid<Geometry::dim>;
    
    struct Boundaries
    {
      using Boundary_Periodic
        = boundary::Periodic_SymmetryPlanes_WithOutsideInfo<geometry::SymmetryPlanes_Bcc>;
      using Boundary_Reflecting_Periodic
        = boundary::ReflectingBeads_Periodic<BeadPack, Boundary_Periodic>;
      
      Boundaries(Geometry const& geometry, BeadPack const& bead_pack)
      : boundary_periodic{ geometry.radius }
      , boundary_reflecting_periodic{ bead_pack, boundary_periodic }
      {}
      
      Boundary_Periodic boundary_periodic;
      Boundary_Reflecting_Periodic boundary_reflecting_periodic;
    };
    
    auto make_bead_pack
    (std::string const& input_dir)
    {
      return beadpack::make_bead_pack<BeadPack, Tag::Radius>(
        Geometry::dim, input_dir + "/" + "beads.dat", 1);
    }
    
    auto make_velocity_field
    (std::string const& input_dir, BeadPack const& bead_pack,
     Boundaries::Boundary_Periodic const& boundary_periodic)
    {
      return beadpack::make_velocity_field
        <VelocityField, 1, 1, 1, 0, 1, Tag::Cylinder>(
        bead_pack,
        input_dir + "/" + "velocities.dat",
        input_dir + "/" + "contacts.dat",
        boundary_periodic, 1e-3, 6, 1);
    }
  }
}


#endif /* BeadPack_Models_h */
