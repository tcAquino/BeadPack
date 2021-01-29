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

  namespace model_bcc_cartesian
  {
    using Geometry = Geometry_Bcc<0>;
    using BeadPack = beadpack::BeadPack<Geometry::dim>;
    using Bead = BeadPack::Bead;
    using VelocityField =
      field::VectorField_LinearInterpolation_UnstructuredGrid<Geometry::dim>;
    
    BeadPack make_bead_pack
    (std::string const& input_dir)
    {
      std::string bead_filename = input_dir + "/" + "beads.dat";
      return beadpack::get_beads_centers_radius<Bead>(
        Geometry::dim, bead_filename, 1);
    }
    
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
    
    VelocityField make_velocity_field
    (std::string const& input_dir, std::string const& output_dir,
     BeadPack const& bead_pack)
    {
      std::string contact_filename = input_dir + "/" + "contacts.dat";
      auto contacts = beadpack::get_contacts(
        Geometry::dim, contact_filename, 1);
      
      std::string velocity_filename = input_dir + "/" + "velocities.dat";
      auto points_velocities = beadpack::get_points_velocities_velocity_point(
        Geometry::dim, velocity_filename, 1);

      for (auto const& point : contacts)
      {
        points_velocities.first.push_back(point);
        points_velocities.second.emplace_back(Geometry::dim, 0.);
      }
      for (auto const& bead : bead_pack.beads())
      {
        points_velocities.first.push_back(bead.center);
        points_velocities.second.emplace_back(Geometry::dim, 0.);
      }
      
      return { points_velocities.first, points_velocities.second };
    }
    
    VelocityField make_velocity_field
    (std::string const& input_dir, std::string const& output_dir,
     BeadPack const& bead_pack,
     Boundaries::Boundary_Periodic const& boundary_periodic)
    {
      return make_velocity_field(input_dir, output_dir, bead_pack);
    }
  }

  namespace model_bcc_symmetryplanes
  {
    using Geometry = Geometry_Bcc<1>;
    
    using BeadPack = beadpack::BeadPack<Geometry::dim>;
    using Bead = BeadPack::Bead;
    using VelocityField =
      field::VectorField_LinearInterpolation_UnstructuredGrid<Geometry::dim, 0>;
    
    BeadPack make_bead_pack
    (std::string const& input_dir)
    {
      std::string bead_filename = input_dir + "/" + "beads.dat";
      return beadpack::get_beads_centers_radius<Bead>(
        Geometry::dim, bead_filename, 1);
    }
    
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
    
    VelocityField make_velocity_field
    (std::string const& input_dir, std::string const& output_dir,
     BeadPack const& bead_pack,
     Boundaries::Boundary_Periodic const& boundary_periodic)
    {
      std::string contact_filename = input_dir + "/" + "cylinders.dat";
      auto contacts = beadpack::get_contacts_from_cylinders(
        Geometry::dim, contact_filename, 1);
      
      std::string velocity_filename = input_dir + "/" + "velocities.dat";
      auto points_velocities = beadpack::get_points_velocities_point_velocity_gradient(
        Geometry::dim, velocity_filename, 6);
      
      for (auto const& point : contacts)
      {
        points_velocities.first.push_back(point);
        points_velocities.second.emplace_back(Geometry::dim, 0.);
      }
      for (auto const& bead : bead_pack.beads())
      {
        points_velocities.first.push_back(bead.center);
        points_velocities.second.emplace_back(Geometry::dim, 0.);
      }
      
      double fraction = 1e-3;
      beadpack::add_periodic_images(points_velocities.first,
                                    points_velocities.second,
                                    fraction, boundary_periodic);
      
      return { points_velocities.first, points_velocities.second };
    }
  }
}


#endif /* BeadPack_Models_h */
