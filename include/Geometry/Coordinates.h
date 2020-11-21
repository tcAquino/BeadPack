//
//  Coordinates.h
//  BeadPack_Reactive
//
//  Created by Tomás Aquino on 17/11/2020.
//  Copyright © 2020 Tomás Aquino. All rights reserved.
//

#ifndef Coordinates_h
#define Coordinates_h

#include <cmath>
#include <vector>
#include "general/Operations.h"

namespace geometry
{
  template <typename Container = std::vector<double>>
  Container cartesian2spherical(Container const& cartesian)
  {
    Container spherical(3);
    
    spherical[0] = operation::abs(cartesian);
    spherical[1] = std::atan2(cartesian[1], cartesian[0]);
    spherical[2] = std::acos(cartesian[2]/spherical[0]);
    
    return spherical;
  }

  template <typename Container = std::vector<double>>
  Container spherical2cartesian(std::vector<double> const& spherical)
  {
    Container cartesian(3);
    
    double sintheta = std::sin(spherical[2]);
    cartesian[0] = spherical[0]*sintheta*std::cos(spherical[1]);
    cartesian[1] = spherical[0]*sintheta*std::sin(spherical[1]);
    cartesian[2] = spherical[0]*std::cos(spherical[2]);
    
    return cartesian;
  }

  template <typename Container = std::vector<double>>
  Container cartesian2polar(std::vector<double> const& cartesian)
  {
    Container polar(2);
    
    polar[0] = operation::abs(cartesian);
    polar[1] = std::atan2(cartesian[1], cartesian[0]);
    
    return polar;
  }

  template <typename Container = std::vector<double>>
  Container polar2cartesian(std::vector<double> const& polar)
  {
    Container cartesian(2);
    
    cartesian[0] = polar[0]*std::cos(polar[1]);
    cartesian[1] = polar[0]*std::sin(polar[1]);
    
    return cartesian;
  }
}

#endif /* Coordinates_h */
