//
// Ranges.h
// general
//
// Created by Tomas Aquino on 8/6/19.
// Copyright © 2019 Tomas Aquino. All rights reserved.
//

// Initialize and set containers to ranges

#ifndef Ranges_h
#define Ranges_h

#include <cmath>
#include <algorithm>

namespace range
{
  // Linear partition of [x_L, x_R]
  // with nr_points points
  template <typename Container = std::vector<double>>
  Container linspace(double x_L, double x_R, std::size_t nr_points)
  {
    double inc = (x_R - x_L) / (nr_points - 1);
    
    Container output(nr_points);
    output[0] = x_L;
    for (std::size_t ii = 1; ii < nr_points - 1; ++ii)
      output[ii] = x_L + ii*inc;
    output[nr_points-1] = x_R;

    return output;
  }

  // log partition of [x_L, x_R]
  // with nr_points points in output
  template <typename Container = std::vector<double>>
  Container logspace(double x_L, double x_R, std::size_t nr_points)
  {
    double inc = std::pow(x_R / x_L, 1. / (nr_points - 1));

    Container output(nr_points);
    output[0] = x_L;
    for (std::size_t ii = 1; ii < nr_points - 1; ++ii)
      output[ii] = output[ii-1] * inc;
    output[nr_points-1] = x_R;

    return output;
  }

  // Linear partition of [x_L, x_R]
  // with spacing given by increment
  template <typename Container = std::vector<double>>
  Container range(double x_L, double increment, double x_R)
  {
    if (x_R < x_L) increment = -std::abs(increment);
    else increment = std::abs(increment);

    // Make sure of size due to roundoffs
    std::size_t size = 0;
    for (double xx = x_L; xx < x_R + increment / 2.; xx += increment)
      ++size;
    
    if (size == 0)
      return { x_L };

    Container output(size);
    double xx = x_L;
    for (std::size_t ii = 0; ii < size; ++ii)
    {
      output[ii] = xx;
      xx += increment;
    }

    return output;
  }

  // Populate container with sequence first first+1 ... first+size-1
  template <typename OutputIterator, typename Size, typename Assignable>
  void iota_n(OutputIterator first, Size size, Assignable value)
  {
    std::generate_n(first, size, [&value]() { return value++; });
  }

  // Use with std::generate to populate vector with the sequence (seed, seed + 1, seed + 2, ...)
  template <typename Type>
  struct gen
  {
    Type val;
    gen(Type seed)
    : val(seed)
    {}

    Type operator()()
    { return val++; }
  };
}

#endif /* Ranges_h */
