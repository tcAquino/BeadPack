//
// Modular.h
// general
//
// Created by Tomas Aquino on 8/6/19.
// Copyright © 2019 Tomas Aquino. All rights reserved.
//

// Modular arithmetic, basis change, and similar operations

#ifndef Modular_h
#define Modular_h

#include <cstdlib>

namespace modular
{
  // Return number + 1, wrapped around with mod = 0
  std::size_t circPlus(std::size_t number, std::size_t mod)
  { return ++number == mod ? 0 : number; }

  // Return number - 1, wrapped around with -1 = mod - 1
  std::size_t circMinus(std::size_t number, std::size_t mod)
  { return number == 0 ? mod - 1 : --number; }

  std::size_t circPlus(std::size_t number, std::size_t cc, std::size_t mod)
  { return (number + cc) % mod; }

  std::size_t circMinus(std::size_t number, std::size_t cc, std::size_t mod)
  { return number < cc ? (mod - (cc - number) % mod) % mod : number - cc; }

  // Convert base basis to base 10
  template <typename Container>
  auto convert(Container const& to_convert, std::size_t basis)
  {
    typename Container::value_type converted = 0;
    typename Container::value_type power = 1;
    for (std::size_t pos = to_convert.size(); pos --> 0;)
    {
      converted += to_convert[pos]*power;
      power *= basis;
    }

    return converted;
  }

  // Convert the base-10 integer <to_convert>
  // to base <basis> returned as a string of digits in <converted>
  // Note: Converted must be passed with the correct size.
  template <typename Container, typename Integer_Type>
  void convert(Integer_Type to_convert, Integer_Type basis, Container& converted)
  {
    Integer_Type to_convert_old;
    for (std::size_t pos = converted.size(); pos --> 0;)
    {
      to_convert_old = to_convert;
      to_convert /= basis;
      converted[pos] = to_convert_old - to_convert*basis;
    }
  }

  // Given 1d index idx, convert to position (x_idx,y_idx,...)
  // The 1d index goes through every x for each value of y, every y for each value of z, etc, in ascending order
  // Note: position must be passed with the correct size.
  template <typename Container, typename IntegerType = std::size_t>
  void indexToPosition(IntegerType idx, Container const& nr_points, Container& position)
  {
    for (size_t dd = 0; dd < nr_points.size(); ++dd)
    {
      position[dd] = idx%nr_points[dd];
      idx /= nr_points[dd];
    }
  }

  // Given position (x_idx,y_idx,...), convert to 1d index
  // The 1d index goes through every x for each value of y, every y for each value of z, etc, in ascending order
  template <typename Container>
  auto positionToIndex(Container const& position, Container const& nr_points)
  {
    typename Container::value_type idx = 0;
    typename Container::value_type offset = 1;
    for (size_t dd = 0; dd < nr_points.size(); ++dd)
    {
      idx += offset * position[dd];
      offset *= nr_points[dd];
    }

    return idx;
  }

  // num are the digits of a basis <basis> number (units digit on num[0])
  // Increment it by 1 (max value wraps to 00....0)
  template <typename Container, typename IntegerType>
  void increment(Container& num, IntegerType basis)
  {
    bool carry;
    for (size_t dig = 0; dig < num.size(); ++dig)
    {
      carry = ++num[dig]/basis;
      num[dig] %= basis;
      if (!carry) break;
    }
  }

  // num are the digits of a number such that each digit ranges from 0 to basis[dig] - 1 (units digit on num[0])
  // Increment it by 1 (max value wraps to 00....0)
  template <typename Container>
  void increment(Container& num, Container const& basis)
  {
    bool carry;
    for (size_t dig = 0; dig < num.size(); ++dig)
    {
      carry = ++num[dig]/basis[dig];
      num[dig] %= basis[dig];
      if (!carry) break;
    }
  }
}

#endif /* Modular_h */
