//
// Roots.h
// general
//
// Created by Tomas Aquino on 8/6/19.
// Copyright © 2019 Tomas Aquino. All rights reserved.
//

// Root-finding algorithms

#ifndef Roots_h
#define Roots_h

namespace roots
{
  class NewtonRaphson
  {
  private:

    bool converged{ 0 };
    double func_val{ 0. };
    double derivative_val{ 0. };

  public:

    template <typename Func, typename Derivative>
    double operator()(Func func, Derivative derivative, double guess, double tolerance, std::size_t max_count)
    {
      std::size_t count = 0;
      func_val = func(guess);

      do
      {
        ++count;
        guess = guess - func_val / derivative(guess);
        func_val = func(guess);
      }
      while (count < max_count && std::abs(func_val) > tolerance);

      converged = (count < max_count);
      return guess;
    }

    bool Converged()
    { return converged; }

    double Func()
    { return func_val; }

    double Derivative()
    { return derivative_val; }
  };

  class Bisection
  {
  public:

    template <typename Func>
    double operator()(Func func, double lower_bracket, double upper_bracket, double tolerance)
    {
      int nr_iterations = std::log((upper_bracket - lower_bracket) / tolerance) / std::log(2.) - 1;

      double mid_point = (lower_bracket + upper_bracket) / 2.;

      for (int count = 0; count < nr_iterations; ++count)
      {
        if (func(lower_bracket) * func(mid_point) < 0.) upper_bracket = mid_point;
        else lower_bracket = mid_point;
        mid_point = (lower_bracket + upper_bracket) / 2.;
      }

      return mid_point;
    }
  };
}

#endif /* Roots_h */
