//
// Operations.h
// general
//
// Created by Tomas Aquino on 8/6/19.
// Copyright © 2019 Tomas Aquino. All rights reserved.
//

// Miscelaneous operations on containers
// Notes:
//    Many methods assume containers with consistent sizes are passed in
//    Many methods require random access and operator []
//    Some operations involve casting.
//    The latter are spelled out explicitly for clarity and to avoid warnings
//    In most cases, the return value type is the type of the first container

#ifndef Operations_h
#define Operations_h

#include <cmath>
#include <functional>
#include "general/useful.h"

namespace operation
{
  // Sum of elements
  template <typename Container>
  auto sum(Container const& input)
  {
    typename Container::value_type output{};
    for (auto const& val : input)
      output += val;
    return output;
  }

  // Product of elements
  template <typename Container>
  auto prod(Container const& input)
  {
    typename Container::value_type output{ 1. };
    for (auto const& val : input)
      output *= val;
    return output;
  }

  // Element-wise sum of scalar
  template <typename Container, typename Scalar, typename Container_out>
  void plus_scalar(Container const& input, Scalar cc, Container_out& output)
  {
    if constexpr (useful::has_plus_v<Container, Scalar>)
      output = input + cc;
    else
    {
      using value_type = typename Container_out::value_type;
      for (std::size_t ii = 0; ii < output.size(); ++ii)
        output[ii] = value_type(input[ii])+value_type(cc);
    }
  }

  template <typename Container, typename Scalar>
  auto plus_scalar(Container const& input, Scalar cc)
  {
    if constexpr (useful::has_plus_v<Container, Scalar>)
      return input + cc;
    else
    {
      Container output(input.size());
      plus_scalar(input, cc, output);
      return output;
    }
  }

  template <typename Container, typename Scalar>
  void plus_scalar_InPlace(Container& input, Scalar cc)
  { plus_scalar(input, cc, input); }

  // Element-wise sum
  template <typename Container_1, typename Container_2, typename Container_out>
  void plus(Container_1 const& input_1, Container_2 const& input_2, Container_out& output)
  {
    if constexpr (useful::has_plus_v<Container_1, Container_2>)
      output = input_1 + input_2;
    else
    {
      using value_type = typename Container_out::value_type;
      for (std::size_t ii = 0; ii < output.size(); ++ii)
        output[ii] = value_type(input_1[ii]) + value_type(input_2[ii]);
    }
  }

  template <typename Container_1, typename Container_2>
  auto plus(Container_1 const& input_1, Container_2 const& input_2)
  {
    if constexpr (useful::has_plus_v<Container_1, Container_2>)
      return input_1 + input_2;
    else
    {
      Container_1 output(input_1.size());
      plus(input_1, input_2, output);
      return output;
    }
  }

  template <typename Container_1, typename Container_2>
  void plus_InPlace(Container_1& input_1, Container_2 const& input_2)
  { plus(input_1, input_2, input_1); }

  // Element-wise subtraction of scalar
  template <typename Container, typename Scalar, typename Container_out>
  void minus_scalar(Container const& input, Scalar cc, Container_out& output)
  {
    if constexpr (useful::has_minus_v<Container, Scalar>)
      output = input - cc;
    else
    {
      using value_type = typename Container_out::value_type;
      for (std::size_t ii = 0; ii < output.size(); ++ii)
        output[ii] = value_type(input[ii]) - value_type(cc);
    }
  }

  template <typename Container, typename Scalar>
  auto minus_scalar(Container const& input, Scalar cc)
  {
    if constexpr (useful::has_minus_v<Container, Scalar>)
      return input - cc;
    else
    {
      Container output(input.size());
      minus_scalar(input, cc, output);
      return output;
    }
  }

  template <typename Container, typename Scalar>
  void minus_scalar_InPlace(Container& input, Scalar cc)
  { minus_scalar(input, cc, input); }
  
  // Element-wise subtraction of scalar
  template <typename Container, typename Scalar, typename Container_out>
  void scalar_minus(Scalar cc, Container const& input, Container_out& output)
  {
    if constexpr (useful::has_minus_v<Container, Scalar>)
      output = cc - input;
    else
    {
      using value_type = typename Container_out::value_type;
      for (std::size_t ii = 0; ii < output.size(); ++ii)
        output[ii] = value_type(cc)-value_type(input[ii]);
    }
  }

  template <typename Container, typename Scalar>
  auto scalar_minus(Scalar cc, Container const& input)
  {
    if constexpr (useful::has_minus_v<Container, Scalar>)
      return cc - input;
    else
    {
      Container output(input.size());
      scalar_minus(cc, input, input);
      return output;
    }
  }

  template <typename Container, typename Scalar>
  void scalar_minus_InPlace(Scalar cc, Container& input)
  { scalar_minus(cc, input, input); }

  // Element-wise subtraction
  template <typename Container_1, typename Container_2, typename Container_out>
  void minus(Container_1 const& input_1, Container_2 const& input_2, Container_out& output)
  {
    if constexpr (useful::has_minus_v<Container_1, Container_2>)
      output = input_1 - input_2;
    else
    {
      using value_type = typename Container_out::value_type;
      for (std::size_t ii = 0; ii < output.size(); ++ii)
        output[ii] = value_type(input_1[ii]) - value_type(input_2[ii]);
    }
  }

  template <typename Container_1, typename Container_2>
  auto minus(Container_1 const& input_1, Container_2 const& input_2)
  {
    if constexpr (useful::has_minus_v<Container_1, Container_2>)
      return input_1 - input_2;
    else
    {
      Container_1 output(input_1.size());
      minus(input_1, input_2, output);
      return output;
    }
  }

  template <typename Container_1, typename Container_2>
  void minus_InPlace(Container_1& input_1, Container_2 const& input_2)
  { minus(input_1, input_2, input_1); }

  // Element-wise multiplication by scalar
  template <typename Container, typename Scalar, typename Container_out>
  void times_scalar(Scalar lambda, Container const& input, Container_out& output)
  {
    if constexpr (useful::has_multiplies_v<Scalar, Container>)
      output = lambda*input;
    else
    {
      using value_type = typename Container_out::value_type;
      for (std::size_t ii = 0; ii < input.size(); ++ii)
        output[ii] = value_type(lambda)*value_type(input[ii]);
    }
  }

  template <typename Container, typename Scalar>
  auto times_scalar(Scalar lambda, Container const& input)
  {
    if constexpr (useful::has_multiplies_v<Scalar, Container>)
      return lambda*input;
    else
    {
      Container output(input.size());
      times_scalar(lambda, input, output);
      return output;
    }
  }

  template <typename Container, typename Scalar>
  void times_scalar_InPlace(Scalar lambda, Container& input)
  { times_scalar(lambda, input, input); }

  // Element-wise multiplication
  template <typename Container_1, typename Container_2, typename Container_out>
  void times(Container_1 const& input_1, Container_2 const& input_2, Container_out& output)
  {
    if constexpr (useful::has_multiplies_v<Container_1, Container_2>)
      output = input_1*input_2;
    else
    {
      using value_type = typename Container_out::value_type;
      for (std::size_t ii = 0; ii < output.size(); ++ii)
        output[ii] = value_type(input_1[ii])*value_type(input_2[ii]);
    }
  }

  template <typename Container_1, typename Container_2>
  auto times(Container_1 const& input_1, Container_2 const& input_2)
  {
    if constexpr (useful::has_multiplies_v<Container_1, Container_2>)
      return input_1*input_2;
    else
    {
      Container_1 output(input_1.size());
      times(input_1, input_2, output);

      return output;
    }
  }

  template <typename Container_1, typename Container_2>
  void times_InPlace(Container_1& input_1, Container_2 const& input_2)
  { times(input_1, input_2, input_1); }

  // Element-wise division by scalar
  template <typename Container, typename Scalar, typename Container_out>
  void div_scalar(Container const& input, Scalar lambda, Container_out& output)
  {
    if constexpr (useful::has_divides_v<Container, Scalar>)
      output = input/lambda;
    else
    {
      using value_type = typename Container_out::value_type;
      for (std::size_t ii = 0; ii < input.size(); ++ii)
        output[ii] = value_type(input[ii]/lambda);
    }
  }

  template <typename Container, typename Scalar>
  auto div_scalar(Container const& input, Scalar lambda)
  {
    if constexpr (useful::has_divides_v<Container, Scalar>)
      return input/lambda;
    else
    {
      Container output(input.size());
      div_scalar(input, lambda, output);

      return output;
    }
  }

  template <typename Container, typename Scalar>
  void div_scalar_InPlace(Container& input, Scalar lambda)
  { div_scalar(input, lambda, input); }

  // Element-wise division
  template <typename Container_1, typename Container_2, typename Container_out>
  void div(Container_1 const& input_1, Container_2 const& input_2, Container_out& output)
  {
    if constexpr (useful::has_divides_v<Container_1, Container_2>)
      output = input_1/input_2;
    else
    {
      using value_type = typename Container_out::value_type;
      for (std::size_t ii = 0; ii < output.size(); ++ii)
        output[ii] = value_type(input_1[ii]/input_2[ii]);
    }
  }

  template <typename Container_1, typename Container_2>
  auto div(Container_1 const& input_1, Container_2 const& input_2)
  {
    if constexpr (useful::has_divides_v<Container_1, Container_2>)
      return input_1/input_2;
    else
    {
      Container_1 output(input_1.size());
      div(input_1, input_2, output);

      return output;
    }
  }

  template <typename Container_1, typename Container_2>
  void div_InPlace(Container_1& input_1, Container_2 const& input_2)
  { div(input_1, input_2, input_1); }

  // lambda_1*input_1 + lambda_2*input2
  template <typename Container_1, typename Type_1,
  typename Container_2, typename Type_2,
  typename Container_out>
  void linearOp
  (Type_1 lambda_1, Container_1 const& input_1,
   Type_2 lambda_2, Container_2 const& input_2,
   Container_out& output)
  {
    if constexpr (useful::has_multiplies_v<Type_1, Container_1>
                  && useful::has_multiplies_v<Type_2, Container_2>
                  && useful::has_plus_v<Container_1, Container_2>)
      output = lambda_1*input_1 + lambda_2*input_2;
    else
    {
      using value_type = typename Container_out::value_type;
      for(size_t ii = 0; ii < output.size(); ++ii)
        output[ii] = value_type(lambda_1) * input_1[ii]
        + value_type(lambda_2) * input_2[ii];
    }
  }

  template <typename Container_1, typename Type_1,
  typename Container_2, typename Type_2>
  auto linearOp
  (Type_1 lambda_1, Container_1 const& input_1,
   Type_2 lambda_2, Container_2 const& input_2)
  {
    if constexpr (useful::has_multiplies_v<Type_1, Container_1>
                && useful::has_multiplies_v<Type_2, Container_2>
                && useful::has_plus_v<Container_1, Container_2>)
      return lambda_1*input_1 + lambda_2*input_2;
    else
    {
      Container_1 output(input_1.size());
      linearOp(lambda_1, input_1, lambda_2, input_2, output);
      return output;
    }
  }

  template <typename Container_1, typename Type_1,
  typename Container_2, typename Type_2>
  void linearOp_InPlace
  (Type_1 lambda_1, Container_1& input_1,
   Type_2 lambda_2, Container_2 const& input_2)
  { linearOp(lambda_1, input_1, lambda_2, input_2, input_1); }

  // lambda*input_1 + input2
  template <typename Container_1, typename Type,
  typename Container_2, typename Container_out>
  void linearOp
  (Type lambda, Container_1 const& input_1,
   Container_2 const& input_2, Container_out& output)
  {
    if constexpr (useful::has_multiplies_v<Type, Container_1>
                  && useful::has_plus_v<Container_1, Container_2>)
      output = lambda*input_1 + input_2;
    else
    {
      using value_type = typename Container_out::value_type;
      for(size_t ii = 0; ii < output.size(); ++ii)
        output[ii] = value_type(lambda)*input_1[ii] +
          input_2[ii];
    }
  }

  template <typename Container_1, typename Type, typename Container_2>
  auto linearOp
  (Type lambda, Container_1 const& input_1, Container_2 const& input_2)
  {
    if constexpr (useful::has_multiplies_v<Type, Container_1>
                  && useful::has_plus_v<Container_1, Container_2>)
      return lambda*input_1 + input_2;
    else
    {
      Container_1 output;
      LinearOp(lambda, input_1, input_2, output);

      return output;
    }
  }

  template <typename Container_1, typename Type,
  typename Container_2>
  void linearOp_InPlace
  (Type lambda, Container_1& input_1, Container_2 const& input_2)
  { LinearOp(lambda, input_1, input_2, input_1); }

  // Element-wise square
  template <typename Container, typename Container_out>
  void square(Container const& input, Container_out& output)
  {
    if constexpr (useful::has_multiplies_v<Container, Container>)
      output = input*input;
    else
    {
      using value_type = typename Container_out::value_type;
      for (std::size_t ii = 0; ii < output.size(); ++ii)
        output[ii] = value_type(input[ii])*value_type(input[ii]);

      return output;
    }
  }

  template <typename Container>
  auto square(Container const& input)
  {
    if constexpr (useful::has_multiplies_v<Container, Container>)
      input*input;
    else
    {
      Container output(input.size());
      square(input,output);

      return output;
    }
  }

  template <typename Container>
  void square_InPlace(Container& input)
  { square(input, input); }

  // Element-wise square root
  template <typename Container, typename Container_out>
  void sqrt(Container const& input, Container_out& output)
  {
    if constexpr (useful::can_call_sqrt_v<Container>)
      output = std::sqrt(input);
    else
    {
      using value_type = typename Container_out::value_type;
      for (std::size_t ii = 0; ii < output.size(); ++ii)
        output[ii] = std::sqrt(value_type(input[ii]));
    }
  }

  template <typename Container>
  auto sqrt(Container const& input)
  {
    if constexpr (useful::can_call_sqrt_v<Container>)
      return std::sqrt(input);
    else
    {
      Container output(input.size());
      sqrt(input, output);

      return output;
    }
  }

  template <typename Container>
  void sqrt_InPlace(Container& input)
  { sqrt(input, input); }

  // Element-wise mean of two vectors
  template <typename Container_1, typename Container_2,
  typename Container_out>
  void mean
  (Container_1 const& input_1, Container_2 const& input_2,
   Container_out& output)
  {
    if constexpr (useful::has_multiplies_v<Container_1, Container_2>
                  && useful::has_multiplies_v<Container_1>)
      output = 0.5*(input_1+input_2);
    else
    {
      using value_type = typename Container_out::value_type;
      for (std::size_t ii = 0; ii <input_1.size(); ++ii)
        output[ii] = (value_type(input_1[ii]) + value_type(input_2[ii]))/2.;
    }
  }

  template <typename Container_1, typename Container_2>
  auto mean(Container_1 const& input_1, Container_2 const& input_2)
  {
    if constexpr (useful::has_multiplies_v<Container_1, Container_2>
                  && useful::has_multiplies_v<Container_1>)
      return 0.5*(input_1+input_2);
    else
    {
      Container_1 output(input_1.size());
      mean(input_1, input_2, output);
      return output;
    }
  }

  template <typename Container>
  void mean_InPlace(Container const& input_1, Container const& input_2)
  { mean(input_1, input_2, input_1); }

  // Euclidean norm squared
  template <typename Container>
  auto abs_sq(Container const& input)
  {
    if constexpr (useful::can_call_abs_v<Container>)
    {
      auto abs = std::abs(input);
      return abs*abs;
    }
    else
    {
      typename Container::value_type abs_sq{ 0. };
      for (auto const& val : input)
        abs_sq += val*val;
      return abs_sq;
    }
  }

  // Euclidean norm
  template <typename Container>
  auto abs(Container const& input)
  {
    if constexpr (useful::can_call_abs_v<Container>)
      return std::abs(input);
    else
      return std::sqrt(abs_sq(input));
  }

  // Like std::adjacent_difference but without the first element
  template
  <class InputIterator, class OutputIterator,
  class BinaryOperation>
  OutputIterator adjacent_difference
  (InputIterator first, InputIterator last,
   OutputIterator result,
   BinaryOperation binary_op = std::minus<typename InputIterator::value_type>{})
  {
    if (first != last)
    {
      InputIterator prev = first++;
      while (first != last)
      {
        InputIterator val = first++;
        *result++ = binary_op(*val, *prev);
        prev = val;
      }
    }
    return result;
  }

  // Dot product
  template <typename Container>
  auto dot(Container const& input_1, Container const& input_2)
  {
    typename Container::value_type result{};
    for(size_t ii = 0; ii < input_1.size(); ++ii)
      result += input_1[ii] * input_2[ii];
    return result;
  }

  // Container of containers dotted into container, = sum_j (input_1)_{ij} (input_2)_j
  template <typename Container_outer, typename Container_inner>
  void dot
  (Container_outer const& input_1, Container_inner const& input_2, Container_inner& output)
  {
    std::size_t counter = 0;
    for (auto const& vec : input_1)
      output[counter++] = Dot(vec, input_2);
  }

  template <typename Container_outer, typename Container_inner>
  auto dot(Container_outer const& input_1, Container_inner const& input_2)
  {
    Container_inner output(input_2.size());
    Dot(input_1, input_2, output);

    return output;
  }

  // Total number of elements in a container of containers
  template< template<class> class Container_outer, typename Container_inner >
  std::size_t nr_elements(Container_outer<Container_inner> const& vv)
  {
    std::size_t nr_elements = 0;
    for( auto const& cont : vv )
      nr_elements += cont.size();
    return nr_elements;
  }

  // Convolution sum
  template <typename Container>
  typename Container::value_type convolution
  (Container const& cont_1, Container const& cont_2,
   std::size_t idx_start, std::size_t idx_end)
  {
    typename Container::value_type res = 0.;
    for (std::size_t ii = idx_start; ii < idx_end; ++ii)
      res += cont_1[ii]*cont_2[idx_end-1-ii];
    return res;
  }

  // Convolution integral using trapezoidal rule
  template <typename Container>
  typename Container::value_type convolution_trap
  (Container const& cont_1, Container const& cont_2,
   std::size_t idx_start, std::size_t idx_end)
  {
    typename Container::value_type res = 0.;
    res += 0.5*(cont_1[idx_start]*cont_2[idx_end] + cont_1[idx_end]*cont_2[idx_start]);
    for (std::size_t ii = idx_start + 1; ii < idx_end; ++ii)
      res += cont_1[ii]*cont_2[idx_end-ii];
    return res;
  }

  // Hamming distance
  template <typename Container>
  auto hamming(Container const& vec, Container const& other_vec)
  {
    typename Container::value_type hamming = 0;
    auto it = vec.begin();
    auto other_it = other_vec.begin();
    for (; it != vec.end(); it++, other_it++)
      hamming += *it > *other_it ? *it - *other_it : *other_it - *it;
    return hamming;
  }
  
  template <std::size_t dd, typename Container>
  auto project(Container const& container)
  { return container[dd]; }
  
  template <>
  auto project<0, double>(double const& val)
  { return val; }
  
  template <>
  auto project<0, int>(int const& val)
  { return val; }
  
  template <>
  auto project<0, std::size_t>(std::size_t const& val)
  { return val; }
  
  std::size_t factorial(std::size_t nn)
  {
    if (nn == 0)
      return 1;
    return nn*factorial(nn-1);
  }

  std::size_t factorial_incomplete(std::size_t nn, std::size_t mm)
  {
    std::size_t result = 1.;
    ++nn;
    for (std::size_t ii = 0; ii < mm && nn --> 0; ++ii)
      result *= nn;
    return result;
  }
}


#endif /* Operations_h */
