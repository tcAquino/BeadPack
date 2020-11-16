//
// Algebra.h
// general
//
// Created by Tomas Aquino on 8/6/19.
// Copyright © 2019 Tomas Aquino. All rights reserved.
//

// Miscellaneous matrix algorithms

#ifndef Algebra_h
#define Algebra_h

//#include <cmath>
#include <vector>
#include "general/Operations.h"

namespace algebra
{
  std::vector<std::vector<double>> gram_schmidt
  (std::vector<std::vector<double>> const& input)
  {
    std::size_t nr_vectors = input.size();
    
    std::vector<std::vector<double>> output;
    output.reserve(nr_vectors);
    output.push_back(operation::div_scalar(input[0], operation::abs(input[0])));
    for (std::size_t ii = 1; ii < nr_vectors; ++ii)
    {
      output.push_back(input[ii]);
      for (std::size_t jj = 0; jj < ii; ++jj)
        operation::minus_InPlace(output[ii],
                                 operation::times_scalar(operation::dot(output[jj],
                                                                        output[ii])
                                                         /operation::abs_sq(output[jj]),
                                                         output[jj]));
      operation::div_scalar_InPlace(output[ii], operation::abs(output[ii]));
    }
    
    return output;
  }
  
  std::vector<std::vector<double>> gram_schmidt(std::vector<double> const& input)
  {
    std::size_t dim = input.size();
    
    std::vector<std::vector<double>> input_vecs;
    input_vecs.reserve(dim);
    input_vecs.push_back(operation::div_scalar(input, operation::abs(input)));
     
    std::vector<double> basis_vector(dim);
    for (std::size_t dd = 0; dd < dim && input_vecs.size() < dim; ++dd)
    {
      basis_vector[dd] = 1.;
      if (input != basis_vector)
        input_vecs.push_back(basis_vector);
      basis_vector[dd] = 0.;
    }
    
    return gram_schmidt(input_vecs);
  }
  
  // Thomas Algorithm for tridiagonal systems
  void Solve_Tridiag
  (std::vector<double> const& l_diag,
   std::vector<double> const& diag,
   std::vector<double> const& u_diag,
   std::vector<double> const& ind,
   std::vector<double>& sol)
  {
    std::size_t dim = diag.size();
    sol.resize(dim);

    // Modified coefficients
    std::vector<double> u_diag_new(dim);
    std::vector<double> ind_new(dim);
    u_diag_new[0] = u_diag[0]/diag[0];
    ind_new[0] = ind[0]/diag[0];
    double temp;
    for (size_t ii = 1; ii < dim; ++ii)
    {
      temp = diag[ii] - l_diag[ii]*u_diag_new[ii-1];
      u_diag_new[ii] = u_diag[ii]/temp;
      ind_new[ii] = (ind[ii] - l_diag[ii]*ind_new[ii-1])/temp;
    }

    // Back-substitution
    sol.back() = ind_new.back();

    for (size_t ii = dim-2; ii != 0; ii--)
      sol[ii] = ind_new[ii] - u_diag_new[ii]*sol[ii+1];
  }

  // Uses Sherman-Morrison and the Thomas Algorithm
  // (see Solve_Tridiag)
  // to solve a cyclic tridiagonal system Ax = ind
  // alpha = A_{ n-1, 0 } , beta = A_{ 0, n-1 }
  void Solve_Cyclic_Tridiag
  (double alpha, double beta,
   std::vector<double> const& l_diag,
   std::vector<double> const& diag,
   std::vector<double> const& u_diag,
   std::vector<double> const& ind,
   std::vector<double>& sol)
  {
    std::size_t nr_points = diag.size();
    sol.resize(nr_points);
    double gamma = - diag[0];

    // Sherman-Morrison v = (1 0 0 ... 0 beta / gamma)
    std::vector<double> zz, uu(nr_points, 0.), diag_mod;
    diag_mod = diag;
    diag_mod[0] -= gamma;
    diag_mod.back() -= alpha*beta/gamma;
    uu[0] = gamma;
    uu.back() = alpha;

    Solve_Tridiag(l_diag, diag_mod, u_diag, ind, sol);
    Solve_Tridiag(l_diag, diag_mod, u_diag, uu, zz);
    double factor = (sol[0] + beta*sol.back()/gamma)/
    (1 + zz[0] +
     beta*zz.back()/gamma);

    for(size_t ii = 0; ii < nr_points; ++ii)
      sol[ii] -= factor * zz[ii];
  }

  // Find the determinant of a Tridiagonal matrix
  // using the continuant recurrence relation
  double Det_Tridiag(std::vector<double> const& l_diag, std::vector<double> const& diag, std::vector<double> const& u_diag)
  {
    double cont_0 = 1.;
    double cont_1 = diag[0];
    double cont_2 = 0.;

    for	(size_t ii = 2; ii <= diag.size(); ++ii)
    {
      cont_2 = diag[ii-1]*cont_1 - u_diag[ii-2]*l_diag[ii-1]*cont_0;
      cont_0 = cont_1;
      cont_1 = cont_2;
    }

    return cont_2;
  }
}


#endif /* Algebra_h */
