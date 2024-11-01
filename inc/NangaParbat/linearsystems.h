//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include <apfel/matrix.h>

namespace NangaParbat
{
  /**
   * @brief Cholesky decomposition of the covariance matrix.
   * @param V: the covariance matrix
   * @return The Cholesky decomposition matrix L such that L<SUP>T</SUP>L = V
   */
  apfel::matrix<double> CholeskyDecomposition(apfel::matrix<double> const V);

  /**
   * @brief Solve lower-diagonal system of equations by forward substitution
   * @param L: lower-diagonal matrix
   * @param y: vector of constants
   * @return the solution vector x
   */
  std::vector<double> SolveLowerSystem(apfel::matrix<double> L, std::vector<double> y);

  /**
   * @brief Solve upper-diagonal system of equations by backward substitution
   * @param U: upper-diagonal matrix
   * @param y: vector of constants
   * @return the solution vector
   */
  std::vector<double> SolveUpperSystem(apfel::matrix<double> U, std::vector<double> y);

  /**
   * @brief Solve symmetric system of equations
   * @param A: symmetric matrix
   * @param rho: vector of constants
   * @return the solution vector
   */
  std::vector<double> SolveSymmetricSystem(apfel::matrix<double> A, std::vector<double> rho);
}
