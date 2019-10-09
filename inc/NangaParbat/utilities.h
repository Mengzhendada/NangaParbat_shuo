//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include <yaml-cpp/yaml.h>
#include <apfel/matrix.h>

namespace NangaParbat
{
  /**
   * @brief Utility function to generate an interpolation grid
   * @param n: the number of nodes of the grid
   * @param min: the lower bound
   * @param max: the upper bound
   * @param ext: the number of extra nodes (default: 0)
   * @param lgt: whether the grid is logarithmically spaced (default: false)
   * @return a linearly-spaced grid with "n" nodes between "min" and
   * "max" plus "ext" extra nodes on the right.
   */
  std::vector<double> GenerateGrid(int const& n, double const& min, double const& max, int const& ext = 0, bool const& lgt = false);

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

  /**
   * @brief b* prescription a la Collins
   * @param b: the impact parameter
   * @param Q: the hard scale
   * @return b*
   */
  double bstar(double const& b, double const& Q);

  /**
   * @brief b* prescription with bmin
   * @param b: the impact parameter
   * @param Q: the hard scale
   * @return b*
   */
  double bstarmin(double const& b, double const& Q);

  /**
   * @brief Function that lists elements in a directory
   * @param path: path to the directory
   * @return a vector of elements
   */
  std::vector<std::string> list_dir(std::string const& path);
}
