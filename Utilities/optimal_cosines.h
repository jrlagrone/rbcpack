/*
  Copyright 2015 John LaGrone

  This file is part of the Yee CRBC Library from RBCPACK.

  The Yee CRBC Library is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by 
  the Free Software Foundation, either version 3 of the License, or (at your 
  option) any later version.

  The Yee CRBC Library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of 
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with the Yee CRBC Library.  If not, see <http://www.gnu.org/licenses/>.
*/

/** \file
 * \brief Functions to compute the optimal cosine parameters
 *
 * This declares the prototypes for functions to compute the optimal cosines
 * {\f$ a_j \f$, j=1,...,2p} which minimize the maximum norm on 0 < x < 1 of:
 *
 *  \f$ e(x) = \exp(-\eta/x)*((1-x)/(1+x))*\prod_{j=1}^{2p} (a_j-x)/(a_j+x) \f$
 *
 * There are two options available. The first is to directly computed the optimal
 * cosines for a set of given parameters. The second is to try to compute a set
 * of of optimal cosines using the fewest terms necessary to meet a given error
 * tolerance.
 *  
 * These are parameters used to contruct optimal local radiation boundary
 * conditions for the wave equation and related systems. Details can be 
 * found in [1]. 
 *
 * As suggested in [1], we will use the Remez algorithm to approximate the
 * solution to the minimax problem. The implementation here is based upon the 
 * second presentation of the Remez algorithm in [2]. A more procedural and 
 * example based description of the algorithm is available at [3].
 *
 * The remez algorithm requires a root finding algorithm. For this, we use
 * Brent's Method (see, e.g., [4,5]).
 *
 * ------------------------------------------------------------------------------
 * References:
 * ------------------------------------------------------------------------------
 * \todo change to BiBTeX style citations
 *
 * [1] T. Hagstrom and T. Warburton. Complete radiation boundary conditions: 
 *     Minimizing the long time error growth of local methods. SIAM Journal on 
 *     Numerical Analysis, 47(5):3678–3704, 2009. 
 *
 *     (A preprint can be found at http://faculty.smu.edu/thagstrom/rbcpac.html
 *      along with a Matlab implementation)
 *
 * [2] P. Petrushev and V. Popov, Rational Approximation of Real Functions, vol.
 *     28 of Encyclopedia of Mathematics, Cambridge University Press, Cambridge, 
 *     1987.
 *
 * [3] The Remez Method. 
 *     http://www.boost.org/doc/libs/1_56_0/libs/math/doc/html/math_toolkit/remez.html
 *
 * [4] Brent, R.P. 1973, Algorithms for Minimization without Derivatives 
 *     (Englewood Cliffs, NJ: Prentice-Hall), Chapter 5.
 *
 * [5] Numerical Recipes in C, The Art of Scientific Computing, Second Edition, 
 *     William H. Press, Saul A. Teukolsky, William T. Vetterling, and 
 *     Brian P. Flannery. Cambridge University Press. 1988, 1992. 
 *
 * ------------------------------------------------------------------------------
 * Required External Libraries
 * ------------------------------------------------------------------------------
 * [a] LAPACK. http://www.netlib.org/lapack/
 *
 * [b] BLAS. http://www.netlib.org/blas/
 *
 */

#ifndef OPTIMAL_COSINES_H
#define OPTIMAL_COSINES_H
# ifdef __cplusplus
extern "C"
{
# endif

/** \brief Function to attempt to generate the minimum number of cosines 
 *         provided an error tolerance.
 *
 *  input:
 *    \param[in] eta  - a dimensionless parameter given by
 *      
 *                   eta = delta / (cT),
 *
 *           where delta is the seperation of the boundary from sources, 
 *           scatterers, and other inhomogeneities, c is the characterist wave 
 *           speed, and T is the final time. Essentially, this is a measure of
 *           how difficult the problem is to simulate.
 *
 *           We enforce the restriction 1e-7 <= eta <= 0.1.
 *
 *    \param[in] pmax - Maximum number of recursions.
 *                      We enforce the resctriction pmax <= 40
 *
 *    \param[in] tol  - Maximum acceptable value for maximum norm of e(x)
 *    
 *
 *  output:
 *    \param[out] a       - array of length 2P containing the cosines in descending 
 *                          order.
 *                          \warning the array 'a' must be allocated to be at 
 *                          least length 2*P on calling.
 * 
 *    \param[out] P       - The number of recursions.
 *
 *    \param[out] emax    - max_{0<x<1} |e(x)|
 *
 *    \return  This is an error flag
 *
 *                0  -- success
 *               -1  -- tolerance not met, returns cosines for largest P<=pmax
 *                      that converged
 *                1  -- eta out of range (eta < 1e-7 or eta > 0.1)
 *                2  -- pmax out of range (pmax < 1 or pmax > 40)
 *                3  -- tol is too small (tol < 1e-8)
 *                4  -- Remez/Newton failed to converge
 *                6 -- Linear solve failed
 */
int optimal_cosines(double eta, int pmax, double tol, double *a, int *P, double *emax);

/** \brief  Function to generated the optimal cosines for a given number of terms P. 
 *          This will generate 2P cosines. 
 *
 *  input:
 *   \param[in] eta - a dimensionless parameter given by
 *      
 *                   eta = delta / (cT),
 *
 *          where delta is the seperation of the boundary from sources, 
 *          scatterers, and other inhomogeneities, c is the characterist wave 
 *          speed, and T is the final time. Essentially, this is a measure of
 *          how difficult the problem is to simulate.
 *
 *          We enforce the restriction 1e-7 <= eta <= 0.1.
 *
 *    \param[in] P   - Number of recursions.
 *                     We enforce the resctriction P <= 40
 *    
 *
 *  output:
 *    \param[out] a       - array of 2P cosines in descending order.
 *                          \warning the array 'a' must be allocated to be at 
 *                          least length 2*pmax on calling.
 *
 *    \param[out] emax    - max_{0<x<1} |e(x)|
 *
 *    \return  - This is an error flag
 *
 *                0 -- success
 *                1 -- eta out of range (eta < 1e-7 or eta > 0.1)
 *                2 -- P out of range (P < 1 or P > 40)
 *                4 -- Remez failed to converge
 *                5 -- Newton failed to converge
 *                6 -- Linear solve failed
 */
int optimal_cosinesP(double eta, int P, double *a, double *emax);

# ifdef __cplusplus
}
# endif
#endif // OPTIMAL_COSINES_H
