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

/*
  This file provides implementations for the functions declared in 
  optimal_cosines.h. Refer to the header file for a description of the two 
  main functions and their respective usage.
*/

// optimal cosine function declarations
#include "optimal_cosines.h"

// math operations (fabs)
#include <math.h>

// malloc, etc.
#include <stdlib.h>

// prototype for LAPACK solve
/* we're not using this for portability reasons
extern void dgesv_(
      int *N,         // order of matrix
      int *NRHS,      // number of right hand sides
      double *A,      // Matrix A 
      int *LDA,       // leading dim of A
      int *IPIV,      // array for pivoting indices (dim = N)
      double *B,      // RHS / solution
      int *LDB,       // leading dim of B  
      int *INFO       // Success flag
      );  
*/

// declarations of some helpful routines

// function we want to minimize
void f(double eta, double *a, int p, double *x, int nx, double *sol);

// derivative of function we want to minimize
void df(double eta, double *a, int p, double *x, int nx, double *sol);

// routine to find zeros (of df) using Brent's Method
double find_zeros(double eta, double *a, int p, double xl, double xr, double tol);

// routine to compute L-infinity norm of a vector
void linfty(double *v, int n, double *sol);

// routine to find the minimum abs value and corresponding index in a vector
void min(double *v, int n, double *minimum, int *ind);

// routine to sort a vector into descending order
void descend_sort(double *v, int n);

// routine to solve Ax=b
int solve(double *A, double *b, int n);

/* Function to try to find the minimal number of cosines needed to meet a 
   specified tolerance;
*/
int optimal_cosines(double eta, int pmax, double tol, double *a, int *p, double *emax)
{

  int flag, p_accept;
  
  p_accept = -1;

  // check inputs
  if ((eta < 1e-7) || (eta > 0.1))
    return 1;
  if ((pmax < 1) || (pmax > 40))
    return 2;
  if (tol < 1e-8)
    return 3;

  *emax = tol + 1;

  // increment P until we get something that works
  for ((*p)=1; (*p) <= pmax; (*p)++) {

    flag = optimal_cosinesP(eta, (*p), a, emax);
    
    if (flag == 0) 
      p_accept = *p;

    if ((*emax) < tol)
      return 0;
  }

  // return the largest value of p that converged
  if (p_accept != -1) {
    *p = p_accept;
    flag = optimal_cosinesP(eta, (*p), a, emax);
    return -1;
  } else {
    return 4; // optimal_cosinesP failed on all attempts
  }

}

/* Function to generated the optimal cosines for a given number of terms P. This
   will generate 2P cosines. 
 
   This is the Remez Alg.
*/
int optimal_cosinesP(double eta, int p, double *a, double *emax)
{

  int i, j, remez_iter, newton_iter, remez_max, newton_max, sgn, *sgns, n;
  double xmin, remez_tol, newton_tol, lambda, damping, *jacobian, *residual, \
         *midpoints, *xm, *fm, *a0, atry, conv_check;

  // check inputs
  if ((eta < 1e-7) || (eta > 0.1))
    return 1;
  if ((p < 1) || (p > 40))
    return 2;

  // size of matrix
  n = 2*p+1;

  // set tolerances
  remez_tol = 1e-12;
  newton_tol = 1e-10;

  // set max iterations
  remez_max = 100;
  newton_max = 60;

  // for x<xmin, |e| < rtol so we don't need to check smaller values
  xmin = eta / log(1.0 / remez_tol);

  // allocate memory
  a0 = malloc(2*p*sizeof(double));
  jacobian = malloc(n*n*sizeof(double));
  sgns = malloc(n*sizeof(int));
  residual = malloc(n*sizeof(double));
  midpoints = malloc(n*sizeof(double));
  xm = malloc(n*sizeof(double));
  fm = malloc(n*sizeof(double));

  // fill the jacobian matrix with 0
  for (i=0; i<n*n; i++)
    jacobian[i] = 0.0;

  // create an initial guess using a geometric distribution
  for (i=0; i<2*p; i++) 
    a[i] = pow(eta, (i+1.0)/(2.0*p));

  // Start Remez
  for (remez_iter=0; remez_iter<remez_max; remez_iter++) {

    // find locations and signs of extrema
    sgn = 1;
    xm[0] = find_zeros(eta, a, p, a[0], 1.0, 1e-14);
    sgns[0] = 1;
    sgn = -sgn;
    for (i=1; i<2*p; i++) {
      xm[i] = find_zeros(eta, a, p, a[i], a[i-1], 1e-14);
      sgns[i] = sgn;
      sgn = -sgn;
    }
    xm[2*p] = find_zeros(eta, a, p, xmin, a[2*p-1], 1e-14);
    sgns[2*p] = sgn;
    sgn = -sgn;

    // get the function values at the extrema and find the minimum
    f(eta, a, p, xm, n, fm);
    min(fm, n, &lambda, &i);

    // save the values for a so we can check convergence
    for (i=0; i<2*p; i++)
      a0[i] = a[i];

    // Now we're going to try to get a better approximation by solving a
    // nonlinear system with the constraint that the extrema of the error
    // function have the same magnitude and alternate signs. We use a damped
    // Newton's method to do this.
    for (newton_iter = 0; newton_iter < newton_max; newton_iter++) {

      // calculate the residual
      f(eta, a, p, xm, n, fm);
      for (i=0; i<n; i++)
        residual[i] = lambda*sgns[i] - fm[i];

      // build the jacobian matrix
      // upper left 2p x 2p block
      for (i=0; i<2*p; i++)  // loop over columns
        for (j=0; j<n; j++)  // loop over rows
          jacobian[j+i*n] = -2.0*fm[j]*xm[j]/((xm[j]-a[i])*(xm[j]+a[i]));

      // last column
      for (j=0; j<n; j++)  // loop over rows
        jacobian[j+(2*p)*n] = -sgns[j];

      // apply the inverse of the jacobian to the residual
      // dgesv_(&n, &nrhs, jacobian, &n, ipiv, residual, &n, &info);
      if (0 != solve(jacobian, residual, n))
        return 6; // solve failed

      // calculate a damping parameter to ensure that we don't move a control
      // point too far.
      damping = 1.0;
     
      // calculate midpoints
      midpoints[0] = (1.0 + a[0]) / 2.0;
      for (i=1; i<2*p; i++)
        midpoints[i] = (a[i-1] + a[i]) / 2.0;
      midpoints[n-1] = (xmin + a[2*p-1]) / 2.0;
      
      // see if the new values for a will work
      for (i=0; i<2*p; i++) {
        atry = a[i] + damping*residual[i];
        damping = (atry > midpoints[i]) ? \
                  (midpoints[i]-a[i])/residual[i] : damping;
        damping = (atry < midpoints[i+1]) ? \
                  (midpoints[i+1]-a[i])/residual[i] : damping;
      }
      damping = ((lambda+damping*residual[i]) < lambda/2.0) ? \
                -lambda / (2.0*residual[i]) : damping;

      // caluclate new values
      for (i=0; i<n; i++) 
        residual[i] = damping*residual[i];
      for (i=0; i<2*p; i++)
        a[i] += residual[i];
      lambda += residual[i];

      // check for convergence
      linfty(residual, n, &conv_check);
      if (conv_check < newton_tol)
        break;

      if (newton_iter == newton_max -1) {

        free(jacobian);
        free(a0);
        free(sgns);
        free(residual);
        free(midpoints);
        free(xm);
        free(fm);
        return 5; // Newton failed to converge
      } 
        
    } // end newton
    
    // sort in desceding order
    descend_sort(a, 2*p);
    
    // calculate difference of iterates
    for (i=0; i<2*p; i++)
      a0[i] -= a[i];

    // calculate the maximum
    linfty(a0, 2*p, &conv_check);
    if (conv_check < remez_tol)
        break;      

    if (remez_iter == remez_max -1) {

      free(jacobian);
      free(a0);
      free(sgns);
      free(residual);
      free(midpoints);
      free(xm);
      free(fm);
      return 4; // Remez failed to converge
    }

  } // end Remez

  // finally calculate an error approximation
  xm[0] = find_zeros(eta, a, p, a[0], 1.0, 1e-14);
  for (i=1; i<2*p; i++) {
    xm[i] = find_zeros(eta, a, p, a[i], a[i-1], 1e-14);
  }
  xm[2*p] = find_zeros(eta, a, p, xmin, a[2*p-1], 1e-14);

  // get the function values at the extrema and find the max
  f(eta, a, p, xm, n, fm);
  linfty(fm, n, emax);
  
  // free memory
  free(jacobian);
  free(a0);
  free(sgns);
  free(residual);
  free(midpoints);
  free(xm);
  free(fm);

  return 0;

}


/* ----------------------------------------------------------------------------- 
                  Definitions of the "useful" routines
------------------------------------------------------------------------------*/

void f(double eta, double *a, int p, double *x, int nx, double *sol)
{

  int i,k;

  for (i=0; i<nx; i++) {
    if (fabs(x[i]) < 1e-14) {
      sol[i] = 0.0;
    } else {
      sol[i] = exp(-eta/x[i])*(1.0-x[i])/(1.0+x[i]);
      for (k=0; k<2*p; k++)
        sol[i] *= (x[i]-a[k])/(x[i]+a[k]);
    }
  }
}

/*
  derivative of the function we want to minimize. We will use this to find local
  extrema in the Remez algorithm
*/               
void df(double eta, double *a, int p, double *x, int nx, double *sol)
{

  int i, j, k, kmin;
  double min_val, w, dw, *work_vec;

  work_vec = malloc(2*p*sizeof(double));

  // zero out solution vector
  for (j=0; j<nx; j++)
    sol[j] = 0.0;

  for (i=0; i<nx; i++) {

    // compute (a-x[i])
    for (j=0; j<2*p; j++)
      work_vec[j] = a[j] - x[i];

    // find closest value to 0 and its index
    min(work_vec, 2*p, &min_val, &kmin);

    // compute derivative
    if (fabs(x[i]) < 1e-12) { // left endpoint
     
      sol[i] = -1e-50;

    } else if (fabs(x[i] - 1.0) < 1e-12) { // right endpoint

      sol[i] = -exp(-eta) / 2.0;
      for (k=0; k<2*p; k++)
        sol[i] *= (1.0 - a[k]) / (1.0 + a[k]);

    } else if (min_val < 1e-12) { // deal with a possible div. by 0

      sol[i] = exp(-eta/x[i])*(1.0-x[i])/((1.0+x[i])*(x[i] + a[kmin]));
      for (k=0; k<2*p; k++)
        sol[i] *= (k != kmin) ? (x[i] - a[k]) / (x[i] + a[k]) : 1.0;

    } else { // rest of the interval

      w = exp(-eta/x[i])*(1.0-x[i])/(1.0+x[i]);
      dw = w*(eta/(x[i]*x[i]) - 2.0/(1.0-x[i]*x[i]));
      for (k=0; k<2*p; k++) {
        dw *= (x[i]-a[k])/(x[i]+a[k]);
        w *= (x[i]-a[k])/(x[i]+a[k]);
      }
      for (k=0; k<2*p; k++)
        dw += 2.0*w*a[k]/(x[i]*x[i] - a[k]*a[k]);
      
      sol[i] = dw;
    }
  }

  free(work_vec);
}

/* 
  This defines a function for Brent's method. This is a root finding
  algorithm that is a combination of the bisection method, secant method, and
  inverse quadratic interpolation.

  The basic idea is to try to get efficiency out of the inverse quadratic
  interpolation and secant method while retaining the robustness of the 
  bisection method.

  This is only implemented for scalar valued functions with a single input.
  Additionally, we require that the solution be bracketed.

  Wikipedia provides a rather good description of the method. As of Aug 15, 2014
  this was available at (http://en.wikipedia.org/wiki/Brent's_method)

  A description is also available in the following:

  Brent, R.P. 1973, Algorithms for Minimization without Derivatives 
   (Englewood Cliffs, NJ: Prentice-Hall), Chapter 5.

  Numerical Recipes in C, The Art of Scientific Computing, Second Edition, 
    William H. Press, Saul A. Teukolsky, William T. Vetterling, and 
    Brian P. Flannery. Cambridge University Press. 1988, 1992. 
*/
double find_zeros(double eta, double *aloc, int p, double xl, double xr, double tol)
{

  double a, b, c, d, s, fa, fb, fs, fc, tmp, x;
  int flag;
  int i, max_it = 2000;

  // save inputs
  a = xl;
  b = xr;

  // get function value at endpoints
  df(eta, aloc, p, &a, 1, &fa);
  df(eta, aloc, p, &b, 1, &fb);

  // make sure the zero is bracketed
  if (fa * fb >= 0.0) {
    return -1;
  }

  // make the side farther from 0 be on the left
  if (fabs(fa) < fabs(fb)) {
    tmp = a;
    a = b;
    b = tmp;
    tmp = fa;
    fa = fb;
    fb = tmp;
  }

  // initialize
  c = a;
  fc = fa;
  flag = 1;
  d = 0;

  // start iterating
  for (i=0; i<max_it; i++) {

    // do inverse quadratic interpolation if we can
    if ((fabs(fa-fc) > 1e-12) && (fabs(fb-fc) > 1e-12)) {
 
      s = (a*fb*fc) / ((fa-fb)*(fa-fc)) \
        + (b*fa*fc) / ((fb-fa)*(fb-fc)) \
        + (c*fa*fb) / ((fc-fa)*(fc-fb));
   
    // otherwise use the secant method
    } else {

      s = b - fb*(b-a) / (fb - fa);    

    }

    // make sure that the value caluclated by either the secant or inv. quad int.
    // is acceptable. Otherwise do bisection, which will always work.
    // Consult a reference for the explanation of the conditions ...
    if ( ((s < (3.*a+b)/4.0) || (s > b))
       || ((flag) && (fabs(s-b) >= fabs(b-c)/2.0))
       || ((!flag) && (fabs(s-b) >= fabs(c-d)/2.0))
       || ((flag) && (fabs(b-c) < tol))
       || ((flag) && (fabs(c-d) < tol)))
    {
     
      // bisection method
      s = (a+b) / 2.0;
      flag = 1;

    } else {
   
      flag = 0;

    }

    // update values
    df(eta, aloc, p, &s, 1, &fs);
    d = c;
    c = b;
    fc = fb;
 
    if (fa*fs < 0.0) {
      b = s;
      fb = fs;
    } else {
      a = s;
      fa = fs;
    }

    // make the side farther from 0 be on the left
    if (fabs(fa) < fabs(fb)) {
      tmp = a;
      a = b;
      b = tmp;
      tmp = fa;
      fa = fb;
      fb = tmp;
    }

    if (fabs(fa) < tol) {
      x = a;
      return x;
    }

    if (fabs(fb) < tol) {
      x = b;
      return x;
    }

    if (fabs(a-b) < tol) {
      x = (a + b) / 2.0;
      return x;
    }
  
  } // end iteration loop

  // if we get here, we reached the max iterations, so return failure
  return -1;
}

// routine to compute L-infinity norm of a vector
void linfty(double *v, int n, double *sol) 
{
  int i;
  double test_max;

  test_max = fabs(v[0]);

  for (i=1; i<n; i++)
    test_max = (fabs(v[i]) > test_max) ? fabs(v[i]) : test_max;

  *sol = test_max;

}

// routine to find the minimum abs value and corresponding index in a vector
void min(double *v, int n, double *minimum, int *ind)
{
  int i, test_ind;
  double test_min;

  test_min = fabs(v[0]);
  test_ind = 0;

  for (i=1; i<n; i++) {
    if (fabs(v[i]) < test_min) {
      test_min = fabs(v[i]);
      test_ind = i;
    }
  }

  *minimum = test_min;
  *ind = test_ind;
}

// routine to sort a vector into descending order
void descend_sort(double *v, int n)
{
  int i, j;
  double max;

  for (i=0; i<n; i++) {
    max = v[i];
    for (j=i+1; j<n; j++) {
      if (v[j] > max) {
        max = v[j];
        v[j] = v[i];
        v[i] = max;
      }
    }
  }
}

// routine to solve Ax = b
int solve(double *A, double *b, int n)
{

  int *p, i, j, k, kp, temp;
  double *s, *y, r, t;

  // initialize pivoting vector
  p = malloc(n*sizeof(int));
  for (i=0; i<n; i++)
    p[i] = i;

  // initilize vector to store row scales
  s = malloc(n*sizeof(double));

  // initilize work vector
  y = malloc(n*sizeof(double));

  // compute the row scales
  for (j=0; j<n; j++) { 
    s[j] = fabs(A[j]);
    for (k=1; k<n; k++) { 
      s[j] = (s[j] > fabs(A[j + k*n])) ? s[j] : fabs(A[j + k*n]);
    }
  }

  // gaussian elimination with scaled partial pivotting. Overwrites A.
  for (k=0; k<n-1; k++) {

    r = fabs(A[p[k]+ k*n] / s[p[k]]);
    kp = k;

    for (i=k+1; i<n; i++) {
       
      t = fabs(A[p[i] + k*n] / s[p[i]]);
      if (t > r) {
        r = t;
        kp = i;
      }
    }
    temp = p[kp];
    p[kp] = p[k];
    p[k] = temp;
    for (i=k+1; i<n; i++) {
      A[p[i]+ k*n] = A[p[i] + k*n] / A[p[k]+ k*n];
      for (j=k+1; j<n; j++) {
        A[p[i]+ j*n] = A[p[i]+ j*n] - A[p[i]+ k*n]*A[ p[k]+j*n];
      } 
    }
  }

  // forward substitution to solve Ly = b. 
  for (i=0; i<n; i++)
    y[i] = 0.0;

  y[0] = b[p[0]];
  for (i=1; i<n; i++) {
    y[i] = b[p[i]];
    for (j=0; j<i; j++) {
      y[i] = y[i] - A[p[i]+ j*n]*y[j];
    }
  }

  // back subsitution to solve Ux = y.
  b[n-1] = y[n-1]/A[p[n-1] + (n-1)*n];
  for (i=n-2; i>=0; i--) {
    b[i] = y[i];
    for (j=i+1; j<n; j++) {
      b[i] = b[i] - A[p[i] + j*n]*b[j];
    }
    b[i] = b[i]/A[p[i]+ i*n];
  }
  
  free(s);
  free(p);
  free(y);

  return 0;

}
