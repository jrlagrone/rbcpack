/*
  Copyright 2016 John LaGrone

  This file is part of RBCPACK.

  The Yee CRBC Library is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by 
  the Free Software Foundation, either version 3 of the License, or (at your 
  option) any later version.

  The RBCPACK Library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of 
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with the RBCPACK Library.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
   This file implements a few solution routines for Maxwell's equations

*/ 

#include <cmath>
#include "solutions.hpp"

namespace maxwell_solutions {

/*******************************************************************************

                            Free Space Solution

*******************************************************************************/
// constructor
MW_FreeSpace::MW_FreeSpace(const double &gamma, 
                                   const double &tau, 
                                   const double &eps,
                                   const double &mu,
                                   const double src_loc[3])
{

  // copy input values
  this->gamma = gamma;
  this->tau = tau;
  this->eps = eps;
  this->mu = mu;
  for (int i=0; i<3; ++i)
    this->src_loc[i] = src_loc[i];

  // calculate wave speed
  c = 1.0 / std::sqrt(eps*mu);

  // set default derivative method to be exact
  dm = EXACT;

  h_set = false;

}

/*******************************************************************************
                        routines to get solutions
*******************************************************************************/
void MW_FreeSpace::get_solution(const double x[3], const double &t, double sol[6]) const
{
  compute_solution(x,t,sol);
}

double MW_FreeSpace::get_Ex_solution(const double x[3], const double &t) const
{
  return compute_Ex_solution(x,t);
}

double MW_FreeSpace::get_Ey_solution(const double x[3], const double &t) const
{
  return compute_Ey_solution(x,t);
}

double MW_FreeSpace::get_Ez_solution(const double x[3], const double &t) const
{
  return compute_Ez_solution(x,t);
}

double MW_FreeSpace::get_Hx_solution(const double x[3], const double &t) const
{
  return compute_Hx_solution(x,t);
}

double MW_FreeSpace::get_Hy_solution(const double x[3], const double &t) const
{
  return compute_Hy_solution(x,t);
}

double MW_FreeSpace::get_Hz_solution(const double x[3], const double &t) const
{
  return compute_Hz_solution(x,t);
}

/*******************************************************************************
                       routines to calculate solutions
*******************************************************************************/
void MW_FreeSpace::compute_solution(const double x[3], const double &t, double sol[6]) const
{

  double xdist, ydist, zdist, xdist2, ydist2, zdist2;

  using std::sqrt;
  using std::exp;

  // calculate commonly used terms
  xdist = x[0] - src_loc[0];
  xdist2 = xdist * xdist;
  ydist = x[1] - src_loc[1];
  ydist2 = ydist * ydist;
  zdist = x[2] - src_loc[2];
  zdist2 = zdist * zdist;
    
  // calculate derivatives
  if (h_set && (dm == FD_YEE)) { // calculate derivatives using FD curl operator

    double arg, dxdt, dydt, dzdt, dxdx, dydy, \
         dzdz, dxdy[3], dxdz[3], dydz[3], rp, rm, argp, argm, expm, expp, rp3, rm3;
    
    arg = c*(t+tau);

    // There is some apparent repetition here, the reason is that we are only
    // computing the outer curl using the finite difference operator and 
    // computing the remaining expressions exactly. This means that dxdy and 
    // dydx are different because one uses FD for dx and the other for dy.

    // x FD
    rp = sqrt(ydist*ydist + (xdist + h[0]/2.0)*(xdist + h[0]/2.0) + zdist*zdist);
    rm = sqrt(ydist*ydist + (xdist - h[0]/2.0)*(xdist - h[0]/2.0) + zdist*zdist);
    rp3 = rp*rp*rp;
    rm3 = rm*rm*rm;
    argp = arg - rp;
    argm = arg - rm;
    expm = exp(-gamma*argm*argm);
    expp = exp(-gamma*argp*argp);
    dxdt = -2*c*gamma*(argp*expp/rp - argm*expm/rm) / h[0];
    dxdx = ((xdist + h[0]/2.0)*(2*gamma*rp*argp-1.0)*expp/rp3 \
                - (xdist - h[0]/2.0)*(2*gamma*rm*argm-1.0)*expm/rm3) / h[0];
    dxdy[1] = (ydist*(2*gamma*rp*argp-1.0)*expp/rp3 \
                - ydist*(2*gamma*rm*argm-1.0)*expm/rm3) / h[0];
    dxdz[2] = (zdist*(2*gamma*rp*argp-1.0)*expp/rp3 \
                - zdist*(2*gamma*rm*argm-1.0)*expm/rm3) / h[0];

    // y FD
    rp = sqrt(xdist*xdist + (ydist + h[1]/2.0)*(ydist + h[1]/2.0) + zdist*zdist);
    rm = sqrt(xdist*xdist + (ydist - h[1]/2.0)*(ydist - h[1]/2.0) + zdist*zdist);
    rp3 = rp*rp*rp;
    rm3 = rm*rm*rm;
    argp = arg - rp;
    argm = arg - rm;
    expm = exp(-gamma*argm*argm);
    expp = exp(-gamma*argp*argp);
    dydt = -2*c*gamma*(argp*expp/rp - argm*expm/rm) / h[1];
    dydy = ((ydist + h[1]/2.0)*(2*gamma*rp*argp-1.0)*expp/rp3 \
                - (ydist - h[1]/2.0)*(2*gamma*rm*argm-1.0)*expm/rm3) / h[1];
    dxdy[0] = (xdist*(2*gamma*rp*argp-1.0)*expp/rp3 \
                - xdist*(2*gamma*rm*argm-1.0)*expm/rm3) / h[1];
    dydz[2] = (zdist*(2*gamma*rp*argp-1.0)*expp/rp3 \
                - zdist*(2*gamma*rm*argm-1.0)*expm/rm3) / h[1];

    // z FD
    rp = sqrt(xdist*xdist + (zdist + h[2]/2.0)*(zdist + h[2]/2.0) + ydist*ydist);
    rm = sqrt(xdist*xdist + (zdist - h[2]/2.0)*(zdist - h[2]/2.0) + ydist*ydist);
    rp3 = rp*rp*rp;
    rm3 = rm*rm*rm;
    argp = arg - rp;
    argm = arg - rm;
    expm = exp(-gamma*argm*argm);
    expp = exp(-gamma*argp*argp);
    dzdt = -2*c*gamma*(argp*expp/rp - argm*expm/rm) / h[2];
    dzdz = ((zdist + h[2]/2.0)*(2*gamma*rp*argp-1.0)*expp/rp3 \
                - (zdist - h[2]/2.0)*(2*gamma*rm*argm-1.0)*expm/rm3) / h[2];
    dxdz[0] = (xdist*(2*gamma*rp*argp-1.0)*expp/rp3 \
                - xdist*(2*gamma*rm*argm-1.0)*expm/rm3) / h[2];
    dydz[1] = (ydist*(2*gamma*rp*argp-1.0)*expp/rp3 \
                - ydist*(2*gamma*rm*argm-1.0)*expm/rm3) / h[2];

    // put everything together
    // Ex = -mu( dw_3/dydt - dw_2/dzdt)
    sol[0] = -mu*(dydt - dzdt);
    
    // Ey = -mu( dw_1/dzdt - dw_3/dxdt)
    sol[1] = -mu*(dzdt - dxdt);

    // Ez = -mu(dw_2/dxdt - dw_1/dydt)
    sol[2] = -mu*(dxdt - dydt);

    // Hx = dw_2/dx/dy - dw_1/dy/dy + dw_3/dx/dz - dw_1/dzdz
    sol[3] = -(dydy + dzdz) + dxdy[1] + dxdz[2];

    // Hy = dw_3/dydz - dw_2/dz/dz - dw_2/dx/dx + dw_1/dxdy
    sol[4] = -(dzdz + dxdx) + dxdy[0] + dydz[2];

    // Hz = dw_1/dx/dz - dw_3/dxdx + dw_2/dy/dz - dw_3/dydy
    sol[5] = -(dydy + dxdx) + dxdz[0] + dydz[1];

  } else { // calculate exact derivatives

    double arg, r, r2, r3, r5, expon, temp, dxdt, dydt, dzdt, dxdx, dydy, dzdz, dxdy;
    double dxdz, dydz;

    r = sqrt(xdist2 + ydist2 + zdist2);
    r2 = r*r;
    r3 = r2*r;
    r5 = r2*r3;

    arg = c*(t+tau) - r;

    temp = gamma * arg * arg;
    expon = exp(-temp);

    dxdt = -2.0 * c * gamma * expon * (r*(2.0 * temp - 1.0) - arg) / r3;
    dydt = ydist * dxdt;
    dzdt = zdist * dxdt;
    dxdt = xdist * dxdt;

    // 2nd partials
    dxdx = expon * (2.0*arg*gamma*r3-r2) / r5;
    temp = expon * (-6.0*r*arg*gamma + r2*(4.0*gamma * temp - 2.0*gamma) + 3.0) / r5;
    dydy = dxdx + ydist2 * temp;
    dzdz = dxdx + zdist2 * temp;
    dxdx = dxdx + xdist2 * temp;

    // mixed partials
    dxdy = xdist * ydist * temp;
    dxdz = xdist * zdist * temp;
    dydz = ydist * zdist * temp;

    // put everything together
    // Ex = -mu( dw/dydt - dw/dzdt)
    sol[0] = -mu*(dydt - dzdt);
    
    // Ey = -mu( dw/dzdt - dw/dxdt)
    sol[1] = -mu*(dzdt - dxdt);

    // Ez = -mu(dw/dxdt - dw/dydt)
    sol[2] = -mu*(dxdt - dydt);

    // Hx = dw/dx/dy - dw/dy/dy + dw/dx/dz - dw/dzdz
    sol[3] = -dydy - dzdz + dxdy + dxdz;

    // Hy = dw/dydz - dw/dz/dz - dw/dx/dx + dw/dxdy
    sol[4] = -dzdz - dxdx + dxdy + dydz;

    // Hz = dw/dx/dz - dw/dxdx + dw/dy/dz - dw/dydy
    sol[5] = -dydy - dxdx + dxdz + dydz;  

  } // end if
} // end compute_solution

// now do just the computations for individual components
double MW_FreeSpace::compute_Ex_solution(const double x[3], const double &t) const
{

  double xdist, ydist, zdist, xdist2, ydist2, zdist2;

  using std::sqrt;
  using std::exp;

  // calculate commonly used terms
  xdist = x[0] - src_loc[0];
  xdist2 = xdist * xdist;
  ydist = x[1] - src_loc[1];
  ydist2 = ydist * ydist;
  zdist = x[2] - src_loc[2];
  zdist2 = zdist * zdist;
    
  // calculate derivatives
  if (h_set && (dm == FD_YEE)) { // calculate derivatives using FD curl operator

    double arg, dydt, dzdt, rp, rm, argp, argm, expm, expp;
    
    arg = c*(t+tau);

    // derivatives we need of w1
    rp = sqrt(xdist*xdist + (ydist + h[1]/2.0)*(ydist + h[1]/2.0) + zdist*zdist);
    rm = sqrt(xdist*xdist + (ydist - h[1]/2.0)*(ydist - h[1]/2.0) + zdist*zdist);
    argp = arg - rp;
    argm = arg - rm;
    expm = exp(-gamma*argm*argm);
    expp = exp(-gamma*argp*argp);
    dydt = -2*c*gamma*(argp*expp/rp - argm*expm/rm) / h[1];

    rp = sqrt(xdist*xdist + (zdist + h[2]/2.0)*(zdist + h[2]/2.0) + ydist*ydist);
    rm = sqrt(xdist*xdist + (zdist - h[2]/2.0)*(zdist - h[2]/2.0) + ydist*ydist);
    argp = arg - rp;
    argm = arg - rm;
    expm = exp(-gamma*argm*argm);
    expp = exp(-gamma*argp*argp);
    dzdt = -2*c*gamma*(argp*expp/rp - argm*expm/rm) / h[2];

    // put everything together
    // Ex = -mu( dw_3/dydt - dw_2/dzdt)
    return -mu*(dydt - dzdt);
    
  } else { // calculate exact derivatives

    double arg, r, r2, r3, expon, temp, dxdt, dydt, dzdt;

    r = sqrt(xdist2 + ydist2 + zdist2);
    r2 = r*r;
    r3 = r2*r;

    arg = c*(t+tau) - r;

    temp = gamma * arg * arg;
    expon = exp(-temp);

    dxdt = -2.0 * c * gamma * expon * (r*(2.0 * temp - 1.0) - arg) / r3;
    dydt = ydist * dxdt;
    dzdt = zdist * dxdt;

    // put everything together
    // Ex = -mu( dw/dydt - dw/dzdt)
    return -mu*(dydt - dzdt);

  } // end if
} // end compute_Ex_solution

double MW_FreeSpace::compute_Ey_solution(const double x[3], const double &t) const
{

  double xdist, ydist, zdist, xdist2, ydist2, zdist2;

  using std::sqrt;
  using std::exp;

  // calculate commonly used terms
  xdist = x[0] - src_loc[0];
  xdist2 = xdist * xdist;
  ydist = x[1] - src_loc[1];
  ydist2 = ydist * ydist;
  zdist = x[2] - src_loc[2];
  zdist2 = zdist * zdist;
    
  // calculate derivatives
  if (h_set && (dm == FD_YEE)) { // calculate derivatives using FD curl operator

    double arg, dxdt, dzdt, rp, rm, argp, argm, expm, expp;
    
    arg = c*(t+tau);

    // derivatives we need of w1
    rp = sqrt(xdist*xdist + (zdist + h[2]/2.0)*(zdist + h[2]/2.0) + ydist*ydist);
    rm = sqrt(xdist*xdist + (zdist - h[2]/2.0)*(zdist - h[2]/2.0) + ydist*ydist);
    argp = arg - rp;
    argm = arg - rm;
    expm = exp(-gamma*argm*argm);
    expp = exp(-gamma*argp*argp);
    dzdt = -2*c*gamma*(argp*expp/rp - argm*expm/rm) / h[2];


     // derivatives we need of w2 for E
    rp = sqrt(ydist*ydist + (xdist + h[0]/2.0)*(xdist + h[0]/2.0) + zdist*zdist);
    rm = sqrt(ydist*ydist + (xdist - h[0]/2.0)*(xdist - h[0]/2.0) + zdist*zdist);
    argp = arg - rp;
    argm = arg - rm;
    expm = exp(-gamma*argm*argm);
    expp = exp(-gamma*argp*argp);
    dxdt = -2*c*gamma*(argp*expp/rp - argm*expm/rm) / h[0];


    // put everything together
    // Ey = -mu( dw_1/dzdt - dw_3/dxdt)
    return -mu*(dzdt - dxdt);

  } else { // calculate exact derivatives

    double arg, r, r2, r3, expon, temp, dxdt, dzdt;

    r = sqrt(xdist2 + ydist2 + zdist2);
    r2 = r*r;
    r3 = r2*r;

    arg = c*(t+tau) - r;

    temp = gamma * arg * arg;
    expon = exp(-temp);

    dxdt = -2.0 * c * gamma * expon * (r*(2.0 * temp - 1.0) - arg) / r3;
    dzdt = zdist * dxdt;
    dxdt = xdist * dxdt;

    // put everything together
    // Ey = -mu( dw/dzdt - dw/dxdt)
    return -mu*(dzdt - dxdt);


  } // end if
} // end compute_Ey_solution

double MW_FreeSpace::compute_Ez_solution(const double x[3], const double &t) const
{

  double xdist, ydist, zdist, xdist2, ydist2, zdist2;

  using std::sqrt;
  using std::exp;

  // calculate commonly used terms
  xdist = x[0] - src_loc[0];
  xdist2 = xdist * xdist;
  ydist = x[1] - src_loc[1];
  ydist2 = ydist * ydist;
  zdist = x[2] - src_loc[2];
  zdist2 = zdist * zdist;
    
  // calculate derivatives
  if (h_set && (dm == FD_YEE)) { // calculate derivatives using FD curl operator

    double arg, dxdt, dydt, rp, rm, argp, argm, expm, expp;
    
    arg = c*(t+tau);

    // derivatives we need of w1
    rp = sqrt(xdist*xdist + (ydist + h[1]/2.0)*(ydist + h[1]/2.0) + zdist*zdist);
    rm = sqrt(xdist*xdist + (ydist - h[1]/2.0)*(ydist - h[1]/2.0) + zdist*zdist);
    argp = arg - rp;
    argm = arg - rm;
    expm = exp(-gamma*argm*argm);
    expp = exp(-gamma*argp*argp);
    dydt = -2*c*gamma*(argp*expp/rp - argm*expm/rm) / h[1];

     // derivatives we need of w2 for E
    rp = sqrt(ydist*ydist + (xdist + h[0]/2.0)*(xdist + h[0]/2.0) + zdist*zdist);
    rm = sqrt(ydist*ydist + (xdist - h[0]/2.0)*(xdist - h[0]/2.0) + zdist*zdist);
    argp = arg - rp;
    argm = arg - rm;
    expm = exp(-gamma*argm*argm);
    expp = exp(-gamma*argp*argp);
    dxdt = -2*c*gamma*(argp*expp/rp - argm*expm/rm) / h[0];


    // put everything together
    // Ez = -mu(dw_2/dxdt - dw_1/dydt)
    return -mu*(dxdt - dydt);

  } else { // calculate exact derivatives

    double arg, r, r2, r3, expon, temp, dxdt, dydt;

    r = sqrt(xdist2 + ydist2 + zdist2);
    r2 = r*r;
    r3 = r2*r;

    arg = c*(t+tau) - r;

    temp = gamma * arg * arg;
    expon = exp(-temp);

    dxdt = -2.0 * c * gamma * expon * (r*(2.0 * temp - 1.0) - arg) / r3;
    dydt = ydist * dxdt;
    dxdt = xdist * dxdt;

    // Ez = -mu(dw/dxdt - dw/dydt)
    return -mu*(dxdt - dydt);

  } // end if
} // end compute_Ez_solution

double MW_FreeSpace::compute_Hx_solution(const double x[3], const double &t) const
{

  double xdist, ydist, zdist, xdist2, ydist2, zdist2;

  using std::sqrt;
  using std::exp;

  // calculate commonly used terms
  xdist = x[0] - src_loc[0];
  xdist2 = xdist * xdist;
  ydist = x[1] - src_loc[1];
  ydist2 = ydist * ydist;
  zdist = x[2] - src_loc[2];
  zdist2 = zdist * zdist;
    
  // calculate derivatives
  if (h_set && (dm == FD_YEE)) { // calculate derivatives using FD curl operator

    double arg, dydy, dzdz, dxdy[3], dxdz[3], rp, rm, argp, argm, expm, expp, rp3, rm3;
    
    arg = c*(t+tau);

    // derivatives we need of w1
    rp = sqrt(xdist*xdist + (ydist + h[1]/2.0)*(ydist + h[1]/2.0) + zdist*zdist);
    rm = sqrt(xdist*xdist + (ydist - h[1]/2.0)*(ydist - h[1]/2.0) + zdist*zdist);
    rp3 = rp*rp*rp;
    rm3 = rm*rm*rm;
    argp = arg - rp;
    argm = arg - rm;
    expm = exp(-gamma*argm*argm);
    expp = exp(-gamma*argp*argp);
    dydy = ((ydist + h[1]/2.0)*(2*gamma*rp*argp-1.0)*expp/rp3 \
                - (ydist - h[1]/2.0)*(2*gamma*rm*argm-1.0)*expm/rm3) / h[1];

    rp = sqrt(xdist*xdist + (zdist + h[2]/2.0)*(zdist + h[2]/2.0) + ydist*ydist);
    rm = sqrt(xdist*xdist + (zdist - h[2]/2.0)*(zdist - h[2]/2.0) + ydist*ydist);
    rp3 = rp*rp*rp;
    rm3 = rm*rm*rm;
    argp = arg - rp;
    argm = arg - rm;
    expm = exp(-gamma*argm*argm);
    expp = exp(-gamma*argp*argp);
    dzdz = ((zdist + h[2]/2.0)*(2*gamma*rp*argp-1.0)*expp/rp3 \
                - (zdist - h[2]/2.0)*(2*gamma*rm*argm-1.0)*expm/rm3) / h[2];

     // derivatives we need of w2 for E
    rp = sqrt(ydist*ydist + (xdist + h[0]/2.0)*(xdist + h[0]/2.0) + zdist*zdist);
    rm = sqrt(ydist*ydist + (xdist - h[0]/2.0)*(xdist - h[0]/2.0) + zdist*zdist);
    rp3 = rp*rp*rp;
    rm3 = rm*rm*rm;
    argp = arg - rp;
    argm = arg - rm;
    expm = exp(-gamma*argm*argm);
    expp = exp(-gamma*argp*argp);
    dxdy[1] = (ydist*(2*gamma*rp*argp-1.0)*expp/rp3 \
                - ydist*(2*gamma*rm*argm-1.0)*expm/rm3) / h[0];

    // derivatives we need of w3 for E
    rp = sqrt(ydist*ydist + (xdist + h[0]/2.0)*(xdist + h[0]/2.0) + zdist*zdist);
    rm = sqrt(ydist*ydist + (xdist - h[0]/2.0)*(xdist - h[0]/2.0) + zdist*zdist);
    rp3 = rp*rp*rp;
    rm3 = rm*rm*rm;
    argp = arg - rp;
    argm = arg - rm;
    expm = exp(-gamma*argm*argm);
    expp = exp(-gamma*argp*argp);
    dxdz[2] = (zdist*(2*gamma*rp*argp-1.0)*expp/rp3 \
                - zdist*(2*gamma*rm*argm-1.0)*expm/rm3) / h[0];

    // put everything together
    // Hx = dw_2/dx/dy - dw_1/dy/dy + dw_3/dx/dz - dw_1/dzdz
    return -(dydy + dzdz) + dxdy[1] + dxdz[2];

  } else { // calculate exact derivatives

    double arg, r, r2, r3, r5, expon, temp, dxdx, dydy, dzdz, dxdy;
    double dxdz;

    r = sqrt(xdist2 + ydist2 + zdist2);
    r2 = r*r;
    r3 = r2*r;
    r5 = r2*r3;

    arg = c*(t+tau) - r;

    temp = gamma * arg * arg;
    expon = exp(-temp);

    // 2nd partials
    dxdx = expon * (2.0*arg*gamma*r3-r2) / r5;
    temp = expon * (-6.0*r*arg*gamma + r2*(4.0*gamma * temp - 2.0*gamma) + 3.0) / r5;
    dydy = dxdx + ydist2 * temp;
    dzdz = dxdx + zdist2 * temp;

    // mixed partials
    dxdy = xdist * ydist * temp;
    dxdz = xdist * zdist * temp;

    // put everything together
    // Hx = dw/dx/dy - dw/dy/dy + dw/dx/dz - dw/dzdz
    return -dydy - dzdz + dxdy + dxdz;

  } // end if
} // end compute_Hx_solution

double MW_FreeSpace::compute_Hy_solution(const double x[3], const double &t) const
{

  double xdist, ydist, zdist, xdist2, ydist2, zdist2;

  using std::sqrt;
  using std::exp;

  // calculate commonly used terms
  xdist = x[0] - src_loc[0];
  xdist2 = xdist * xdist;
  ydist = x[1] - src_loc[1];
  ydist2 = ydist * ydist;
  zdist = x[2] - src_loc[2];
  zdist2 = zdist * zdist;
    
  // calculate derivatives
  if (h_set && (dm == FD_YEE)) { // calculate derivatives using FD curl operator

    double arg, dxdx, dzdz, dxdy[3], dydz[3], rp, rm, argp, argm, expm, expp, rp3, rm3;
    
    arg = c*(t+tau);

    // There is some apparent repetition here, the reason is that we are only
    // computing the outer curl using the finite difference operator and 
    // computing the remaining expressions exactly. This means that dxdy and 
    // dydx are different because one uses FD for dx and the other for dy.

    // x FD
    rp = sqrt(ydist*ydist + (xdist + h[0]/2.0)*(xdist + h[0]/2.0) + zdist*zdist);
    rm = sqrt(ydist*ydist + (xdist - h[0]/2.0)*(xdist - h[0]/2.0) + zdist*zdist);
    rp3 = rp*rp*rp;
    rm3 = rm*rm*rm;
    argp = arg - rp;
    argm = arg - rm;
    expm = exp(-gamma*argm*argm);
    expp = exp(-gamma*argp*argp);
    dxdx = ((xdist + h[0]/2.0)*(2*gamma*rp*argp-1.0)*expp/rp3 \
                - (xdist - h[0]/2.0)*(2*gamma*rm*argm-1.0)*expm/rm3) / h[0];

    // y FD
    rp = sqrt(xdist*xdist + (ydist + h[1]/2.0)*(ydist + h[1]/2.0) + zdist*zdist);
    rm = sqrt(xdist*xdist + (ydist - h[1]/2.0)*(ydist - h[1]/2.0) + zdist*zdist);
    rp3 = rp*rp*rp;
    rm3 = rm*rm*rm;
    argp = arg - rp;
    argm = arg - rm;
    expm = exp(-gamma*argm*argm);
    expp = exp(-gamma*argp*argp);
    dydz[2] = (zdist*(2*gamma*rp*argp-1.0)*expp/rp3 \
                - zdist*(2*gamma*rm*argm-1.0)*expm/rm3) / h[1];
    dxdy[0] = (xdist*(2*gamma*rp*argp-1.0)*expp/rp3 \
                - xdist*(2*gamma*rm*argm-1.0)*expm/rm3) / h[1];

    // z FD
    rp = sqrt(xdist*xdist + (zdist + h[2]/2.0)*(zdist + h[2]/2.0) + ydist*ydist);
    rm = sqrt(xdist*xdist + (zdist - h[2]/2.0)*(zdist - h[2]/2.0) + ydist*ydist);
    rp3 = rp*rp*rp;
    rm3 = rm*rm*rm;
    argp = arg - rp;
    argm = arg - rm;
    expm = exp(-gamma*argm*argm);
    expp = exp(-gamma*argp*argp);
    dzdz = ((zdist + h[2]/2.0)*(2*gamma*rp*argp-1.0)*expp/rp3 \
                - (zdist - h[2]/2.0)*(2*gamma*rm*argm-1.0)*expm/rm3) / h[2];

    // put everything together
    // Hy = dw_3/dydz - dw_2/dz/dz - dw_2/dx/dx + dw_1/dxdy
    return -(dzdz + dxdx) + dxdy[0] + dydz[2];

  } else { // calculate exact derivatives

    double arg, r, r2, r3, r5, expon, temp, dxdx,  dzdz, dxdy, dydz;

    r = sqrt(xdist2 + ydist2 + zdist2);
    r2 = r*r;
    r3 = r2*r;
    r5 = r2*r3;

    arg = c*(t+tau) - r;

    temp = gamma * arg * arg;
    expon = exp(-temp);

    // 2nd partials
    dxdx = expon * (2.0*arg*gamma*r3-r2) / r5;
    temp = expon * (-6.0*r*arg*gamma + r2*(4.0*gamma * temp - 2.0*gamma) + 3.0) / r5;
    dzdz = dxdx + zdist2 * temp;
    dxdx = dxdx + xdist2 * temp;

    // mixed partials
    dxdy = xdist * ydist * temp;
    dydz = ydist * zdist * temp;

    // put everything together
    // Hy = dw/dydz - dw/dz/dz - dw/dx/dx + dw/dxdy
    return -dzdz - dxdx + dxdy + dydz;

  } // end if
} // end compute_Hy_solution

double MW_FreeSpace::compute_Hz_solution(const double x[3], const double &t) const
{

  double xdist, ydist, zdist, xdist2, ydist2, zdist2;

  using std::sqrt;
  using std::exp;

  // calculate commonly used terms
  xdist = x[0] - src_loc[0];
  xdist2 = xdist * xdist;
  ydist = x[1] - src_loc[1];
  ydist2 = ydist * ydist;
  zdist = x[2] - src_loc[2];
  zdist2 = zdist * zdist;
    
  // calculate derivatives
  if (h_set && (dm == FD_YEE)) { // calculate derivatives using FD curl operator

    double arg, dxdx, dydy, dxdz[3], dydz[3], rp, rm, argp, argm, expm, expp, rp3, rm3;
    
    arg = c*(t+tau);

    // There is some apparent repetition here, the reason is that we are only
    // computing the outer curl using the finite difference operator and 
    // computing the remaining expressions exactly. This means that dxdy and 
    // dydx are different because one uses FD for dx and the other for dy.

    // x FD
    rp = sqrt(ydist*ydist + (xdist + h[0]/2.0)*(xdist + h[0]/2.0) + zdist*zdist);
    rm = sqrt(ydist*ydist + (xdist - h[0]/2.0)*(xdist - h[0]/2.0) + zdist*zdist);
    rp3 = rp*rp*rp;
    rm3 = rm*rm*rm;
    argp = arg - rp;
    argm = arg - rm;
    expm = exp(-gamma*argm*argm);
    expp = exp(-gamma*argp*argp);
    dxdx = ((xdist + h[0]/2.0)*(2*gamma*rp*argp-1.0)*expp/rp3 \
                - (xdist - h[0]/2.0)*(2*gamma*rm*argm-1.0)*expm/rm3) / h[0];

    // y FD
    rp = sqrt(xdist*xdist + (ydist + h[1]/2.0)*(ydist + h[1]/2.0) + zdist*zdist);
    rm = sqrt(xdist*xdist + (ydist - h[1]/2.0)*(ydist - h[1]/2.0) + zdist*zdist);
    rp3 = rp*rp*rp;
    rm3 = rm*rm*rm;
    argp = arg - rp;
    argm = arg - rm;
    expm = exp(-gamma*argm*argm);
    expp = exp(-gamma*argp*argp);
    dydy = ((ydist + h[1]/2.0)*(2*gamma*rp*argp-1.0)*expp/rp3 \
                - (ydist - h[1]/2.0)*(2*gamma*rm*argm-1.0)*expm/rm3) / h[1];

    // z FD
    rp = sqrt(xdist*xdist + (zdist + h[2]/2.0)*(zdist + h[2]/2.0) + ydist*ydist);
    rm = sqrt(xdist*xdist + (zdist - h[2]/2.0)*(zdist - h[2]/2.0) + ydist*ydist);
    rp3 = rp*rp*rp;
    rm3 = rm*rm*rm;
    argp = arg - rp;
    argm = arg - rm;
    expm = exp(-gamma*argm*argm);
    expp = exp(-gamma*argp*argp);
    dxdz[0] = (xdist*(2*gamma*rp*argp-1.0)*expp/rp3 \
                - xdist*(2*gamma*rm*argm-1.0)*expm/rm3) / h[2];
    dydz[1] = (ydist*(2*gamma*rp*argp-1.0)*expp/rp3 \
                - ydist*(2*gamma*rm*argm-1.0)*expm/rm3) / h[2];

    // put everything together
    // Hz = dw_1/dx/dz - dw_3/dxdx + dw_2/dy/dz - dw_3/dydy
    return -(dydy + dxdx) + dxdz[0] + dydz[1];

  } else { // calculate exact derivatives

    double arg, r, r2, r3, r5, expon, temp, dxdx, dydy, dxdz, dydz;

    r = sqrt(xdist2 + ydist2 + zdist2);
    r2 = r*r;
    r3 = r2*r;
    r5 = r2*r3;

    arg = c*(t+tau) - r;

    temp = gamma * arg * arg;
    expon = exp(-temp);

    // 2nd partials
    dxdx = expon * (2.0*arg*gamma*r3-r2) / r5;
    temp = expon * (-6.0*r*arg*gamma + r2*(4.0*gamma * temp - 2.0*gamma) + 3.0) / r5;
    dydy = dxdx + ydist2 * temp;
    dxdx = dxdx + xdist2 * temp;

    // mixed partials
    dxdz = xdist * zdist * temp;
    dydz = ydist * zdist * temp;

    // put everything together
    // Hz = dw/dx/dz - dw/dxdx + dw/dy/dz - dw/dydy
    return -dydy - dxdx + dxdz + dydz;  

  } // end if
} // end compute_Hz_solution

} // namespace
