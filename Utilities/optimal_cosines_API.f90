! This module implements a fortran interface to the optimal cosines routines
!
!  Last Updated:
!  Sept 16, 2014
!  John LaGrone
!  Southern Methodist University

module optimal_cosines_module

  use, intrinsic :: ISO_C_Binding, only: C_INT, C_DOUBLE, C_PTR
  implicit none
  private

  interface 
    ! declarations for the C functions from optimal_cosines
    function c_optimal_cosines(eta, pmax, tol, a, P, emax) result(flag) bind(C, name="optimal_cosines")
      import
      REAL (C_DOUBLE), value :: eta, tol
      INTEGER (C_INT), value :: pmax
      INTEGER (C_INT) :: P
      INTEGER (C_INT) :: flag
      REAL (C_DOUBLE) :: emax
      REAL (C_DOUBLE), dimension(*) :: a
    end function c_optimal_cosines

    function c_optimal_cosinesP(eta, P, a, emax) result(flag) bind(C, name="optimal_cosinesP")
      import
      REAL (C_DOUBLE), value :: eta
      INTEGER (C_INT), value :: P
      INTEGER (C_INT) :: flag
      REAL (C_DOUBLE) :: emax
      REAL (C_DOUBLE), dimension(*) :: a
    end function c_optimal_cosinesP

  end interface

  ! declare the fortran interfaces
  interface get_optimal_cosines
    module procedure get_optimal_cosines
  end interface get_optimal_cosines

  interface get_optimal_cosinesP
    module procedure get_optimal_cosinesP
  end interface get_optimal_cosinesP

  public get_optimal_cosines, get_optimal_cosinesP

  ! now actually give the definitions for the interface
  contains

    subroutine get_optimal_cosines(eta, pmax, tol, a, P, emax, flag)
      double precision, intent(in) :: eta, tol
      integer, intent(in) :: pmax
      double precision, dimension(*), intent(out) :: a
      integer, intent(out) :: P, flag
      double precision, intent(out) :: emax

      flag = c_optimal_cosines(eta, pmax, tol, a, P, emax)

    end subroutine get_optimal_cosines

    subroutine get_optimal_cosinesP(eta, P, a, emax, flag)
      double precision, intent(in) :: eta
      double precision, dimension(*), intent(out) :: a
      integer, intent(in) :: P
      double precision, intent(out) :: emax
      integer, intent(out) :: flag


      flag = c_optimal_cosinesP(eta, P, a, emax)

    end subroutine get_optimal_cosinesP


end module optimal_cosines_module
