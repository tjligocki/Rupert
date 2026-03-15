!*****************************************************************************

!  MPFUN-MPFR: An MPFR-based arbitrary precision computation package
!  Language interface module (module MPFUNG)
!  Variant 1: Precision level specifications are *optional*; no real*16 support.
!  Search for !> for variant differences.

!  Revision date:  27 Apr 2016

!  AUTHOR:
!     David H. Bailey
!     Lawrence Berkeley National Lab (retired) and University of California, Davis
!     Email: dhbailey@lbl.gov

!  COPYRIGHT AND DISCLAIMER:
!  All software in this package (c) 2016 David H. Bailey.
!  By downloading or using this software you agree to the copyright, disclaimer
!  and license agreement in the accompanying file DISCLAIMER.txt.

!  PURPOSE OF PACKAGE:
!    This package permits one to perform floating-point computations (real and
!    complex) to arbitrarily high numeric precision, by making only relatively
!    minor changes to existing Fortran-90 programs.  All basic arithmetic
!    operations and transcendental functions are supported, together with several
!    special functions.

!    This version differs from the MPFUN-Fort version by the same author in that
!    it is based on MPFR, which presently is the fastest available low-level
!    package for high-precision floating-point computation.  Thus most user
!    applications typically run 3X faster.  In addition, the developers of the
!    MPFR package have taken considerable pains to ensure that the many functions
!    return correctly rounded (to the last bit) results for each input.  At the
!    Fortran user level, application codes written for MPFUN-Fort may be compiled
!    and executed with MPFUN-MPFR --i.e., MPFUN-MPFR is "plug compatible" with
!    MPFUN-Fort.

!  DOCUMENTATION:
!    A detailed description of this package, and instructions for compiling
!    and testing this program on various specific systems are included in the
!    README file accompanying this package, and, in more detail, in the
!    following technical paper:

!    David H. Bailey, "MPFUN2015: A thread-safe arbitrary precision package," 
!    available at http://www.davidhbailey.com/dhbpapers/mpfun2015.pdf.
  
!  DESCRIPTION OF THIS MODULE (MPFUNG):
!    This module contains all high-level Fortran-90 language interfaces.
!    There are two variants of this module, which take different approaches
!    to dynamically changing working precision within an application.  Variant
!    1 allows mixed-mode assignments, and does not require precision level
!    specifications in certain functions, whereas variant 2 does not permit
!    mixed-mode assignments, and requires precision level specifications in
!    certain functions.  See documentation for details.
 
module mptjl
use mpfuna
use mpfunf
use mpfung

implicit none
  
!  Algebraic, transcendental and type conversion functions:
  
private &
  mp_fixlocr,mp_fform

!  Operator extension interface blocks:

!  MP generic function interface blogs, listed alphabetically by interface name:

interface mpfform
  module procedure mp_fform
end interface

contains
  subroutine mp_fixlocr (ra)
    implicit none
    type (mp_real):: ra
    ra%mpr(4) = loc (ra%mpr(4)) + 8
    return
  end subroutine

  subroutine mp_fform (ra, nb, nd, b)
    implicit none
    type (mp_real), intent (in):: ra
    integer, intent (in):: nb, nd
    character(1), intent (out):: b(nb)
    character(1) b1(nb+8)
    integer ix, i, j
    integer(8) iexp
    character(16) str16

    call mp_fixlocr (ra)

!  Call mpfrgetstr to convert numbmer.

    call mpfrgetstr (ra%mpr(1), b1, nb, iexp, mprnd)

!  Check for overflow of field.
    if ( &
      b1(1) /= '-' .and. iexp > 0 .and. nb < iexp + nd + 1 .or. &
      b1(1) /= '-' .and. iexp <= 0 .and. nb < nd + 2 .or. &
      b1(1) == '-' .and. iexp > 0 .and. nb < iexp + nd + 2 .or. &
      b1(1) == '-' .and. iexp <= 0 .and. nb < nd + 3) then
      do i = 1, nb
        b(i) = '*'
      enddo
    endif        

    if (b1(1) == '0' .or. iexp <= - nd) then

!  Output is zero.

      do i = 1, nb - nd - 2
        b(i) = ' '
      enddo
      ix = nb - nd - 2
      b(ix+1) = '0'
      b(ix+2) = '.'
      ix = ix + 2
      do i = 1, nd
        b(ix+i) = '0'
      enddo
    elseif (b1(1) /= '-' .and. iexp > 0) then

!  Value is positive and exponent is positive.

      do i = 1, nb - iexp - nd - 1
        b(i) = ' '
      enddo
      ix = nb - iexp - nd - 1
      do i = 1, iexp
        b(ix+i) = b1(i)
      enddo
      ix = ix + iexp
      b(ix+1) = '.'
      ix = ix + 1
      do i = 1, nd
        b(ix+i) = b1(i+iexp)
      enddo
    elseif (b1(1) /= '-' .and. iexp <= 0) then

!  Value is positive and exponent is negative or zero.

      do i = 1, nb - nd - 2
        b(i) = ' '
      enddo
      ix = nb - nd - 2
      b(ix+1) = '0'
      b(ix+2) = '.'
      ix = ix + 2
      do i = 1, abs (iexp)
        b(ix+i) = '0'
      enddo
      ix = ix + abs (iexp)
      do i = 1, nd - abs (iexp)
        b(ix+i) = b1(i)
      enddo
    elseif (b1(1) == '-' .and. iexp > 0) then

!  Value is negative and exponent is positive.

      do i = 1, nb - iexp - nd - 2
        b(i) = ' '
      enddo
      ix = nb - iexp - nd - 2
      do i = 1, iexp + 1
        b(ix+i) = b1(i)
      enddo
      ix = ix + iexp + 1
      b(ix+1) = '.'
      ix = ix + 1
      do i = 1, nd
        b(ix+i) = b1(i+iexp+1)
      enddo
    elseif (b1(1) /= '-' .and. iexp <= 0) then

!  Value is negative and exponent is negative or zero.

      do i = 1, nb - nd - 3
        b(i) = ' '
      enddo
      ix = nb - nd - 3
      b(ix+1) = '-'
      b(ix+2) = '0'
      b(ix+3) = '.'
      ix = ix + 3
      do i = 1, abs (iexp)
        b(ix+i) = '0'
      enddo
      ix = ix + abs (iexp)
      do i = 1, nd - abs (iexp)
        b(ix+i) = b1(i)
      enddo
    endif
    return
  end subroutine

end module mptjl

