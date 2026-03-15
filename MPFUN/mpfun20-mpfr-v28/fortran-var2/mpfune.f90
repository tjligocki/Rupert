!*****************************************************************************

!  MPFUN20-MPFR  A thread-safe arbitrary precision computation package
!  Special functions module (module MPFUNE)

!  Revision date: 2 Feb 2023

!  AUTHOR:
!    David H. Bailey
!    Lawrence Berkeley National Lab (retired)
!    Email: dhbailey@lbl.gov

!  COPYRIGHT AND DISCLAIMER:
!    All software in this package (c) 2023 David H. Bailey.
!    By downloading or using this software you agree to the copyright, disclaimer
!    and license agreement in the accompanying file DISCLAIMER.txt.

!  PURPOSE OF PACKAGE:
!    This package permits one to perform floating-point computations (real and
!    complex) to arbitrarily high numeric precision, by making only relatively
!    minor changes to existing Fortran-90 programs. All basic arithmetic
!    operations and transcendental functions are supported, together with numerous
!    special functions.

!    In addition to fast execution times, one key feature of this package is a
!    100% THREAD-SAFE design, which means that user-level applications can be
!    easily converted for parallel execution, say using a threaded parallel
!    environment such as OpenMP.

!  DOCUMENTATION:
!    A detailed description of this package, and instructions for compiling
!    and testing this program on various specific systems are included in the
!    README file accompanying this package, and, in more detail, in the
!    following technical paper:

!    David H. Bailey, "MPFUN2020: A new thread-safe arbitrary precision package,"
!    available at http://www.davidhbailey.com/dhbpapers/mpfun2020.pdf.

!  DESCRIPTION OF THIS MODULE (MPFUNE):
!    This module contains subroutines to perform special functions. Additional
!    functions will be added as they are completed.

!  NOTE ON PROGRAMMING CONVENTION FOR THIS MODULE:
!    This module is designed to facilitate easy translation (using the program
!    convmpfr.f90), for use in the MPFUN-MPFR package, and also for translation
!    (using the program convmpdq.f90), for use in the DQFUN double-quad package.

module mpfune
use mpfuna
! use mpfunb
! use mpfunc
! use mpfund

contains

!  These routines perform simple operations on the MP data structure. Those
!  listed between !> and !>> are for MPFUN20-Fort; those between !>> and !>>>
!  are for MPFUN20-MPFR. Those between !>>> and !>>>> are for DQFUN. The
!  translator selects the proper set.
!>

! subroutine mpinitwds (ra, mpnw)
! implicit none
! integer (mpiknd), intent(out):: ra(0:)
! integer, intent(in):: mpnw

! ra(0) = mpnw + 6
! ra(1) = mpnw
! ra(2) = 0; ra(3) = 0; ra(4) = 0
! return
! end subroutine mpinitwds

! integer function mpwprecr (ra)
! implicit none
! integer (mpiknd), intent(in):: ra(0:)

! mpwprecr = ra(1)
! return
! end function mpwprecr

! integer function mpspacer (ra)
! implicit none
! integer (mpiknd), intent(in):: ra(0:)

! mpspacer = ra(0)
! return
! end function mpspacer

!>>

subroutine mpabrt (ier)
implicit none
integer, intent(in):: ier
write (mpldb, 1) ier
1 format ('*** MPABRT: Execution terminated, error code =',i4)
stop
end subroutine mpabrt

subroutine mpinitwds (ra, mpnw)
implicit none
integer (mpiknd), intent(out):: ra(0:)
integer, intent(in):: mpnw

ra(0) = mpnw + 6
ra(1) = mpnw * mpnbt
ra(2) = 1
ra(3) = mpnan
ra(4) = loc (ra(4)) + 8
ra(mpnw+5) = 0
return
end subroutine mpinitwds

subroutine mpfixlocr (ia)
implicit none
integer (mpiknd), intent(out):: ia(0:)

ia(4) = loc (ia(4)) + 8
return
end subroutine

integer function mpwprecr (ra)
implicit none
integer (mpiknd), intent(in):: ra(0:)

mpwprecr = ra(1) / mpnbt
return
end function mpwprecr

integer function mpspacer (ra)
implicit none
integer (mpiknd), intent(in):: ra(0:)

mpspacer = ra(0)
return
end function mpspacer

!>>>

! integer function dqwprecr (ra)
! implicit none
! real (dqknd), intent(in):: ra(2)
!
! dqwprecr = 2
! return
! end function dqwprecr

! integer function dqspacer (ra)
! implicit none
! real (dqknd), intent(in):: ra(2)
!
! dqspacer = 2
! return
! end function dqspacer

!>>>>

subroutine mpberner (nb1, nb2, berne, mpnw)

!   This returns the array berne, containing Bernoulli numbers indexed 2*k for
!   k = 1 to n, to mpnw words precision. This is done by first computing
!   zeta(2*k), based on the following known formulas:

!   coth (pi*x) = cosh (pi*x) / sinh (pi*x)

!            1      1 + (pi*x)^2/2! + (pi*x)^4/4! + ...
!        =  ---- * -------------------------------------
!           pi*x    1 + (pi*x)^2/3! + (pi*x)^4/5! + ...

!        = 1/(pi*x) * (1 + (pi*x)^2/3 - (pi*x)^4/45 + 2*(pi*x)^6/945 - ...)

!        = 2/(pi*x) * Sum_{k >= 1} (-1)^(k+1) * zeta(2*k) * x^{2*k}

!   The strategy is to calculate the coefficients of the series by polynomial
!   operations. Polynomial division is performed by computing the reciprocal
!   of the denominator polynomial, by a polynomial Newton iteration, as follows.
!   Let N(x) be the polynomial approximation to the numerator series; let D(x) be
!   a polynomial approximation to the numerator numerator series; and let Q_k(x)
!   be polynomial approximations to R(x) = 1/D(x). Then iterate:

!   Q_{k+1} = Q_k(x) + [1 - D(x)*Q_k(x)]*Q_k(x)

!   In these iterations, both the degree of the polynomial Q_k(x) and the
!   precision level in words are initially set to 4. When convergence is
!   achieved at this precision level, the degree is doubled, and iterations are
!   continued, etc., until the final desired degree is achieved. Then the
!   precision level is doubled and iterations are performed in a similar way,
!   until the final desired precision level is achieved. The reciprocal polynomial
!   R(x) produced by this process is then multiplied by the numerator polynomial
!   N(x) to yield an approximation to the quotient series. The even zeta values
!   are then the coefficients of this series, scaled according to the formula above.

!   Once the even integer zeta values have been computed in this way, the even
!   Bernoulli numbers are computed via the formula (for n > 0):

!   B(2*n) = (-1)^(n-1) * 2 * (2*n)! * zeta(2*n) / (2*pi)^(2*n)

!   Note: The notation in the description above is not the same as in the code below.

implicit none
integer, intent(in):: nb1, nb2, mpnw
integer (mpiknd), intent(out):: berne(0:nb1+5,nb2)
integer, parameter:: ibz = 66, idb = 0
real (mprknd), parameter:: alog102 = 0.30102999566398119d0, pi = 3.1415926535897932385d0
integer i, i1, ic1, j, kn, mpnw1, mpnw2, n, n1, nn1
real (mprknd) d1, dd1, dd2, dd3
integer (mpiknd) c1(0:mpnw+6,0:nb2), cp2(0:mpnw+6), p1(0:mpnw+6,0:nb2), &
  p2(0:mpnw+6,0:nb2), q(0:mpnw+6,0:nb2), q1(0:mpnw+6), &
  r(0:mpnw+6,0:nb2), s(0:mpnw+6,0:nb2), t1(0:mpnw+6), t2(0:mpnw+6), &
  t3(0:mpnw+6), t4(0:mpnw+6), eps(0:mpnw+6)

integer, external:: mpfrsgn
integer, external:: mpfrcmp
real (mprknd), external:: mpfrgetd2exp
integer (mpiknd) ix8
 
!  End of declaration

i1 = 1000000000

do i = 1, nb2
  i1 = min (i1, mpspacer (berne(0:nb1+5,i)))
enddo

if (mpnw < 4 .or. i1 < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPBERNER: Uninitialized or inadequately sized arrays')
  call mpabrt ( 501)
endif

!  Check if PI has been precomputed.

!+ call mpmdc (mppicon, d1, n1, mpnw+1
call mpfixlocr (mppicon)
d1 = 2.d0 * mpfrgetd2exp (ix8, mppicon(1:), mprnd)
if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
if (n1 /= 1 .or. abs (d1 * 2.d0**n1 - pi) > mprdfz &
  .or. mpwprecr (mppicon) < mpnw+1) then
  write (mpldb, 2) mpnw+1
2 format ('*** MPBERNER: Pi must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 502)
endif

!  End of initial input array check

n = nb2
mpnw1 = min (mpnw + 1, mpnwx)
mpnw2 = mpnw1

if (idb > 0) write (mpldb, 3) n, mpnw
3 format ('Even Bernoulli number calculation; n, mpnw =',2i6)

call mpinitwds (cp2, mpnw2)
call mpinitwds (q1, mpnw2)
call mpinitwds (t1, mpnw2)
call mpinitwds (t2, mpnw2)
call mpinitwds (t3, mpnw2)
call mpinitwds (t4, mpnw2)
call mpinitwds (eps, mpnw2)

do i = 0, nb2
  call mpinitwds (c1(0:mpnw1+5,i), mpnw2)
  call mpinitwds (p1(0:mpnw1+5,i), mpnw2)
  call mpinitwds (p2(0:mpnw1+5,i), mpnw2)
  call mpinitwds (q(0:mpnw1+5,i), mpnw2)
  call mpinitwds (r(0:mpnw1+5,i), mpnw2)
  call mpinitwds (s(0:mpnw1+5,i), mpnw2)
enddo

!+ call mpmul (mppicon, mppicon, cp2, mpnw2
call mpfixlocr (mppicon)
call mpfrsetprec (cp2(1:), mpnw2)
call mpfrmul (cp2(1:), mppicon(1:), mppicon(1:), mprnd)
 
!+ call mpdmc (1.d0, 0, c1(0:mpnw1+5,0), mpnw2
call mpfrsetprec (c1(1:mpnw1+5,0), mpnw2)
call mpfrsetd (c1(1:mpnw1+5,0), 1.d0, mprnd)
 
!+ call mpdmc (1.d0, 0, p1(0:mpnw1+5,0), mpnw2
call mpfrsetprec (p1(1:mpnw1+5,0), mpnw2)
call mpfrsetd (p1(1:mpnw1+5,0), 1.d0, mprnd)
 
!+ call mpdmc (1.d0, 0, p2(0:mpnw1+5,0), mpnw2
call mpfrsetprec (p2(1:mpnw1+5,0), mpnw2)
call mpfrsetd (p2(1:mpnw1+5,0), 1.d0, mprnd)
 
!+ call mpdmc (1.d0, 0, q(0:mpnw1+5,0), mpnw2
call mpfrsetprec (q(1:mpnw1+5,0), mpnw2)
call mpfrsetd (q(1:mpnw1+5,0), 1.d0, mprnd)
 

!   Construct numerator and denominator polynomials.

do i = 1, n
!+   call mpdmc (0.d0, 0, c1(0:mpnw1+5,i), mpnw2
  call mpfrsetprec (c1(1:mpnw1+5,i), mpnw2)
  call mpfrsetd (c1(1:mpnw1+5,i), 0.d0, mprnd)
 
  dd1 = 2.d0 * (i + 1) - 3.d0
  dd2 = dd1 + 1.d0
  dd3 = dd2 + 1.d0
!+   call mpmul (cp2, p1(0:mpnw1+5,i-1), t1, mpnw2
  call mpfrsetprec (t1(1:), mpnw2)
  call mpfrmul (t1(1:), cp2(1:), p1(1:mpnw1+5,i-1), mprnd)
 
!+   call mpdivd (t1, dd1 * dd2, p1(0:mpnw1+5,i), mpnw2
  call mpfrsetprec (p1(1:mpnw1+5,i), mpnw2)
  call mpfrdivd (p1(1:mpnw1+5,i), t1(1:), dd1 * dd2, mprnd)
 
!+   call mpmul (cp2, p2(0:mpnw1+5,i-1), t1, mpnw2
  call mpfrsetprec (t1(1:), mpnw2)
  call mpfrmul (t1(1:), cp2(1:), p2(1:mpnw1+5,i-1), mprnd)
 
!+   call mpdivd (t1, dd2 * dd3, p2(0:mpnw1+5,i), mpnw2
  call mpfrsetprec (p2(1:mpnw1+5,i), mpnw2)
  call mpfrdivd (p2(1:mpnw1+5,i), t1(1:), dd2 * dd3, mprnd)
 
!+   call mpdmc (0.d0, 0, q(0:mpnw1+5,i), mpnw2
  call mpfrsetprec (q(1:mpnw1+5,i), mpnw2)
  call mpfrsetd (q(1:mpnw1+5,i), 0.d0, mprnd)
 
enddo

kn = 4
mpnw2 = min (4, mpnwx)

!+ call mpdmc (2.d0, 0, t1, mpnw2
call mpfrsetprec (t1(1:), mpnw2)
call mpfrsetd (t1(1:), 2.d0, mprnd)
 
!+ call mpnpwr (t1, ibz - mpnw2 * mpnbt, eps, mpnw2
call mpfrsetprec (eps(1:), mpnw2)
call mpfrpowsi (eps(1:), t1(1:), ibz - mpnw2 * mpnbt, mprnd)
 
if (idb > 0) then
!+   call mpmdc (eps, dd1, nn1, mpnw2
  dd1 = 2.d0 * mpfrgetd2exp (ix8, eps(1:), mprnd)
  if (dd1 == 0) then; nn1 = 0; else; nn1 = ix8 - 1; endif
 
  write (mpldb, 4) mpnw2, nint (alog102*nn1)
4 format ('mpnw2, log10eps =',2i6)
endif
!+ call mpdmc (0.d0, 0, q1, mpnw2
call mpfrsetprec (q1(1:), mpnw2)
call mpfrsetd (q1(1:), 0.d0, mprnd)
 

!   Perform Newton iterations with dynamic precision levels, using an
!   iteration formula similar to that used to evaluate reciprocals.

do j = 1, 10000
  if (idb > 0) write (mpldb, 5) j, kn, mpnw2
5 format ('j, kn, mpnw2 =',3i6)

  call mppolymul (mpnw1, kn, p2, q, r, mpnw2)
  call mppolysub (mpnw1, kn, c1, r, s, mpnw2)
  call mppolymul (mpnw1, kn, s, q, r, mpnw2)
  call mppolyadd (mpnw1, kn, q, r, q, mpnw2)
!+   call mpsub (q(0:mpnw1+5,kn), q1, t1, mpnw2
  call mpfrsetprec (t1(1:), mpnw2)
  call mpfrsub (t1(1:), q(1:mpnw1+5,kn), q1(1:), mprnd)
 

  if (idb > 0) then
!+     call mpmdc (t1, dd1, nn1, mpnw2
    dd1 = 2.d0 * mpfrgetd2exp (ix8, t1(1:), mprnd)
    if (dd1 == 0) then; nn1 = 0; else; nn1 = ix8 - 1; endif
 
    if (dd1 == 0.d0) then
      write (mpldb, 6)
6     format ('Newton error = 0')
    else
      write (mpldb, 7) nint (alog102*nn1)
7     format ('Newton error = 10^',i6)
    endif
  endif

!+   call mpabs (t1, t2, mpnw2
  call mpfrsetprec (t2(1:), mpnw2)
  call mpfrabs (t2(1:), t1(1:), mprnd)
 
!+   call mpcpr (t2, eps, ic1, mpnw2
  ic1 = mpfrcmp (t2(1:), eps(1:))
 
  if (ic1 < 0) then
    if (kn == n .and. mpnw2 == mpnw1) goto 100
    if (kn < n) then
      kn = min (2 * kn, n)
!+       call mpdmc (0.d0, 0, q1, mpnw2
      call mpfrsetprec (q1(1:), mpnw2)
      call mpfrsetd (q1(1:), 0.d0, mprnd)
 
    elseif (mpnw2 < mpnw1) then
      mpnw2 = min (2 * mpnw2, mpnw1)
!+       call mpdmc (2.d0, 0, t1, mpnw2
      call mpfrsetprec (t1(1:), mpnw2)
      call mpfrsetd (t1(1:), 2.d0, mprnd)
 
!+       call mpnpwr (t1, ibz - mpnw2 * mpnbt, eps, mpnw2
      call mpfrsetprec (eps(1:), mpnw2)
      call mpfrpowsi (eps(1:), t1(1:), ibz - mpnw2 * mpnbt, mprnd)
 
!+       call mpdmc (0.d0, 0, q1, mpnw2
      call mpfrsetprec (q1(1:), mpnw2)
      call mpfrsetd (q1(1:), 0.d0, mprnd)
 
      if (idb > 0) then
!+         call mpmdc (eps, dd1, nn1, mpnw2
        dd1 = 2.d0 * mpfrgetd2exp (ix8, eps(1:), mprnd)
        if (dd1 == 0) then; nn1 = 0; else; nn1 = ix8 - 1; endif
 
        write (mpldb, 4) mpnw2, nint (alog102*nn1)
      endif
    endif
  else
!+     call mpeq (q(0:mpnw1+5,kn), q1, mpnw2
    call mpfrsetprec (q1(1:), mpnw2)
    call mpfrset (q1(1:), q(1:mpnw1+5,kn), mprnd)
 
  endif
enddo

write (mpldb, 8)
8 format ('*** MPBERNER: Loop end error')
call mpabrt ( 503)

100 continue

if (idb > 0) write (mpldb, 9)
9 format ('Even zeta computation complete')

!   Multiply numerator polynomial by reciprocal of denominator polynomial.

call mppolymul (mpnw1, n, p1, q, r, mpnw2)

!   Apply formula to produce Bernoulli numbers.

!+ call mpdmc (-2.d0, 0, t1, mpnw2
call mpfrsetprec (t1(1:), mpnw2)
call mpfrsetd (t1(1:), -2.d0, mprnd)
 
!+ call mpdmc (1.d0, 0, t2, mpnw2
call mpfrsetprec (t2(1:), mpnw2)
call mpfrsetd (t2(1:), 1.d0, mprnd)
 

do i = 1, n
  d1 = - dble (2*i-1) * dble (2*i)
!+   call mpmuld (t1, d1, t3, mpnw2
  call mpfrsetprec (t3(1:), mpnw2)
  call mpfrmuld (t3(1:), t1(1:), d1, mprnd)
 
!+   call mpeq (t3, t1, mpnw2
  call mpfrsetprec (t1(1:), mpnw2)
  call mpfrset (t1(1:), t3(1:), mprnd)
 
!+   call mpmuld (cp2, 4.d0, t3, mpnw2
  call mpfrsetprec (t3(1:), mpnw2)
  call mpfrmuld (t3(1:), cp2(1:), 4.d0, mprnd)
 
!+   call mpmul (t3, t2, t4, mpnw2
  call mpfrsetprec (t4(1:), mpnw2)
  call mpfrmul (t4(1:), t3(1:), t2(1:), mprnd)
 
!+   call mpeq (t4, t2, mpnw2
  call mpfrsetprec (t2(1:), mpnw2)
  call mpfrset (t2(1:), t4(1:), mprnd)
 
!+   call mpmuld (t1, 0.5d0, t3, mpnw2
  call mpfrsetprec (t3(1:), mpnw2)
  call mpfrmuld (t3(1:), t1(1:), 0.5d0, mprnd)
 
!+   call mpdiv (t3, t2, t4, mpnw2
  call mpfrsetprec (t4(1:), mpnw2)
  call mpfrdiv (t4(1:), t3(1:), t2(1:), mprnd)
 
!+   call mpabs (r(0:mpnw1+5,i), t3, mpnw2
  call mpfrsetprec (t3(1:), mpnw2)
  call mpfrabs (t3(1:), r(1:mpnw1+5,i), mprnd)
 
!+   call mpmul (t4, t3, berne(0:nb1+5,i), mpnw
  call mpfrsetprec (berne(1:nb1+5,i), mpnw)
  call mpfrmul (berne(1:nb1+5,i), t4(1:), t3(1:), mprnd)
 
enddo

if (idb > 0) write (mpldb, 10)
10 format ('Bernoulli number computation complete')
return
end subroutine mpberner

subroutine mppolyadd (nb1, n, a, b, c, mpnw)

!   This adds two polynomials, as is required by mpberne.
!   The output array C may be the same as A or B.

implicit none
integer, intent(in):: nb1, n, mpnw
integer (mpiknd), intent(in):: a(0:nb1+5,0:n), b(0:nb1+5,0:n)
integer (mpiknd), intent(out):: c(0:nb1+5,0:n)
integer k
integer (mpiknd) t1(0:mpnw+5), t2(0:mpnw+5)

integer, external:: mpfrsgn
integer, external:: mpfrcmp
real (mprknd), external:: mpfrgetd2exp
integer (mpiknd) ix8
 
!  End of declaration
!  End of initial input array check

call mpinitwds (t1, mpnw)
call mpinitwds (t2, mpnw)

do k = 0, n
!+   call mpeq (a(0:nb1+5,k), t1, mpnw
  call mpfrsetprec (t1(1:), mpnw)
  call mpfrset (t1(1:), a(1:nb1+5,k), mprnd)
 
!+   call mpeq (b(0:nb1+5,k), t2, mpnw
  call mpfrsetprec (t2(1:), mpnw)
  call mpfrset (t2(1:), b(1:nb1+5,k), mprnd)
 
!+   call mpadd (t1, t2, c(0:nb1+5,k), mpnw
  call mpfrsetprec (c(1:nb1+5,k), mpnw)
  call mpfradd (c(1:nb1+5,k), t1(1:), t2(1:), mprnd)
 
enddo

return
end subroutine mppolyadd

subroutine mppolysub (nb1, n, a, b, c, mpnw)

!   This adds two polynomials, as is required by mpberne.
!   The output array C may be the same as A or B.

implicit none
integer, intent(in):: nb1, n, mpnw
integer (mpiknd), intent(in):: a(0:nb1+5,0:n), b(0:nb1+5,0:n)
integer (mpiknd), intent(out):: c(0:nb1+5,0:n)
integer k
integer (mpiknd) t1(0:mpnw+5), t2(0:mpnw+5)

integer, external:: mpfrsgn
integer, external:: mpfrcmp
real (mprknd), external:: mpfrgetd2exp
integer (mpiknd) ix8
 
!  End of declaration
!  End of initial input array check

call mpinitwds (t1, mpnw)
call mpinitwds (t2, mpnw)

do k = 0, n
!+   call mpeq (a(0:nb1+5,k), t1, mpnw
  call mpfrsetprec (t1(1:), mpnw)
  call mpfrset (t1(1:), a(1:nb1+5,k), mprnd)
 
!+   call mpeq (b(0:nb1+5,k), t2, mpnw
  call mpfrsetprec (t2(1:), mpnw)
  call mpfrset (t2(1:), b(1:nb1+5,k), mprnd)
 
!+   call mpsub (t1, t2, c(0:nb1+5,k), mpnw
  call mpfrsetprec (c(1:nb1+5,k), mpnw)
  call mpfrsub (c(1:nb1+5,k), t1(1:), t2(1:), mprnd)
 
enddo

return
end subroutine mppolysub

subroutine mppolymul (nb1, n, a, b, c, mpnw)

!   This adds two polynomials (ignoring high-order terms), as is required
!   by mpberne. The output array C may not be the same as A or B.

implicit none
integer, intent(in):: nb1, n, mpnw
integer (mpiknd), intent(in):: a(0:nb1+5,0:n), b(0:nb1+5,0:n)
integer (mpiknd), intent(out):: c(0:nb1+5,0:n)
integer j, k
integer (mpiknd) t0(0:mpnw+5), t1(0:mpnw+5), t2(0:mpnw+5), t3(0:mpnw+5)

integer, external:: mpfrsgn
integer, external:: mpfrcmp
real (mprknd), external:: mpfrgetd2exp
integer (mpiknd) ix8
 
!  End of declaration
!  End of initial input array check

call mpinitwds (t0, mpnw)
call mpinitwds (t1, mpnw)
call mpinitwds (t2, mpnw)
call mpinitwds (t3, mpnw)

do k = 0, n
!+   call mpdmc (0.d0, 0, t0, mpnw
  call mpfrsetprec (t0(1:), mpnw)
  call mpfrsetd (t0(1:), 0.d0, mprnd)
 

  do j = 0, k
!+     call mpeq (a(0:nb1+5,j), t1, mpnw
    call mpfrsetprec (t1(1:), mpnw)
    call mpfrset (t1(1:), a(1:nb1+5,j), mprnd)
 
!+     call mpeq (b(0:nb1+5,k-j), t2, mpnw
    call mpfrsetprec (t2(1:), mpnw)
    call mpfrset (t2(1:), b(1:nb1+5,k-j), mprnd)
 
!+     call mpmul (t1, t2, t3, mpnw
    call mpfrsetprec (t3(1:), mpnw)
    call mpfrmul (t3(1:), t1(1:), t2(1:), mprnd)
 
!+     call mpadd (t0, t3, t2, mpnw
    call mpfrsetprec (t2(1:), mpnw)
    call mpfradd (t2(1:), t0(1:), t3(1:), mprnd)
 
!+     call mpeq (t2, t0, mpnw
    call mpfrsetprec (t0(1:), mpnw)
    call mpfrset (t0(1:), t2(1:), mprnd)
 
  enddo

!+   call mpeq (t0, c(0:nb1+5,k), mpnw
  call mpfrsetprec (c(1:nb1+5,k), mpnw)
  call mpfrset (c(1:nb1+5,k), t0(1:), mprnd)
 
enddo

return
end subroutine mppolymul

subroutine mpbesselinr (nu, rr, ss, mpnw)

!   This evaluates the modified Bessel function BesselI (NU,RR).
!   NU is an integer. The algorithm is DLMF formula 10.25.2 for modest RR,
!   and DLMF 10.40.1 for large RR, relative to precision.

implicit none
integer, intent(in):: nu, mpnw
integer (mpiknd), intent(in):: rr(0:)
integer (mpiknd), intent(out):: ss(0:)
integer, parameter:: itrmax = 1000000
real (mprknd), parameter:: dfrac = 1.5d0, pi = 3.1415926535897932385d0
integer ic1, k, mpnw1, nua, n1
real (mprknd) d1
integer (mpiknd) f1(0:mpnw+6), f2(0:mpnw+6), sum(0:mpnw+6), td(0:mpnw+6), &
  tn(0:mpnw+6), t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), &
  rra(0:mpnw+6), tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)

integer, external:: mpfrsgn
integer, external:: mpfrcmp
real (mprknd), external:: mpfrgetd2exp
integer (mpiknd) ix8
 
!  End of declaration

if (mpnw < 4 .or. mpspacer (rr) < mpnw + 4 .or. mpspacer (ss) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPBESSELINR: Uninitialized or inadequately sized arrays')
  call mpabrt ( 504)
endif

!  Check if PI has been precomputed.

!+ call mpmdc (mppicon, d1, n1, mpnw+1
call mpfixlocr (mppicon)
d1 = 2.d0 * mpfrgetd2exp (ix8, mppicon(1:), mprnd)
if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
if (n1 /= 1 .or. abs (d1 * 2.d0**n1 - pi) > mprdfz &
  .or. mpwprecr (mppicon) < mpnw+1) then
  write (mpldb, 3) mpnw+1
3 format ('*** MPBESSELINR: Pi must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 506)
endif

!  End of initial input array check

!   Check for RR = 0.

!+ if (mpsgn (rr) == 0) the
if (mpfrsgn (rr(1:)) == 0) then
 
  write (mpldb, 2)
2 format ('*** MPBESSELINR: Second argument is zero')
  call mpabrt ( 505)
endif

mpnw1 = min (mpnw + 1, mpnwx)
nua = abs (nu)
call mpinitwds (f1, mpnw1)
call mpinitwds (f2, mpnw1)
call mpinitwds (sum, mpnw1)
call mpinitwds (td, mpnw1)
call mpinitwds (tn, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (rra, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
!+ call mpdmc (2.d0, 0, tc1, 4
call mpfrsetprec (tc1(1:), 4)
call mpfrsetd (tc1(1:), 2.d0, mprnd)
 
!+ call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4
call mpfrsetprec (eps(1:), 4)
call mpfrpowsi (eps(1:), tc1(1:), -mpnw1*mpnbt, mprnd)
 
!+ call mpabs (rr, rra, mpnw1
call mpfrsetprec (rra(1:), mpnw1)
call mpfrabs (rra(1:), rr(1:), mprnd)
 
!+ call mpmdc (rra, d1, n1, 4
d1 = 2.d0 * mpfrgetd2exp (ix8, rra(1:), mprnd)
if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
d1 = d1 * 2.d0 ** n1

if (d1 < dfrac * mpnw1 * mpdpw) then
!+   call mpdmc (1.d0, 0, tn, mpnw1
  call mpfrsetprec (tn(1:), mpnw1)
  call mpfrsetd (tn(1:), 1.d0, mprnd)
 
!+   call mpdmc (1.d0, 0, f1, mpnw1
  call mpfrsetprec (f1(1:), mpnw1)
  call mpfrsetd (f1(1:), 1.d0, mprnd)
 
!+   call mpdmc (1.d0, 0, f2, mpnw1
  call mpfrsetprec (f2(1:), mpnw1)
  call mpfrsetd (f2(1:), 1.d0, mprnd)
 
!+   call mpmul (rra, rra, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmul (t2(1:), rra(1:), rra(1:), mprnd)
 
!+   call mpmuld (t2, 0.25d0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmuld (t1(1:), t2(1:), 0.25d0, mprnd)
 

  do k = 1, nua
!+     call mpmuld (f2, dble (k), t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmuld (t2(1:), f2(1:), dble (k), mprnd)
 
!+     call mpeq (t2, f2, mpnw1
    call mpfrsetprec (f2(1:), mpnw1)
    call mpfrset (f2(1:), t2(1:), mprnd)
 
  enddo

!+   call mpmul (f1, f2, td, mpnw1
  call mpfrsetprec (td(1:), mpnw1)
  call mpfrmul (td(1:), f1(1:), f2(1:), mprnd)
 
!+   call mpdiv (tn, td, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrdiv (t2(1:), tn(1:), td(1:), mprnd)
 
!+   call mpeq (t2, sum, mpnw1
  call mpfrsetprec (sum(1:), mpnw1)
  call mpfrset (sum(1:), t2(1:), mprnd)
 

  do k = 1, itrmax
!+     call mpmuld (f1, dble (k), t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmuld (t2(1:), f1(1:), dble (k), mprnd)
 
!+     call mpeq (t2, f1, mpnw1
    call mpfrsetprec (f1(1:), mpnw1)
    call mpfrset (f1(1:), t2(1:), mprnd)
 
!+     call mpmuld (f2, dble (k + nua), t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmuld (t2(1:), f2(1:), dble (k + nua), mprnd)
 
!+     call mpeq (t2, f2, mpnw1
    call mpfrsetprec (f2(1:), mpnw1)
    call mpfrset (f2(1:), t2(1:), mprnd)
 
!+     call mpmul (t1, tn, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmul (t2(1:), t1(1:), tn(1:), mprnd)
 
!+     call mpeq (t2, tn, mpnw1
    call mpfrsetprec (tn(1:), mpnw1)
    call mpfrset (tn(1:), t2(1:), mprnd)
 
!+     call mpmul (f1, f2, td, mpnw1
    call mpfrsetprec (td(1:), mpnw1)
    call mpfrmul (td(1:), f1(1:), f2(1:), mprnd)
 
!+     call mpdiv (tn, td, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrdiv (t2(1:), tn(1:), td(1:), mprnd)
 
!+     call mpadd (sum, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfradd (t3(1:), sum(1:), t2(1:), mprnd)
 
!+     call mpeq (t3, sum, mpnw1
    call mpfrsetprec (sum(1:), mpnw1)
    call mpfrset (sum(1:), t3(1:), mprnd)
 

!+     call mpabs (t2, tc1, 4
    call mpfrsetprec (tc1(1:), 4)
    call mpfrabs (tc1(1:), t2(1:), mprnd)
 
!+     call mpmul (eps, sum, tc3, 4
    call mpfrsetprec (tc3(1:), 4)
    call mpfrmul (tc3(1:), eps(1:), sum(1:), mprnd)
 
!+     call mpabs (tc3, tc2, 4
    call mpfrsetprec (tc2(1:), 4)
    call mpfrabs (tc2(1:), tc3(1:), mprnd)
 
!+     call mpcpr (tc1, tc2, ic1, 4
    ic1 = mpfrcmp (tc1(1:), tc2(1:))
 
    if (ic1 <= 0) goto 100
  enddo

  write (mpldb, 4)
  4 format ('*** MPBESSELINR: Loop end error 1')
  call mpabrt ( 507)

100 continue

!+   call mpmuld (rra, 0.5d0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmuld (t1(1:), rra(1:), 0.5d0, mprnd)
 
!+   call mpnpwr (t1, nua, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrpowsi (t2(1:), t1(1:), nua, mprnd)
 
!+   call mpmul (sum, t2, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrmul (t3(1:), sum(1:), t2(1:), mprnd)
 
else
!  sum1 = mpreal (1.d0, mpnw)
!  t1 = mpreal (4.d0 * dble (nua)**2, mpnw)
!  t2 = mpreal (1.d0, mpnw)
!  t3 = mpreal (1.d0, mpnw)

!+   call mpdmc (1.d0, 0, sum, mpnw1
  call mpfrsetprec (sum(1:), mpnw1)
  call mpfrsetd (sum(1:), 1.d0, mprnd)
 
  d1 = 4.d0 * dble (nua)**2
!+   call mpdmc (d1, 0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrsetd (t1(1:), d1, mprnd)
 
!+   call mpdmc (1.d0, 0, tn, mpnw1
  call mpfrsetprec (tn(1:), mpnw1)
  call mpfrsetd (tn(1:), 1.d0, mprnd)
 
!+   call mpdmc (1.d0, 0, td, mpnw1
  call mpfrsetprec (td(1:), mpnw1)
  call mpfrsetd (td(1:), 1.d0, mprnd)
 

  do k = 1, itrmax
!  t2 = -t2 * (t1 - (2.d0*k - 1.d0)**2)
!  t3 = t3 * dble (k) * 8.d0 * xa
!  t4 = t2 / t3
!  sum1 = sum1 + t4

    d1 = 2.d0 * k - 1.d0
!+     call mpdmc (d1, 0, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrsetd (t2(1:), d1, mprnd)
 
!+     call mpmul (t2, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrmul (t3(1:), t2(1:), t2(1:), mprnd)
 
!+     call mpsub (t1, t3, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrsub (t2(1:), t1(1:), t3(1:), mprnd)
 
!+     call mpmul (tn, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrmul (t3(1:), tn(1:), t2(1:), mprnd)
 
!+     call mpneg (t3, tn, mpnw1
    call mpfrsetprec (tn(1:), mpnw1)
    call mpfrneg (tn(1:), t3(1:), mprnd)
 
!+     call mpmuld (rra, 8.d0 * k, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmuld (t2(1:), rra(1:), 8.d0 * k, mprnd)
 
!+     call mpmul (td, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrmul (t3(1:), td(1:), t2(1:), mprnd)
 
!+     call mpeq (t3, td, mpnw1
    call mpfrsetprec (td(1:), mpnw1)
    call mpfrset (td(1:), t3(1:), mprnd)
 
!+     call mpdiv (tn, td, t4, mpnw1
    call mpfrsetprec (t4(1:), mpnw1)
    call mpfrdiv (t4(1:), tn(1:), td(1:), mprnd)
 
!+     call mpadd (sum, t4, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfradd (t3(1:), sum(1:), t4(1:), mprnd)
 
!+     call mpeq (t3, sum, mpnw1
    call mpfrsetprec (sum(1:), mpnw1)
    call mpfrset (sum(1:), t3(1:), mprnd)
 

!   if (abs (t4) / abs (sum1) < eps) goto 110

!+     call mpabs (t4, tc1, 4
    call mpfrsetprec (tc1(1:), 4)
    call mpfrabs (tc1(1:), t4(1:), mprnd)
 
!+     call mpmul (eps, sum, tc3, 4
    call mpfrsetprec (tc3(1:), 4)
    call mpfrmul (tc3(1:), eps(1:), sum(1:), mprnd)
 
!+     call mpabs (tc3, tc2, 4
    call mpfrsetprec (tc2(1:), 4)
    call mpfrabs (tc2(1:), tc3(1:), mprnd)
 
!+     call mpcpr (tc1, tc2, ic1, 4
    ic1 = mpfrcmp (tc1(1:), tc2(1:))
 
    if (ic1 <= 0) goto 110
  enddo

write (mpldb, 5)
5 format ('*** MPBESSELINR: Loop end error 2')
call mpabrt ( 508)

110 continue

! t1 = exp (xa) / sqrt (2.d0 * mppi (mpnw) * xa)
! besseli = t1 * sum1

!+   call mpexp (rra, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrexp (t1(1:), rra(1:), mprnd)
 
!+   call mpmuld (mppicon, 2.d0, t2, mpnw1
call mpfixlocr (mppicon)
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmuld (t2(1:), mppicon(1:), 2.d0, mprnd)
 
!+   call mpmul (t2, rra, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrmul (t3(1:), t2(1:), rra(1:), mprnd)
 
!+   call mpsqrt (t3, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfrsqrt (t4(1:), t3(1:), mprnd)
 
!+   call mpdiv (t1, t4, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrdiv (t2(1:), t1(1:), t4(1:), mprnd)
 
!+   call mpmul (t2, sum, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrmul (t3(1:), t2(1:), sum(1:), mprnd)
 
endif

! if (x < 0.d0 .and. mod (nu, 2) /= 0) besseli = - besseli

!+ if (mpsgn (rr) < 0 .and. mod (nu, 2) /= 0) the
if (mpfrsgn (rr(1:)) < 0 .and. mod (nu, 2) /= 0) then
 
!+   call mpneg (t3, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfrneg (t4(1:), t3(1:), mprnd)
 
!+   call mpeq (t4, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrset (t3(1:), t4(1:), mprnd)
 
endif

!+ call mproun (t3, mpnw
call mpfrprecround (t3(1:), mpnw, mprnd)
 
!+ call mpeq (t3, ss, mpnw
call mpfrsetprec (ss(1:), mpnw)
call mpfrset (ss(1:), t3(1:), mprnd)
 

return
end subroutine mpbesselinr

subroutine mpbesselir (qq, rr, ss, mpnw)

!   This evaluates the modified Bessel function BesselI (QQ,RR) for QQ and RR
!   both MPR. The algorithm is DLMF formula 10.25.2 for modest RR, and
!   DLMF 10.40.1 for large RR, relative to precision.

implicit none
integer (mpiknd), intent(in):: qq(0:), rr(0:)
integer (mpiknd), intent(out):: ss(0:)
integer, intent(in):: mpnw
integer ic1, i1, i2, k, mpnw1, n1
integer, parameter:: itrmax = 1000000
real (mprknd), parameter:: dfrac = 1.5d0, pi = 3.1415926535897932385d0
real (mprknd) d1
integer (mpiknd) f1(0:mpnw+6), f2(0:mpnw+6), sum(0:mpnw+6), td(0:mpnw+6), &
  tn(0:mpnw+6), t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), &
  rra(0:mpnw+6), tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)

integer, external:: mpfrsgn
integer, external:: mpfrcmp
real (mprknd), external:: mpfrgetd2exp
integer (mpiknd) ix8
 
!  End of declaration

if (mpnw < 4 .or. mpspacer (rr) < mpnw + 4 .or. mpspacer (qq) < mpnw + 4 &
  .or. mpspacer (ss) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPBESSELIR: Uninitialized or inadequately sized arrays')
  call mpabrt ( 509)
endif

!  Check if PI has been precomputed.

!+ call mpmdc (mppicon, d1, n1, mpnw+1
call mpfixlocr (mppicon)
d1 = 2.d0 * mpfrgetd2exp (ix8, mppicon(1:), mprnd)
if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
if (n1 /= 1 .or. abs (d1 * 2.d0**n1 - pi) > mprdfz &
  .or. mpwprecr (mppicon) < mpnw+1) then
  write (mpldb, 2) mpnw+1
2 format ('*** MPBESSELIR: Pi must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 510)
endif

!  End of initial input array check

mpnw1 = min (mpnw + 1, mpnwx)
call mpinitwds (f1, mpnw1)
call mpinitwds (f2, mpnw1)
call mpinitwds (sum, mpnw1)
call mpinitwds (td, mpnw1)
call mpinitwds (tn, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (rra, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
!+ call mpdmc (2.d0, 0, tc1, 4
call mpfrsetprec (tc1(1:), 4)
call mpfrsetd (tc1(1:), 2.d0, mprnd)
 
!+ call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4
call mpfrsetprec (eps(1:), 4)
call mpfrpowsi (eps(1:), tc1(1:), -mpnw1*mpnbt, mprnd)
 

!   If QQ is integer, call mpbesselinr; if qq < 0 and rr <= 0, then error.

!+ call mpinfr (qq, t1, t2, mpnw
call mpfrsetprec (t1(1:), mpnw)
call mpfrsetprec (t2(1:), mpnw)
call mpfrtrunc (t1(1:), qq(1:))
call mpfrsub (t2(1:), qq(1:), t1(1:), mprnd)
 
!+ i1 = mpsgn (qq
i1 = mpfrsgn (qq(1:))
 
!+ i2 = mpsgn (rr
i2 = mpfrsgn (rr(1:))
 
!+ if (mpsgn (t2) == 0) the
if (mpfrsgn (t2(1:)) == 0) then
 
!+   call mpmdc (qq, d1, n1, mpnw
  d1 = 2.d0 * mpfrgetd2exp (ix8, qq(1:), mprnd)
  if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
  d1 = d1 * 2.d0**n1
  n1 = nint (d1)
  call mpbesselinr (n1, rr, t3, mpnw)
  goto 120
elseif (i1 < 0 .and. i2 <= 0) then
  write (mpldb, 3)
3 format ('*** MPBESSELIR: First argument < 0 and second argument <= 0')
  call mpabrt ( 511)
endif

!+ call mpabs (rr, rra, mpnw1
call mpfrsetprec (rra(1:), mpnw1)
call mpfrabs (rra(1:), rr(1:), mprnd)
 
!+ call mpmdc (rra, d1, n1, 4
d1 = 2.d0 * mpfrgetd2exp (ix8, rra(1:), mprnd)
if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
d1 = d1 * 2.d0 ** n1

if (d1 < dfrac * mpnw1 * mpdpw) then
!+   call mpdmc (1.d0, 0, tn, mpnw1
  call mpfrsetprec (tn(1:), mpnw1)
  call mpfrsetd (tn(1:), 1.d0, mprnd)
 
!+   call mpdmc (1.d0, 0, f1, mpnw1
  call mpfrsetprec (f1(1:), mpnw1)
  call mpfrsetd (f1(1:), 1.d0, mprnd)
 
!+   call mpadd (qq, f1, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfradd (t1(1:), qq(1:), f1(1:), mprnd)
 
!+   call mpgammar (t1, f2, mpnw1
  call mpfrsetprec (f2(1:), mpnw1)
  call mpfrgamma (f2(1:), t1(1:), mprnd)
 
!+   call mpmul (rra, rra, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmul (t2(1:), rra(1:), rra(1:), mprnd)
 
!+   call mpmuld (t2, 0.25d0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmuld (t1(1:), t2(1:), 0.25d0, mprnd)
 

!+   call mpmul (f1, f2, td, mpnw1
  call mpfrsetprec (td(1:), mpnw1)
  call mpfrmul (td(1:), f1(1:), f2(1:), mprnd)
 
!+   call mpdiv (tn, td, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrdiv (t2(1:), tn(1:), td(1:), mprnd)
 
!+   call mpeq (t2, sum, mpnw1
  call mpfrsetprec (sum(1:), mpnw1)
  call mpfrset (sum(1:), t2(1:), mprnd)
 

  do k = 1, itrmax
!+     call mpmuld (f1, dble (k), t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmuld (t2(1:), f1(1:), dble (k), mprnd)
 
!+     call mpeq (t2, f1, mpnw1
    call mpfrsetprec (f1(1:), mpnw1)
    call mpfrset (f1(1:), t2(1:), mprnd)
 
!+     call mpdmc (dble (k), 0, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrsetd (t3(1:), dble (k), mprnd)
 
!+     call mpadd (qq, t3, t4, mpnw1
    call mpfrsetprec (t4(1:), mpnw1)
    call mpfradd (t4(1:), qq(1:), t3(1:), mprnd)
 
!+     call mpmul (f2, t4, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrmul (t3(1:), f2(1:), t4(1:), mprnd)
 
!+     call mpeq (t3, f2, mpnw1
    call mpfrsetprec (f2(1:), mpnw1)
    call mpfrset (f2(1:), t3(1:), mprnd)
 
!+     call mpmul (t1, tn, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmul (t2(1:), t1(1:), tn(1:), mprnd)
 
!+     call mpeq (t2, tn, mpnw1
    call mpfrsetprec (tn(1:), mpnw1)
    call mpfrset (tn(1:), t2(1:), mprnd)
 
!+     call mpmul (f1, f2, td, mpnw1
    call mpfrsetprec (td(1:), mpnw1)
    call mpfrmul (td(1:), f1(1:), f2(1:), mprnd)
 
!+     call mpdiv (tn, td, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrdiv (t2(1:), tn(1:), td(1:), mprnd)
 
!+     call mpadd (sum, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfradd (t3(1:), sum(1:), t2(1:), mprnd)
 
!+     call mpeq (t3, sum, mpnw1
    call mpfrsetprec (sum(1:), mpnw1)
    call mpfrset (sum(1:), t3(1:), mprnd)
 

!+     call mpabs (t2, tc1, 4
    call mpfrsetprec (tc1(1:), 4)
    call mpfrabs (tc1(1:), t2(1:), mprnd)
 
!+     call mpmul (eps, sum, tc3, 4
    call mpfrsetprec (tc3(1:), 4)
    call mpfrmul (tc3(1:), eps(1:), sum(1:), mprnd)
 
!+     call mpabs (tc3, tc2, 4
    call mpfrsetprec (tc2(1:), 4)
    call mpfrabs (tc2(1:), tc3(1:), mprnd)
 
!+     call mpcpr (tc1, tc2, ic1, 4
    ic1 = mpfrcmp (tc1(1:), tc2(1:))
 
    if (ic1 <= 0) goto 100
  enddo

  write (mpldb, 4)
  4 format ('*** MPBESSELIR: Loop end error 1')
  call mpabrt ( 512)

100 continue

!+   call mpmuld (rra, 0.5d0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmuld (t1(1:), rra(1:), 0.5d0, mprnd)
 
!+   call mppower (t1, qq, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrpow (t2(1:), t1(1:), qq(1:), mprnd)
 
!+   call mpmul (sum, t2, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrmul (t3(1:), sum(1:), t2(1:), mprnd)
 
else
!  sum1 = mpreal (1.d0, mpnw)
!  t1 = mpreal (4.d0 * dble (nua)**2, mpnw)
!  t2 = mpreal (1.d0, mpnw)
!  t3 = mpreal (1.d0, mpnw)

!+   call mpdmc (1.d0, 0, sum, mpnw1
  call mpfrsetprec (sum(1:), mpnw1)
  call mpfrsetd (sum(1:), 1.d0, mprnd)
 
!+   call mpmul (qq, qq, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmul (t2(1:), qq(1:), qq(1:), mprnd)
 
!+   call mpmuld (t2, 4.d0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmuld (t1(1:), t2(1:), 4.d0, mprnd)
 
!+   call mpdmc (1.d0, 0, tn, mpnw1
  call mpfrsetprec (tn(1:), mpnw1)
  call mpfrsetd (tn(1:), 1.d0, mprnd)
 
!+   call mpdmc (1.d0, 0, td, mpnw1
  call mpfrsetprec (td(1:), mpnw1)
  call mpfrsetd (td(1:), 1.d0, mprnd)
 

  do k = 1, itrmax
!  t2 = -t2 * (t1 - (2.d0*k - 1.d0)**2)
!  t3 = t3 * dble (k) * 8.d0 * xa
!  t4 = t2 / t3
!  sum1 = sum1 + t4

    d1 = 2.d0 * k - 1.d0
!+     call mpdmc (d1, 0, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrsetd (t2(1:), d1, mprnd)
 
!+     call mpmul (t2, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrmul (t3(1:), t2(1:), t2(1:), mprnd)
 
!+     call mpsub (t1, t3, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrsub (t2(1:), t1(1:), t3(1:), mprnd)
 
!+     call mpmul (tn, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrmul (t3(1:), tn(1:), t2(1:), mprnd)
 
!+     call mpneg (t3, tn, mpnw1
    call mpfrsetprec (tn(1:), mpnw1)
    call mpfrneg (tn(1:), t3(1:), mprnd)
 
!+     call mpmuld (rra, 8.d0 * k, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmuld (t2(1:), rra(1:), 8.d0 * k, mprnd)
 
!+     call mpmul (td, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrmul (t3(1:), td(1:), t2(1:), mprnd)
 
!+     call mpeq (t3, td, mpnw1
    call mpfrsetprec (td(1:), mpnw1)
    call mpfrset (td(1:), t3(1:), mprnd)
 
!+     call mpdiv (tn, td, t4, mpnw1
    call mpfrsetprec (t4(1:), mpnw1)
    call mpfrdiv (t4(1:), tn(1:), td(1:), mprnd)
 
!+     call mpadd (sum, t4, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfradd (t3(1:), sum(1:), t4(1:), mprnd)
 
!+     call mpeq (t3, sum, mpnw1
    call mpfrsetprec (sum(1:), mpnw1)
    call mpfrset (sum(1:), t3(1:), mprnd)
 

!   if (abs (t4) / abs (sum1) < eps) goto 110

!+     call mpabs (t4, tc1, 4
    call mpfrsetprec (tc1(1:), 4)
    call mpfrabs (tc1(1:), t4(1:), mprnd)
 
!+     call mpmul (eps, sum, tc3, 4
    call mpfrsetprec (tc3(1:), 4)
    call mpfrmul (tc3(1:), eps(1:), sum(1:), mprnd)
 
!+     call mpabs (tc3, tc2, 4
    call mpfrsetprec (tc2(1:), 4)
    call mpfrabs (tc2(1:), tc3(1:), mprnd)
 
!+     call mpcpr (tc1, tc2, ic1, 4
    ic1 = mpfrcmp (tc1(1:), tc2(1:))
 
    if (ic1 <= 0) goto 110
  enddo

write (mpldb, 5)
5 format ('*** MPBESSELIR: Loop end error 2')
call mpabrt ( 513)

110 continue

! t1 = exp (xa) / sqrt (2.d0 * mppi (mpnw) * xa)
! besseli = t1 * sum1

!+   call mpexp (rra, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrexp (t1(1:), rra(1:), mprnd)
 
!+   call mpmuld (mppicon, 2.d0, t2, mpnw1
call mpfixlocr (mppicon)
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmuld (t2(1:), mppicon(1:), 2.d0, mprnd)
 
!+   call mpmul (t2, rra, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrmul (t3(1:), t2(1:), rra(1:), mprnd)
 
!+   call mpsqrt (t3, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfrsqrt (t4(1:), t3(1:), mprnd)
 
!+   call mpdiv (t1, t4, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrdiv (t2(1:), t1(1:), t4(1:), mprnd)
 
!+   call mpmul (t2, sum, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrmul (t3(1:), t2(1:), sum(1:), mprnd)
 
endif

120 continue

!+ call mproun (t3, mpnw
call mpfrprecround (t3(1:), mpnw, mprnd)
 
!+ call mpeq (t3, ss, mpnw
call mpfrsetprec (ss(1:), mpnw)
call mpfrset (ss(1:), t3(1:), mprnd)
 

return
end subroutine mpbesselir

subroutine mpbesseljnr (nu, rr, ss, mpnw)

!   This evaluates the modified Bessel function BesselJ (NU,RR).
!   NU is an integer. The algorithm is DLMF formula 10.2.2 for modest RR,
!   and DLMF 10.17.3 for large RR, relative to precision.

implicit none
integer, intent(in):: nu, mpnw
integer (mpiknd), intent(in):: rr(0:)
integer (mpiknd), intent(out):: ss(0:)
integer, parameter:: itrmax = 1000000
real (mprknd), parameter:: dfrac = 1.5d0, pi = 3.1415926535897932385d0
integer ic1, ic2, k, mpnw1, nua, n1
real (mprknd) d1, d2
integer (mpiknd) f1(0:2*mpnw+6), f2(0:2*mpnw+6), sum1(0:2*mpnw+6), &
  sum2(0:2*mpnw+6), td1(0:2*mpnw+6), td2(0:2*mpnw+6), tn1(0:2*mpnw+6), &
  tn2(0:2*mpnw+6), t1(0:2*mpnw+6), t2(0:2*mpnw+6), t3(0:2*mpnw+6), &
  t41(0:2*mpnw+6), t42(0:2*mpnw+6), t5(0:2*mpnw+6), rra(0:2*mpnw+6), &
  rr2(0:2*mpnw+6), tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)

integer, external:: mpfrsgn
integer, external:: mpfrcmp
real (mprknd), external:: mpfrgetd2exp
integer (mpiknd) ix8
 
!  End of declaration

if (mpnw < 4 .or. mpspacer (rr) < mpnw + 4 .or. mpspacer (ss) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPBESSELJNR: Uninitialized or inadequately sized arrays')
  call mpabrt ( 514)
endif

!  Check if PI has been precomputed.

!+ call mpmdc (mppicon, d1, n1, mpnw+1
call mpfixlocr (mppicon)
d1 = 2.d0 * mpfrgetd2exp (ix8, mppicon(1:), mprnd)
if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
if (n1 /= 1 .or. abs (d1 * 2.d0**n1 - pi) > mprdfz &
  .or. mpwprecr (mppicon) < mpnw+1) then
  write (mpldb, 3) mpnw+1
3 format ('*** MPBESSELJNR: Pi must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 516)
endif

!  End of initial input array check

!   Check for RR = 0.

!+ if (mpsgn (rr) == 0) the
if (mpfrsgn (rr(1:)) == 0) then
 
  write (mpldb, 2)
2 format ('*** MPBESSELJNR: Second argument is zero')
  call mpabrt ( 515)
endif

mpnw1 = min (2 * mpnw + 1, mpnwx)
nua = abs (nu)
call mpinitwds (f1, mpnw1)
call mpinitwds (f2, mpnw1)
call mpinitwds (sum1, mpnw1)
call mpinitwds (sum2, mpnw1)
call mpinitwds (td1, mpnw1)
call mpinitwds (tn1, mpnw1)
call mpinitwds (td2, mpnw1)
call mpinitwds (tn2, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t41, mpnw1)
call mpinitwds (t42, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (rra, mpnw1)
call mpinitwds (rr2, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
!+ call mpdmc (2.d0, 0, tc1, 4
call mpfrsetprec (tc1(1:), 4)
call mpfrsetd (tc1(1:), 2.d0, mprnd)
 
!+ call mpnpwr (tc1, -mpnw*mpnbt, eps, 4
call mpfrsetprec (eps(1:), 4)
call mpfrpowsi (eps(1:), tc1(1:), -mpnw*mpnbt, mprnd)
 
!+ call mpmdc (rr, d1, n1, 4
d1 = 2.d0 * mpfrgetd2exp (ix8, rr(1:), mprnd)
if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
d1 = abs (d1) * 2.d0 ** n1

if (d1 < dfrac * mpnw1 * mpdpw) then
  mpnw1 = min (mpnw + nint (d1 / (dfrac * mpdpw)), 2 * mpnw + 1, mpnwx)
!+   call mpabs (rr, rra, mpnw1
  call mpfrsetprec (rra(1:), mpnw1)
  call mpfrabs (rra(1:), rr(1:), mprnd)
 
!+   call mpdmc (1.d0, 0, tn1, mpnw1
  call mpfrsetprec (tn1(1:), mpnw1)
  call mpfrsetd (tn1(1:), 1.d0, mprnd)
 
!+   call mpdmc (1.d0, 0, f1, mpnw1
  call mpfrsetprec (f1(1:), mpnw1)
  call mpfrsetd (f1(1:), 1.d0, mprnd)
 
!+   call mpdmc (1.d0, 0, f2, mpnw1
  call mpfrsetprec (f2(1:), mpnw1)
  call mpfrsetd (f2(1:), 1.d0, mprnd)
 
!+   call mpmul (rra, rra, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmul (t2(1:), rra(1:), rra(1:), mprnd)
 
!+   call mpmuld (t2, 0.25d0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmuld (t1(1:), t2(1:), 0.25d0, mprnd)
 

  do k = 1, nua
!+     call mpmuld (f2, dble (k), t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmuld (t2(1:), f2(1:), dble (k), mprnd)
 
!+     call mpeq (t2, f2, mpnw1
    call mpfrsetprec (f2(1:), mpnw1)
    call mpfrset (f2(1:), t2(1:), mprnd)
 
  enddo

!+   call mpmul (f1, f2, td1, mpnw1
  call mpfrsetprec (td1(1:), mpnw1)
  call mpfrmul (td1(1:), f1(1:), f2(1:), mprnd)
 
!+   call mpdiv (tn1, td1, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrdiv (t2(1:), tn1(1:), td1(1:), mprnd)
 
!+   call mpeq (t2, sum1, mpnw1
  call mpfrsetprec (sum1(1:), mpnw1)
  call mpfrset (sum1(1:), t2(1:), mprnd)
 

  do k = 1, itrmax
!+     call mpmuld (f1, dble (k), t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmuld (t2(1:), f1(1:), dble (k), mprnd)
 
!+     call mpeq (t2, f1, mpnw1
    call mpfrsetprec (f1(1:), mpnw1)
    call mpfrset (f1(1:), t2(1:), mprnd)
 
!+     call mpmuld (f2, dble (k + nua), t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmuld (t2(1:), f2(1:), dble (k + nua), mprnd)
 
!+     call mpeq (t2, f2, mpnw1
    call mpfrsetprec (f2(1:), mpnw1)
    call mpfrset (f2(1:), t2(1:), mprnd)
 
!+     call mpmul (t1, tn1, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmul (t2(1:), t1(1:), tn1(1:), mprnd)
 
!+     call mpneg (t2, tn1, mpnw1
    call mpfrsetprec (tn1(1:), mpnw1)
    call mpfrneg (tn1(1:), t2(1:), mprnd)
 
!+     call mpmul (f1, f2, td1, mpnw1
    call mpfrsetprec (td1(1:), mpnw1)
    call mpfrmul (td1(1:), f1(1:), f2(1:), mprnd)
 
!+     call mpdiv (tn1, td1, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrdiv (t2(1:), tn1(1:), td1(1:), mprnd)
 
!+     call mpadd (sum1, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfradd (t3(1:), sum1(1:), t2(1:), mprnd)
 
!+     call mpeq (t3, sum1, mpnw1
    call mpfrsetprec (sum1(1:), mpnw1)
    call mpfrset (sum1(1:), t3(1:), mprnd)
 

!+     call mpabs (t2, tc1, 4
    call mpfrsetprec (tc1(1:), 4)
    call mpfrabs (tc1(1:), t2(1:), mprnd)
 
!+     call mpmul (eps, sum1, tc3, 4
    call mpfrsetprec (tc3(1:), 4)
    call mpfrmul (tc3(1:), eps(1:), sum1(1:), mprnd)
 
!+     call mpabs (tc3, tc2, 4
    call mpfrsetprec (tc2(1:), 4)
    call mpfrabs (tc2(1:), tc3(1:), mprnd)
 
!+     call mpcpr (tc1, tc2, ic1, 4
    ic1 = mpfrcmp (tc1(1:), tc2(1:))
 
    if (ic1 <= 0) goto 100
  enddo

  write (mpldb, 4)
4 format ('*** MPBESSELJNR: Loop end error 1')
  call mpabrt ( 517)

100 continue

!+   call mpmuld (rra, 0.5d0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmuld (t1(1:), rra(1:), 0.5d0, mprnd)
 
!+   call mpnpwr (t1, nua, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrpowsi (t2(1:), t1(1:), nua, mprnd)
 
!+   call mpmul (sum1, t2, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrmul (t3(1:), sum1(1:), t2(1:), mprnd)
 
else
! xa2 = xa**2
! t1 = mpreal (4.d0 * dble (nua)**2, mpnw)
! tn1 = mpreal (1.d0, mpnw)
! tn2 = (t1 - 1.d0) / 8.d0
! td1 = mpreal (1.d0, mpnw)
! td2 = xa
! sum1 = tn1 / td1
! sum2 = tn2 / t32

  mpnw1 = min (mpnw + 1, mpnwx)
!+   call mpabs (rr, rra, mpnw1
  call mpfrsetprec (rra(1:), mpnw1)
  call mpfrabs (rra(1:), rr(1:), mprnd)
 
!+   call mpmul (rra, rra, rr2, mpnw1
  call mpfrsetprec (rr2(1:), mpnw1)
  call mpfrmul (rr2(1:), rra(1:), rra(1:), mprnd)
 
  d1 = 4.d0 * dble (nua)**2
!+   call mpdmc (d1, 0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrsetd (t1(1:), d1, mprnd)
 
!+   call mpdmc (1.d0, 0, tn1, mpnw1
  call mpfrsetprec (tn1(1:), mpnw1)
  call mpfrsetd (tn1(1:), 1.d0, mprnd)
 
!+   call mpsub (t1, tn1, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrsub (t2(1:), t1(1:), tn1(1:), mprnd)
 
!+   call mpdivd (t2, 8.d0, tn2, mpnw1
  call mpfrsetprec (tn2(1:), mpnw1)
  call mpfrdivd (tn2(1:), t2(1:), 8.d0, mprnd)
 
!+   call mpdmc (1.d0, 0, td1, mpnw1
  call mpfrsetprec (td1(1:), mpnw1)
  call mpfrsetd (td1(1:), 1.d0, mprnd)
 
!+   call mpeq (rra, td2, mpnw1
  call mpfrsetprec (td2(1:), mpnw1)
  call mpfrset (td2(1:), rra(1:), mprnd)
 
!+   call mpdiv (tn1, td1, sum1, mpnw1
  call mpfrsetprec (sum1(1:), mpnw1)
  call mpfrdiv (sum1(1:), tn1(1:), td1(1:), mprnd)
 
!+   call mpdiv (tn2, td2, sum2, mpnw1
  call mpfrsetprec (sum2(1:), mpnw1)
  call mpfrdiv (sum2(1:), tn2(1:), td2(1:), mprnd)
 

  do k = 1, itrmax
!   tn1 = -tn1 * (t1 - (2.d0*(2.d0*k-1.d0) - 1.d0)**2) * (t1 - (2.d0*(2.d0*k) - 1.d0)**2)

    d1 = (2.d0*(2.d0*k-1.d0) - 1.d0)**2
    d2 = (2.d0*(2.d0*k) - 1.d0)**2
!+     call mpdmc (d1, 0, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrsetd (t2(1:), d1, mprnd)
 
!+     call mpsub (t1, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrsub (t3(1:), t1(1:), t2(1:), mprnd)
 
!+     call mpdmc (d2, 0, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrsetd (t2(1:), d2, mprnd)
 
!+     call mpsub (t1, t2, t5, mpnw1
    call mpfrsetprec (t5(1:), mpnw1)
    call mpfrsub (t5(1:), t1(1:), t2(1:), mprnd)
 
!+     call mpmul (t3, t5, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmul (t2(1:), t3(1:), t5(1:), mprnd)
 
!+     call mpmul (tn1, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrmul (t3(1:), tn1(1:), t2(1:), mprnd)
 
!+     call mpneg (t3, tn1, mpnw1
    call mpfrsetprec (tn1(1:), mpnw1)
    call mpfrneg (tn1(1:), t3(1:), mprnd)
 

!   td1 = td1 * dble (2*k-1) * dble (2*k) * 64.d0 * xa2

    d1 = dble (2*k-1) * dble (2*k) * 64.d0
!+     call mpmuld (td1, d1, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmuld (t2(1:), td1(1:), d1, mprnd)
 
!+     call mpmul (t2, rr2, td1, mpnw1
    call mpfrsetprec (td1(1:), mpnw1)
    call mpfrmul (td1(1:), t2(1:), rr2(1:), mprnd)
 

!   t41 = tn1 / td1
!   sum1 = sum1 + t41

!+     call mpdiv (tn1, td1, t41, mpnw1
    call mpfrsetprec (t41(1:), mpnw1)
    call mpfrdiv (t41(1:), tn1(1:), td1(1:), mprnd)
 
!+     call mpadd (sum1, t41, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfradd (t2(1:), sum1(1:), t41(1:), mprnd)
 
!+     call mpeq (t2, sum1, mpnw1
    call mpfrsetprec (sum1(1:), mpnw1)
    call mpfrset (sum1(1:), t2(1:), mprnd)
 

!   tn2 = -tn2 * (t1 - (2.d0*(2.d0*k) - 1.d0)**2) * (t1 - (2.d0*(2.d0*k+1.d0) - 1.d0)**2)

    d1 = (2.d0*(2.d0*k) - 1.d0)**2
    d2 = (2.d0*(2.d0*k+1.d0) - 1.d0)**2
!+     call mpdmc (d1, 0, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrsetd (t2(1:), d1, mprnd)
 
!+     call mpsub (t1, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrsub (t3(1:), t1(1:), t2(1:), mprnd)
 
!+     call mpdmc (d2, 0, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrsetd (t2(1:), d2, mprnd)
 
!+     call mpsub (t1, t2, t5, mpnw1
    call mpfrsetprec (t5(1:), mpnw1)
    call mpfrsub (t5(1:), t1(1:), t2(1:), mprnd)
 
!+     call mpmul (t3, t5, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmul (t2(1:), t3(1:), t5(1:), mprnd)
 
!+     call mpmul (tn2, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrmul (t3(1:), tn2(1:), t2(1:), mprnd)
 
!+     call mpneg (t3, tn2, mpnw1
    call mpfrsetprec (tn2(1:), mpnw1)
    call mpfrneg (tn2(1:), t3(1:), mprnd)
 

!   td2 = td2 * dble (2*k) * dble (2*k+1) * 64.d0 * xa2

    d1 = dble (2*k) * dble (2*k+1) * 64.d0
!+     call mpmuld (td2, d1, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmuld (t2(1:), td2(1:), d1, mprnd)
 
!+     call mpmul (t2, rr2, td2, mpnw1
    call mpfrsetprec (td2(1:), mpnw1)
    call mpfrmul (td2(1:), t2(1:), rr2(1:), mprnd)
 

!  t42 = tn2 / td2
!  sum2 = sum2 + t42

!+     call mpdiv (tn2, td2, t42, mpnw1
    call mpfrsetprec (t42(1:), mpnw1)
    call mpfrdiv (t42(1:), tn2(1:), td2(1:), mprnd)
 
!+     call mpadd (sum2, t42, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfradd (t2(1:), sum2(1:), t42(1:), mprnd)
 
!+     call mpeq (t2, sum2, mpnw1
    call mpfrsetprec (sum2(1:), mpnw1)
    call mpfrset (sum2(1:), t2(1:), mprnd)
 

!  if (abs (t41) / abs (sum1) < eps .and. abs (t42) / abs (sum2) < eps ) goto 110

!+     call mpabs (t41, tc1, 4
    call mpfrsetprec (tc1(1:), 4)
    call mpfrabs (tc1(1:), t41(1:), mprnd)
 
!+     call mpmul (eps, sum1, tc3, 4
    call mpfrsetprec (tc3(1:), 4)
    call mpfrmul (tc3(1:), eps(1:), sum1(1:), mprnd)
 
!+     call mpabs (tc3, tc2, 4
    call mpfrsetprec (tc2(1:), 4)
    call mpfrabs (tc2(1:), tc3(1:), mprnd)
 
!+     call mpcpr (tc1, tc2, ic1, 4
    ic1 = mpfrcmp (tc1(1:), tc2(1:))
 
!+     call mpabs (t42, tc1, 4
    call mpfrsetprec (tc1(1:), 4)
    call mpfrabs (tc1(1:), t42(1:), mprnd)
 
!+     call mpmul (eps, sum2, tc3, 4
    call mpfrsetprec (tc3(1:), 4)
    call mpfrmul (tc3(1:), eps(1:), sum2(1:), mprnd)
 
!+     call mpabs (tc3, tc2, 4
    call mpfrsetprec (tc2(1:), 4)
    call mpfrabs (tc2(1:), tc3(1:), mprnd)
 
!+     call mpcpr (tc1, tc2, ic2, 4
    ic2 = mpfrcmp (tc1(1:), tc2(1:))
 
    if (ic1 <= 0 .and. ic2 <= 0) goto 110
  enddo

  write (mpldb, 5)
5 format ('*** MPBESSELJNR: Loop end error 2')
  call mpabrt ( 518)

110 continue

! t1 = xa - 0.5d0 * nua * pi - 0.25d0 * pi
! besselj = sqrt (2.d0 / (pi * xa)) * (cos (t1) * sum1 - sin (t1) * sum2)

!+   call mpmuld (mppicon, 0.5d0 * nua, t1, mpnw1
call mpfixlocr (mppicon)
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmuld (t1(1:), mppicon(1:), 0.5d0 * nua, mprnd)
 
!+   call mpsub (rra, t1, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrsub (t2(1:), rra(1:), t1(1:), mprnd)
 
!+   call mpmuld (mppicon, 0.25d0, t1, mpnw1
call mpfixlocr (mppicon)
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmuld (t1(1:), mppicon(1:), 0.25d0, mprnd)
 
!+   call mpsub (t2, t1, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrsub (t3(1:), t2(1:), t1(1:), mprnd)
 
!+   call mpcssnr (t3, t41, t42, mpnw1
  call mpfrsetprec (t42(1:), mpnw1)
  call mpfrsincos (t42(1:), t41(1:), t3(1:), mprnd)
 
!+   call mpmul (t41, sum1, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmul (t1(1:), t41(1:), sum1(1:), mprnd)
 
!+   call mpmul (t42, sum2, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmul (t2(1:), t42(1:), sum2(1:), mprnd)
 
!+   call mpsub (t1, t2, t5, mpnw1
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfrsub (t5(1:), t1(1:), t2(1:), mprnd)
 
!+   call mpmul (mppicon, rra, t1, mpnw1
call mpfixlocr (mppicon)
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmul (t1(1:), mppicon(1:), rra(1:), mprnd)
 
!+   call mpdmc (2.d0, 0, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrsetd (t2(1:), 2.d0, mprnd)
 
!+   call mpdiv (t2, t1, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrdiv (t3(1:), t2(1:), t1(1:), mprnd)
 
!+   call mpsqrt (t3, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrsqrt (t1(1:), t3(1:), mprnd)
 
!+   call mpmul (t1, t5, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrmul (t3(1:), t1(1:), t5(1:), mprnd)
 
endif

if (mod (nu, 2) /= 0) then
!  if (nu < 0 .and. x > 0.d0 .or. nu > 0 .and. x < 0.d0) besselj = - besselj

!+   ic1 = mpsgn (rr
  ic1 = mpfrsgn (rr(1:))
 
  if (nu < 0 .and. ic1 > 0 .or. nu > 0 .and. ic1 < 0) then
!+     call mpneg (t3, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrneg (t2(1:), t3(1:), mprnd)
 
!+     call mpeq (t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrset (t3(1:), t2(1:), mprnd)
 
  endif
endif

!+ call mproun (t3, mpnw
call mpfrprecround (t3(1:), mpnw, mprnd)
 
!+ call mpeq (t3, ss, mpnw
call mpfrsetprec (ss(1:), mpnw)
call mpfrset (ss(1:), t3(1:), mprnd)
 

return
end subroutine mpbesseljnr

subroutine mpbesseljr (qq, rr, ss, mpnw)

!   This evaluates the modified Bessel function BesselJ (QQ,RR) for QQ and RR
!   both MPR. The algorithm is DLMF formula 10.2.2 for modest RR,
!   and DLMF 10.17.3 for large RR, relative to precision.

implicit none
integer (mpiknd), intent(in):: qq(0:), rr(0:)
integer (mpiknd), intent(out):: ss(0:)
integer, intent(in):: mpnw
integer ic1, ic2, i1, i2, k, mpnw1, n1
real (mprknd) d1, d2
integer, parameter:: itrmax = 1000000
real (mprknd), parameter:: dfrac = 1.5d0, pi = 3.1415926535897932385d0
integer (mpiknd) f1(0:2*mpnw+6), f2(0:2*mpnw+6), sum1(0:2*mpnw+6), &
  sum2(0:2*mpnw+6), td1(0:2*mpnw+6), td2(0:2*mpnw+6), tn1(0:2*mpnw+6), &
  tn2(0:2*mpnw+6), t1(0:2*mpnw+6), t2(0:2*mpnw+6), t3(0:2*mpnw+6), &
  t4(0:2*mpnw+6), t41(0:2*mpnw+6), t42(0:2*mpnw+6), t5(0:2*mpnw+6), &
  rra(0:2*mpnw+6), rr2(0:2*mpnw+6), tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)

integer, external:: mpfrsgn
integer, external:: mpfrcmp
real (mprknd), external:: mpfrgetd2exp
integer (mpiknd) ix8
 
!  End of declaration

if (mpnw < 4 .or. mpspacer (rr) < mpnw + 4 .or. mpspacer (ss) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPBESSELJR: Uninitialized or inadequately sized arrays')
  call mpabrt ( 519)
endif

!  Check if PI has been precomputed.

!+ call mpmdc (mppicon, d1, n1, mpnw+1
call mpfixlocr (mppicon)
d1 = 2.d0 * mpfrgetd2exp (ix8, mppicon(1:), mprnd)
if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
if (n1 /= 1 .or. abs (d1 * 2.d0**n1 - pi) > mprdfz &
  .or. mpwprecr (mppicon) < mpnw+1) then
  write (mpldb, 2) mpnw+1
2 format ('*** MPBESSELJR: Pi must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 520)
endif

!  End of initial input array check

mpnw1 = min (2 * mpnw + 1, mpnwx)
call mpinitwds (f1, mpnw1)
call mpinitwds (f2, mpnw1)
call mpinitwds (sum1, mpnw1)
call mpinitwds (sum2, mpnw1)
call mpinitwds (td1, mpnw1)
call mpinitwds (tn1, mpnw1)
call mpinitwds (td2, mpnw1)
call mpinitwds (tn2, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (t41, mpnw1)
call mpinitwds (t42, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (rra, mpnw1)
call mpinitwds (rr2, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
!+ call mpdmc (2.d0, 0, tc1, 4
call mpfrsetprec (tc1(1:), 4)
call mpfrsetd (tc1(1:), 2.d0, mprnd)
 
!+ call mpnpwr (tc1, -mpnw*mpnbt, eps, 4
call mpfrsetprec (eps(1:), 4)
call mpfrpowsi (eps(1:), tc1(1:), -mpnw*mpnbt, mprnd)
 

!   If QQ is integer, call mpbesseljnr; if RR <= 0, then error.

!+ call mpinfr (qq, t1, t2, mpnw
call mpfrsetprec (t1(1:), mpnw)
call mpfrsetprec (t2(1:), mpnw)
call mpfrtrunc (t1(1:), qq(1:))
call mpfrsub (t2(1:), qq(1:), t1(1:), mprnd)
 
!+ i1 = mpsgn (qq
i1 = mpfrsgn (qq(1:))
 
!+ i2 = mpsgn (rr
i2 = mpfrsgn (rr(1:))
 
!+ if (mpsgn (t2) == 0) the
if (mpfrsgn (t2(1:)) == 0) then
 
!+   call mpmdc (qq, d1, n1, mpnw
  d1 = 2.d0 * mpfrgetd2exp (ix8, qq(1:), mprnd)
  if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
  d1 = d1 * 2.d0**n1
  n1 = nint (d1)
  call mpbesseljnr (n1, rr, t3, mpnw)
  goto 120
elseif (i2 <= 0) then
  write (mpldb, 3)
3 format ('*** MPBESSELJR: Second argument <= 0')
  call mpabrt ( 521)
endif

!+ call mpmdc (rr, d1, n1, 4
d1 = 2.d0 * mpfrgetd2exp (ix8, rr(1:), mprnd)
if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
d1 = abs (d1) * 2.d0 ** n1

if (d1 < dfrac * mpnw1 * mpdpw) then
  mpnw1 = min (mpnw + nint (d1 / (dfrac * mpdpw)), 2 * mpnw + 1, mpnwx)
!+   call mpabs (rr, rra, mpnw1
  call mpfrsetprec (rra(1:), mpnw1)
  call mpfrabs (rra(1:), rr(1:), mprnd)
 
!+   call mpdmc (1.d0, 0, tn1, mpnw1
  call mpfrsetprec (tn1(1:), mpnw1)
  call mpfrsetd (tn1(1:), 1.d0, mprnd)
 
!+   call mpdmc (1.d0, 0, f1, mpnw1
  call mpfrsetprec (f1(1:), mpnw1)
  call mpfrsetd (f1(1:), 1.d0, mprnd)
 
!+   call mpadd (qq, f1, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfradd (t1(1:), qq(1:), f1(1:), mprnd)
 
!+   call mpgammar (t1, f2, mpnw1
  call mpfrsetprec (f2(1:), mpnw1)
  call mpfrgamma (f2(1:), t1(1:), mprnd)
 
!+   call mpmul (rra, rra, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmul (t2(1:), rra(1:), rra(1:), mprnd)
 
!+   call mpmuld (t2, 0.25d0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmuld (t1(1:), t2(1:), 0.25d0, mprnd)
 

!+   call mpmul (f1, f2, td1, mpnw1
  call mpfrsetprec (td1(1:), mpnw1)
  call mpfrmul (td1(1:), f1(1:), f2(1:), mprnd)
 
!+   call mpdiv (tn1, td1, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrdiv (t2(1:), tn1(1:), td1(1:), mprnd)
 
!+   call mpeq (t2, sum1, mpnw1
  call mpfrsetprec (sum1(1:), mpnw1)
  call mpfrset (sum1(1:), t2(1:), mprnd)
 

  do k = 1, itrmax
!+     call mpmuld (f1, dble (k), t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmuld (t2(1:), f1(1:), dble (k), mprnd)
 
!+     call mpeq (t2, f1, mpnw1
    call mpfrsetprec (f1(1:), mpnw1)
    call mpfrset (f1(1:), t2(1:), mprnd)
 
!+     call mpdmc (dble (k), 0, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrsetd (t3(1:), dble (k), mprnd)
 
!+     call mpadd (qq, t3, t4, mpnw1
    call mpfrsetprec (t4(1:), mpnw1)
    call mpfradd (t4(1:), qq(1:), t3(1:), mprnd)
 
!+     call mpmul (f2, t4, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrmul (t3(1:), f2(1:), t4(1:), mprnd)
 
!+     call mpeq (t3, f2, mpnw1
    call mpfrsetprec (f2(1:), mpnw1)
    call mpfrset (f2(1:), t3(1:), mprnd)
 
!+     call mpmul (t1, tn1, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmul (t2(1:), t1(1:), tn1(1:), mprnd)
 
!+     call mpneg (t2, tn1, mpnw1
    call mpfrsetprec (tn1(1:), mpnw1)
    call mpfrneg (tn1(1:), t2(1:), mprnd)
 
!+     call mpmul (f1, f2, td1, mpnw1
    call mpfrsetprec (td1(1:), mpnw1)
    call mpfrmul (td1(1:), f1(1:), f2(1:), mprnd)
 
!+     call mpdiv (tn1, td1, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrdiv (t2(1:), tn1(1:), td1(1:), mprnd)
 
!+     call mpadd (sum1, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfradd (t3(1:), sum1(1:), t2(1:), mprnd)
 
!+     call mpeq (t3, sum1, mpnw1
    call mpfrsetprec (sum1(1:), mpnw1)
    call mpfrset (sum1(1:), t3(1:), mprnd)
 

!+     call mpabs (t2, tc1, 4
    call mpfrsetprec (tc1(1:), 4)
    call mpfrabs (tc1(1:), t2(1:), mprnd)
 
!+     call mpmul (eps, sum1, tc3, 4
    call mpfrsetprec (tc3(1:), 4)
    call mpfrmul (tc3(1:), eps(1:), sum1(1:), mprnd)
 
!+     call mpabs (tc3, tc2, 4
    call mpfrsetprec (tc2(1:), 4)
    call mpfrabs (tc2(1:), tc3(1:), mprnd)
 
!+     call mpcpr (tc1, tc2, ic1, 4
    ic1 = mpfrcmp (tc1(1:), tc2(1:))
 
    if (ic1 <= 0) goto 100
  enddo

  write (mpldb, 4)
4 format ('*** MPBESSELJR: Loop end error 1')
  call mpabrt ( 522)

100 continue

!+   call mpmuld (rr, 0.5d0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmuld (t1(1:), rr(1:), 0.5d0, mprnd)
 
!+   call mppower (t1, qq, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrpow (t2(1:), t1(1:), qq(1:), mprnd)
 
!+   call mpmul (sum1, t2, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrmul (t3(1:), sum1(1:), t2(1:), mprnd)
 
else
! xa2 = xa**2
! t1 = mpreal (4.d0 * dble (nua)**2, mpnw)
! tn1 = mpreal (1.d0, mpnw)
! tn2 = (t1 - 1.d0) / 8.d0
! td1 = mpreal (1.d0, mpnw)
! td2 = xa
! sum1 = tn1 / td1
! sum2 = tn2 / t32

  mpnw1 = min (mpnw + 1, mpnwx)
!+   call mpabs (rr, rra, mpnw1
  call mpfrsetprec (rra(1:), mpnw1)
  call mpfrabs (rra(1:), rr(1:), mprnd)
 
!+   call mpmul (rra, rra, rr2, mpnw1
  call mpfrsetprec (rr2(1:), mpnw1)
  call mpfrmul (rr2(1:), rra(1:), rra(1:), mprnd)
 
!+   call mpmul (qq, qq, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmul (t2(1:), qq(1:), qq(1:), mprnd)
 
!+   call mpmuld (t2, 4.d0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmuld (t1(1:), t2(1:), 4.d0, mprnd)
 
!+   call mpdmc (1.d0, 0, tn1, mpnw1
  call mpfrsetprec (tn1(1:), mpnw1)
  call mpfrsetd (tn1(1:), 1.d0, mprnd)
 
!+   call mpsub (t1, tn1, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrsub (t2(1:), t1(1:), tn1(1:), mprnd)
 
!+   call mpdivd (t2, 8.d0, tn2, mpnw1
  call mpfrsetprec (tn2(1:), mpnw1)
  call mpfrdivd (tn2(1:), t2(1:), 8.d0, mprnd)
 
!+   call mpdmc (1.d0, 0, td1, mpnw1
  call mpfrsetprec (td1(1:), mpnw1)
  call mpfrsetd (td1(1:), 1.d0, mprnd)
 
!+   call mpeq (rra, td2, mpnw1
  call mpfrsetprec (td2(1:), mpnw1)
  call mpfrset (td2(1:), rra(1:), mprnd)
 
!+   call mpdiv (tn1, td1, sum1, mpnw1
  call mpfrsetprec (sum1(1:), mpnw1)
  call mpfrdiv (sum1(1:), tn1(1:), td1(1:), mprnd)
 
!+   call mpdiv (tn2, td2, sum2, mpnw1
  call mpfrsetprec (sum2(1:), mpnw1)
  call mpfrdiv (sum2(1:), tn2(1:), td2(1:), mprnd)
 

  do k = 1, itrmax
!   tn1 = -tn1 * (t1 - (2.d0*(2.d0*k-1.d0) - 1.d0)**2) * (t1 - (2.d0*(2.d0*k) - 1.d0)**2)

    d1 = (2.d0*(2.d0*k-1.d0) - 1.d0)**2
    d2 = (2.d0*(2.d0*k) - 1.d0)**2
!+     call mpdmc (d1, 0, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrsetd (t2(1:), d1, mprnd)
 
!+     call mpsub (t1, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrsub (t3(1:), t1(1:), t2(1:), mprnd)
 
!+     call mpdmc (d2, 0, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrsetd (t2(1:), d2, mprnd)
 
!+     call mpsub (t1, t2, t5, mpnw1
    call mpfrsetprec (t5(1:), mpnw1)
    call mpfrsub (t5(1:), t1(1:), t2(1:), mprnd)
 
!+     call mpmul (t3, t5, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmul (t2(1:), t3(1:), t5(1:), mprnd)
 
!+     call mpmul (tn1, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrmul (t3(1:), tn1(1:), t2(1:), mprnd)
 
!+     call mpneg (t3, tn1, mpnw1
    call mpfrsetprec (tn1(1:), mpnw1)
    call mpfrneg (tn1(1:), t3(1:), mprnd)
 

!   td1 = td1 * dble (2*k-1) * dble (2*k) * 64.d0 * xa2

    d1 = dble (2*k-1) * dble (2*k) * 64.d0
!+     call mpmuld (td1, d1, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmuld (t2(1:), td1(1:), d1, mprnd)
 
!+     call mpmul (t2, rr2, td1, mpnw1
    call mpfrsetprec (td1(1:), mpnw1)
    call mpfrmul (td1(1:), t2(1:), rr2(1:), mprnd)
 

!   t41 = tn1 / td1
!   sum1 = sum1 + t41

!+     call mpdiv (tn1, td1, t41, mpnw1
    call mpfrsetprec (t41(1:), mpnw1)
    call mpfrdiv (t41(1:), tn1(1:), td1(1:), mprnd)
 
!+     call mpadd (sum1, t41, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfradd (t2(1:), sum1(1:), t41(1:), mprnd)
 
!+     call mpeq (t2, sum1, mpnw1
    call mpfrsetprec (sum1(1:), mpnw1)
    call mpfrset (sum1(1:), t2(1:), mprnd)
 

!   tn2 = -tn2 * (t1 - (2.d0*(2.d0*k) - 1.d0)**2) * (t1 - (2.d0*(2.d0*k+1.d0) - 1.d0)**2)

    d1 = (2.d0*(2.d0*k) - 1.d0)**2
    d2 = (2.d0*(2.d0*k+1.d0) - 1.d0)**2
!+     call mpdmc (d1, 0, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrsetd (t2(1:), d1, mprnd)
 
!+     call mpsub (t1, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrsub (t3(1:), t1(1:), t2(1:), mprnd)
 
!+     call mpdmc (d2, 0, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrsetd (t2(1:), d2, mprnd)
 
!+     call mpsub (t1, t2, t5, mpnw1
    call mpfrsetprec (t5(1:), mpnw1)
    call mpfrsub (t5(1:), t1(1:), t2(1:), mprnd)
 
!+     call mpmul (t3, t5, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmul (t2(1:), t3(1:), t5(1:), mprnd)
 
!+     call mpmul (tn2, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrmul (t3(1:), tn2(1:), t2(1:), mprnd)
 
!+     call mpneg (t3, tn2, mpnw1
    call mpfrsetprec (tn2(1:), mpnw1)
    call mpfrneg (tn2(1:), t3(1:), mprnd)
 

!   td2 = td2 * dble (2*k) * dble (2*k+1) * 64.d0 * xa2

    d1 = dble (2*k) * dble (2*k+1) * 64.d0
!+     call mpmuld (td2, d1, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmuld (t2(1:), td2(1:), d1, mprnd)
 
!+     call mpmul (t2, rr2, td2, mpnw1
    call mpfrsetprec (td2(1:), mpnw1)
    call mpfrmul (td2(1:), t2(1:), rr2(1:), mprnd)
 

!  t42 = tn2 / td2
!  sum2 = sum2 + t42

!+     call mpdiv (tn2, td2, t42, mpnw1
    call mpfrsetprec (t42(1:), mpnw1)
    call mpfrdiv (t42(1:), tn2(1:), td2(1:), mprnd)
 
!+     call mpadd (sum2, t42, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfradd (t2(1:), sum2(1:), t42(1:), mprnd)
 
!+     call mpeq (t2, sum2, mpnw1
    call mpfrsetprec (sum2(1:), mpnw1)
    call mpfrset (sum2(1:), t2(1:), mprnd)
 

!  if (abs (t41) / abs (sum1) < eps .and. abs (t42) / abs (sum2) < eps ) goto 110

!+     call mpabs (t41, tc1, 4
    call mpfrsetprec (tc1(1:), 4)
    call mpfrabs (tc1(1:), t41(1:), mprnd)
 
!+     call mpmul (eps, sum1, tc3, 4
    call mpfrsetprec (tc3(1:), 4)
    call mpfrmul (tc3(1:), eps(1:), sum1(1:), mprnd)
 
!+     call mpabs (tc3, tc2, 4
    call mpfrsetprec (tc2(1:), 4)
    call mpfrabs (tc2(1:), tc3(1:), mprnd)
 
!+     call mpcpr (tc1, tc2, ic1, 4
    ic1 = mpfrcmp (tc1(1:), tc2(1:))
 
!+     call mpabs (t42, tc1, 4
    call mpfrsetprec (tc1(1:), 4)
    call mpfrabs (tc1(1:), t42(1:), mprnd)
 
!+     call mpmul (eps, sum2, tc3, 4
    call mpfrsetprec (tc3(1:), 4)
    call mpfrmul (tc3(1:), eps(1:), sum2(1:), mprnd)
 
!+     call mpabs (tc3, tc2, 4
    call mpfrsetprec (tc2(1:), 4)
    call mpfrabs (tc2(1:), tc3(1:), mprnd)
 
!+     call mpcpr (tc1, tc2, ic2, 4
    ic2 = mpfrcmp (tc1(1:), tc2(1:))
 
    if (ic1 <= 0 .and. ic2 <= 0) goto 110
  enddo

  write (mpldb, 5)
5 format ('*** MPBESSELJR: Loop end error 2')
  call mpabrt ( 523)

110 continue

! t1 = xa - 0.5d0 * nua * pi - 0.25d0 * pi
! besselj = sqrt (2.d0 / (pi * xa)) * (cos (t1) * sum1 - sin (t1) * sum2)

!   call mpmuld (mppicon, 0.5d0 * nua, t1, mpnw1)

!+   call mpmul (mppicon, qq, t2, mpnw1
call mpfixlocr (mppicon)
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmul (t2(1:), mppicon(1:), qq(1:), mprnd)
 
!+   call mpmuld (t2, 0.5d0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmuld (t1(1:), t2(1:), 0.5d0, mprnd)
 
!+   call mpsub (rra, t1, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrsub (t2(1:), rra(1:), t1(1:), mprnd)
 
!+   call mpmuld (mppicon, 0.25d0, t1, mpnw1
call mpfixlocr (mppicon)
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmuld (t1(1:), mppicon(1:), 0.25d0, mprnd)
 
!+   call mpsub (t2, t1, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrsub (t3(1:), t2(1:), t1(1:), mprnd)
 
!+   call mpcssnr (t3, t41, t42, mpnw1
  call mpfrsetprec (t42(1:), mpnw1)
  call mpfrsincos (t42(1:), t41(1:), t3(1:), mprnd)
 
!+   call mpmul (t41, sum1, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmul (t1(1:), t41(1:), sum1(1:), mprnd)
 
!+   call mpmul (t42, sum2, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmul (t2(1:), t42(1:), sum2(1:), mprnd)
 
!+   call mpsub (t1, t2, t5, mpnw1
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfrsub (t5(1:), t1(1:), t2(1:), mprnd)
 
!+   call mpmul (mppicon, rra, t1, mpnw1
call mpfixlocr (mppicon)
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmul (t1(1:), mppicon(1:), rra(1:), mprnd)
 
!+   call mpdmc (2.d0, 0, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrsetd (t2(1:), 2.d0, mprnd)
 
!+   call mpdiv (t2, t1, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrdiv (t3(1:), t2(1:), t1(1:), mprnd)
 
!+   call mpsqrt (t3, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrsqrt (t1(1:), t3(1:), mprnd)
 
!+   call mpmul (t1, t5, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrmul (t3(1:), t1(1:), t5(1:), mprnd)
 
endif

! if (mod (nu, 2) /= 0) then
!  if (nu < 0 .and. x > 0.d0 .or. nu > 0 .and. x < 0.d0) besselj = - besselj

!   ic1 = mpsgn (rr)
!   if (nu < 0 .and. ic1 > 0 .or. nu > 0 .and. ic1 < 0) then
!     call mpneg (t3, t2, mpnw1)
!     call mpeq (t2, t3, mpnw1)
!   endif
! endif

120 continue

!+ call mproun (t3, mpnw
call mpfrprecround (t3(1:), mpnw, mprnd)
 
!+ call mpeq (t3, ss, mpnw
call mpfrsetprec (ss(1:), mpnw)
call mpfrset (ss(1:), t3(1:), mprnd)
 

return
end subroutine mpbesseljr

subroutine mpbesselknr (nu, rr, ss, mpnw)

!   This evaluates the modified Bessel function BesselK (NU,RR).
!   NU is an integer. The algorithm is DLMF formula 10.31.1 for modest RR,
!   and DLMF 10.40.2 for large RR, relative to precision.

implicit none
integer, intent(in):: nu, mpnw
integer (mpiknd), intent(in):: rr(0:)
integer (mpiknd), intent(out):: ss(0:)
integer, parameter:: itrmax = 1000000
integer ic1, k, mpnw1, nua, n1
real (mprknd) d1
real (mprknd), parameter:: dfrac = 1.5d0, egam = 0.5772156649015328606d0, &
  pi = 3.1415926535897932385d0
integer (mpiknd) f1(0:mpnw+6), f2(0:mpnw+6), f3(0:mpnw+6), f4(0:mpnw+6), &
  f5(0:mpnw+6), sum1(0:mpnw+6), sum2(0:mpnw+6), sum3(0:mpnw+6), td(0:mpnw+6), &
  tn(0:mpnw+6), t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), &
  rra(0:mpnw+6), tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)

integer, external:: mpfrsgn
integer, external:: mpfrcmp
real (mprknd), external:: mpfrgetd2exp
integer (mpiknd) ix8
 
!  End of declaration

if (mpnw < 4 .or. mpspacer (rr) < mpnw + 4 .or. mpspacer (ss) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPBESSELKNR: Uninitialized or inadequately sized arrays')
  call mpabrt ( 524)
endif

!   Check if EGAMMA and PI have been precomputed.

!+ call mpmdc (mpegammacon, d1, n1, mpnw+1
call mpfixlocr (mpegammacon)
d1 = 2.d0 * mpfrgetd2exp (ix8, mpegammacon(1:), mprnd)
if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
if (n1 /= -1 .or. abs (d1 * 2.d0**n1 - egam) > mprdfz &
  .or. mpwprecr (mpegammacon) < mpnw+1) then
  write (mpldb, 3) mpnw+1
3 format ('*** MPBESSELKNR: EGAMMA must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 526)
endif
!+ call mpmdc (mppicon, d1, n1, mpnw+1
call mpfixlocr (mppicon)
d1 = 2.d0 * mpfrgetd2exp (ix8, mppicon(1:), mprnd)
if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
if (n1 /= 1 .or. abs (d1 * 2.d0**n1 - pi) > mprdfz &
  .or. mpwprecr (mppicon) < mpnw+1) then
  write (mpldb, 4) mpnw+1
4 format ('*** MPBESSELKNR: Pi must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 527)
endif

!  End of initial input array check

!   Check for RR = 0.

!+ if (mpsgn (rr) == 0) the
if (mpfrsgn (rr(1:)) == 0) then
 
  write (mpldb, 2)
2 format ('*** MPBESSELKNR: Second argument is zero')
  call mpabrt ( 525)
endif

mpnw1 = min (mpnw + 1, mpnwx)
nua = abs (nu)
mpnw1 = min (mpnw + 1, mpnwx)
call mpinitwds (f1, mpnw1)
call mpinitwds (f2, mpnw1)
call mpinitwds (f3, mpnw1)
call mpinitwds (f4, mpnw1)
call mpinitwds (f5, mpnw1)
call mpinitwds (sum1, mpnw1)
call mpinitwds (sum2, mpnw1)
call mpinitwds (sum3, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (rra, mpnw1)
call mpinitwds (td, mpnw1)
call mpinitwds (tn, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
!+ call mpdmc (2.d0, 0, tc1, 4
call mpfrsetprec (tc1(1:), 4)
call mpfrsetd (tc1(1:), 2.d0, mprnd)
 
!+ call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4
call mpfrsetprec (eps(1:), 4)
call mpfrpowsi (eps(1:), tc1(1:), -mpnw1*mpnbt, mprnd)
 
!+ call mpabs (rr, rra, mpnw1
call mpfrsetprec (rra(1:), mpnw1)
call mpfrabs (rra(1:), rr(1:), mprnd)
 
!+ call mpmdc (rra, d1, n1, 4
d1 = 2.d0 * mpfrgetd2exp (ix8, rra(1:), mprnd)
if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
d1 = d1 * 2.d0 ** n1

if (d1 < dfrac * mpnw1 * mpdpw) then
!+   call mpmul (rra, rra, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmul (t2(1:), rra(1:), rra(1:), mprnd)
 
!+   call mpmuld (t2, 0.25d0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmuld (t1(1:), t2(1:), 0.25d0, mprnd)
 
!+   call mpdmc (1.d0, 0, f1, mpnw1
  call mpfrsetprec (f1(1:), mpnw1)
  call mpfrsetd (f1(1:), 1.d0, mprnd)
 
!+   call mpdmc (1.d0, 0, f2, mpnw1
  call mpfrsetprec (f2(1:), mpnw1)
  call mpfrsetd (f2(1:), 1.d0, mprnd)
 
!+   call mpdmc (1.d0, 0, f3, mpnw1
  call mpfrsetprec (f3(1:), mpnw1)
  call mpfrsetd (f3(1:), 1.d0, mprnd)
 
!+   call mpdmc (0.d0, 0,  sum1, mpnw1
  call mpfrsetprec (sum1(1:), mpnw1)
  call mpfrsetd (sum1(1:), 0.d0, mprnd)
 

  do k = 1, nua - 1
!+     call mpmuld (f1, dble (k), t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmuld (t2(1:), f1(1:), dble (k), mprnd)
 
!+     call mpeq (t2, f1, mpnw1
    call mpfrsetprec (f1(1:), mpnw1)
    call mpfrset (f1(1:), t2(1:), mprnd)
 
  enddo

  do k = 0, nua - 1
    if (k > 0) then
!+       call mpdivd (f1, dble (nua - k), t2, mpnw1
      call mpfrsetprec (t2(1:), mpnw1)
      call mpfrdivd (t2(1:), f1(1:), dble (nua - k), mprnd)
 
!+       call mpeq (t2, f1, mpnw1
      call mpfrsetprec (f1(1:), mpnw1)
      call mpfrset (f1(1:), t2(1:), mprnd)
 
!+       call mpmul (t1, f2, t2, mpnw1
      call mpfrsetprec (t2(1:), mpnw1)
      call mpfrmul (t2(1:), t1(1:), f2(1:), mprnd)
 
!+       call mpneg (t2, f2, mpnw1
      call mpfrsetprec (f2(1:), mpnw1)
      call mpfrneg (f2(1:), t2(1:), mprnd)
 
!+       call mpmuld (f3, dble (k), t2, mpnw1
      call mpfrsetprec (t2(1:), mpnw1)
      call mpfrmuld (t2(1:), f3(1:), dble (k), mprnd)
 
!+       call mpeq (t2, f3, mpnw1
      call mpfrsetprec (f3(1:), mpnw1)
      call mpfrset (f3(1:), t2(1:), mprnd)
 
    endif
!+     call mpmul (f1, f2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrmul (t3(1:), f1(1:), f2(1:), mprnd)
 
!+     call mpdiv (t3, f3, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrdiv (t2(1:), t3(1:), f3(1:), mprnd)
 
!+     call mpadd (sum1, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfradd (t3(1:), sum1(1:), t2(1:), mprnd)
 
!+     call mpeq (t3, sum1, mpnw1
    call mpfrsetprec (sum1(1:), mpnw1)
    call mpfrset (sum1(1:), t3(1:), mprnd)
 
  enddo

!+   call mpmuld (sum1, 0.5d0, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmuld (t2(1:), sum1(1:), 0.5d0, mprnd)
 
!+   call mpmuld (rra, 0.5d0, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrmuld (t3(1:), rra(1:), 0.5d0, mprnd)
 
!+   call mpnpwr (t3, nua, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfrpowsi (t4(1:), t3(1:), nua, mprnd)
 
!+   call mpdiv (t2, t4, sum1, mpnw1
  call mpfrsetprec (sum1(1:), mpnw1)
  call mpfrdiv (sum1(1:), t2(1:), t4(1:), mprnd)
 

!+   call mpmuld (rra, 0.5d0, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmuld (t2(1:), rra(1:), 0.5d0, mprnd)
 
!+   call mplog (t2, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrlog (t3(1:), t2(1:), mprnd)
 
  d1 = (-1.d0) ** (nua + 1)
!+   call mpmuld (t3, d1, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmuld (t2(1:), t3(1:), d1, mprnd)
 
  call mpbesselinr (nua, rra, t3, mpnw1)
!+   call mpmul (t2, t3, sum2, mpnw1
  call mpfrsetprec (sum2(1:), mpnw1)
  call mpfrmul (sum2(1:), t2(1:), t3(1:), mprnd)
 

!+   call mpneg (mpegammacon, f1, mpnw1
call mpfixlocr (mpegammacon)
  call mpfrsetprec (f1(1:), mpnw1)
  call mpfrneg (f1(1:), mpegammacon(1:), mprnd)
 
!+   call mpeq (f1, f2, mpnw1
  call mpfrsetprec (f2(1:), mpnw1)
  call mpfrset (f2(1:), f1(1:), mprnd)
 
!+   call mpdmc (1.d0, 0, f3, mpnw1
  call mpfrsetprec (f3(1:), mpnw1)
  call mpfrsetd (f3(1:), 1.d0, mprnd)
 
!+   call mpdmc (1.d0, 0, f4, mpnw1
  call mpfrsetprec (f4(1:), mpnw1)
  call mpfrsetd (f4(1:), 1.d0, mprnd)
 
!+   call mpdmc (1.d0, 0, f5, mpnw1
  call mpfrsetprec (f5(1:), mpnw1)
  call mpfrsetd (f5(1:), 1.d0, mprnd)
 

  do k = 1, nua
!+     call mpdmc (1.d0, 0, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrsetd (t2(1:), 1.d0, mprnd)
 
!+     call mpdivd (t2, dble (k), t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrdivd (t3(1:), t2(1:), dble (k), mprnd)
 
!+     call mpadd (f2, t3, t4, mpnw1
    call mpfrsetprec (t4(1:), mpnw1)
    call mpfradd (t4(1:), f2(1:), t3(1:), mprnd)
 
!+     call mpeq (t4, f2, mpnw1
    call mpfrsetprec (f2(1:), mpnw1)
    call mpfrset (f2(1:), t4(1:), mprnd)
 
!+     call mpmuld (f5, dble (k), t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmuld (t2(1:), f5(1:), dble (k), mprnd)
 
!+     call mpeq (t2, f5, mpnw1
    call mpfrsetprec (f5(1:), mpnw1)
    call mpfrset (f5(1:), t2(1:), mprnd)
 
  enddo

!+   call mpadd (f1, f2, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfradd (t2(1:), f1(1:), f2(1:), mprnd)
 
!+   call mpmul (t2, f3, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrmul (t3(1:), t2(1:), f3(1:), mprnd)
 
!+   call mpmul (f4, f5, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfrmul (t4(1:), f4(1:), f5(1:), mprnd)
 
!+   call mpdiv (t3, t4, sum3, mpnw1
  call mpfrsetprec (sum3(1:), mpnw1)
  call mpfrdiv (sum3(1:), t3(1:), t4(1:), mprnd)
 

  do k = 1, itrmax
!+     call mpdmc (1.d0, 0, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrsetd (t2(1:), 1.d0, mprnd)
 
!+     call mpdivd (t2, dble (k), t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrdivd (t3(1:), t2(1:), dble (k), mprnd)
 
!+     call mpadd (f1, t3, t4, mpnw1
    call mpfrsetprec (t4(1:), mpnw1)
    call mpfradd (t4(1:), f1(1:), t3(1:), mprnd)
 
!+     call mpeq (t4, f1, mpnw1
    call mpfrsetprec (f1(1:), mpnw1)
    call mpfrset (f1(1:), t4(1:), mprnd)
 
!+     call mpdivd (t2, dble (nua + k), t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrdivd (t3(1:), t2(1:), dble (nua + k), mprnd)
 
!+     call mpadd (f2, t3, t4, mpnw1
    call mpfrsetprec (t4(1:), mpnw1)
    call mpfradd (t4(1:), f2(1:), t3(1:), mprnd)
 
!+     call mpeq (t4, f2, mpnw1
    call mpfrsetprec (f2(1:), mpnw1)
    call mpfrset (f2(1:), t4(1:), mprnd)
 
!+     call mpmul (t1, f3, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmul (t2(1:), t1(1:), f3(1:), mprnd)
 
!+     call mpeq (t2, f3, mpnw1
    call mpfrsetprec (f3(1:), mpnw1)
    call mpfrset (f3(1:), t2(1:), mprnd)
 
!+     call mpmuld (f4, dble (k), t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmuld (t2(1:), f4(1:), dble (k), mprnd)
 
!+     call mpeq (t2, f4, mpnw1
    call mpfrsetprec (f4(1:), mpnw1)
    call mpfrset (f4(1:), t2(1:), mprnd)
 
!+     call mpmuld (f5, dble (nua + k), t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmuld (t2(1:), f5(1:), dble (nua + k), mprnd)
 
!+     call mpeq (t2, f5, mpnw1
    call mpfrsetprec (f5(1:), mpnw1)
    call mpfrset (f5(1:), t2(1:), mprnd)
 
!+     call mpadd (f1, f2, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfradd (t2(1:), f1(1:), f2(1:), mprnd)
 
!+     call mpmul (t2, f3, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrmul (t3(1:), t2(1:), f3(1:), mprnd)
 
!+     call mpmul (f4, f5, t4, mpnw1
    call mpfrsetprec (t4(1:), mpnw1)
    call mpfrmul (t4(1:), f4(1:), f5(1:), mprnd)
 
!+     call mpdiv (t3, t4, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrdiv (t2(1:), t3(1:), t4(1:), mprnd)
 
!+     call mpadd (sum3, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfradd (t3(1:), sum3(1:), t2(1:), mprnd)
 
!+     call mpeq (t3, sum3, mpnw1
    call mpfrsetprec (sum3(1:), mpnw1)
    call mpfrset (sum3(1:), t3(1:), mprnd)
 

!+     call mpabs (t2, tc1, 4
    call mpfrsetprec (tc1(1:), 4)
    call mpfrabs (tc1(1:), t2(1:), mprnd)
 
!+     call mpmul (eps, sum3, tc3, 4
    call mpfrsetprec (tc3(1:), 4)
    call mpfrmul (tc3(1:), eps(1:), sum3(1:), mprnd)
 
!+     call mpabs (tc3, tc2, 4
    call mpfrsetprec (tc2(1:), 4)
    call mpfrabs (tc2(1:), tc3(1:), mprnd)
 
!+     call mpcpr (tc1, tc2, ic1, 4
    ic1 = mpfrcmp (tc1(1:), tc2(1:))
 
    if (ic1 <= 0) goto 100
  enddo

  write (mpldb, 5)
5 format ('*** MPBESSELKNR: Loop end error 1')
  call mpabrt ( 528)

100 continue

!+   call mpmuld (rra, 0.5d0, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmuld (t2(1:), rra(1:), 0.5d0, mprnd)
 
!+   call mpnpwr (t2, nua, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrpowsi (t3(1:), t2(1:), nua, mprnd)
 
  d1 = (-1.d0)**nua * 0.5d0
!+   call mpmuld (t3, d1, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfrmuld (t4(1:), t3(1:), d1, mprnd)
 
!+   call mpmul (t4, sum3, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmul (t2(1:), t4(1:), sum3(1:), mprnd)
 
!+   call mpeq (t2, sum3, mpnw1
  call mpfrsetprec (sum3(1:), mpnw1)
  call mpfrset (sum3(1:), t2(1:), mprnd)
 
!+   call mpadd (sum1, sum2, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfradd (t2(1:), sum1(1:), sum2(1:), mprnd)
 
!+   call mpadd (t2, sum3, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfradd (t3(1:), t2(1:), sum3(1:), mprnd)
 
else
!  sum1 = mpreal (1.d0, mpnw)
!  t1 = mpreal (4.d0 * dble (nua)**2, mpnw)
!  t2 = mpreal (1.d0, mpnw)
!  t3 = mpreal (1.d0, mpnw)

!+   call mpdmc (1.d0, 0, sum1, mpnw1
  call mpfrsetprec (sum1(1:), mpnw1)
  call mpfrsetd (sum1(1:), 1.d0, mprnd)
 
  d1 = 4.d0 * dble (nua)**2
!+   call mpdmc (d1, 0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrsetd (t1(1:), d1, mprnd)
 
!+   call mpdmc (1.d0, 0, tn, mpnw1
  call mpfrsetprec (tn(1:), mpnw1)
  call mpfrsetd (tn(1:), 1.d0, mprnd)
 
!+   call mpdmc (1.d0, 0, td, mpnw1
  call mpfrsetprec (td(1:), mpnw1)
  call mpfrsetd (td(1:), 1.d0, mprnd)
 

  do k = 1, itrmax
!  t2 = t2 * (t1 - (2.d0*k - 1.d0)**2)
!  t3 = t3 * dble (k) * 8.d0 * xa
!  t4 = t2 / t3
!  sum1 = sum1 + t4

    d1 = 2.d0 * k - 1.d0
!+     call mpdmc (d1, 0, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrsetd (t2(1:), d1, mprnd)
 
!+     call mpmul (t2, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrmul (t3(1:), t2(1:), t2(1:), mprnd)
 
!+     call mpsub (t1, t3, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrsub (t2(1:), t1(1:), t3(1:), mprnd)
 
!+     call mpmul (tn, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrmul (t3(1:), tn(1:), t2(1:), mprnd)
 
!+     call mpeq (t3, tn, mpnw1
    call mpfrsetprec (tn(1:), mpnw1)
    call mpfrset (tn(1:), t3(1:), mprnd)
 
!+     call mpmuld (rra, 8.d0 * k, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmuld (t2(1:), rra(1:), 8.d0 * k, mprnd)
 
!+     call mpmul (td, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrmul (t3(1:), td(1:), t2(1:), mprnd)
 
!+     call mpeq (t3, td, mpnw1
    call mpfrsetprec (td(1:), mpnw1)
    call mpfrset (td(1:), t3(1:), mprnd)
 
!+     call mpdiv (tn, td, t4, mpnw1
    call mpfrsetprec (t4(1:), mpnw1)
    call mpfrdiv (t4(1:), tn(1:), td(1:), mprnd)
 
!+     call mpadd (sum1, t4, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfradd (t3(1:), sum1(1:), t4(1:), mprnd)
 
!+     call mpeq (t3, sum1, mpnw1
    call mpfrsetprec (sum1(1:), mpnw1)
    call mpfrset (sum1(1:), t3(1:), mprnd)
 

!   if (abs (t4) / abs (sum1) < eps) goto 110

!+     call mpabs (t4, tc1, 4
    call mpfrsetprec (tc1(1:), 4)
    call mpfrabs (tc1(1:), t4(1:), mprnd)
 
!+     call mpmul (eps, sum1, tc3, 4
    call mpfrsetprec (tc3(1:), 4)
    call mpfrmul (tc3(1:), eps(1:), sum1(1:), mprnd)
 
!+     call mpabs (tc3, tc2, 4
    call mpfrsetprec (tc2(1:), 4)
    call mpfrabs (tc2(1:), tc3(1:), mprnd)
 
!+     call mpcpr (tc1, tc2, ic1, 4
    ic1 = mpfrcmp (tc1(1:), tc2(1:))
 
    if (ic1 <= 0) goto 110
  enddo

write (mpldb, 6)
6 format ('*** MPBESSELKNR: Loop end error 2')
call mpabrt ( 529)

110 continue

! t1 = sqrt (mppi (mpnw) / (2.d0 * xa)) / exp (xa)
! besseli = t1 * sum1

!+   call mpexp (rra, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrexp (t1(1:), rra(1:), mprnd)
 
!+   call mpmuld (rra, 2.d0, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmuld (t2(1:), rra(1:), 2.d0, mprnd)
 
!+   call mpdiv (mppicon, t2, t3, mpnw1
call mpfixlocr (mppicon)
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrdiv (t3(1:), mppicon(1:), t2(1:), mprnd)
 
!+   call mpsqrt (t3, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfrsqrt (t4(1:), t3(1:), mprnd)
 
!+   call mpdiv (t4, t1, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrdiv (t2(1:), t4(1:), t1(1:), mprnd)
 
!+   call mpmul (t2, sum1, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrmul (t3(1:), t2(1:), sum1(1:), mprnd)
 
endif

! if (x < 0.d0 .and. mod (nu, 2) /= 0) besselk = - besselk

!+ if (mpsgn (rr) < 0 .and. mod (nu, 2) /= 0) the
if (mpfrsgn (rr(1:)) < 0 .and. mod (nu, 2) /= 0) then
 
!+   call mpneg (t3, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfrneg (t4(1:), t3(1:), mprnd)
 
!+   call mpeq (t4, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrset (t3(1:), t4(1:), mprnd)
 
endif
!+ call mproun (t3, mpnw
call mpfrprecround (t3(1:), mpnw, mprnd)
 
!+ call mpeq (t3, ss, mpnw
call mpfrsetprec (ss(1:), mpnw)
call mpfrset (ss(1:), t3(1:), mprnd)
 
return
end subroutine mpbesselknr

subroutine mpbesselkr (qq, rr, ss, mpnw)

!   This evaluates the Bessel function BesselK (QQ,RR) for QQ and RR
!   both MPR. This uses DLMF formula 10.27.4.

implicit none
integer (mpiknd), intent(in):: qq(0:), rr(0:)
integer (mpiknd), intent(out):: ss(0:)
integer, intent(in):: mpnw
integer i1, i2, mpnw1, n1
real (mprknd) d1
integer (mpiknd) t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6)

integer, external:: mpfrsgn
integer, external:: mpfrcmp
real (mprknd), external:: mpfrgetd2exp
integer (mpiknd) ix8
 
!  End of declaration

if (mpnw < 4 .or. mpspacer (qq) < mpnw + 4 .or. mpspacer (rr) < mpnw + 4 &
  .or. mpspacer (ss) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPBESSELKR: Uninitialized or inadequately sized arrays')
  call mpabrt ( 530)
endif

!  End of initial input array check

mpnw1 = min (mpnw + 1, mpnwx)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)

!   If QQ is integer, call mpbesselknr; if qq < 0 and rr <= 0, then error.

!+ call mpinfr (qq, t1, t2, mpnw
call mpfrsetprec (t1(1:), mpnw)
call mpfrsetprec (t2(1:), mpnw)
call mpfrtrunc (t1(1:), qq(1:))
call mpfrsub (t2(1:), qq(1:), t1(1:), mprnd)
 
!+ i1 = mpsgn (qq
i1 = mpfrsgn (qq(1:))
 
!+ i2 = mpsgn (rr
i2 = mpfrsgn (rr(1:))
 
!+ if (mpsgn (t2) == 0) the
if (mpfrsgn (t2(1:)) == 0) then
 
!+   call mpmdc (qq, d1, n1, mpnw
  d1 = 2.d0 * mpfrgetd2exp (ix8, qq(1:), mprnd)
  if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
  d1 = d1 * 2.d0**n1
  n1 = nint (d1)
  call mpbesselknr (n1, rr, t1, mpnw)
  goto 120
elseif (i1 < 0 .and. i2 <= 0) then
  write (mpldb, 2)
2 format ('*** MPBESSELKR: First argument < 0 and second argument <= 0')
  call mpabrt ( 531)
endif

!+ call mpneg (qq, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrneg (t1(1:), qq(1:), mprnd)
 
call mpbesselir (t1, rr, t2, mpnw1)
call mpbesselir (qq, rr, t3, mpnw1)
!+ call mpsub (t2, t3, t4, mpnw1
call mpfrsetprec (t4(1:), mpnw1)
call mpfrsub (t4(1:), t2(1:), t3(1:), mprnd)
 
!+ call mpmul (qq, mppicon, t1, mpnw1
call mpfixlocr (mppicon)
call mpfrsetprec (t1(1:), mpnw1)
call mpfrmul (t1(1:), qq(1:), mppicon(1:), mprnd)
 
!+ call mpcssnr (t1, t2, t3, mpnw1
call mpfrsetprec (t3(1:), mpnw1)
call mpfrsincos (t3(1:), t2(1:), t1(1:), mprnd)
 
!+ call mpdiv (t4, t3, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfrdiv (t2(1:), t4(1:), t3(1:), mprnd)
 
!+ call mpmul (mppicon, t2, t3, mpnw1
call mpfixlocr (mppicon)
call mpfrsetprec (t3(1:), mpnw1)
call mpfrmul (t3(1:), mppicon(1:), t2(1:), mprnd)
 
!+ call mpmuld (t3, 0.5d0, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrmuld (t1(1:), t3(1:), 0.5d0, mprnd)
 

120 continue

!+ call mproun (t1, mpnw
call mpfrprecround (t1(1:), mpnw, mprnd)
 
!+ call mpeq (t1, ss, mpnw
call mpfrsetprec (ss(1:), mpnw)
call mpfrset (ss(1:), t1(1:), mprnd)
 
return
end subroutine mpbesselkr

subroutine mpbesselynr (nu, rr, ss, mpnw)

!   This evaluates the modified Bessel function BesselY (NU,RR).
!   NU is an integer. The algorithm is DLMF formula 10.8.1 for modest RR,
!   and DLMF 10.17.4 for large RR, relative to precision.

implicit none
integer, intent(in):: nu, mpnw
integer (mpiknd), intent(in):: rr(0:)
integer (mpiknd), intent(out):: ss(0:)
integer, parameter:: itrmax = 1000000
real (mprknd), parameter:: dfrac = 1.5d0, egam = 0.5772156649015328606d0, &
  pi = 3.1415926535897932385d0
integer ic1, ic2, k, mpnw1, nua, n1
real (mprknd) d1, d2
integer (mpiknd) f1(0:2*mpnw+6), f2(0:2*mpnw+6), f3(0:2*mpnw+6), f4(0:2*mpnw+6), &
  f5(0:2*mpnw+6), rra(0:2*mpnw+6), rr2(0:2*mpnw+6), sum1(0:2*mpnw+6), &
  sum2(0:2*mpnw+6), sum3(0:2*mpnw+6), td1(0:2*mpnw+6), td2(0:2*mpnw+6), &
  tn1(0:2*mpnw+6), tn2(0:2*mpnw+6), t1(0:2*mpnw+6), t2(0:2*mpnw+6), &
  t3(0:2*mpnw+6), t4(0:2*mpnw+6), t41(0:2*mpnw+6), t42(0:2*mpnw+6), &
  t5(0:2*mpnw+6), tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)

integer, external:: mpfrsgn
integer, external:: mpfrcmp
real (mprknd), external:: mpfrgetd2exp
integer (mpiknd) ix8
 
!  End of declaration

if (mpnw < 4 .or. mpspacer (rr) < mpnw + 4 .or. mpspacer (ss) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPBESSELYNR: Uninitialized or inadequately sized arrays')
  call mpabrt ( 532)
endif

!   Check if EGAMMA and PI have been precomputed.

!+ call mpmdc (mpegammacon, d1, n1, mpnw+1
call mpfixlocr (mpegammacon)
d1 = 2.d0 * mpfrgetd2exp (ix8, mpegammacon(1:), mprnd)
if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
if (n1 /= -1 .or. abs (d1 * 2.d0**n1 - egam) > mprdfz &
  .or. mpwprecr (mpegammacon) < mpnw+1) then
  write (mpldb, 3) mpnw+1
3 format ('*** MPBESSELYNR: EGAMMA must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 534)
endif
!+ call mpmdc (mppicon, d1, n1, mpnw+1
call mpfixlocr (mppicon)
d1 = 2.d0 * mpfrgetd2exp (ix8, mppicon(1:), mprnd)
if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
if (n1 /= 1 .or. abs (d1 * 2.d0**n1 - pi) > mprdfz &
  .or. mpwprecr (mppicon) < mpnw+1) then
  write (mpldb, 4) mpnw+1
4 format ('*** MPBESSELYNR: Pi must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 535)
endif

!  End of initial input array check

!   Check for RR = 0.

!+ if (mpsgn (rr) == 0) the
if (mpfrsgn (rr(1:)) == 0) then
 
  write (mpldb, 2)
2 format ('*** MPBESSELYNR: argument is negative or too large')
  call mpabrt ( 533)
endif

mpnw1 = min (2 * mpnw + 1, mpnwx)
nua = abs (nu)
call mpinitwds (f1, mpnw1)
call mpinitwds (f2, mpnw1)
call mpinitwds (f3, mpnw1)
call mpinitwds (f4, mpnw1)
call mpinitwds (f5, mpnw1)
call mpinitwds (rra, mpnw1)
call mpinitwds (rr2, mpnw1)
call mpinitwds (sum1, mpnw1)
call mpinitwds (sum2, mpnw1)
call mpinitwds (sum3, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (t41, mpnw1)
call mpinitwds (t42, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (td1, mpnw1)
call mpinitwds (td2, mpnw1)
call mpinitwds (tn1, mpnw1)
call mpinitwds (tn2, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
!+ call mpdmc (2.d0, 0, tc1, 4
call mpfrsetprec (tc1(1:), 4)
call mpfrsetd (tc1(1:), 2.d0, mprnd)
 
!+ call mpnpwr (tc1, -mpnw*mpnbt, eps, 4
call mpfrsetprec (eps(1:), 4)
call mpfrpowsi (eps(1:), tc1(1:), -mpnw*mpnbt, mprnd)
 
!+ call mpmdc (rr, d1, n1, 4
d1 = 2.d0 * mpfrgetd2exp (ix8, rr(1:), mprnd)
if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
d1 = abs (d1) * 2.d0 ** n1

if (d1 < dfrac * mpnw1 * mpdpw) then
  mpnw1 = min (mpnw + nint (d1 / (dfrac * mpdpw)), 2 * mpnw + 1, mpnwx)
!+   call mpabs (rr, rra, mpnw1
  call mpfrsetprec (rra(1:), mpnw1)
  call mpfrabs (rra(1:), rr(1:), mprnd)
 
!+   call mpmul (rra, rra, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmul (t2(1:), rra(1:), rra(1:), mprnd)
 
!+   call mpmuld (t2, 0.25d0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmuld (t1(1:), t2(1:), 0.25d0, mprnd)
 
!+   call mpdmc (1.d0, 0, f1, mpnw1
  call mpfrsetprec (f1(1:), mpnw1)
  call mpfrsetd (f1(1:), 1.d0, mprnd)
 
!+   call mpdmc (1.d0, 0, f2, mpnw1
  call mpfrsetprec (f2(1:), mpnw1)
  call mpfrsetd (f2(1:), 1.d0, mprnd)
 
!+   call mpdmc (1.d0, 0, f3, mpnw1
  call mpfrsetprec (f3(1:), mpnw1)
  call mpfrsetd (f3(1:), 1.d0, mprnd)
 
!+   call mpdmc (0.d0, 0,  sum1, mpnw1
  call mpfrsetprec (sum1(1:), mpnw1)
  call mpfrsetd (sum1(1:), 0.d0, mprnd)
 

  do k = 1, nua - 1
!+     call mpmuld (f1, dble (k), t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmuld (t2(1:), f1(1:), dble (k), mprnd)
 
!+     call mpeq (t2, f1, mpnw1
    call mpfrsetprec (f1(1:), mpnw1)
    call mpfrset (f1(1:), t2(1:), mprnd)
 
  enddo

  do k = 0, nua - 1
    if (k > 0) then
!+       call mpdivd (f1, dble (nua - k), t2, mpnw1
      call mpfrsetprec (t2(1:), mpnw1)
      call mpfrdivd (t2(1:), f1(1:), dble (nua - k), mprnd)
 
!+       call mpeq (t2, f1, mpnw1
      call mpfrsetprec (f1(1:), mpnw1)
      call mpfrset (f1(1:), t2(1:), mprnd)
 
!+       call mpmul (t1, f2, t2, mpnw1
      call mpfrsetprec (t2(1:), mpnw1)
      call mpfrmul (t2(1:), t1(1:), f2(1:), mprnd)
 
!+       call mpeq (t2, f2, mpnw1
      call mpfrsetprec (f2(1:), mpnw1)
      call mpfrset (f2(1:), t2(1:), mprnd)
 
!+       call mpmuld (f3, dble (k), t2, mpnw1
      call mpfrsetprec (t2(1:), mpnw1)
      call mpfrmuld (t2(1:), f3(1:), dble (k), mprnd)
 
!+       call mpeq (t2, f3, mpnw1
      call mpfrsetprec (f3(1:), mpnw1)
      call mpfrset (f3(1:), t2(1:), mprnd)
 
    endif
!+     call mpmul (f1, f2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrmul (t3(1:), f1(1:), f2(1:), mprnd)
 
!+     call mpdiv (t3, f3, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrdiv (t2(1:), t3(1:), f3(1:), mprnd)
 
!+     call mpadd (sum1, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfradd (t3(1:), sum1(1:), t2(1:), mprnd)
 
!+     call mpeq (t3, sum1, mpnw1
    call mpfrsetprec (sum1(1:), mpnw1)
    call mpfrset (sum1(1:), t3(1:), mprnd)
 
  enddo

!+   call mpmuld (rra, 0.5d0, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrmuld (t3(1:), rra(1:), 0.5d0, mprnd)
 
!+   call mpnpwr (t3, nua, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfrpowsi (t4(1:), t3(1:), nua, mprnd)
 
!+   call mpdiv (sum1, t4, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrdiv (t3(1:), sum1(1:), t4(1:), mprnd)
 
!+   call mpneg (t3, sum1, mpnw1
  call mpfrsetprec (sum1(1:), mpnw1)
  call mpfrneg (sum1(1:), t3(1:), mprnd)
 

!+   call mpmuld (rra, 0.5d0, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmuld (t2(1:), rra(1:), 0.5d0, mprnd)
 
!+   call mplog (t2, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrlog (t3(1:), t2(1:), mprnd)
 
!+   call mpmuld (t3, 2.d0, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmuld (t2(1:), t3(1:), 2.d0, mprnd)
 
  call mpbesseljnr (nua, rra, t3, mpnw1)
!+   call mpmul (t2, t3, sum2, mpnw1
  call mpfrsetprec (sum2(1:), mpnw1)
  call mpfrmul (sum2(1:), t2(1:), t3(1:), mprnd)
 

!+   call mpneg (mpegammacon, f1, mpnw1
call mpfixlocr (mpegammacon)
  call mpfrsetprec (f1(1:), mpnw1)
  call mpfrneg (f1(1:), mpegammacon(1:), mprnd)
 
!+   call mpeq (f1, f2, mpnw1
  call mpfrsetprec (f2(1:), mpnw1)
  call mpfrset (f2(1:), f1(1:), mprnd)
 
!+   call mpdmc (1.d0, 0, f3, mpnw1
  call mpfrsetprec (f3(1:), mpnw1)
  call mpfrsetd (f3(1:), 1.d0, mprnd)
 
!+   call mpdmc (1.d0, 0, f4, mpnw1
  call mpfrsetprec (f4(1:), mpnw1)
  call mpfrsetd (f4(1:), 1.d0, mprnd)
 
!+   call mpdmc (1.d0, 0, f5, mpnw1
  call mpfrsetprec (f5(1:), mpnw1)
  call mpfrsetd (f5(1:), 1.d0, mprnd)
 

  do k = 1, nua
!+     call mpdmc (1.d0, 0, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrsetd (t2(1:), 1.d0, mprnd)
 
!+     call mpdivd (t2, dble (k), t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrdivd (t3(1:), t2(1:), dble (k), mprnd)
 
!+     call mpadd (f2, t3, t4, mpnw1
    call mpfrsetprec (t4(1:), mpnw1)
    call mpfradd (t4(1:), f2(1:), t3(1:), mprnd)
 
!+     call mpeq (t4, f2, mpnw1
    call mpfrsetprec (f2(1:), mpnw1)
    call mpfrset (f2(1:), t4(1:), mprnd)
 
!+     call mpmuld (f5, dble (k), t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmuld (t2(1:), f5(1:), dble (k), mprnd)
 
!+     call mpeq (t2, f5, mpnw1
    call mpfrsetprec (f5(1:), mpnw1)
    call mpfrset (f5(1:), t2(1:), mprnd)
 
  enddo

!+   call mpadd (f1, f2, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfradd (t2(1:), f1(1:), f2(1:), mprnd)
 
!+   call mpmul (t2, f3, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrmul (t3(1:), t2(1:), f3(1:), mprnd)
 
!+   call mpmul (f4, f5, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfrmul (t4(1:), f4(1:), f5(1:), mprnd)
 
!+   call mpdiv (t3, t4, sum3, mpnw1
  call mpfrsetprec (sum3(1:), mpnw1)
  call mpfrdiv (sum3(1:), t3(1:), t4(1:), mprnd)
 

  do k = 1, itrmax
!+     call mpdmc (1.d0, 0, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrsetd (t2(1:), 1.d0, mprnd)
 
!+     call mpdivd (t2, dble (k), t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrdivd (t3(1:), t2(1:), dble (k), mprnd)
 
!+     call mpadd (f1, t3, t4, mpnw1
    call mpfrsetprec (t4(1:), mpnw1)
    call mpfradd (t4(1:), f1(1:), t3(1:), mprnd)
 
!+     call mpeq (t4, f1, mpnw1
    call mpfrsetprec (f1(1:), mpnw1)
    call mpfrset (f1(1:), t4(1:), mprnd)
 
!+     call mpdivd (t2, dble (nua + k), t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrdivd (t3(1:), t2(1:), dble (nua + k), mprnd)
 
!+     call mpadd (f2, t3, t4, mpnw1
    call mpfrsetprec (t4(1:), mpnw1)
    call mpfradd (t4(1:), f2(1:), t3(1:), mprnd)
 
!+     call mpeq (t4, f2, mpnw1
    call mpfrsetprec (f2(1:), mpnw1)
    call mpfrset (f2(1:), t4(1:), mprnd)
 
!+     call mpmul (t1, f3, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmul (t2(1:), t1(1:), f3(1:), mprnd)
 
!+     call mpneg (t2, f3, mpnw1
    call mpfrsetprec (f3(1:), mpnw1)
    call mpfrneg (f3(1:), t2(1:), mprnd)
 
!+     call mpmuld (f4, dble (k), t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmuld (t2(1:), f4(1:), dble (k), mprnd)
 
!+     call mpeq (t2, f4, mpnw1
    call mpfrsetprec (f4(1:), mpnw1)
    call mpfrset (f4(1:), t2(1:), mprnd)
 
!+     call mpmuld (f5, dble (nua + k), t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmuld (t2(1:), f5(1:), dble (nua + k), mprnd)
 
!+     call mpeq (t2, f5, mpnw1
    call mpfrsetprec (f5(1:), mpnw1)
    call mpfrset (f5(1:), t2(1:), mprnd)
 
!+     call mpadd (f1, f2, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfradd (t2(1:), f1(1:), f2(1:), mprnd)
 
!+     call mpmul (t2, f3, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrmul (t3(1:), t2(1:), f3(1:), mprnd)
 
!+     call mpmul (f4, f5, t4, mpnw1
    call mpfrsetprec (t4(1:), mpnw1)
    call mpfrmul (t4(1:), f4(1:), f5(1:), mprnd)
 
!+     call mpdiv (t3, t4, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrdiv (t2(1:), t3(1:), t4(1:), mprnd)
 
!+     call mpadd (sum3, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfradd (t3(1:), sum3(1:), t2(1:), mprnd)
 
!+     call mpeq (t3, sum3, mpnw1
    call mpfrsetprec (sum3(1:), mpnw1)
    call mpfrset (sum3(1:), t3(1:), mprnd)
 

!+     call mpabs (t2, tc1, 4
    call mpfrsetprec (tc1(1:), 4)
    call mpfrabs (tc1(1:), t2(1:), mprnd)
 
!+     call mpmul (eps, sum3, tc3, 4
    call mpfrsetprec (tc3(1:), 4)
    call mpfrmul (tc3(1:), eps(1:), sum3(1:), mprnd)
 
!+     call mpabs (tc3, tc2, 4
    call mpfrsetprec (tc2(1:), 4)
    call mpfrabs (tc2(1:), tc3(1:), mprnd)
 
!+     call mpcpr (tc1, tc2, ic1, 4
    ic1 = mpfrcmp (tc1(1:), tc2(1:))
 
    if (ic1 <= 0) goto 100
  enddo

  write (mpldb, 6)
6 format ('*** MPBESSELYNR: Loop end error 1')
  call mpabrt ( 536)

100 continue

!+   call mpmuld (rra, 0.5d0, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmuld (t2(1:), rra(1:), 0.5d0, mprnd)
 
!+   call mpnpwr (t2, nua, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrpowsi (t3(1:), t2(1:), nua, mprnd)
 
!+   call mpmul (t3, sum3, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmul (t2(1:), t3(1:), sum3(1:), mprnd)
 
!+   call mpneg (t2, sum3, mpnw1
  call mpfrsetprec (sum3(1:), mpnw1)
  call mpfrneg (sum3(1:), t2(1:), mprnd)
 

!+   call mpadd (sum1, sum2, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfradd (t2(1:), sum1(1:), sum2(1:), mprnd)
 
!+   call mpadd (t2, sum3, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfradd (t4(1:), t2(1:), sum3(1:), mprnd)
 
!+   call mpeq (mppicon, t2, mpnw1
call mpfixlocr (mppicon)
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrset (t2(1:), mppicon(1:), mprnd)
 
!+   call mpdiv (t4, t2, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrdiv (t3(1:), t4(1:), t2(1:), mprnd)
 
else

! xa2 = xa**2
! t1 = mpreal (4.d0 * dble (nua)**2, mpnw)
! tn1 = mpreal (1.d0, mpnw)
! tn2 = (t1 - 1.d0) / 8.d0
! td1 = mpreal (1.d0, mpnw)
! td2 = xa
! sum1 = tn1 / td1
! sum2 = tn2 / t32

  mpnw1 = min (mpnw + 1, mpnwx)
!+   call mpabs (rr, rra, mpnw1
  call mpfrsetprec (rra(1:), mpnw1)
  call mpfrabs (rra(1:), rr(1:), mprnd)
 
!+   call mpmul (rra, rra, rr2, mpnw1
  call mpfrsetprec (rr2(1:), mpnw1)
  call mpfrmul (rr2(1:), rra(1:), rra(1:), mprnd)
 
  d1 = 4.d0 * dble (nua)**2
!+   call mpdmc (d1, 0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrsetd (t1(1:), d1, mprnd)
 
!+   call mpdmc (1.d0, 0, tn1, mpnw1
  call mpfrsetprec (tn1(1:), mpnw1)
  call mpfrsetd (tn1(1:), 1.d0, mprnd)
 
!+   call mpsub (t1, tn1, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrsub (t2(1:), t1(1:), tn1(1:), mprnd)
 
!+   call mpdivd (t2, 8.d0, tn2, mpnw1
  call mpfrsetprec (tn2(1:), mpnw1)
  call mpfrdivd (tn2(1:), t2(1:), 8.d0, mprnd)
 
!+   call mpdmc (1.d0, 0, td1, mpnw1
  call mpfrsetprec (td1(1:), mpnw1)
  call mpfrsetd (td1(1:), 1.d0, mprnd)
 
!+   call mpeq (rra, td2, mpnw1
  call mpfrsetprec (td2(1:), mpnw1)
  call mpfrset (td2(1:), rra(1:), mprnd)
 
!+   call mpdiv (tn1, td1, sum1, mpnw1
  call mpfrsetprec (sum1(1:), mpnw1)
  call mpfrdiv (sum1(1:), tn1(1:), td1(1:), mprnd)
 
!+   call mpdiv (tn2, td2, sum2, mpnw1
  call mpfrsetprec (sum2(1:), mpnw1)
  call mpfrdiv (sum2(1:), tn2(1:), td2(1:), mprnd)
 

  do k = 1, itrmax
!   tn1 = -tn1 * (t1 - (2.d0*(2.d0*k-1.d0) - 1.d0)**2) * (t1 - (2.d0*(2.d0*k) - 1.d0)**2)

    d1 = (2.d0*(2.d0*k-1.d0) - 1.d0)**2
    d2 = (2.d0*(2.d0*k) - 1.d0)**2
!+     call mpdmc (d1, 0, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrsetd (t2(1:), d1, mprnd)
 
!+     call mpsub (t1, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrsub (t3(1:), t1(1:), t2(1:), mprnd)
 
!+     call mpdmc (d2, 0, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrsetd (t2(1:), d2, mprnd)
 
!+     call mpsub (t1, t2, t5, mpnw1
    call mpfrsetprec (t5(1:), mpnw1)
    call mpfrsub (t5(1:), t1(1:), t2(1:), mprnd)
 
!+     call mpmul (t3, t5, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmul (t2(1:), t3(1:), t5(1:), mprnd)
 
!+     call mpmul (tn1, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrmul (t3(1:), tn1(1:), t2(1:), mprnd)
 
!+     call mpneg (t3, tn1, mpnw1
    call mpfrsetprec (tn1(1:), mpnw1)
    call mpfrneg (tn1(1:), t3(1:), mprnd)
 

!   td1 = td1 * dble (2*k-1) * dble (2*k) * 64.d0 * xa2

    d1 = dble (2*k-1) * dble (2*k) * 64.d0
!+     call mpmuld (td1, d1, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmuld (t2(1:), td1(1:), d1, mprnd)
 
!+     call mpmul (t2, rr2, td1, mpnw1
    call mpfrsetprec (td1(1:), mpnw1)
    call mpfrmul (td1(1:), t2(1:), rr2(1:), mprnd)
 

!   t41 = tn1 / td1
!   sum1 = sum1 + t41

!+     call mpdiv (tn1, td1, t41, mpnw1
    call mpfrsetprec (t41(1:), mpnw1)
    call mpfrdiv (t41(1:), tn1(1:), td1(1:), mprnd)
 
!+     call mpadd (sum1, t41, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfradd (t2(1:), sum1(1:), t41(1:), mprnd)
 
!+     call mpeq (t2, sum1, mpnw1
    call mpfrsetprec (sum1(1:), mpnw1)
    call mpfrset (sum1(1:), t2(1:), mprnd)
 

!   tn2 = -tn2 * (t1 - (2.d0*(2.d0*k) - 1.d0)**2) * (t1 - (2.d0*(2.d0*k+1.d0) - 1.d0)**2)

    d1 = (2.d0*(2.d0*k) - 1.d0)**2
    d2 = (2.d0*(2.d0*k+1.d0) - 1.d0)**2
!+     call mpdmc (d1, 0, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrsetd (t2(1:), d1, mprnd)
 
!+     call mpsub (t1, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrsub (t3(1:), t1(1:), t2(1:), mprnd)
 
!+     call mpdmc (d2, 0, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrsetd (t2(1:), d2, mprnd)
 
!+     call mpsub (t1, t2, t5, mpnw1
    call mpfrsetprec (t5(1:), mpnw1)
    call mpfrsub (t5(1:), t1(1:), t2(1:), mprnd)
 
!+     call mpmul (t3, t5, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmul (t2(1:), t3(1:), t5(1:), mprnd)
 
!+     call mpmul (tn2, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrmul (t3(1:), tn2(1:), t2(1:), mprnd)
 
!+     call mpneg (t3, tn2, mpnw1
    call mpfrsetprec (tn2(1:), mpnw1)
    call mpfrneg (tn2(1:), t3(1:), mprnd)
 

!   td2 = td2 * dble (2*k) * dble (2*k+1) * 64.d0 * xa2

    d1 = dble (2*k) * dble (2*k+1) * 64.d0
!+     call mpmuld (td2, d1, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmuld (t2(1:), td2(1:), d1, mprnd)
 
!+     call mpmul (t2, rr2, td2, mpnw1
    call mpfrsetprec (td2(1:), mpnw1)
    call mpfrmul (td2(1:), t2(1:), rr2(1:), mprnd)
 

!  t42 = tn2 / td2
!  sum2 = sum2 + t42

!+     call mpdiv (tn2, td2, t42, mpnw1
    call mpfrsetprec (t42(1:), mpnw1)
    call mpfrdiv (t42(1:), tn2(1:), td2(1:), mprnd)
 
!+     call mpadd (sum2, t42, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfradd (t2(1:), sum2(1:), t42(1:), mprnd)
 
!+     call mpeq (t2, sum2, mpnw1
    call mpfrsetprec (sum2(1:), mpnw1)
    call mpfrset (sum2(1:), t2(1:), mprnd)
 

!  if (abs (t41) / abs (sum1) < eps .and. abs (t42) / abs (sum2) < eps ) goto 110

!+     call mpabs (t41, tc1, 4
    call mpfrsetprec (tc1(1:), 4)
    call mpfrabs (tc1(1:), t41(1:), mprnd)
 
!+     call mpmul (eps, sum1, tc3, 4
    call mpfrsetprec (tc3(1:), 4)
    call mpfrmul (tc3(1:), eps(1:), sum1(1:), mprnd)
 
!+     call mpabs (tc3, tc2, 4
    call mpfrsetprec (tc2(1:), 4)
    call mpfrabs (tc2(1:), tc3(1:), mprnd)
 
!+     call mpcpr (tc1, tc2, ic1, 4
    ic1 = mpfrcmp (tc1(1:), tc2(1:))
 
!+     call mpabs (t42, tc1, 4
    call mpfrsetprec (tc1(1:), 4)
    call mpfrabs (tc1(1:), t42(1:), mprnd)
 
!+     call mpmul (eps, sum2, tc3, 4
    call mpfrsetprec (tc3(1:), 4)
    call mpfrmul (tc3(1:), eps(1:), sum2(1:), mprnd)
 
!+     call mpabs (tc3, tc2, 4
    call mpfrsetprec (tc2(1:), 4)
    call mpfrabs (tc2(1:), tc3(1:), mprnd)
 
!+     call mpcpr (tc1, tc2, ic2, 4
    ic2 = mpfrcmp (tc1(1:), tc2(1:))
 
    if (ic1 <= 0 .and. ic2 <= 0) goto 110
  enddo

  write (mpldb, 5)
5 format ('*** MPBESSELYNR: Loop end error 2')
  call mpabrt ( 537)

110 continue

! t1 = xa - 0.5d0 * nua * pi - 0.25d0 * pi
! besselj = sqrt (2.d0 / (pi * xa)) * (cos (t1) * sum1 - sin (t1) * sum2)

!+   call mpmuld (mppicon, 0.5d0 * nua, t1, mpnw1
call mpfixlocr (mppicon)
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmuld (t1(1:), mppicon(1:), 0.5d0 * nua, mprnd)
 
!+   call mpsub (rra, t1, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrsub (t2(1:), rra(1:), t1(1:), mprnd)
 
!+   call mpmuld (mppicon, 0.25d0, t1, mpnw1
call mpfixlocr (mppicon)
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmuld (t1(1:), mppicon(1:), 0.25d0, mprnd)
 
!+   call mpsub (t2, t1, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrsub (t3(1:), t2(1:), t1(1:), mprnd)
 
!+   call mpcssnr (t3, t41, t42, mpnw1
  call mpfrsetprec (t42(1:), mpnw1)
  call mpfrsincos (t42(1:), t41(1:), t3(1:), mprnd)
 
!+   call mpmul (t42, sum1, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmul (t1(1:), t42(1:), sum1(1:), mprnd)
 
!+   call mpmul (t41, sum2, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmul (t2(1:), t41(1:), sum2(1:), mprnd)
 
!+   call mpadd (t1, t2, t5, mpnw1
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfradd (t5(1:), t1(1:), t2(1:), mprnd)
 
!+   call mpmul (mppicon, rra, t1, mpnw1
call mpfixlocr (mppicon)
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmul (t1(1:), mppicon(1:), rra(1:), mprnd)
 
!+   call mpdmc (2.d0, 0, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrsetd (t2(1:), 2.d0, mprnd)
 
!+   call mpdiv (t2, t1, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrdiv (t3(1:), t2(1:), t1(1:), mprnd)
 
!+   call mpsqrt (t3, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrsqrt (t1(1:), t3(1:), mprnd)
 
!+   call mpmul (t1, t5, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrmul (t3(1:), t1(1:), t5(1:), mprnd)
 
endif

if (mod (nu, 2) /= 0) then
!   if (nu < 0 .and. x > 0.d0 .or. nu > 0 .and. x < 0.d0) bessely = - bessely

!+   ic1 = mpsgn (rr
  ic1 = mpfrsgn (rr(1:))
 
  if (nu < 0 .and. ic1 > 0 .or. nu > 0 .and. ic1 < 0) then
!+     call mpneg (t3, t4, mpnw1
    call mpfrsetprec (t4(1:), mpnw1)
    call mpfrneg (t4(1:), t3(1:), mprnd)
 
!+     call mpeq (t4, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrset (t3(1:), t4(1:), mprnd)
 
  endif
endif

!+ call mproun (t3, mpnw
call mpfrprecround (t3(1:), mpnw, mprnd)
 
!+ call mpeq (t3, ss, mpnw
call mpfrsetprec (ss(1:), mpnw)
call mpfrset (ss(1:), t3(1:), mprnd)
 
return
end subroutine mpbesselynr

subroutine mpbesselyr (qq, rr, ss, mpnw)

!   This evaluates the modified Bessel function BesselY (QQ,RR).
!   NU is an integer. The algorithm is DLMF formula 10.2.2.

implicit none
integer (mpiknd), intent(in):: qq(0:), rr(0:)
integer (mpiknd), intent(out):: ss(0:)
integer, intent(in):: mpnw
integer i1, i2, mpnw1, n1
real (mprknd) d1
integer (mpiknd) t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6)

integer, external:: mpfrsgn
integer, external:: mpfrcmp
real (mprknd), external:: mpfrgetd2exp
integer (mpiknd) ix8
 
!  End of declaration

if (mpnw < 4 .or. mpspacer (qq) < mpnw + 4 .or. mpspacer (rr) < mpnw + 4 &
  .or. mpspacer (ss) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPBESSELYR: Uninitialized or inadequately sized arrays')
  call mpabrt ( 538)
endif

!  End of initial input array check

mpnw1 = min (mpnw + 1, mpnwx)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)

!   If QQ is integer, call mpbesselynr; if qq < 0 and rr <= 0, then error.

!+ call mpinfr (qq, t1, t2, mpnw
call mpfrsetprec (t1(1:), mpnw)
call mpfrsetprec (t2(1:), mpnw)
call mpfrtrunc (t1(1:), qq(1:))
call mpfrsub (t2(1:), qq(1:), t1(1:), mprnd)
 
!+ i1 = mpsgn (qq
i1 = mpfrsgn (qq(1:))
 
!+ i2 = mpsgn (rr
i2 = mpfrsgn (rr(1:))
 
!+ if (mpsgn (t2) == 0) the
if (mpfrsgn (t2(1:)) == 0) then
 
!+   call mpmdc (qq, d1, n1, mpnw
  d1 = 2.d0 * mpfrgetd2exp (ix8, qq(1:), mprnd)
  if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
  d1 = d1 * 2.d0**n1
  n1 = nint (d1)
  call mpbesselynr (n1, rr, t1, mpnw)
  goto 120
elseif (i1 < 0 .and. i2 <= 0) then
  write (mpldb, 2)
2 format ('*** MPBESSELYR: First argument < 0 and second argument <= 0')
  call mpabrt ( 539)
endif

!+ call mpmul (qq, mppicon, t1, mpnw1
call mpfixlocr (mppicon)
call mpfrsetprec (t1(1:), mpnw1)
call mpfrmul (t1(1:), qq(1:), mppicon(1:), mprnd)
 
!+ call mpcssnr (t1, t2, t3, mpnw1
call mpfrsetprec (t3(1:), mpnw1)
call mpfrsincos (t3(1:), t2(1:), t1(1:), mprnd)
 
call mpbesseljr (qq, rr, t4, mpnw1)
!+ call mpmul (t4, t2, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrmul (t1(1:), t4(1:), t2(1:), mprnd)
 
!+ call mpneg (qq, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfrneg (t2(1:), qq(1:), mprnd)
 
call mpbesseljr (t2, rr, t4, mpnw1)
!+ call mpsub (t1, t4, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfrsub (t2(1:), t1(1:), t4(1:), mprnd)
 
!+ call mpdiv (t2, t3, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrdiv (t1(1:), t2(1:), t3(1:), mprnd)
 

120 continue

!+ call mproun (t1, mpnw
call mpfrprecround (t1(1:), mpnw, mprnd)
 
!+ call mpeq (t1, ss, mpnw
call mpfrsetprec (ss(1:), mpnw)
call mpfrset (ss(1:), t1(1:), mprnd)
 
return
end subroutine mpbesselyr

subroutine mpdigammabe (nb1, nb2, berne, x, y, mpnw)

!  This evaluates the digamma function, using asymptotic formula DLMF 5.11.2:
!  dig(x) ~ log(x) - 1/(2*x) - Sum_{k=1}^inf B[2k] / (2*k*x^(2*k)).
!  Before using this formula, the recursion dig(x+1) = dig(x) + 1/x is used
!  to shift the argument up by IQQ, where IQQ is set based on MPNW below.
!  The array berne contains precomputed even Bernoulli numbers (see MPBERNER
!  above). Its dimensions must be as shown below. NB2 must be greater than
!  1.4 x precision in decimal digits.

implicit none
integer, intent (in):: nb1, nb2, mpnw
integer (mpiknd), intent(in):: berne(0:nb1+5,nb2), x(0:)
integer (mpiknd), intent(out):: y(0:)
real (mprknd), parameter:: dber = 1.4d0, dfrac = 0.4d0
integer k, i1, i2, ic1, iqq, mpnw1, n1
real (mprknd) d1
integer (mpiknd) f1(0:mpnw+6), sum1(0:mpnw+6), sum2(0:mpnw+6), &
  t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), &
  t5(0:mpnw+6), xq(0:mpnw+6), tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)

integer, external:: mpfrsgn
integer, external:: mpfrcmp
real (mprknd), external:: mpfrgetd2exp
integer (mpiknd) ix8
 
!  End of declaration

if (mpnw < 4 .or. mpspacer (x) < mpnw + 4 .or. mpspacer (y) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPDIGAMMABE: Uninitialized or inadequately sized arrays')
  call mpabrt ( 540)
endif

!  End of initial input array check

mpnw1 = min (mpnw + 1, mpnwx)
iqq = dfrac * mpnw1 * mpdpw
call mpinitwds (f1, mpnw1)
call mpinitwds (sum1, mpnw1)
call mpinitwds (sum2, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (xq, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
!+ call mpdmc (2.d0, 0, tc1, 4
call mpfrsetprec (tc1(1:), 4)
call mpfrsetd (tc1(1:), 2.d0, mprnd)
 
!+ call mpnpwr (tc1, -mpnw*mpnbt, eps, 4
call mpfrsetprec (eps(1:), 4)
call mpfrpowsi (eps(1:), tc1(1:), -mpnw*mpnbt, mprnd)
 
!+ call mpdmc (1.d0, 0, f1, mpnw1
call mpfrsetprec (f1(1:), mpnw1)
call mpfrsetd (f1(1:), 1.d0, mprnd)
 

!   Check if argument is less than or equal to 0 -- undefined.

!+ if (mpsgn (x) <= 0) the
if (mpfrsgn (x(1:)) <= 0) then
 
  write (mpldb, 2)
2 format ('*** MPDIGAMMABE: Argument is less than or equal to 0')
  call mpabrt ( 541)
endif

!   Check if berne array has been initialized.

!+ call mpmdc (berne(0:nb1+5,1), d1, n1, mpnw
d1 = 2.d0 * mpfrgetd2exp (ix8, berne(1:nb1+5,1), mprnd)
if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
d1 = d1 * 2.d0 ** n1
if (mpwprecr (berne(0:nb1+5,1)) < mpnw .or. &
  abs (d1 - 1.d0 / 6.d0) > mprdfz .or. nb2 < int (dber * mpdpw * mpnw)) then
  write (mpldb, 3) int (dber * mpdpw * mpnw)
3 format ('*** MPDIGAMMABE: Array of even Bernoulli coefficients must be initialized'/ &
   'with at least',i8,' entries using MPBERNE or MPBERNER.')
  call mpabrt ( 542)
endif

! sum1 = mpreal (0.d0, nwds)
! sum2 = mpreal (0.d0, nwds)
! xq = x + dble (iqq)

!+ call mpdmc (0.d0, 0, sum1, mpnw1
call mpfrsetprec (sum1(1:), mpnw1)
call mpfrsetd (sum1(1:), 0.d0, mprnd)
 
!+ call mpdmc (0.d0, 0, sum2, mpnw1
call mpfrsetprec (sum2(1:), mpnw1)
call mpfrsetd (sum2(1:), 0.d0, mprnd)
 
!+ call mpdmc (dble (iqq), 0, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrsetd (t1(1:), dble (iqq), mprnd)
 
!+ call mpadd (x, t1, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfradd (t2(1:), x(1:), t1(1:), mprnd)
 
!+ call mpeq (t2, xq, mpnw1
call mpfrsetprec (xq(1:), mpnw1)
call mpfrset (xq(1:), t2(1:), mprnd)
 

do k = 0, iqq - 1
!   sum1 = sum1 + f1 / (x + dble (k))

!+   call mpdmc (dble (k), 0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrsetd (t1(1:), dble (k), mprnd)
 
!+   call mpadd (x, t1, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfradd (t2(1:), x(1:), t1(1:), mprnd)
 
!+   call mpdiv (f1, t2, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrdiv (t3(1:), f1(1:), t2(1:), mprnd)
 
!+   call mpadd (sum1, t3, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfradd (t1(1:), sum1(1:), t3(1:), mprnd)
 
!+   call mpeq (t1, sum1, mpnw1
  call mpfrsetprec (sum1(1:), mpnw1)
  call mpfrset (sum1(1:), t1(1:), mprnd)
 
enddo

! t1 = mpreal (1.d0, nwds)
! t2 = xq ** 2

!+ call mpdmc (1.d0, 0, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrsetd (t1(1:), 1.d0, mprnd)
 
!+ call mpmul (xq, xq, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfrmul (t2(1:), xq(1:), xq(1:), mprnd)
 

do k = 1, nb2
!  t1 = t1 * t2
!  t3 = bb(k) / (2.d0 * dble (k) * t1)
!  sum2 = sum2 + t3

!+   call mpmul (t1, t2, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrmul (t3(1:), t1(1:), t2(1:), mprnd)
 
!+   call mpeq (t3, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrset (t1(1:), t3(1:), mprnd)
 
!+   call mpmuld (t1, 2.d0 * dble (k), t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfrmuld (t4(1:), t1(1:), 2.d0 * dble (k), mprnd)
 
!+   call mpdiv (berne(0:nb1+5,k), t4, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrdiv (t3(1:), berne(1:nb1+5,k), t4(1:), mprnd)
 
!+   call mpadd (sum2, t3, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfradd (t4(1:), sum2(1:), t3(1:), mprnd)
 
!+   call mpeq (t4, sum2, mpnw1
  call mpfrsetprec (sum2(1:), mpnw1)
  call mpfrset (sum2(1:), t4(1:), mprnd)
 

!  if (abs (t3 / sum2) < eps) goto 100

!+   call mpabs (t3, tc1, 4
  call mpfrsetprec (tc1(1:), 4)
  call mpfrabs (tc1(1:), t3(1:), mprnd)
 
!+   call mpmul (eps, sum2, tc3, 4
  call mpfrsetprec (tc3(1:), 4)
  call mpfrmul (tc3(1:), eps(1:), sum2(1:), mprnd)
 
!+   call mpabs (tc3, tc2, 4
  call mpfrsetprec (tc2(1:), 4)
  call mpfrabs (tc2(1:), tc3(1:), mprnd)
 
!+   call mpcpr (tc1, tc2, ic1, 4
  ic1 = mpfrcmp (tc1(1:), tc2(1:))
 
  if (ic1 <= 0) goto 110
enddo

write (mpldb, 4)
4 format ('*** MPDIGAMMABE: Loop end error: Increase NB2')
call mpabrt ( 543)

110 continue

! digammax = -sum1 + log (xq) - 1.d0 / (2.d0 * xq) - sum2

!+ call mpneg (sum1, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrneg (t1(1:), sum1(1:), mprnd)
 
!+ call mplog (xq, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfrlog (t2(1:), xq(1:), mprnd)
 
!+ call mpadd (t1, t2, t3, mpnw1
call mpfrsetprec (t3(1:), mpnw1)
call mpfradd (t3(1:), t1(1:), t2(1:), mprnd)
 
!+ call mpmuld (xq, 2.d0, t4, mpnw1
call mpfrsetprec (t4(1:), mpnw1)
call mpfrmuld (t4(1:), xq(1:), 2.d0, mprnd)
 
!+ call mpdiv (f1, t4, t5, mpnw1
call mpfrsetprec (t5(1:), mpnw1)
call mpfrdiv (t5(1:), f1(1:), t4(1:), mprnd)
 
!+ call mpsub (t3, t5, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfrsub (t2(1:), t3(1:), t5(1:), mprnd)
 
!+ call mpsub (t2, sum2, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrsub (t1(1:), t2(1:), sum2(1:), mprnd)
 
!+ call mproun (t1, mpnw
call mpfrprecround (t1(1:), mpnw, mprnd)
 
!+ call mpeq (t1, y, mpnw
call mpfrsetprec (y(1:), mpnw)
call mpfrset (y(1:), t1(1:), mprnd)
 
return
end subroutine mpdigammabe

subroutine mperfr (z, terf, mpnw)

!   This evaluates the erf function, using a combination of two series.
!   In particular, the algorithm is (where B = (mpnw + 1) * mpnbt, and
!   dcon is a constant defined below):

!   if (t == 0) then
!     erf = 0
!   elseif (z > sqrt(B*log(2))) then
!     erf = 1
!   elseif (z < -sqrt(B*log(2))) then
!     erf = -1
!   elseif (abs(z) < B/dcon + 8) then
!     erf = 2 / (sqrt(pi)*exp(z^2)) * Sum_{k>=0} 2^k * z^(2*k+1)
!             / (1.3....(2*k+1))
!   else
!     erf = 1 - 1 / (sqrt(pi)*exp(z^2))
!             * Sum_{k>=0} (-1)^k * (1.3...(2*k-1)) / (2^k * z^(2*k+1))
!   endif

implicit none
integer (mpiknd), intent(in):: z(0:)
integer (mpiknd), intent(out):: terf(0:)
integer, intent(in):: mpnw
integer, parameter:: itrmx = 100000
real (mprknd), parameter:: dcon = 100.d0, pi = 3.1415926535897932385d0
integer ic1, ic2, ic3, k, mpnw1, nbt, n1
real (mprknd) d1, d2

integer (mpiknd) t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), &
  t5(0:mpnw+6), t6(0:mpnw+6), t7(0:mpnw+6), z2(0:mpnw+6), tc1(0:9), &
  tc2(0:9), tc3(0:9), eps(0:9)

integer, external:: mpfrsgn
integer, external:: mpfrcmp
real (mprknd), external:: mpfrgetd2exp
integer (mpiknd) ix8
 
!  End of declaration

if (mpnw < 4 .or. mpspacer (z) < mpnw + 4 .or. mpspacer (terf) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPERFR: Uninitialized or inadequately sized arrays')
  call mpabrt ( 544)
endif

!  Check if PI has been precomputed.

!+ call mpmdc (mppicon, d1, n1, mpnw+1
call mpfixlocr (mppicon)
d1 = 2.d0 * mpfrgetd2exp (ix8, mppicon(1:), mprnd)
if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
if (n1 /= 1 .or. abs (d1 * 2.d0**n1 - pi) > mprdfz &
  .or. mpwprecr (mppicon) < mpnw+1) then
  write (mpldb, 4) mpnw+1
4 format ('*** MPERFR: Pi must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 545)
endif

!  End of initial input array check

mpnw1 = min (mpnw + 1, mpnwx)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (t6, mpnw1)
call mpinitwds (t7, mpnw1)
call mpinitwds (z2, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
!+ call mpdmc (2.d0, 0, tc1, 4
call mpfrsetprec (tc1(1:), 4)
call mpfrsetd (tc1(1:), 2.d0, mprnd)
 
!+ call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4
call mpfrsetprec (eps(1:), 4)
call mpfrpowsi (eps(1:), tc1(1:), -mpnw1*mpnbt, mprnd)
 

nbt = mpnw * mpnbt
d1 = aint (1.d0 + sqrt (nbt * log (2.d0)))
d2 = aint (nbt / dcon + 8.d0)
!+ call mpdmc (d1, 0, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrsetd (t1(1:), d1, mprnd)
 
!+ call mpdmc (d2, 0, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfrsetd (t2(1:), d2, mprnd)
 
!+ call mpcpr (z, t1, ic1, mpnw1
ic1 = mpfrcmp (z(1:), t1(1:))
 
! t1(2) = - t1(2)
!+ call mpneg (t1, t3, mpnw1
call mpfrsetprec (t3(1:), mpnw1)
call mpfrneg (t3(1:), t1(1:), mprnd)
 
!+ call mpeq (t3, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrset (t1(1:), t3(1:), mprnd)
 
!+ call mpcpr (z, t1, ic2, mpnw1
ic2 = mpfrcmp (z(1:), t1(1:))
 
!+ call mpcpr (z, t2, ic3, mpnw1
ic3 = mpfrcmp (z(1:), t2(1:))
 

!+ if (mpsgn (z) == 0) the
if (mpfrsgn (z(1:)) == 0) then
 
!+   call mpdmc (0.d0, 0, terf, mpnw
  call mpfrsetprec (terf(1:), mpnw)
  call mpfrsetd (terf(1:), 0.d0, mprnd)
 
elseif (ic1 > 0) then
!+   call mpdmc (1.d0, 0, terf, mpnw
  call mpfrsetprec (terf(1:), mpnw)
  call mpfrsetd (terf(1:), 1.d0, mprnd)
 
elseif (ic2 < 0) then
!+   call mpdmc (-1.d0, 0, terf, mpnw
  call mpfrsetprec (terf(1:), mpnw)
  call mpfrsetd (terf(1:), -1.d0, mprnd)
 
elseif (ic3 < 0) then
!+   call mpmul (z, z, z2, mpnw1
  call mpfrsetprec (z2(1:), mpnw1)
  call mpfrmul (z2(1:), z(1:), z(1:), mprnd)
 
!+   call mpdmc (0.d0, 0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrsetd (t1(1:), 0.d0, mprnd)
 
!+   call mpeq (z, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrset (t2(1:), z(1:), mprnd)
 
!+   call mpdmc (1.d0, 0, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrsetd (t3(1:), 1.d0, mprnd)
 
!+   call mpdmc (1.d10, 0, t5, 4
  call mpfrsetprec (t5(1:), 4)
  call mpfrsetd (t5(1:), 1.d10, mprnd)
 

  do k = 0, itrmx
    if (k > 0) then
!+       call mpmuld (z2, 2.d0, t6, mpnw1
      call mpfrsetprec (t6(1:), mpnw1)
      call mpfrmuld (t6(1:), z2(1:), 2.d0, mprnd)
 
!+       call mpmul (t6, t2, t7, mpnw1
      call mpfrsetprec (t7(1:), mpnw1)
      call mpfrmul (t7(1:), t6(1:), t2(1:), mprnd)
 
!+       call mpeq (t7, t2, mpnw1
      call mpfrsetprec (t2(1:), mpnw1)
      call mpfrset (t2(1:), t7(1:), mprnd)
 
      d1 = 2.d0 * k + 1.d0
!+       call mpmuld (t3, d1, t6, mpnw1
      call mpfrsetprec (t6(1:), mpnw1)
      call mpfrmuld (t6(1:), t3(1:), d1, mprnd)
 
!+       call mpeq (t6, t3, mpnw1
      call mpfrsetprec (t3(1:), mpnw1)
      call mpfrset (t3(1:), t6(1:), mprnd)
 
    endif

!+     call mpdiv (t2, t3, t4, mpnw1
    call mpfrsetprec (t4(1:), mpnw1)
    call mpfrdiv (t4(1:), t2(1:), t3(1:), mprnd)
 
!+     call mpadd (t1, t4, t6, mpnw1
    call mpfrsetprec (t6(1:), mpnw1)
    call mpfradd (t6(1:), t1(1:), t4(1:), mprnd)
 
!+     call mpeq (t6, t1, mpnw1
    call mpfrsetprec (t1(1:), mpnw1)
    call mpfrset (t1(1:), t6(1:), mprnd)
 
!+     call mpdiv (t4, t1, t6, 4
    call mpfrsetprec (t6(1:), 4)
    call mpfrdiv (t6(1:), t4(1:), t1(1:), mprnd)
 
!+     call mpcpr (t6, eps, ic1, 4
    ic1 = mpfrcmp (t6(1:), eps(1:))
 
!+     call mpcpr (t6, t5, ic2, 4
    ic2 = mpfrcmp (t6(1:), t5(1:))
 
    if (ic1 <= 0 .or. ic2 >= 0) goto 120
!+     call mpeq (t6, t5, 4
    call mpfrsetprec (t5(1:), 4)
    call mpfrset (t5(1:), t6(1:), mprnd)
 
  enddo

write (mpldb, 3) 1, itrmx
3 format ('*** MPERFR: iteration limit exceeded',2i10)
call mpabrt ( 546)

120 continue

!+   call mpmuld (t1, 2.d0, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrmuld (t3(1:), t1(1:), 2.d0, mprnd)
 
!+   call mpsqrt (mppicon, t4, mpnw1
call mpfixlocr (mppicon)
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfrsqrt (t4(1:), mppicon(1:), mprnd)
 
!+   call mpexp (z2, t5, mpnw1
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfrexp (t5(1:), z2(1:), mprnd)
 
!+   call mpmul (t4, t5, t6, mpnw1
  call mpfrsetprec (t6(1:), mpnw1)
  call mpfrmul (t6(1:), t4(1:), t5(1:), mprnd)
 
!+   call mpdiv (t3, t6, t7, mpnw1
  call mpfrsetprec (t7(1:), mpnw1)
  call mpfrdiv (t7(1:), t3(1:), t6(1:), mprnd)
 
!+   call mproun (t7, mpnw
  call mpfrprecround (t7(1:), mpnw, mprnd)
 
!+   call mpeq (t7, terf, mpnw
  call mpfrsetprec (terf(1:), mpnw)
  call mpfrset (terf(1:), t7(1:), mprnd)
 
else
!+   call mpmul (z, z, z2, mpnw1
  call mpfrsetprec (z2(1:), mpnw1)
  call mpfrmul (z2(1:), z(1:), z(1:), mprnd)
 
!+   call mpdmc (0.d0, 0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrsetd (t1(1:), 0.d0, mprnd)
 
!+   call mpdmc (1.d0, 0, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrsetd (t2(1:), 1.d0, mprnd)
 
!+   call mpabs (z, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrabs (t3(1:), z(1:), mprnd)
 
!+   call mpdmc (1.d10, 0, t5, 4
  call mpfrsetprec (t5(1:), 4)
  call mpfrsetd (t5(1:), 1.d10, mprnd)
 

  do k = 0, itrmx
    if (k > 0) then
      d1 = -(2.d0 * k - 1.d0)
!+       call mpmuld (t2, d1, t6, mpnw1
      call mpfrsetprec (t6(1:), mpnw1)
      call mpfrmuld (t6(1:), t2(1:), d1, mprnd)
 
!+       call mpeq (t6, t2, mpnw1
      call mpfrsetprec (t2(1:), mpnw1)
      call mpfrset (t2(1:), t6(1:), mprnd)
 
!+       call mpmul (t2, t3, t6, mpnw1
      call mpfrsetprec (t6(1:), mpnw1)
      call mpfrmul (t6(1:), t2(1:), t3(1:), mprnd)
 
!+       call mpeq (t6, t3, mpnw1
      call mpfrsetprec (t3(1:), mpnw1)
      call mpfrset (t3(1:), t6(1:), mprnd)
 
    endif

!+     call mpdiv (t2, t3, t4, mpnw1
    call mpfrsetprec (t4(1:), mpnw1)
    call mpfrdiv (t4(1:), t2(1:), t3(1:), mprnd)
 
!+     call mpadd (t1, t4, t6, mpnw1
    call mpfrsetprec (t6(1:), mpnw1)
    call mpfradd (t6(1:), t1(1:), t4(1:), mprnd)
 
!+     call mpeq (t6, t1, mpnw1
    call mpfrsetprec (t1(1:), mpnw1)
    call mpfrset (t1(1:), t6(1:), mprnd)
 
!+     call mpdiv (t4, t1, t6, 4
    call mpfrsetprec (t6(1:), 4)
    call mpfrdiv (t6(1:), t4(1:), t1(1:), mprnd)
 
!+     call mpcpr (t6, eps, ic1, 4
    ic1 = mpfrcmp (t6(1:), eps(1:))
 
!+     call mpcpr (t6, t5, ic2, 4
    ic2 = mpfrcmp (t6(1:), t5(1:))
 
    if (ic1 <= 0 .or. ic2 >= 0) goto 130
!+     call mpeq (t6, t5, 4
    call mpfrsetprec (t5(1:), 4)
    call mpfrset (t5(1:), t6(1:), mprnd)
 
  enddo

write (mpldb, 3) 2, itrmx
call mpabrt ( 547)

130 continue

!+   call mpdmc (1.d0, 0, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrsetd (t2(1:), 1.d0, mprnd)
 
!+   call mpsqrt (mppicon, t3, mpnw1
call mpfixlocr (mppicon)
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrsqrt (t3(1:), mppicon(1:), mprnd)
 
!+   call mpexp (z2, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfrexp (t4(1:), z2(1:), mprnd)
 
!+   call mpmul (t3, t4, t5, mpnw1
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfrmul (t5(1:), t3(1:), t4(1:), mprnd)
 
!+   call mpdiv (t1, t5, t6, mpnw1
  call mpfrsetprec (t6(1:), mpnw1)
  call mpfrdiv (t6(1:), t1(1:), t5(1:), mprnd)
 
!+   call mpsub (t2, t6, t7, mpnw1
  call mpfrsetprec (t7(1:), mpnw1)
  call mpfrsub (t7(1:), t2(1:), t6(1:), mprnd)
 
!+   call mproun (t7, mpnw
  call mpfrprecround (t7(1:), mpnw, mprnd)
 
!+   call mpeq (t7, terf, mpnw
  call mpfrsetprec (terf(1:), mpnw)
  call mpfrset (terf(1:), t7(1:), mprnd)
 
!+   if (mpsgn (z) < 0) the
  if (mpfrsgn (z(1:)) < 0) then
 
!+     call mpneg (terf, t6, mpnw
    call mpfrsetprec (t6(1:), mpnw)
    call mpfrneg (t6(1:), terf(1:), mprnd)
 
!+     call mpeq (t6, terf, mpnw
    call mpfrsetprec (terf(1:), mpnw)
    call mpfrset (terf(1:), t6(1:), mprnd)
 
  endif
endif

return
end subroutine mperfr

subroutine mperfcr (z, terfc, mpnw)

!   This evaluates the erfc function, using a combination of two series.
!   In particular, the algorithm is (where B = (mpnw + 1) * mpnbt, and
!   dcon is a constant defined below):

!   if (t == 0) then
!     erfc = 1
!   elseif (z > sqrt(B*log(2))) then
!     erfc = 0
!   elseif (z < -sqrt(B*log(2))) then
!     erfc = 2
!   elseif (abs(z) < B/dcon + 8) then
!     erfc = 1 - 2 / (sqrt(pi)*exp(z^2)) * Sum_{k>=0} 2^k * z^(2*k+1)
!               / (1.3....(2*k+1))
!   else
!     erfc = 1 / (sqrt(pi)*exp(z^2))
!             * Sum_{k>=0} (-1)^k * (1.3...(2*k-1)) / (2^k * z^(2*k+1))
!   endif

implicit none
integer (mpiknd), intent(in):: z(0:)
integer (mpiknd), intent(out):: terfc(0:)
integer, intent(in):: mpnw
integer, parameter:: itrmx = 100000
real (mprknd), parameter:: dcon = 100.d0, pi = 3.1415926535897932385d0
integer ic1, ic2, ic3, k, mpnw1, nbt, n1
real (mprknd) d1, d2
integer (mpiknd) t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), &
  t5(0:mpnw+6), t6(0:mpnw+6), t7(0:mpnw+6), z2(0:mpnw+6), tc1(0:9), &
  tc2(0:9), tc3(0:9), eps(0:9)

integer, external:: mpfrsgn
integer, external:: mpfrcmp
real (mprknd), external:: mpfrgetd2exp
integer (mpiknd) ix8
 
!  End of declaration

if (mpnw < 4 .or. mpspacer (z) < mpnw + 4 .or. mpspacer (terfc) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPERFCR: Uninitialized or inadequately sized arrays')
  call mpabrt ( 548)
endif

!  Check if PI has been precomputed.

!+ call mpmdc (mppicon, d1, n1, mpnw+1
call mpfixlocr (mppicon)
d1 = 2.d0 * mpfrgetd2exp (ix8, mppicon(1:), mprnd)
if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
if (n1 /= 1 .or. abs (d1 * 2.d0**n1 - pi) > mprdfz &
  .or. mpwprecr (mppicon) < mpnw+1) then
  write (mpldb, 4) mpnw+1
4 format ('*** MPERFCR: Pi must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 549)
endif

!  End of initial input array check

mpnw1 = min (mpnw + 1, mpnwx)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (t6, mpnw1)
call mpinitwds (t7, mpnw1)
call mpinitwds (z2, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
!+ call mpdmc (2.d0, 0, tc1, 4
call mpfrsetprec (tc1(1:), 4)
call mpfrsetd (tc1(1:), 2.d0, mprnd)
 
!+ call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4
call mpfrsetprec (eps(1:), 4)
call mpfrpowsi (eps(1:), tc1(1:), -mpnw1*mpnbt, mprnd)
 

nbt = mpnw * mpnbt
d1 = aint (1.d0 + sqrt (nbt * log (2.d0)))
d2 = aint (nbt / dcon + 8.d0)
!+ call mpdmc (d1, 0, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrsetd (t1(1:), d1, mprnd)
 
!+ call mpdmc (d2, 0, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfrsetd (t2(1:), d2, mprnd)
 
!+ call mpcpr (z, t1, ic1, mpnw1
ic1 = mpfrcmp (z(1:), t1(1:))
 
!+ call mpneg (t1, t3, mpnw1
call mpfrsetprec (t3(1:), mpnw1)
call mpfrneg (t3(1:), t1(1:), mprnd)
 
!+ call mpeq (t3, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrset (t1(1:), t3(1:), mprnd)
 
!+ call mpcpr (z, t1, ic2, mpnw1
ic2 = mpfrcmp (z(1:), t1(1:))
 
!+ call mpcpr (z, t2, ic3, mpnw1
ic3 = mpfrcmp (z(1:), t2(1:))
 

!+ if (mpsgn (z) == 0) the
if (mpfrsgn (z(1:)) == 0) then
 
!+   call mpdmc (1.d0, 0, terfc, mpnw
  call mpfrsetprec (terfc(1:), mpnw)
  call mpfrsetd (terfc(1:), 1.d0, mprnd)
 
elseif (ic1 > 0) then
!+   call mpdmc (0.d0, 0, terfc, mpnw
  call mpfrsetprec (terfc(1:), mpnw)
  call mpfrsetd (terfc(1:), 0.d0, mprnd)
 
elseif (ic2 < 0) then
!+   call mpdmc (2.d0, 0, terfc, mpnw
  call mpfrsetprec (terfc(1:), mpnw)
  call mpfrsetd (terfc(1:), 2.d0, mprnd)
 
elseif (ic3 < 0) then
!+   call mpmul (z, z, z2, mpnw1
  call mpfrsetprec (z2(1:), mpnw1)
  call mpfrmul (z2(1:), z(1:), z(1:), mprnd)
 
!+   call mpdmc (0.d0, 0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrsetd (t1(1:), 0.d0, mprnd)
 
!+   call mpeq (z, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrset (t2(1:), z(1:), mprnd)
 
!+   call mpdmc (1.d0, 0, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrsetd (t3(1:), 1.d0, mprnd)
 
!+   call mpdmc (1.d10, 0, t5, 4
  call mpfrsetprec (t5(1:), 4)
  call mpfrsetd (t5(1:), 1.d10, mprnd)
 

  do k = 0, itrmx
    if (k > 0) then
!+       call mpmuld (z2, 2.d0, t6, mpnw1
      call mpfrsetprec (t6(1:), mpnw1)
      call mpfrmuld (t6(1:), z2(1:), 2.d0, mprnd)
 
!+       call mpmul (t6, t2, t7, mpnw1
      call mpfrsetprec (t7(1:), mpnw1)
      call mpfrmul (t7(1:), t6(1:), t2(1:), mprnd)
 
!+       call mpeq (t7, t2, mpnw1
      call mpfrsetprec (t2(1:), mpnw1)
      call mpfrset (t2(1:), t7(1:), mprnd)
 
      d1 = 2.d0 * k + 1.d0
!+       call mpmuld (t3, d1, t6, mpnw1
      call mpfrsetprec (t6(1:), mpnw1)
      call mpfrmuld (t6(1:), t3(1:), d1, mprnd)
 
!+       call mpeq (t6, t3, mpnw1
      call mpfrsetprec (t3(1:), mpnw1)
      call mpfrset (t3(1:), t6(1:), mprnd)
 
    endif

!+     call mpdiv (t2, t3, t4, mpnw1
    call mpfrsetprec (t4(1:), mpnw1)
    call mpfrdiv (t4(1:), t2(1:), t3(1:), mprnd)
 
!+     call mpadd (t1, t4, t6, mpnw1
    call mpfrsetprec (t6(1:), mpnw1)
    call mpfradd (t6(1:), t1(1:), t4(1:), mprnd)
 
!+     call mpeq (t6, t1, mpnw1
    call mpfrsetprec (t1(1:), mpnw1)
    call mpfrset (t1(1:), t6(1:), mprnd)
 
!+     call mpdiv (t4, t1, t6, 4
    call mpfrsetprec (t6(1:), 4)
    call mpfrdiv (t6(1:), t4(1:), t1(1:), mprnd)
 
!+     call mpcpr (t6, eps, ic1, 4
    ic1 = mpfrcmp (t6(1:), eps(1:))
 
!+     call mpcpr (t6, t5, ic2, 4
    ic2 = mpfrcmp (t6(1:), t5(1:))
 
    if (ic1 <= 0 .or. ic2 >= 0) goto 120
!+     call mpeq (t6, t5, 4
    call mpfrsetprec (t5(1:), 4)
    call mpfrset (t5(1:), t6(1:), mprnd)
 
  enddo

write (mpldb, 3) 1, itrmx
3 format ('*** MPERFCR: iteration limit exceeded',2i10)
call mpabrt ( 550)

120 continue

!+   call mpdmc (1.d0, 0, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrsetd (t2(1:), 1.d0, mprnd)
 
!+   call mpmuld (t1, 2.d0, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrmuld (t3(1:), t1(1:), 2.d0, mprnd)
 
!+   call mpsqrt (mppicon, t4, mpnw1
call mpfixlocr (mppicon)
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfrsqrt (t4(1:), mppicon(1:), mprnd)
 
!+   call mpexp (z2, t5, mpnw1
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfrexp (t5(1:), z2(1:), mprnd)
 
!+   call mpmul (t4, t5, t6, mpnw1
  call mpfrsetprec (t6(1:), mpnw1)
  call mpfrmul (t6(1:), t4(1:), t5(1:), mprnd)
 
!+   call mpdiv (t3, t6, t7, mpnw1
  call mpfrsetprec (t7(1:), mpnw1)
  call mpfrdiv (t7(1:), t3(1:), t6(1:), mprnd)
 
!+   call mpsub (t2, t7, t6, mpnw1
  call mpfrsetprec (t6(1:), mpnw1)
  call mpfrsub (t6(1:), t2(1:), t7(1:), mprnd)
 
!+   call mproun (t6, mpnw
  call mpfrprecround (t6(1:), mpnw, mprnd)
 
!+   call mpeq (t6, terfc, mpnw
  call mpfrsetprec (terfc(1:), mpnw)
  call mpfrset (terfc(1:), t6(1:), mprnd)
 
else
!+   call mpmul (z, z, z2, mpnw1
  call mpfrsetprec (z2(1:), mpnw1)
  call mpfrmul (z2(1:), z(1:), z(1:), mprnd)
 
!+   call mpdmc (0.d0, 0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrsetd (t1(1:), 0.d0, mprnd)
 
!+   call mpdmc (1.d0, 0, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrsetd (t2(1:), 1.d0, mprnd)
 
!+   call mpabs (z, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrabs (t3(1:), z(1:), mprnd)
 
!+   call mpdmc (1.d10, 0, t5, 4
  call mpfrsetprec (t5(1:), 4)
  call mpfrsetd (t5(1:), 1.d10, mprnd)
 

  do k = 0, itrmx
    if (k > 0) then
      d1 = -(2.d0 * k - 1.d0)
!+       call mpmuld (t2, d1, t6, mpnw1
      call mpfrsetprec (t6(1:), mpnw1)
      call mpfrmuld (t6(1:), t2(1:), d1, mprnd)
 
!+       call mpeq (t6, t2, mpnw1
      call mpfrsetprec (t2(1:), mpnw1)
      call mpfrset (t2(1:), t6(1:), mprnd)
 
!+       call mpmul (t2, t3, t6, mpnw1
      call mpfrsetprec (t6(1:), mpnw1)
      call mpfrmul (t6(1:), t2(1:), t3(1:), mprnd)
 
!+       call mpeq (t6, t3, mpnw1
      call mpfrsetprec (t3(1:), mpnw1)
      call mpfrset (t3(1:), t6(1:), mprnd)
 
    endif

!+     call mpdiv (t2, t3, t4, mpnw1
    call mpfrsetprec (t4(1:), mpnw1)
    call mpfrdiv (t4(1:), t2(1:), t3(1:), mprnd)
 
!+     call mpadd (t1, t4, t6, mpnw1
    call mpfrsetprec (t6(1:), mpnw1)
    call mpfradd (t6(1:), t1(1:), t4(1:), mprnd)
 
!+     call mpeq (t6, t1, mpnw1
    call mpfrsetprec (t1(1:), mpnw1)
    call mpfrset (t1(1:), t6(1:), mprnd)
 
!+     call mpdiv (t4, t1, t6, 4
    call mpfrsetprec (t6(1:), 4)
    call mpfrdiv (t6(1:), t4(1:), t1(1:), mprnd)
 
!+     call mpcpr (t6, eps, ic1, 4
    ic1 = mpfrcmp (t6(1:), eps(1:))
 
!+     call mpcpr (t6, t5, ic2, 4
    ic2 = mpfrcmp (t6(1:), t5(1:))
 
    if (ic1 <= 0 .or. ic2 >= 0) goto 130
!+     call mpeq (t6, t5, 4
    call mpfrsetprec (t5(1:), 4)
    call mpfrset (t5(1:), t6(1:), mprnd)
 
  enddo

write (mpldb, 3) 2, itrmx
call mpabrt ( 551)

130 continue

!+   call mpsqrt (mppicon, t3, mpnw1
call mpfixlocr (mppicon)
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrsqrt (t3(1:), mppicon(1:), mprnd)
 
!+   call mpexp (z2, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfrexp (t4(1:), z2(1:), mprnd)
 
!+   call mpmul (t3, t4, t5, mpnw1
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfrmul (t5(1:), t3(1:), t4(1:), mprnd)
 
!+   call mpdiv (t1, t5, t6, mpnw1
  call mpfrsetprec (t6(1:), mpnw1)
  call mpfrdiv (t6(1:), t1(1:), t5(1:), mprnd)
 
!+   if (mpsgn (z) < 0) the
  if (mpfrsgn (z(1:)) < 0) then
 
!+     call mpdmc (2.d0, 0, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrsetd (t2(1:), 2.d0, mprnd)
 
!+     call mpsub (t2, t6, t7, mpnw1
    call mpfrsetprec (t7(1:), mpnw1)
    call mpfrsub (t7(1:), t2(1:), t6(1:), mprnd)
 
!+     call mpeq (t7, t6, mpnw1
    call mpfrsetprec (t6(1:), mpnw1)
    call mpfrset (t6(1:), t7(1:), mprnd)
 
  endif

!+   call mproun (t6, mpnw
  call mpfrprecround (t6(1:), mpnw, mprnd)
 
!+   call mpeq (t6, terfc, mpnw
  call mpfrsetprec (terfc(1:), mpnw)
  call mpfrset (terfc(1:), t6(1:), mprnd)
 
endif

return
end subroutine mperfcr

subroutine mpexpint (x, y, mpnw)

!   This evaluates the exponential integral function Ei(x):
!   Ei(x) = - incgamma (0, -x)

implicit none
integer (mpiknd), intent(in):: x(0:)
integer (mpiknd), intent(out):: y(0:)
integer (mpiknd) t1(0:mpnw+5), t2(0:mpnw+5), t3(0:mpnw+5)
integer, intent(in):: mpnw

integer, external:: mpfrsgn
integer, external:: mpfrcmp
real (mprknd), external:: mpfrgetd2exp
integer (mpiknd) ix8
 
!  End of declaration

if (mpnw < 4 .or. mpspacer (x) < mpnw + 4 .or. mpspacer (y) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPEXPINT: Uninitialized or inadequately sized arrays')
  call mpabrt ( 552)
endif

!  End of initial input array check

!+ if (mpsgn (x) == 0) the
if (mpfrsgn (x(1:)) == 0) then
 
  write (mpldb, 2)
2 format ('*** MPEXPINT: argument is zero')
  call mpabrt ( 553)
endif

call mpinitwds (t1, mpnw)
call mpinitwds (t2, mpnw)
call mpinitwds (t3, mpnw)
!+ call mpdmc (0.d0, 0, t1, mpnw
call mpfrsetprec (t1(1:), mpnw)
call mpfrsetd (t1(1:), 0.d0, mprnd)
 
!+ call mpneg (x, t2, mpnw
call mpfrsetprec (t2(1:), mpnw)
call mpfrneg (t2(1:), x(1:), mprnd)
 
call mpincgammar (t1, t2, t3, mpnw)
!+ call mpneg (t3, y, mpnw
call mpfrsetprec (y(1:), mpnw)
call mpfrneg (y(1:), t3(1:), mprnd)
 
return
end subroutine mpexpint

subroutine mpgammar (t, z, mpnw)

!   This evaluates the gamma function, using an algorithm of R. W. Potter.
!   The argument t must not exceed 10^8 in size (this limit is set below),
!   must not be zero, and if negative must not be integer.

!   In the parameter statement below:
!     itrmx = limit of number of iterations in series; default = 100000.
!     con1 = 1/2 * log (10) to DP accuracy.
!     dmax = maximum size of input argument.

implicit none
integer (mpiknd), intent(in):: t(0:)
integer (mpiknd), intent(out):: z(0:)
integer, intent(in):: mpnw
integer, parameter:: itrmx = 100000
real (mprknd), parameter:: al2 = 0.69314718055994530942d0, dmax = 1d8, &
  pi = 3.1415926535897932385d0
integer i, i1, ic1, j, mpnw1, nt, n1, n2, n3
real (mprknd) alpha, d1, d2, d3
integer (mpiknd) f1(0:mpnw+6), sum1(0:mpnw+6), sum2(0:mpnw+6), tn(0:mpnw+6), &
  t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), t5(0:mpnw+6), &
  t6(0:mpnw+6), tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)

integer, external:: mpfrsgn
integer, external:: mpfrcmp
real (mprknd), external:: mpfrgetd2exp
integer (mpiknd) ix8
 
!  End of declaration

if (mpnw < 4 .or. mpspacer (t) < mpnw + 4 .or. mpspacer (z) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPGAMMAR: Uninitialized or inadequately sized arrays')
  call mpabrt ( 554)
endif

!  Check if PI has been precomputed.

!+ call mpmdc (mppicon, d1, n1, mpnw+1
call mpfixlocr (mppicon)
d1 = 2.d0 * mpfrgetd2exp (ix8, mppicon(1:), mprnd)
if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
if (n1 /= 1 .or. abs (d1 * 2.d0**n1 - pi) > mprdfz &
  .or. mpwprecr (mppicon) < mpnw+1) then
  write (mpldb, 4) mpnw+1
4 format ('*** MPGAMMAR: Pi must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 555)
endif

!  End of initial input array check

mpnw1 = min (mpnw + 1, mpnwx)
call mpinitwds (f1, mpnw1)
call mpinitwds (sum1, mpnw1)
call mpinitwds (sum2, mpnw1)
call mpinitwds (tn, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (t6, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
!+ call mpdmc (2.d0, 0, tc1, 4
call mpfrsetprec (tc1(1:), 4)
call mpfrsetd (tc1(1:), 2.d0, mprnd)
 
!+ call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4
call mpfrsetprec (eps(1:), 4)
call mpfrpowsi (eps(1:), tc1(1:), -mpnw1*mpnbt, mprnd)
 

!+ call mpmdc (t, d1, n1, mpnw
d1 = 2.d0 * mpfrgetd2exp (ix8, t(1:), mprnd)
if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
d1 = d1 * 2.d0**n1
!+ call mpnint (t, t1, mpnw
call mpfrround (t1(1:), t(1:))
 
!+ call mpcpr (t, t1, ic1, mpnw
ic1 = mpfrcmp (t(1:), t1(1:))
 
!+ i1 = mpsgn (t
i1 = mpfrsgn (t(1:))
 
if (i1 == 0 .or. d1 > dmax .or. (i1 < 0 .and. ic1 == 0)) then
  write (mpldb, 2) dmax
2 format ('*** MPGAMMAR: input argument must have absolute value <=',f10.0,','/ &
  'must not be zero, and if negative must not be an integer.')
  call mpabrt ( 556)
endif

!+ call mpdmc (1.d0, 0, f1, mpnw1
call mpfrsetprec (f1(1:), mpnw1)
call mpfrsetd (f1(1:), 1.d0, mprnd)
 

!   Find the integer and fractional parts of t.

!+ call mpinfr (t, t2, t3, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfrsetprec (t3(1:), mpnw1)
call mpfrtrunc (t2(1:), t(1:))
call mpfrsub (t3(1:), t(1:), t2(1:), mprnd)
 

!+ if (mpsgn (t3) == 0) the
if (mpfrsgn (t3(1:)) == 0) then
 

!   If t is a positive integer, then apply the usual factorial recursion.

!+   call mpmdc (t2, d2, n2, mpnw1
  d2 = 2.d0 * mpfrgetd2exp (ix8, t2(1:), mprnd)
  if (d2 == 0) then; n2 = 0; else; n2 = ix8 - 1; endif
 
  nt = d2 * 2.d0 ** n2
!+   call mpeq (f1, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrset (t1(1:), f1(1:), mprnd)
 

  do i = 2, nt - 1
!+     call mpmuld (t1, dble (i), t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmuld (t2(1:), t1(1:), dble (i), mprnd)
 
!+     call mpeq (t2, t1, mpnw1
    call mpfrsetprec (t1(1:), mpnw1)
    call mpfrset (t1(1:), t2(1:), mprnd)
 
  enddo

!+   call mproun (t1, mpnw
  call mpfrprecround (t1(1:), mpnw, mprnd)
 
!+   call mpeq (t1, z, mpnw
  call mpfrsetprec (z(1:), mpnw)
  call mpfrset (z(1:), t1(1:), mprnd)
 
  goto 120
!+ elseif (mpsgn (t) > 0) the
elseif (mpfrsgn (t(1:)) > 0) then
 

!   Apply the identity Gamma[t+1] = t * Gamma[t] to reduce the input argument
!   to the unit interval.

!+   call mpmdc (t2, d2, n2, mpnw1
  d2 = 2.d0 * mpfrgetd2exp (ix8, t2(1:), mprnd)
  if (d2 == 0) then; n2 = 0; else; n2 = ix8 - 1; endif
 
  nt = d2 * 2.d0 ** n2
!+   call mpeq (f1, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrset (t1(1:), f1(1:), mprnd)
 
!+   call mpeq (t3, tn, mpnw1
  call mpfrsetprec (tn(1:), mpnw1)
  call mpfrset (tn(1:), t3(1:), mprnd)
 

  do i = 1, nt
!+     call mpdmc (dble (i), 0, t4, mpnw1
    call mpfrsetprec (t4(1:), mpnw1)
    call mpfrsetd (t4(1:), dble (i), mprnd)
 
!+     call mpsub (t, t4, t5, mpnw1
    call mpfrsetprec (t5(1:), mpnw1)
    call mpfrsub (t5(1:), t(1:), t4(1:), mprnd)
 
!+     call mpmul (t1, t5, t6, mpnw1
    call mpfrsetprec (t6(1:), mpnw1)
    call mpfrmul (t6(1:), t1(1:), t5(1:), mprnd)
 
!+     call mpeq (t6, t1, mpnw1
    call mpfrsetprec (t1(1:), mpnw1)
    call mpfrset (t1(1:), t6(1:), mprnd)
 
  enddo
else

!   Apply the gamma identity to reduce a negative argument to the unit interval.

!+   call mpsub (f1, t, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfrsub (t4(1:), f1(1:), t(1:), mprnd)
 
!+   call mpinfr (t4, t3, t5, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfrtrunc (t3(1:), t4(1:))
  call mpfrsub (t5(1:), t4(1:), t3(1:), mprnd)
 
!+   call mpmdc (t3, d3, n3, mpnw1
  d3 = 2.d0 * mpfrgetd2exp (ix8, t3(1:), mprnd)
  if (d3 == 0) then; n3 = 0; else; n3 = ix8 - 1; endif
 
  nt = d3 * 2.d0 ** n3

!+   call mpeq (f1, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrset (t1(1:), f1(1:), mprnd)
 
!+   call mpsub (f1, t5, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrsub (t2(1:), f1(1:), t5(1:), mprnd)
 
!+   call mpeq (t2, tn, mpnw1
  call mpfrsetprec (tn(1:), mpnw1)
  call mpfrset (tn(1:), t2(1:), mprnd)
 

  do i = 0, nt - 1
!+     call mpdmc (dble (i), 0, t4, mpnw1
    call mpfrsetprec (t4(1:), mpnw1)
    call mpfrsetd (t4(1:), dble (i), mprnd)
 
!+     call mpadd (t, t4, t5, mpnw1
    call mpfrsetprec (t5(1:), mpnw1)
    call mpfradd (t5(1:), t(1:), t4(1:), mprnd)
 
!+     call mpdiv (t1, t5, t6, mpnw1
    call mpfrsetprec (t6(1:), mpnw1)
    call mpfrdiv (t6(1:), t1(1:), t5(1:), mprnd)
 
!+     call mpeq (t6, t1, mpnw1
    call mpfrsetprec (t1(1:), mpnw1)
    call mpfrset (t1(1:), t6(1:), mprnd)
 
  enddo
endif

!   Calculate alpha = bits of precision * log(2) / 2, then take the next even
!   integer value, so that alpha/2 and alpha^2/4 can be calculated exactly in DP.

alpha = 2.d0 * aint (0.25d0 * (mpnw1 + 1) * mpnbt * al2 + 1.d0)
d2 = 0.25d0 * alpha**2
!+ call mpeq (tn, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfrset (t2(1:), tn(1:), mprnd)
 
!+ call mpdiv (f1, t2, t3, mpnw1
call mpfrsetprec (t3(1:), mpnw1)
call mpfrdiv (t3(1:), f1(1:), t2(1:), mprnd)
 
!+ call mpeq (t3, sum1, mpnw1
call mpfrsetprec (sum1(1:), mpnw1)
call mpfrset (sum1(1:), t3(1:), mprnd)
 

!   Evaluate the series with t.

do j = 1, itrmx
!+   call mpdmc (dble (j), 0, t6, mpnw1
  call mpfrsetprec (t6(1:), mpnw1)
  call mpfrsetd (t6(1:), dble (j), mprnd)
 
!+   call mpadd (t2, t6, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfradd (t4(1:), t2(1:), t6(1:), mprnd)
 
!+   call mpmuld (t4, dble (j), t5, mpnw1
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfrmuld (t5(1:), t4(1:), dble (j), mprnd)
 
!+   call mpdiv (t3, t5, t6, mpnw1
  call mpfrsetprec (t6(1:), mpnw1)
  call mpfrdiv (t6(1:), t3(1:), t5(1:), mprnd)
 
!+   call mpmuld (t6, d2, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrmuld (t3(1:), t6(1:), d2, mprnd)
 
!+   call mpadd (sum1, t3, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfradd (t4(1:), sum1(1:), t3(1:), mprnd)
 
!+   call mpeq (t4, sum1, mpnw1
  call mpfrsetprec (sum1(1:), mpnw1)
  call mpfrset (sum1(1:), t4(1:), mprnd)
 

!+   call mpabs (t3, tc1, 4
  call mpfrsetprec (tc1(1:), 4)
  call mpfrabs (tc1(1:), t3(1:), mprnd)
 
!+   call mpmul (eps, sum1, tc3, 4
  call mpfrsetprec (tc3(1:), 4)
  call mpfrmul (tc3(1:), eps(1:), sum1(1:), mprnd)
 
!+   call mpabs (tc3, tc2, 4
  call mpfrsetprec (tc2(1:), 4)
  call mpfrabs (tc2(1:), tc3(1:), mprnd)
 
!+   call mpcpr (tc1, tc2, ic1, 4
  ic1 = mpfrcmp (tc1(1:), tc2(1:))
 
  if (ic1 <= 0) goto 100
enddo

write (mpldb, 3) 1, itrmx
3 format ('*** MPGAMMAR: iteration limit exceeded',2i10)
call mpabrt ( 557)

100 continue

!+ call mpneg (tn, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfrneg (t2(1:), tn(1:), mprnd)
 
!+ call mpdiv (f1, t2, t3, mpnw1
call mpfrsetprec (t3(1:), mpnw1)
call mpfrdiv (t3(1:), f1(1:), t2(1:), mprnd)
 
!+ call mpeq (t3, sum2, mpnw1
call mpfrsetprec (sum2(1:), mpnw1)
call mpfrset (sum2(1:), t3(1:), mprnd)
 

!   Evaluate the same series with -t.

do j = 1, itrmx
!+   call mpdmc (dble (j), 0, t6, mpnw1
  call mpfrsetprec (t6(1:), mpnw1)
  call mpfrsetd (t6(1:), dble (j), mprnd)
 
!+   call mpadd (t2, t6, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfradd (t4(1:), t2(1:), t6(1:), mprnd)
 
!+   call mpmuld (t4, dble (j), t5, mpnw1
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfrmuld (t5(1:), t4(1:), dble (j), mprnd)
 
!+   call mpdiv (t3, t5, t6, mpnw1
  call mpfrsetprec (t6(1:), mpnw1)
  call mpfrdiv (t6(1:), t3(1:), t5(1:), mprnd)
 
!+   call mpmuld (t6, d2, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrmuld (t3(1:), t6(1:), d2, mprnd)
 
!+   call mpadd (sum2, t3, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfradd (t4(1:), sum2(1:), t3(1:), mprnd)
 
!+   call mpeq (t4, sum2, mpnw1
  call mpfrsetprec (sum2(1:), mpnw1)
  call mpfrset (sum2(1:), t4(1:), mprnd)
 

!+   call mpabs (t3, tc1, 4
  call mpfrsetprec (tc1(1:), 4)
  call mpfrabs (tc1(1:), t3(1:), mprnd)
 
!+   call mpmul (eps, sum2, tc3, 4
  call mpfrsetprec (tc3(1:), 4)
  call mpfrmul (tc3(1:), eps(1:), sum2(1:), mprnd)
 
!+   call mpabs (tc3, tc2, 4
  call mpfrsetprec (tc2(1:), 4)
  call mpfrabs (tc2(1:), tc3(1:), mprnd)
 
!+   call mpcpr (tc1, tc2, ic1, 4
  ic1 = mpfrcmp (tc1(1:), tc2(1:))
 
  if (ic1 <= 0) goto 110
enddo

write (mpldb, 3) 2, itrmx
call mpabrt ( 558)

110 continue

!   Compute sqrt (pi * sum1 / (tn * sin (pi * tn) * sum2))
!   and (alpha/2)^tn terms. Also, multiply by the factor t1, from the
!   If block above.

!+ call mpeq (mppicon, t2, mpnw1
call mpfixlocr (mppicon)
call mpfrsetprec (t2(1:), mpnw1)
call mpfrset (t2(1:), mppicon(1:), mprnd)
 
!+ call mpmul (t2, tn, t3, mpnw1
call mpfrsetprec (t3(1:), mpnw1)
call mpfrmul (t3(1:), t2(1:), tn(1:), mprnd)
 
!+ call mpcssnr (t3, t4, t5, mpnw1
call mpfrsetprec (t5(1:), mpnw1)
call mpfrsincos (t5(1:), t4(1:), t3(1:), mprnd)
 
!+ call mpmul (t5, sum2, t6, mpnw1
call mpfrsetprec (t6(1:), mpnw1)
call mpfrmul (t6(1:), t5(1:), sum2(1:), mprnd)
 
!+ call mpmul (tn, t6, t5, mpnw1
call mpfrsetprec (t5(1:), mpnw1)
call mpfrmul (t5(1:), tn(1:), t6(1:), mprnd)
 
!+ call mpmul (t2, sum1, t3, mpnw1
call mpfrsetprec (t3(1:), mpnw1)
call mpfrmul (t3(1:), t2(1:), sum1(1:), mprnd)
 
!+ call mpdiv (t3, t5, t6, mpnw1
call mpfrsetprec (t6(1:), mpnw1)
call mpfrdiv (t6(1:), t3(1:), t5(1:), mprnd)
 
!+ call mpneg (t6, t4, mpnw1
call mpfrsetprec (t4(1:), mpnw1)
call mpfrneg (t4(1:), t6(1:), mprnd)
 
!+ call mpeq (t4, t6, mpnw1
call mpfrsetprec (t6(1:), mpnw1)
call mpfrset (t6(1:), t4(1:), mprnd)
 
!+ call mpsqrt (t6, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfrsqrt (t2(1:), t6(1:), mprnd)
 
!+ call mpdmc (0.5d0 * alpha, 0, t3, mpnw1
call mpfrsetprec (t3(1:), mpnw1)
call mpfrsetd (t3(1:), 0.5d0 * alpha, mprnd)
 
!+ call mplog (t3, t4, mpnw1
call mpfrsetprec (t4(1:), mpnw1)
call mpfrlog (t4(1:), t3(1:), mprnd)
 
!+ call mpmul (tn, t4, t5, mpnw1
call mpfrsetprec (t5(1:), mpnw1)
call mpfrmul (t5(1:), tn(1:), t4(1:), mprnd)
 
!+ call mpexp (t5, t6, mpnw1
call mpfrsetprec (t6(1:), mpnw1)
call mpfrexp (t6(1:), t5(1:), mprnd)
 
!+ call mpmul (t2, t6, t3, mpnw1
call mpfrsetprec (t3(1:), mpnw1)
call mpfrmul (t3(1:), t2(1:), t6(1:), mprnd)
 
!+ call mpmul (t1, t3, t4, mpnw1
call mpfrsetprec (t4(1:), mpnw1)
call mpfrmul (t4(1:), t1(1:), t3(1:), mprnd)
 

!   Round to mpnw words precision.

!+ call mproun (t4, mpnw
call mpfrprecround (t4(1:), mpnw, mprnd)
 
!+ call mpeq (t4, z, mpnw
call mpfrsetprec (z(1:), mpnw)
call mpfrset (z(1:), t4(1:), mprnd)
 

120 continue

return
end subroutine mpgammar

subroutine mphurwitzzetan (is, aa, zz, mpnw)

!   This returns the Hurwitz zeta function of IS and AA, using an algorithm from:
!   David H. Bailey and Jonathan M. Borwein, "Crandall's computation of the
!   incomplete gamma function and the Hurwitz zeta function with applications to
!   Dirichlet L-series," Applied Mathematics and Computation, vol. 268C (Oct 2015),
!   pg. 462-477, preprint at:
!   https://www.davidhbailey.com/dhbpapers/lerch.pdf
!   This is limited to IS >= 2 and 0 < AA < 1.

implicit none
integer, intent(in):: is, mpnw
integer (mpiknd), intent(in):: aa(0:)
integer (mpiknd), intent(out):: zz(0:)
integer, parameter:: itrmax = 1000000
real (mprknd), parameter:: pi = 3.1415926535897932385d0
integer i1, ic1, ic2, ic3, k, n1, mpnw1
real (mprknd) d1, dk
integer (mpiknd) gs1(0:mpnw+6), gs2(0:mpnw+6), ss(0:mpnw+6), &
  sum1(0:mpnw+6), sum2(0:mpnw+6), sum3(0:mpnw+6), ss1(0:mpnw+6), ss2(0:mpnw+6), &
  ss3(0:mpnw+6), ss4(0:mpnw+6), s1(0:mpnw+6), s2(0:mpnw+6), s3(0:mpnw+6), &
  t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), t5(0:mpnw+6), &
  t6(0:mpnw+6), t7(0:mpnw+6), t8(0:mpnw+6), tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)

integer, external:: mpfrsgn
integer, external:: mpfrcmp
real (mprknd), external:: mpfrgetd2exp
integer (mpiknd) ix8
 
!  End of declaration

if (mpnw < 4 .or. mpspacer (aa) < mpnw + 4 .or. mpspacer (zz) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPHURWITZZETAN: Uninitialized or inadequately sized arrays')
  call mpabrt ( 559)
endif

!  Check if PI has been precomputed.

!+ call mpmdc (mppicon, d1, n1, mpnw+1
call mpfixlocr (mppicon)
d1 = 2.d0 * mpfrgetd2exp (ix8, mppicon(1:), mprnd)
if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
if (n1 /= 1 .or. abs (d1 * 2.d0**n1 - pi) > mprdfz &
  .or. mpwprecr (mppicon) < mpnw+1) then
  write (mpldb, 2) mpnw+1
2 format ('*** MPHURWITZZETAN: Pi must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 560)
endif

!  End of initial input array check

mpnw1 = min (mpnw + 1, mpnwx)
call mpinitwds (gs1, mpnw1)
call mpinitwds (gs2, mpnw1)
call mpinitwds (ss, mpnw1)
call mpinitwds (sum1, mpnw1)
call mpinitwds (sum2, mpnw1)
call mpinitwds (sum3, mpnw1)
call mpinitwds (ss1, mpnw1)
call mpinitwds (ss2, mpnw1)
call mpinitwds (ss3, mpnw1)
call mpinitwds (ss4, mpnw1)
call mpinitwds (s1, mpnw1)
call mpinitwds (s2, mpnw1)
call mpinitwds (s3, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (t6, mpnw1)
call mpinitwds (t7, mpnw1)
call mpinitwds (t8, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
!+ call mpdmc (2.d0, 0, tc1, 4
call mpfrsetprec (tc1(1:), 4)
call mpfrsetd (tc1(1:), 2.d0, mprnd)
 
!+ call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4
call mpfrsetprec (eps(1:), 4)
call mpfrpowsi (eps(1:), tc1(1:), -mpnw1*mpnbt, mprnd)
 

if (is <= 0) then
  write (mpldb, 3)
3 format ('*** MPHURWITZZETAN: IS less than or equal to 0:')
  call mpabrt ( 561)
endif

!+ call mpdmc (0.d0, 0, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrsetd (t1(1:), 0.d0, mprnd)
 
!+ call mpdmc (1.d0, 0, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfrsetd (t2(1:), 1.d0, mprnd)
 
!+ call mpcpr (aa, t1, ic1, mpnw1
ic1 = mpfrcmp (aa(1:), t1(1:))
 
!+ call mpcpr (aa, t2, ic2, mpnw1
ic2 = mpfrcmp (aa(1:), t2(1:))
 
if (ic1 <= 0 .or. ic2 >= 0) then
  write (mpldb, 4)
4 format ('*** MPHURWITZZETAN: AA must be in the range (0, 1)')
  call mpabrt ( 562)
endif

! ss = mpreal (dble (is), mpnw)
! ss1 = 0.5d0 * ss
! ss2 = 0.5d0 * (ss + 1.d0)
! ss3 = 0.5d0 * (1.d0 - ss)
! ss4 = 1.d0 - 0.5d0 * ss

!+ call mpdmc (dble(is), 0, ss, mpnw1
call mpfrsetprec (ss(1:), mpnw1)
call mpfrsetd (ss(1:), dble(is), mprnd)
 
!+ call mpmuld (ss, 0.5d0, ss1, mpnw1
call mpfrsetprec (ss1(1:), mpnw1)
call mpfrmuld (ss1(1:), ss(1:), 0.5d0, mprnd)
 
!+ call mpdmc (1.d0, 0, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrsetd (t1(1:), 1.d0, mprnd)
 
!+ call mpadd (t1, ss, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfradd (t2(1:), t1(1:), ss(1:), mprnd)
 
!+ call mpmuld (t2, 0.5d0, ss2, mpnw1
call mpfrsetprec (ss2(1:), mpnw1)
call mpfrmuld (ss2(1:), t2(1:), 0.5d0, mprnd)
 
!+ call mpsub (t1, ss, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfrsub (t2(1:), t1(1:), ss(1:), mprnd)
 
!+ call mpmuld (t2, 0.5d0, ss3, mpnw1
call mpfrsetprec (ss3(1:), mpnw1)
call mpfrmuld (ss3(1:), t2(1:), 0.5d0, mprnd)
 
!+ call mpsub (t1, ss1, ss4, mpnw1
call mpfrsetprec (ss4(1:), mpnw1)
call mpfrsub (ss4(1:), t1(1:), ss1(1:), mprnd)
 

! gs1 = gamma (ss1)
! gs2 = gamma (ss2)
! t1 = pi * aa ** 2

!+ call mpgammar (ss1, gs1, mpnw1
call mpfrsetprec (gs1(1:), mpnw1)
call mpfrgamma (gs1(1:), ss1(1:), mprnd)
 
!+ call mpgammar (ss2, gs2, mpnw1
call mpfrsetprec (gs2(1:), mpnw1)
call mpfrgamma (gs2(1:), ss2(1:), mprnd)
 
!+ call mpmul (aa, aa, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfrmul (t2(1:), aa(1:), aa(1:), mprnd)
 
!+ call mpmul (mppicon, t2, t1, mpnw1
call mpfixlocr (mppicon)
call mpfrsetprec (t1(1:), mpnw1)
call mpfrmul (t1(1:), mppicon(1:), t2(1:), mprnd)
 

! sum1 = (incgamma (ss1, t1) / gs1 + incgamma (ss2, t1) / gs2) / abs (aa)**is
! sum2 = mpreal (0.d0, mpnw)
! sum3 = mpreal (0.d0, mpnw)

call mpincgammar (ss1, t1, t2, mpnw1)
!+ call mpdiv (t2, gs1, t3, mpnw1
call mpfrsetprec (t3(1:), mpnw1)
call mpfrdiv (t3(1:), t2(1:), gs1(1:), mprnd)
 
call mpincgammar (ss2, t1, t2, mpnw1)
!+ call mpdiv (t2, gs2, t4, mpnw1
call mpfrsetprec (t4(1:), mpnw1)
call mpfrdiv (t4(1:), t2(1:), gs2(1:), mprnd)
 
!+ call mpadd (t3, t4, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfradd (t2(1:), t3(1:), t4(1:), mprnd)
 
!+ call mpabs (aa, t3, mpnw1
call mpfrsetprec (t3(1:), mpnw1)
call mpfrabs (t3(1:), aa(1:), mprnd)
 
!+ call mpnpwr (t3, is, t4, mpnw1
call mpfrsetprec (t4(1:), mpnw1)
call mpfrpowsi (t4(1:), t3(1:), is, mprnd)
 
!+ call mpdiv (t2, t4, sum1, mpnw1
call mpfrsetprec (sum1(1:), mpnw1)
call mpfrdiv (sum1(1:), t2(1:), t4(1:), mprnd)
 
!+ call mpdmc (0.d0, 0, sum2, mpnw1
call mpfrsetprec (sum2(1:), mpnw1)
call mpfrsetd (sum2(1:), 0.d0, mprnd)
 
!+ call mpdmc (0.d0, 0, sum3, mpnw1
call mpfrsetprec (sum3(1:), mpnw1)
call mpfrsetd (sum3(1:), 0.d0, mprnd)
 

do k = 1, itrmax
  dk = dble (k)

!  t1 = pi * (dk + aa)**2
!  t2 = pi * (-dk + aa)**2
!  t3 = dk**2 * pi
!  t4 = 2.d0 * pi * dk * aa

!+   call mpdmc (dk, 0, t5, mpnw1
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfrsetd (t5(1:), dk, mprnd)
 
!+   call mpadd (t5, aa, t6, mpnw1
  call mpfrsetprec (t6(1:), mpnw1)
  call mpfradd (t6(1:), t5(1:), aa(1:), mprnd)
 
!+   call mpmul (t6, t6, t5, mpnw1
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfrmul (t5(1:), t6(1:), t6(1:), mprnd)
 
!+   call mpmul (mppicon, t5, t1, mpnw1
call mpfixlocr (mppicon)
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmul (t1(1:), mppicon(1:), t5(1:), mprnd)
 
!+   call mpdmc (-dk, 0, t5, mpnw1
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfrsetd (t5(1:), -dk, mprnd)
 
!+   call mpadd (t5, aa, t6, mpnw1
  call mpfrsetprec (t6(1:), mpnw1)
  call mpfradd (t6(1:), t5(1:), aa(1:), mprnd)
 
!+   call mpmul (t6, t6, t7, mpnw1
  call mpfrsetprec (t7(1:), mpnw1)
  call mpfrmul (t7(1:), t6(1:), t6(1:), mprnd)
 
!+   call mpmul (mppicon, t7, t2, mpnw1
call mpfixlocr (mppicon)
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmul (t2(1:), mppicon(1:), t7(1:), mprnd)
 
!+   call mpmul (t5, t5, t6, mpnw1
  call mpfrsetprec (t6(1:), mpnw1)
  call mpfrmul (t6(1:), t5(1:), t5(1:), mprnd)
 
!+   call mpmul (mppicon, t6, t3, mpnw1
call mpfixlocr (mppicon)
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrmul (t3(1:), mppicon(1:), t6(1:), mprnd)
 
!+   call mpmuld (mppicon, 2.d0 * dk, t5, mpnw1
call mpfixlocr (mppicon)
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfrmuld (t5(1:), mppicon(1:), 2.d0 * dk, mprnd)
 
!+   call mpmul (t5, aa, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfrmul (t4(1:), t5(1:), aa(1:), mprnd)
 

!  s1 = (incgamma (ss1, t1) / gs1 + incgamma (ss2, t1) / gs2) / abs (dk + aa)**is

  call mpincgammar (ss1, t1, t5, mpnw1)
!+   call mpdiv (t5, gs1, t6, mpnw1
  call mpfrsetprec (t6(1:), mpnw1)
  call mpfrdiv (t6(1:), t5(1:), gs1(1:), mprnd)
 
  call mpincgammar (ss2, t1, t5, mpnw1)
!+   call mpdiv (t5, gs2, t7, mpnw1
  call mpfrsetprec (t7(1:), mpnw1)
  call mpfrdiv (t7(1:), t5(1:), gs2(1:), mprnd)
 
!+   call mpadd (t6, t7, t5, mpnw1
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfradd (t5(1:), t6(1:), t7(1:), mprnd)
 
!+   call mpdmc (dk, 0, t6, mpnw1
  call mpfrsetprec (t6(1:), mpnw1)
  call mpfrsetd (t6(1:), dk, mprnd)
 
!+   call mpadd (t6, aa, t7, mpnw1
  call mpfrsetprec (t7(1:), mpnw1)
  call mpfradd (t7(1:), t6(1:), aa(1:), mprnd)
 
!+   call mpabs (t7, t6, mpnw1
  call mpfrsetprec (t6(1:), mpnw1)
  call mpfrabs (t6(1:), t7(1:), mprnd)
 
!+   call mpnpwr (t6, is, t7, mpnw1
  call mpfrsetprec (t7(1:), mpnw1)
  call mpfrpowsi (t7(1:), t6(1:), is, mprnd)
 
!+   call mpdiv (t5, t7, s1, mpnw1
  call mpfrsetprec (s1(1:), mpnw1)
  call mpfrdiv (s1(1:), t5(1:), t7(1:), mprnd)
 

!  s2 = (incgamma (ss1, t2) / gs1 - incgamma (ss2, t2) / gs2) / abs (-dk + aa)**is

  call mpincgammar (ss1, t2, t5, mpnw1)
!+   call mpdiv (t5, gs1, t6, mpnw1
  call mpfrsetprec (t6(1:), mpnw1)
  call mpfrdiv (t6(1:), t5(1:), gs1(1:), mprnd)
 
  call mpincgammar (ss2, t2, t5, mpnw1)
!+   call mpdiv (t5, gs2, t7, mpnw1
  call mpfrsetprec (t7(1:), mpnw1)
  call mpfrdiv (t7(1:), t5(1:), gs2(1:), mprnd)
 
!+   call mpsub (t6, t7, t5, mpnw1
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfrsub (t5(1:), t6(1:), t7(1:), mprnd)
 
!+   call mpdmc (-dk, 0, t6, mpnw1
  call mpfrsetprec (t6(1:), mpnw1)
  call mpfrsetd (t6(1:), -dk, mprnd)
 
!+   call mpadd (t6, aa, t7, mpnw1
  call mpfrsetprec (t7(1:), mpnw1)
  call mpfradd (t7(1:), t6(1:), aa(1:), mprnd)
 
!+   call mpabs (t7, t6, mpnw1
  call mpfrsetprec (t6(1:), mpnw1)
  call mpfrabs (t6(1:), t7(1:), mprnd)
 
!+   call mpnpwr (t6, is, t7, mpnw1
  call mpfrsetprec (t7(1:), mpnw1)
  call mpfrpowsi (t7(1:), t6(1:), is, mprnd)
 
!+   call mpdiv (t5, t7, s2, mpnw1
  call mpfrsetprec (s2(1:), mpnw1)
  call mpfrdiv (s2(1:), t5(1:), t7(1:), mprnd)
 

!  sum1 = sum1 + s1
!  sum2 = sum2 + s2

!+   call mpadd (sum1, s1, t5, mpnw1
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfradd (t5(1:), sum1(1:), s1(1:), mprnd)
 
!+   call mpeq (t5, sum1, mpnw1
  call mpfrsetprec (sum1(1:), mpnw1)
  call mpfrset (sum1(1:), t5(1:), mprnd)
 
!+   call mpadd (sum2, s2, t5, mpnw1
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfradd (t5(1:), sum2(1:), s2(1:), mprnd)
 
!+   call mpeq (t5, sum2, mpnw1
  call mpfrsetprec (sum2(1:), mpnw1)
  call mpfrset (sum2(1:), t5(1:), mprnd)
 

!  s3 = (incgamma (ss3, t3) * cos (t4)/ gs1 + incgamma (ss4, t3) * sin (t4) / gs2) &
!    / mpreal (dk, mpnw)**(1-is)

  call mpincgammar (ss3, t3, t5, mpnw1)
!+   call mpcssnr (t4, t6, t7, mpnw1
  call mpfrsetprec (t7(1:), mpnw1)
  call mpfrsincos (t7(1:), t6(1:), t4(1:), mprnd)
 
!+   call mpmul (t5, t6, t8, mpnw1
  call mpfrsetprec (t8(1:), mpnw1)
  call mpfrmul (t8(1:), t5(1:), t6(1:), mprnd)
 
!+   call mpdiv (t8, gs1, t5, mpnw1
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfrdiv (t5(1:), t8(1:), gs1(1:), mprnd)
 
  call mpincgammar (ss4, t3, t6, mpnw1)
!+   call mpmul (t6, t7, t8, mpnw1
  call mpfrsetprec (t8(1:), mpnw1)
  call mpfrmul (t8(1:), t6(1:), t7(1:), mprnd)
 
!+   call mpdiv (t8, gs2, t6, mpnw1
  call mpfrsetprec (t6(1:), mpnw1)
  call mpfrdiv (t6(1:), t8(1:), gs2(1:), mprnd)
 
!+   call mpadd (t5, t6, t7, mpnw1
  call mpfrsetprec (t7(1:), mpnw1)
  call mpfradd (t7(1:), t5(1:), t6(1:), mprnd)
 
!+   call mpdmc (dk, 0, t5, mpnw1
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfrsetd (t5(1:), dk, mprnd)
 
  i1 = 1 - is
!+   call mpnpwr (t5, i1, t6, mpnw1
  call mpfrsetprec (t6(1:), mpnw1)
  call mpfrpowsi (t6(1:), t5(1:), i1, mprnd)
 
!+   call mpdiv (t7, t6, s3, mpnw1
  call mpfrsetprec (s3(1:), mpnw1)
  call mpfrdiv (s3(1:), t7(1:), t6(1:), mprnd)
 

!  sum3 = sum3 + s3

!+   call mpadd (sum3, s3, t5, mpnw1
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfradd (t5(1:), sum3(1:), s3(1:), mprnd)
 
!+   call mpeq (t5, sum3, mpnw1
  call mpfrsetprec (sum3(1:), mpnw1)
  call mpfrset (sum3(1:), t5(1:), mprnd)
 

!  if (abs (s1) < eps * abs (sum1) .and. abs (s2) < eps * abs (sum2) .and. &
!    abs (s3) < eps * abs (sum3)) goto 100

!+   call mpabs (s1, tc1, 4
  call mpfrsetprec (tc1(1:), 4)
  call mpfrabs (tc1(1:), s1(1:), mprnd)
 
!+   call mpmul (eps, sum1, tc3, 4
  call mpfrsetprec (tc3(1:), 4)
  call mpfrmul (tc3(1:), eps(1:), sum1(1:), mprnd)
 
!+   call mpabs (tc3, tc2, 4
  call mpfrsetprec (tc2(1:), 4)
  call mpfrabs (tc2(1:), tc3(1:), mprnd)
 
!+   call mpcpr (tc1, tc2, ic1, 4
  ic1 = mpfrcmp (tc1(1:), tc2(1:))
 
!+   call mpabs (s2, tc1, 4
  call mpfrsetprec (tc1(1:), 4)
  call mpfrabs (tc1(1:), s2(1:), mprnd)
 
!+   call mpmul (eps, sum2, tc3, 4
  call mpfrsetprec (tc3(1:), 4)
  call mpfrmul (tc3(1:), eps(1:), sum2(1:), mprnd)
 
!+   call mpabs (tc3, tc2, 4
  call mpfrsetprec (tc2(1:), 4)
  call mpfrabs (tc2(1:), tc3(1:), mprnd)
 
!+   call mpcpr (tc1, tc2, ic2, 4
  ic2 = mpfrcmp (tc1(1:), tc2(1:))
 
!+   call mpabs (s3, tc1, 4
  call mpfrsetprec (tc1(1:), 4)
  call mpfrabs (tc1(1:), s3(1:), mprnd)
 
!+   call mpmul (eps, sum3, tc3, 4
  call mpfrsetprec (tc3(1:), 4)
  call mpfrmul (tc3(1:), eps(1:), sum3(1:), mprnd)
 
!+   call mpabs (tc3, tc2, 4
  call mpfrsetprec (tc2(1:), 4)
  call mpfrabs (tc2(1:), tc3(1:), mprnd)
 
!+   call mpcpr (tc1, tc2, ic3, 4
  ic3 = mpfrcmp (tc1(1:), tc2(1:))
 
  if (ic1 <= 0 .and. ic2 <= 0 .and. ic3 <= 0) goto 100
enddo

write (mpldb, 5)
5 format ('*** MPHURWITZZETAN: Loop end error')
call mpabrt ( 563)

100 continue

if (mod (is, 2) == 0) then
!  t1 = pi ** (is / 2) / ((ss - 1.d0) * gamma (ss1))

  i1 = is / 2
!+   call mpnpwr (mppicon, i1, t2, mpnw1
call mpfixlocr (mppicon)
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrpowsi (t2(1:), mppicon(1:), i1, mprnd)
 
!+   call mpdmc (1.d0, 0, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrsetd (t3(1:), 1.d0, mprnd)
 
!+   call mpsub (ss, t3, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfrsub (t4(1:), ss(1:), t3(1:), mprnd)
 
!+   call mpgammar (ss1, t5, mpnw1
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfrgamma (t5(1:), ss1(1:), mprnd)
 
!+   call mpmul (t4, t5, t6, mpnw1
  call mpfrsetprec (t6(1:), mpnw1)
  call mpfrmul (t6(1:), t4(1:), t5(1:), mprnd)
 
!+   call mpdiv (t2, t6, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrdiv (t1(1:), t2(1:), t6(1:), mprnd)
 
else
!  t1 = sqrt (pi) * pi ** ((is - 1) / 2) / ((ss - 1.d0) * gamma (ss1))

  i1 = (is - 1) / 2
!+   call mpnpwr (mppicon, i1, t2, mpnw1
call mpfixlocr (mppicon)
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrpowsi (t2(1:), mppicon(1:), i1, mprnd)
 
!+   call mpsqrt (mppicon, t3, mpnw1
call mpfixlocr (mppicon)
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrsqrt (t3(1:), mppicon(1:), mprnd)
 
!+   call mpmul (t2, t3, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfrmul (t4(1:), t2(1:), t3(1:), mprnd)
 
!+   call mpdmc (1.d0, 0, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrsetd (t2(1:), 1.d0, mprnd)
 
!+   call mpsub (ss, t2, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrsub (t3(1:), ss(1:), t2(1:), mprnd)
 
!+   call mpgammar (ss1, t5, mpnw1
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfrgamma (t5(1:), ss1(1:), mprnd)
 
!+   call mpmul (t3, t5, t6, mpnw1
  call mpfrsetprec (t6(1:), mpnw1)
  call mpfrmul (t6(1:), t3(1:), t5(1:), mprnd)
 
!+   call mpdiv (t4, t6, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrdiv (t1(1:), t4(1:), t6(1:), mprnd)
 
endif

!t2 = pi ** is / sqrt (pi)

!+ call mpnpwr (mppicon, is, t3, mpnw1
call mpfixlocr (mppicon)
call mpfrsetprec (t3(1:), mpnw1)
call mpfrpowsi (t3(1:), mppicon(1:), is, mprnd)
 
!+ call mpsqrt (mppicon, t4, mpnw1
call mpfixlocr (mppicon)
call mpfrsetprec (t4(1:), mpnw1)
call mpfrsqrt (t4(1:), mppicon(1:), mprnd)
 
!+ call mpdiv (t3, t4, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfrdiv (t2(1:), t3(1:), t4(1:), mprnd)
 

! hurwitzzetan = t1 + 0.5d0 * sum1 + 0.5d0 * sum2 + t2 * sum3

!+ call mpmuld (sum1, 0.5d0, t3, mpnw1
call mpfrsetprec (t3(1:), mpnw1)
call mpfrmuld (t3(1:), sum1(1:), 0.5d0, mprnd)
 
!+ call mpmuld (sum2, 0.5d0, t4, mpnw1
call mpfrsetprec (t4(1:), mpnw1)
call mpfrmuld (t4(1:), sum2(1:), 0.5d0, mprnd)
 
!+ call mpmul (sum3, t2, t5, mpnw1
call mpfrsetprec (t5(1:), mpnw1)
call mpfrmul (t5(1:), sum3(1:), t2(1:), mprnd)
 
!+ call mpadd (t1, t3, t6, mpnw1
call mpfrsetprec (t6(1:), mpnw1)
call mpfradd (t6(1:), t1(1:), t3(1:), mprnd)
 
!+ call mpadd (t6, t4, t7, mpnw1
call mpfrsetprec (t7(1:), mpnw1)
call mpfradd (t7(1:), t6(1:), t4(1:), mprnd)
 
!+ call mpadd (t7, t5, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfradd (t1(1:), t7(1:), t5(1:), mprnd)
 
!+ call mproun (t1, mpnw
call mpfrprecround (t1(1:), mpnw, mprnd)
 
!+ call mpeq (t1, zz, mpnw
call mpfrsetprec (zz(1:), mpnw)
call mpfrset (zz(1:), t1(1:), mprnd)
 

return
end subroutine mphurwitzzetan

subroutine mphurwitzzetanbe (nb1, nb2, berne, iss, aa, zz, mpnw)

!  This evaluates the Hurwitz zeta function, using the combination of
!  the definition formula (for large iss), and an Euler-Maclaurin scheme
!  (see formula 25.2.9 of the DLMF). The array berne contains precomputed
!  even Bernoulli numbers (see MPBERNER above). Its dimensions must be as
!  shown below. NB2 must be greater than 1.4 x precision in decimal digits.

implicit none
integer, intent(in):: nb1, nb2, iss
integer (mpiknd), intent(in):: berne(0:nb1+5,nb2), aa(0:)
integer (mpiknd), intent(out):: zz(0:)
integer, parameter:: itrmax = 1000000
real (mprknd), parameter:: dber = 1.4d0, dcon = 0.6d0
integer i1, i2, ic1, iqq, k, n1, mpnw, mpnw1
real (mprknd) d1, dp
integer (mpiknd) aq(0:mpnw+6), aq2(0:mpnw+6), s1(0:mpnw+6), s2(0:mpnw+6), &
  s3(0:mpnw+6), s4(0:mpnw+6), t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), &
  t4(0:mpnw+6), t5(0:mpnw+6), t6(0:mpnw+6), eps(0:9), f1(0:9), tc1(0:9), &
  tc2(0:9), tc3(0:9)

integer, external:: mpfrsgn
integer, external:: mpfrcmp
real (mprknd), external:: mpfrgetd2exp
integer (mpiknd) ix8
 
!  End of declaration

if (mpnw < 4 .or. mpspacer (aa) < mpnw + 4 .or. mpspacer (zz) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPHURWITZZETANBE: Uninitialized or inadequately sized arrays')
  call mpabrt ( 564)
endif

!  End of initial input array check

!   Check if berne array has been initialized.

!+ call mpmdc (berne(0:nb1+5,1), d1, n1, mpnw
d1 = 2.d0 * mpfrgetd2exp (ix8, berne(1:nb1+5,1), mprnd)
if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
d1 = d1 * 2.d0 ** n1
if (mpwprecr (berne(0:nb1+5,1)) < mpnw .or. &
  abs (d1 - 1.d0 / 6.d0) > mprdfz .or. nb2 < int (dber * mpdpw * mpnw)) then
  write (mpldb, 3) int (dber * mpdpw * mpnw)
3 format ('*** MPHURWITZZETANBE: Array of even Bernoulli coefficients must be initialized'/ &
   'with at least',i8,' entries by calling MPBERNE or MPBERNER.')
  call mpabrt ( 565)
endif

mpnw1 = min (mpnw + 1, mpnwx)
call mpinitwds (aq, mpnw1)
call mpinitwds (aq2, mpnw1)
call mpinitwds (s1, mpnw1)
call mpinitwds (s2, mpnw1)
call mpinitwds (s3, mpnw1)
call mpinitwds (s4, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (t6, mpnw1)
call mpinitwds (eps, 4)
call mpinitwds (f1, 4)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
!+ call mpdmc (2.d0, 0, tc1, 4
call mpfrsetprec (tc1(1:), 4)
call mpfrsetd (tc1(1:), 2.d0, mprnd)
 
!+ call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4
call mpfrsetprec (eps(1:), 4)
call mpfrpowsi (eps(1:), tc1(1:), -mpnw1*mpnbt, mprnd)
 
!+ call mpdmc (1.d0, 0, f1, 4
call mpfrsetprec (f1(1:), 4)
call mpfrsetd (f1(1:), 1.d0, mprnd)
 
!+ call mpdmc (0.d0, 0, s1, mpnw1
call mpfrsetprec (s1(1:), mpnw1)
call mpfrsetd (s1(1:), 0.d0, mprnd)
 
!+ call mpdmc (0.d0, 0, s2, mpnw1
call mpfrsetprec (s2(1:), mpnw1)
call mpfrsetd (s2(1:), 0.d0, mprnd)
 
!+ call mpdmc (0.d0, 0, s3, mpnw1
call mpfrsetprec (s3(1:), mpnw1)
call mpfrsetd (s3(1:), 0.d0, mprnd)
 
!+ call mpdmc (0.d0, 0, s4, mpnw1
call mpfrsetprec (s4(1:), mpnw1)
call mpfrsetd (s4(1:), 0.d0, mprnd)
 

if (iss <= 0) then
  write (mpldb, 4)
4 format ('*** MPHURWITZZETANBE: ISS <= 0')
  call mpabrt ( 566)
endif

!+ if (mpsgn (aa) < 0) the
if (mpfrsgn (aa(1:)) < 0) then
 
  write (mpldb, 5)
5 format ('*** MPHURWITZZETANBE: AA < 0')
  call mpabrt ( 567)
endif

dp = anint (mpnw1 * mpdpw)

!   If iss > a certain value, then use definition formula.

if (iss > 2.303d0 * dp / log (2.515d0 * dp)) then
  do k = 0, itrmax
!    t1 = 1.d0 / (aa + dble (k))**iss
!    s1 = s1 + t1

!+     call mpdmc (dble (k), 0, t1, mpnw1
    call mpfrsetprec (t1(1:), mpnw1)
    call mpfrsetd (t1(1:), dble (k), mprnd)
 
!+     call mpadd (aa, t1, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfradd (t2(1:), aa(1:), t1(1:), mprnd)
 
!+     call mpnpwr (t2, iss, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrpowsi (t3(1:), t2(1:), iss, mprnd)
 
!+     call mpdiv (f1, t3, t1, mpnw1
    call mpfrsetprec (t1(1:), mpnw1)
    call mpfrdiv (t1(1:), f1(1:), t3(1:), mprnd)
 
!+     call mpadd (s1, t1, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfradd (t2(1:), s1(1:), t1(1:), mprnd)
 
!+     call mpeq (t2, s1, mpnw1
    call mpfrsetprec (s1(1:), mpnw1)
    call mpfrset (s1(1:), t2(1:), mprnd)
 

!    if (abs (t1 / s1) < eps) goto 110

!+     call mpabs (t1, tc1, 4
    call mpfrsetprec (tc1(1:), 4)
    call mpfrabs (tc1(1:), t1(1:), mprnd)
 
!+     call mpmul (eps, s1, tc3, 4
    call mpfrsetprec (tc3(1:), 4)
    call mpfrmul (tc3(1:), eps(1:), s1(1:), mprnd)
 
!+     call mpabs (tc3, tc2, 4
    call mpfrsetprec (tc2(1:), 4)
    call mpfrabs (tc2(1:), tc3(1:), mprnd)
 
!+     call mpcpr (tc1, tc2, ic1, 4
    ic1 = mpfrcmp (tc1(1:), tc2(1:))
 
    if (ic1 <= 0) goto 110
  enddo

  write (6, 6)
6 format ('*** MPHURWITZZETANBE: Loop end error 1')
  call mpabrt ( 568)
endif

!+ call mpmdc (aa, d1, n1, mpnw1
d1 = 2.d0 * mpfrgetd2exp (ix8, aa(1:), mprnd)
if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
d1 = d1 * 2.d0**n1
iqq = max (dcon * mpnw1 * mpdpw - d1, 0.d0)

do k = 0, iqq - 1
!  s1 = s1 + 1.d0 / (aa + dble (k))**iss

!+   call mpdmc (dble (k), 0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrsetd (t1(1:), dble (k), mprnd)
 
!+   call mpadd (aa, t1, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfradd (t2(1:), aa(1:), t1(1:), mprnd)
 
!+   call mpnpwr (t2, iss, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrpowsi (t3(1:), t2(1:), iss, mprnd)
 
!+   call mpdiv (f1, t3, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrdiv (t1(1:), f1(1:), t3(1:), mprnd)
 
!+   call mpadd (s1, t1, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfradd (t2(1:), s1(1:), t1(1:), mprnd)
 
!+   call mpeq (t2, s1, mpnw1
  call mpfrsetprec (s1(1:), mpnw1)
  call mpfrset (s1(1:), t2(1:), mprnd)
 
enddo

! aq = aa + dble (iqq)

!+ call mpdmc (dble (iqq), 0, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrsetd (t1(1:), dble (iqq), mprnd)
 
!+ call mpadd (aa, t1, aq, mpnw1
call mpfrsetprec (aq(1:), mpnw1)
call mpfradd (aq(1:), aa(1:), t1(1:), mprnd)
 

! s2 = 1.d0 / (dble (iss - 1) * aq**(iss -  1))

!+ call mpdmc (dble (iss - 1), 0, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrsetd (t1(1:), dble (iss - 1), mprnd)
 
!+ call mpnpwr (aq, iss - 1, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfrpowsi (t2(1:), aq(1:), iss - 1, mprnd)
 
!+ call mpmul (t1, t2, t3, mpnw1
call mpfrsetprec (t3(1:), mpnw1)
call mpfrmul (t3(1:), t1(1:), t2(1:), mprnd)
 
!+ call mpdiv (f1, t3, s2, mpnw1
call mpfrsetprec (s2(1:), mpnw1)
call mpfrdiv (s2(1:), f1(1:), t3(1:), mprnd)
 

! s3 = 1.d0 / (2.d0 * aq**iss)

!+ call mpnpwr (aq, iss, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrpowsi (t1(1:), aq(1:), iss, mprnd)
 
!+ call mpmuld (t1, 2.d0, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfrmuld (t2(1:), t1(1:), 2.d0, mprnd)
 
!+ call mpdiv (f1, t2, s3, mpnw1
call mpfrsetprec (s3(1:), mpnw1)
call mpfrdiv (s3(1:), f1(1:), t2(1:), mprnd)
 

! t1 = mpreal (dble (iss), nwds)
! t2 = mpreal (1.d0, nwds)
! t3 = aq**(iss - 1)
! aq2 = aq**2

!+ call mpdmc (dble (iss), 0, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrsetd (t1(1:), dble (iss), mprnd)
 
!+ call mpdmc (1.d0, 0, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfrsetd (t2(1:), 1.d0, mprnd)
 
!+ call mpnpwr (aq, iss - 1, t3, mpnw1
call mpfrsetprec (t3(1:), mpnw1)
call mpfrpowsi (t3(1:), aq(1:), iss - 1, mprnd)
 
!+ call mpmul (aq, aq, aq2, mpnw1
call mpfrsetprec (aq2(1:), mpnw1)
call mpfrmul (aq2(1:), aq(1:), aq(1:), mprnd)
 

do k = 1, nb2
!  if (k > 1) t1 = t1 * dble (iss + 2*k - 3) * dble (iss + 2*k - 2)

  if (k > 1) then
!+     call mpmuld (t1, dble (iss + 2*k - 3), t5, mpnw1
    call mpfrsetprec (t5(1:), mpnw1)
    call mpfrmuld (t5(1:), t1(1:), dble (iss + 2*k - 3), mprnd)
 
!+     call mpmuld (t5, dble (iss + 2*k - 2), t1, mpnw1
    call mpfrsetprec (t1(1:), mpnw1)
    call mpfrmuld (t1(1:), t5(1:), dble (iss + 2*k - 2), mprnd)
 
  endif

!  t2 = t2 * dble (2 * k - 1) * dble (2 * k)

!+   call mpmuld (t2, dble (2 * k - 1), t5, mpnw1
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfrmuld (t5(1:), t2(1:), dble (2 * k - 1), mprnd)
 
!+   call mpmuld (t5, dble (2 * k), t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmuld (t2(1:), t5(1:), dble (2 * k), mprnd)
 

!  t3 = t3 * aq2
!  t4 = rb(k) * t1 / (t2 * t3)
!  s4 = s4 + t4

!+   call mpmul (t3, aq2, t5, mpnw1
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfrmul (t5(1:), t3(1:), aq2(1:), mprnd)
 
!+   call mpeq (t5, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrset (t3(1:), t5(1:), mprnd)
 
!+   call mpmul (berne(0:nb1+5,k), t1, t5, mpnw1
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfrmul (t5(1:), berne(1:nb1+5,k), t1(1:), mprnd)
 
!+   call mpmul (t2, t3, t6, mpnw1
  call mpfrsetprec (t6(1:), mpnw1)
  call mpfrmul (t6(1:), t2(1:), t3(1:), mprnd)
 
!+   call mpdiv (t5, t6, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfrdiv (t4(1:), t5(1:), t6(1:), mprnd)
 
!+   call mpadd (s4, t4, t5, mpnw1
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfradd (t5(1:), s4(1:), t4(1:), mprnd)
 
!+   call mpeq (t5, s4, mpnw1
  call mpfrsetprec (s4(1:), mpnw1)
  call mpfrset (s4(1:), t5(1:), mprnd)
 

!  if (abs (t4) < eps) goto 110

!+   call mpabs (t4, tc1, 4
  call mpfrsetprec (tc1(1:), 4)
  call mpfrabs (tc1(1:), t4(1:), mprnd)
 
!+   call mpmul (eps, s4, tc3, 4
  call mpfrsetprec (tc3(1:), 4)
  call mpfrmul (tc3(1:), eps(1:), s4(1:), mprnd)
 
!+   call mpabs (tc3, tc2, 4
  call mpfrsetprec (tc2(1:), 4)
  call mpfrabs (tc2(1:), tc3(1:), mprnd)
 
!+   call mpcpr (tc1, tc2, ic1, 4
  ic1 = mpfrcmp (tc1(1:), tc2(1:))
 
  if (ic1 <= 0) goto 110
enddo

write (6, 7)
7 format ('*** MPHURWITZZETANBE: End loop error 2; call MPBERNE with larger NB.')
call mpabrt ( 569)

110 continue

! hurwitz_be = s1 + s2 + s3 + s4

!+ call mpadd (s1, s2, t5, mpnw1
call mpfrsetprec (t5(1:), mpnw1)
call mpfradd (t5(1:), s1(1:), s2(1:), mprnd)
 
!+ call mpadd (t5, s3, t6, mpnw1
call mpfrsetprec (t6(1:), mpnw1)
call mpfradd (t6(1:), t5(1:), s3(1:), mprnd)
 
!+ call mpadd (t6, s4, s1, mpnw1
call mpfrsetprec (s1(1:), mpnw1)
call mpfradd (s1(1:), t6(1:), s4(1:), mprnd)
 
!+ call mproun (s1, mpnw
call mpfrprecround (s1(1:), mpnw, mprnd)
 
!+ call mpeq (s1, zz, mpnw
call mpfrsetprec (zz(1:), mpnw)
call mpfrset (zz(1:), s1(1:), mprnd)
 

return
end subroutine mphurwitzzetanbe

subroutine mphypergeompfq (np, nq, nw, aa, bb, xx, yy, mpnw)

!  This returns the HypergeometricPFQ function, namely the sum of the infinite series

!  Sum_0^infinity poch(aa(1),n)*poch(aa(2),n)*...*poch(aa(np),n) /
!      poch(bb(1),n)*poch(bb(2),n)*...*poch(bb(nq),n) * xx^n / n!

!  This subroutine evaluates the HypergeometricPFQ function directly according to
!  this definitional formula. The arrays aa and bb must be dimensioned as shown below.
!  NP and NQ are limited to [1,10].

implicit none
integer, intent(in):: np, nq, nw, mpnw
integer (mpiknd), intent(in):: aa(0:nw+5,np), bb(0:nw+5,nq), xx(0:)
integer (mpiknd), intent(out):: yy(0:)
integer, parameter:: itrmax = 1000000, npq = 10
integer i1, i2, ic1, j, k, mpnw1
integer (mpiknd) sum(0:mpnw+6), td(0:mpnw+6), tn(0:mpnw+6), t1(0:mpnw+6), &
  t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)

integer, external:: mpfrsgn
integer, external:: mpfrcmp
real (mprknd), external:: mpfrgetd2exp
integer (mpiknd) ix8
 
!  End of declaration

i1 = 1000000
i2 = 1000000

do k = 1, np
  i1 = min (i1, mpspacer(aa(0:nw+5,k)))
enddo

do k = 1, nq
  i2 = min (i2, mpspacer(bb(0:nw+5,k)))
enddo

if (mpnw < 4 .or. min (i1, i2) < mpnw + 4 .or. mpspacer (xx) < mpnw + 4 .or. &
  mpspacer (yy) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPHYPERGEOMPFQ: Uninitialized or inadequately sized arrays')
  call mpabrt ( 570)
endif

!  End of initial input array check

mpnw1 = min (mpnw + 1, mpnwx)
call mpinitwds (sum, mpnw1)
call mpinitwds (td, mpnw1)
call mpinitwds (tn, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
!+ call mpdmc (2.d0, 0, tc1, 4
call mpfrsetprec (tc1(1:), 4)
call mpfrsetd (tc1(1:), 2.d0, mprnd)
 
!+ call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4
call mpfrsetprec (eps(1:), 4)
call mpfrpowsi (eps(1:), tc1(1:), -mpnw1*mpnbt, mprnd)
 

if (np < 1 .or. np > npq .or. nq < 1 .or. nq > npq) then
  write (mpldb, 2) npq
2 format ('*** MPHYPERGEOMPFQ: NP and NQ must be between 1 and',i4)
  call mpabrt ( 571)
endif

!+ call mpdmc (1.d0, 0, sum, mpnw1
call mpfrsetprec (sum(1:), mpnw1)
call mpfrsetd (sum(1:), 1.d0, mprnd)
 
!+ call mpdmc (1.d0, 0, td, mpnw1
call mpfrsetprec (td(1:), mpnw1)
call mpfrsetd (td(1:), 1.d0, mprnd)
 
!+ call mpdmc (1.d0, 0, tn, mpnw1
call mpfrsetprec (tn(1:), mpnw1)
call mpfrsetd (tn(1:), 1.d0, mprnd)
 

do k = 1, itrmax
!+   call mpdmc (dble (k - 1), 0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrsetd (t1(1:), dble (k - 1), mprnd)
 

  do j = 1, np
!+     call mpadd (t1, aa(0:nw+5,j), t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfradd (t2(1:), t1(1:), aa(1:nw+5,j), mprnd)
 
!+     call mpmul (tn, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrmul (t3(1:), tn(1:), t2(1:), mprnd)
 
!+     call mpeq (t3, tn, mpnw1
    call mpfrsetprec (tn(1:), mpnw1)
    call mpfrset (tn(1:), t3(1:), mprnd)
 
  enddo

  do j = 1, nq
!+     call mpadd (t1, bb(0:nw+5,j), t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfradd (t2(1:), t1(1:), bb(1:nw+5,j), mprnd)
 
!+     call mpmul (td, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrmul (t3(1:), td(1:), t2(1:), mprnd)
 
!+     call mpeq (t3, td, mpnw1
    call mpfrsetprec (td(1:), mpnw1)
    call mpfrset (td(1:), t3(1:), mprnd)
 
  enddo

!+   call mpmul (tn, xx, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmul (t2(1:), tn(1:), xx(1:), mprnd)
 
!+   call mpeq (t2, tn, mpnw1
  call mpfrsetprec (tn(1:), mpnw1)
  call mpfrset (tn(1:), t2(1:), mprnd)
 
!+   call mpmuld (td, dble (k), t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrmuld (t3(1:), td(1:), dble (k), mprnd)
 
!+   call mpeq (t3, td, mpnw1
  call mpfrsetprec (td(1:), mpnw1)
  call mpfrset (td(1:), t3(1:), mprnd)
 
!+   call mpdiv (tn, td, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrdiv (t1(1:), tn(1:), td(1:), mprnd)
 
!+   call mpadd (sum, t1, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfradd (t2(1:), sum(1:), t1(1:), mprnd)
 
!+   call mpeq (t2, sum, mpnw1
  call mpfrsetprec (sum(1:), mpnw1)
  call mpfrset (sum(1:), t2(1:), mprnd)
 

!+   call mpabs (t1, tc1, 4
  call mpfrsetprec (tc1(1:), 4)
  call mpfrabs (tc1(1:), t1(1:), mprnd)
 
!+   call mpmul (eps, sum, tc3, 4
  call mpfrsetprec (tc3(1:), 4)
  call mpfrmul (tc3(1:), eps(1:), sum(1:), mprnd)
 
!+   call mpabs (tc3, tc2, 4
  call mpfrsetprec (tc2(1:), 4)
  call mpfrabs (tc2(1:), tc3(1:), mprnd)
 
!+   call mpcpr (tc1, tc2, ic1, 4
  ic1 = mpfrcmp (tc1(1:), tc2(1:))
 
  if (ic1 <= 0) goto 100
enddo

    write (mpldb, 3) itrmax
3   format ('*** MPHYPERGEOMPFQ: Loop end error',i10)
    call mpabrt ( 572)

100  continue

!+ call mproun (sum, mpnw
call mpfrprecround (sum(1:), mpnw, mprnd)
 
!+ call mpeq (sum, yy, mpnw
call mpfrsetprec (yy(1:), mpnw)
call mpfrset (yy(1:), sum(1:), mprnd)
 
return
end subroutine mphypergeompfq

subroutine mpincgammar (s, z, g, mpnw)

!  This returns the incomplete gamma function, using a combination of formula
!  8.7.3 of the DLMF (for modest-sized z), formula 8.11.2 (for large z),
!  a formula from the Wikipedia page for the case S = 0, and another formula
!  from the Wikipedia page for the case S = negative integer. The formula
!  for the case S = 0 requires increased working precision, up to 2.5X normal,
!  depending on the size of Z.

implicit none
integer (mpiknd), intent(in):: s(0:), z(0:)
integer (mpiknd), intent(out):: g(0:)
integer, intent(in):: mpnw
integer, parameter:: itrmax = 1000000
real (mprknd), parameter:: dmax = 0.833d0, egam = 0.5772156649015328606d0
integer ic1, k, mpnw1, mpnw2, nn, n1, n2
real (mprknd) d1, d2, bits
integer (mpiknd) t0(0:5*mpnw/2+6), t1(0:5*mpnw/2+6), t2(0:5*mpnw/2+6), &
  t3(0:5*mpnw/2+6), t4(0:5*mpnw/2+6), t5(0:5*mpnw/2+6), f1(0:5*mpnw/2+6), &
  tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)

integer, external:: mpfrsgn
integer, external:: mpfrcmp
real (mprknd), external:: mpfrgetd2exp
integer (mpiknd) ix8
 
!  End of declaration

if (mpnw < 4 .or. mpspacer (s) < mpnw + 4 .or. mpspacer (z) < mpnw + 4 &
  .or. mpspacer (g) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPINCGAMMAR: Uninitialized or inadequately sized arrays')
  call mpabrt ( 573)
endif

!   Check if EGAMMA has been precomputed.

!+ call mpmdc (mpegammacon, d1, n1, mpnw+1
call mpfixlocr (mpegammacon)
d1 = 2.d0 * mpfrgetd2exp (ix8, mpegammacon(1:), mprnd)
if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
if (n1 /= -1 .or. abs (d1 * 2.d0**n1 - egam) > mprdfz &
  .or. mpwprecr (mpegammacon) < mpnw+1) then
  write (mpldb, 3) mpnw+1
3 format ('*** MPINCGAMMAR: EGAMMA must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 575)
endif

!  End of initial input array check

!+ n1 = mpsgn (s
n1 = mpfrsgn (s(1:))
 
!+ n2 = mpsgn (z
n2 = mpfrsgn (z(1:))
 
if (n2 == 0 .or. n1 /= 0 .and. n2 < 0) then
  write (mpldb, 2)
2 format ('*** MPINCGAMMAR: The second argument must not be zero,'/ &
    'and must not be negative unless the first is zero.')
  call mpabrt ( 574)
endif

mpnw1 = min (mpnw + 1, mpnwx)
call mpinitwds (t0, 5*mpnw/2+1)
call mpinitwds (t1, 5*mpnw/2+1)
call mpinitwds (t2, 5*mpnw/2+1)
call mpinitwds (t3, 5*mpnw/2+1)
call mpinitwds (t4, 5*mpnw/2+1)
call mpinitwds (t5, 5*mpnw/2+1)
call mpinitwds (f1, 5*mpnw/2+1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
!+ call mpdmc (2.d0, 0, tc1, 4
call mpfrsetprec (tc1(1:), 4)
call mpfrsetd (tc1(1:), 2.d0, mprnd)
 
!+ call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4
call mpfrsetprec (eps(1:), 4)
call mpfrpowsi (eps(1:), tc1(1:), -mpnw1*mpnbt, mprnd)
 

!+ call mpdmc (1.d0, 0, f1, mpnw1
call mpfrsetprec (f1(1:), mpnw1)
call mpfrsetd (f1(1:), 1.d0, mprnd)
 
!+ call mpmdc (z, d1, n1, mpnw1
d1 = 2.d0 * mpfrgetd2exp (ix8, z(1:), mprnd)
if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
d1 = d1 * 2.d0 ** n1
bits = mpnw1 * mpnbt

if (abs (d1) < dmax * bits) then

!   This is for modest-sized z.

!+   call mpinfr (s, t1, t2, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrtrunc (t1(1:), s(1:))
  call mpfrsub (t2(1:), s(1:), t1(1:), mprnd)
 
!+   call mpcpr (s, t1, ic1, mpnw1
  ic1 = mpfrcmp (s(1:), t1(1:))
 
!+   call mpmdc (s, d2, n2, mpnw1
  d2 = 2.d0 * mpfrgetd2exp (ix8, s(1:), mprnd)
  if (d2 == 0) then; n2 = 0; else; n2 = ix8 - 1; endif
 
  nn = d2 * 2.d0**n2

  if (ic1 == 0 .and. nn == 1) then

!   S = 1; result is exp (-z).

!+     call mpneg (z, t0, mpnw1
    call mpfrsetprec (t0(1:), mpnw1)
    call mpfrneg (t0(1:), z(1:), mprnd)
 
!+     call mpexp (t0, t1, mpnw1
    call mpfrsetprec (t1(1:), mpnw1)
    call mpfrexp (t1(1:), t0(1:), mprnd)
 
    goto 200
  elseif (ic1 == 0 .and. nn <= 0) then

!    S is zero or a negative integer -- use a different algorithm. In
!    either event, first compute incgamma for S = 0. For large Z, the
!    working precision must be increased, up to 2.5X times normal.

!    mpnw2 = min (mpnw1 + 1.5d0 * d1 / (dmax * bits) * mpnw, 5*mpnw/2+1.d0)
    mpnw2 = mpnw1
!+     call mpeq (z, t0, mpnw2
    call mpfrsetprec (t0(1:), mpnw2)
    call mpfrset (t0(1:), z(1:), mprnd)
 
!+     call mpeq (z, t1, mpnw2
    call mpfrsetprec (t1(1:), mpnw2)
    call mpfrset (t1(1:), z(1:), mprnd)
 
!+     call mpdmc (1.d0, 0, t2, mpnw2
    call mpfrsetprec (t2(1:), mpnw2)
    call mpfrsetd (t2(1:), 1.d0, mprnd)
 

    do k = 2, itrmax
      if (mod (k, 2) == 1) then
        d1 = dble (k)
!+         call mpdivd (f1, d1, t3, mpnw2
        call mpfrsetprec (t3(1:), mpnw2)
        call mpfrdivd (t3(1:), f1(1:), d1, mprnd)
 
!+         call mpadd (t2, t3, t4, mpnw2
        call mpfrsetprec (t4(1:), mpnw2)
        call mpfradd (t4(1:), t2(1:), t3(1:), mprnd)
 
!+         call mpeq (t4, t2, mpnw2
        call mpfrsetprec (t2(1:), mpnw2)
        call mpfrset (t2(1:), t4(1:), mprnd)
 
      endif
!+       call mpmul (z, t1, t3, mpnw2
      call mpfrsetprec (t3(1:), mpnw2)
      call mpfrmul (t3(1:), z(1:), t1(1:), mprnd)
 
      d1 = 2.d0 * dble (k)
!+       call mpdivd (t3, d1, t1, mpnw2
      call mpfrsetprec (t1(1:), mpnw2)
      call mpfrdivd (t1(1:), t3(1:), d1, mprnd)
 
!+       call mpmul (t1, t2, t3, mpnw2
      call mpfrsetprec (t3(1:), mpnw2)
      call mpfrmul (t3(1:), t1(1:), t2(1:), mprnd)
 
!+       call mpadd (t0, t3, t4, mpnw2
      call mpfrsetprec (t4(1:), mpnw2)
      call mpfradd (t4(1:), t0(1:), t3(1:), mprnd)
 
!+       call mpeq (t4, t0, mpnw2
      call mpfrsetprec (t0(1:), mpnw2)
      call mpfrset (t0(1:), t4(1:), mprnd)
 

!+       call mpabs (t3, tc1, 4
      call mpfrsetprec (tc1(1:), 4)
      call mpfrabs (tc1(1:), t3(1:), mprnd)
 
!+       call mpmul (eps, t0, tc3, 4
      call mpfrsetprec (tc3(1:), 4)
      call mpfrmul (tc3(1:), eps(1:), t0(1:), mprnd)
 
!+       call mpabs (tc3, tc2, 4
      call mpfrsetprec (tc2(1:), 4)
      call mpfrabs (tc2(1:), tc3(1:), mprnd)
 
!+       call mpcpr (tc1, tc2, ic1, 4
      ic1 = mpfrcmp (tc1(1:), tc2(1:))
 
      if (ic1 <= 0) goto 100
    enddo

    write (mpldb, 4)
4   format ('*** MPINCGAMMAR: Loop end error 1')
    call mpabrt ( 576)

100  continue

!+     call mpneg (mpegammacon, t1, mpnw1
call mpfixlocr (mpegammacon)
    call mpfrsetprec (t1(1:), mpnw1)
    call mpfrneg (t1(1:), mpegammacon(1:), mprnd)
 
!+     call mpabs (z, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrabs (t3(1:), z(1:), mprnd)
 
!+     call mplog (t3, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrlog (t2(1:), t3(1:), mprnd)
 
!+     call mpsub (t1, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrsub (t3(1:), t1(1:), t2(1:), mprnd)
 
!+     call mpmuld (z, -0.5d0, t4, mpnw1
    call mpfrsetprec (t4(1:), mpnw1)
    call mpfrmuld (t4(1:), z(1:), -0.5d0, mprnd)
 
!+     call mpexp (t4, t5, mpnw1
    call mpfrsetprec (t5(1:), mpnw1)
    call mpfrexp (t5(1:), t4(1:), mprnd)
 
!+     call mpmul (t5, t0, t4, mpnw1
    call mpfrsetprec (t4(1:), mpnw1)
    call mpfrmul (t4(1:), t5(1:), t0(1:), mprnd)
 
!+     call mpadd (t3, t4, t1, mpnw1
    call mpfrsetprec (t1(1:), mpnw1)
    call mpfradd (t1(1:), t3(1:), t4(1:), mprnd)
 
    if (nn == 0) goto 200

!   S is negative integer (not zero).

    nn = abs (nn)
!+     call mpdmc (1.d0, 0, t0, mpnw1
    call mpfrsetprec (t0(1:), mpnw1)
    call mpfrsetd (t0(1:), 1.d0, mprnd)
 
!+     call mpeq (t0, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrset (t2(1:), t0(1:), mprnd)
 

    do k = 1, nn - 1
!+       call mpmuld (t0, dble (k), t2, mpnw1
      call mpfrsetprec (t2(1:), mpnw1)
      call mpfrmuld (t2(1:), t0(1:), dble (k), mprnd)
 
!+       call mpeq (t2, t0, mpnw1
      call mpfrsetprec (t0(1:), mpnw1)
      call mpfrset (t0(1:), t2(1:), mprnd)
 
    enddo

!+     call mpmuld (t0, dble (nn), t5, mpnw1
    call mpfrsetprec (t5(1:), mpnw1)
    call mpfrmuld (t5(1:), t0(1:), dble (nn), mprnd)
 

    do k = 1, nn - 1
!+       call mpmul (t2, z, t3, mpnw1
      call mpfrsetprec (t3(1:), mpnw1)
      call mpfrmul (t3(1:), t2(1:), z(1:), mprnd)
 
!+       call mpdivd (t3, dble (nn - k), t4, mpnw1
      call mpfrsetprec (t4(1:), mpnw1)
      call mpfrdivd (t4(1:), t3(1:), dble (nn - k), mprnd)
 
!+       call mpneg (t4, t2, mpnw1
      call mpfrsetprec (t2(1:), mpnw1)
      call mpfrneg (t2(1:), t4(1:), mprnd)
 
!+       call mpadd (t0, t2, t3, mpnw1
      call mpfrsetprec (t3(1:), mpnw1)
      call mpfradd (t3(1:), t0(1:), t2(1:), mprnd)
 
!+       call mpeq (t3, t0, mpnw1
      call mpfrsetprec (t0(1:), mpnw1)
      call mpfrset (t0(1:), t3(1:), mprnd)
 
    enddo

!+     call mpexp (z, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrexp (t2(1:), z(1:), mprnd)
 
!+     call mpdiv (t0, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrdiv (t3(1:), t0(1:), t2(1:), mprnd)
 
!+     call mpnpwr (z, nn, t4, mpnw1
    call mpfrsetprec (t4(1:), mpnw1)
    call mpfrpowsi (t4(1:), z(1:), nn, mprnd)
 
!+     call mpdiv (t3, t4, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrdiv (t2(1:), t3(1:), t4(1:), mprnd)
 

    if (mod (nn, 2) == 0) then
!+       call mpadd (t2, t1, t3, mpnw1
      call mpfrsetprec (t3(1:), mpnw1)
      call mpfradd (t3(1:), t2(1:), t1(1:), mprnd)
 
    else
!+       call mpsub (t2, t1, t3, mpnw1
      call mpfrsetprec (t3(1:), mpnw1)
      call mpfrsub (t3(1:), t2(1:), t1(1:), mprnd)
 
    endif
!+     call mpdiv (t3, t5, t1, mpnw1
    call mpfrsetprec (t1(1:), mpnw1)
    call mpfrdiv (t1(1:), t3(1:), t5(1:), mprnd)
 
    goto 200
  endif

!+   call mpgammar (s, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrgamma (t1(1:), s(1:), mprnd)
 
!+   call mpmul (s, t1, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrmul (t3(1:), s(1:), t1(1:), mprnd)
 
!+   call mpdiv (f1, t3, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrdiv (t2(1:), f1(1:), t3(1:), mprnd)
 
!+   call mpeq (t2, t0, mpnw1
  call mpfrsetprec (t0(1:), mpnw1)
  call mpfrset (t0(1:), t2(1:), mprnd)
 

  do k = 1, itrmax
!+     call mpmul (t2, z, t5, mpnw1
    call mpfrsetprec (t5(1:), mpnw1)
    call mpfrmul (t5(1:), t2(1:), z(1:), mprnd)
 
!+     call mpdmc (dble (k), 0, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrsetd (t3(1:), dble (k), mprnd)
 
!+     call mpadd (s, t3, t4, mpnw1
    call mpfrsetprec (t4(1:), mpnw1)
    call mpfradd (t4(1:), s(1:), t3(1:), mprnd)
 
!+     call mpdiv (t5, t4, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrdiv (t2(1:), t5(1:), t4(1:), mprnd)
 
!+     call mpadd (t0, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfradd (t3(1:), t0(1:), t2(1:), mprnd)
 
!+     call mpeq (t3, t0, mpnw1
    call mpfrsetprec (t0(1:), mpnw1)
    call mpfrset (t0(1:), t3(1:), mprnd)
 

!+     call mpabs (t2, tc1, 4
    call mpfrsetprec (tc1(1:), 4)
    call mpfrabs (tc1(1:), t2(1:), mprnd)
 
!+     call mpmul (eps, t0, tc3, 4
    call mpfrsetprec (tc3(1:), 4)
    call mpfrmul (tc3(1:), eps(1:), t0(1:), mprnd)
 
!+     call mpabs (tc3, tc2, 4
    call mpfrsetprec (tc2(1:), 4)
    call mpfrabs (tc2(1:), tc3(1:), mprnd)
 
!+     call mpcpr (tc1, tc2, ic1, 4
    ic1 = mpfrcmp (tc1(1:), tc2(1:))
 
    if (ic1 <= 0) goto 110
  enddo

  write (mpldb, 5) itrmax
5   format ('*** MPINCGAMMAR: Loop end error 1')
  call mpabrt ( 577)

110 continue

!+   call mppower (z, s, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrpow (t2(1:), z(1:), s(1:), mprnd)
 
!+   call mpexp (z, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrexp (t3(1:), z(1:), mprnd)
 
!+   call mpdiv (t2, t3, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfrdiv (t4(1:), t2(1:), t3(1:), mprnd)
 
!+   call mpmul (t4, t0, t5, mpnw1
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfrmul (t5(1:), t4(1:), t0(1:), mprnd)
 
!+   call mpsub (f1, t5, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrsub (t2(1:), f1(1:), t5(1:), mprnd)
 
!+   call mpmul (t1, t2, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrmul (t3(1:), t1(1:), t2(1:), mprnd)
 
!+   call mpeq (t3, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrset (t1(1:), t3(1:), mprnd)
 
  goto 200
else

!   This is for large z. Note that if S is a positive integer, this loop
!   is finite.

!+   call mpdmc (1.d0, 0, t0, mpnw1
  call mpfrsetprec (t0(1:), mpnw1)
  call mpfrsetd (t0(1:), 1.d0, mprnd)
 
!+   call mpdmc (1.d0, 0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrsetd (t1(1:), 1.d0, mprnd)
 

  do k = 1, itrmax
!+     call mpdmc (dble (k), 0, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrsetd (t2(1:), dble (k), mprnd)
 
!+     call mpsub (s, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrsub (t3(1:), s(1:), t2(1:), mprnd)
 
!+     call mpmul (t1, t3, t4, mpnw1
    call mpfrsetprec (t4(1:), mpnw1)
    call mpfrmul (t4(1:), t1(1:), t3(1:), mprnd)
 
!+     call mpdiv (t4, z, t1, mpnw1
    call mpfrsetprec (t1(1:), mpnw1)
    call mpfrdiv (t1(1:), t4(1:), z(1:), mprnd)
 
!+     call mpadd (t0, t1, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfradd (t2(1:), t0(1:), t1(1:), mprnd)
 
!+     call mpeq (t2, t0, mpnw1
    call mpfrsetprec (t0(1:), mpnw1)
    call mpfrset (t0(1:), t2(1:), mprnd)
 

!+     call mpabs (t1, tc1, 4
    call mpfrsetprec (tc1(1:), 4)
    call mpfrabs (tc1(1:), t1(1:), mprnd)
 
!+     call mpmul (eps, t0, tc3, 4
    call mpfrsetprec (tc3(1:), 4)
    call mpfrmul (tc3(1:), eps(1:), t0(1:), mprnd)
 
!+     call mpabs (tc3, tc2, 4
    call mpfrsetprec (tc2(1:), 4)
    call mpfrabs (tc2(1:), tc3(1:), mprnd)
 
!+     call mpcpr (tc1, tc2, ic1, 4
    ic1 = mpfrcmp (tc1(1:), tc2(1:))
 
    if (ic1 <= 0) goto 120
  enddo

  write (mpldb, 6)
6 format ('*** MPINCGAMMAR: Loop end error 3')
  call mpabrt ( 578)

120 continue

!+    call mpsub (s, f1, t2, mpnw1
   call mpfrsetprec (t2(1:), mpnw1)
   call mpfrsub (t2(1:), s(1:), f1(1:), mprnd)
 
!+    call mppower (z, t2, t3, mpnw1
   call mpfrsetprec (t3(1:), mpnw1)
   call mpfrpow (t3(1:), z(1:), t2(1:), mprnd)
 
!+    call mpexp (z, t4, mpnw1
   call mpfrsetprec (t4(1:), mpnw1)
   call mpfrexp (t4(1:), z(1:), mprnd)
 
!+    call mpdiv (t3, t4, t2, mpnw1
   call mpfrsetprec (t2(1:), mpnw1)
   call mpfrdiv (t2(1:), t3(1:), t4(1:), mprnd)
 
!+    call mpmul (t2, t0, t1, mpnw1
   call mpfrsetprec (t1(1:), mpnw1)
   call mpfrmul (t1(1:), t2(1:), t0(1:), mprnd)
 
   goto 200
endif

200 continue

!+ call mproun (t1, mpnw
call mpfrprecround (t1(1:), mpnw, mprnd)
 
!+ call mpeq (t1, g, mpnw
call mpfrsetprec (g(1:), mpnw)
call mpfrset (g(1:), t1(1:), mprnd)
 

return
end subroutine mpincgammar

subroutine mppolygamma (nn, x, y, mpnw)

!   This returns polygamma (nn, x) for nn >= 0 and 0 < x < 1, by calling
!   mphurwitzzetan.

implicit none
integer, intent(in):: nn, mpnw
integer (mpiknd), intent(in):: x(0:)
integer (mpiknd), intent(out):: y(0:)
integer ic1, ic2, k, mpnw1
integer (mpiknd) t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6)

integer, external:: mpfrsgn
integer, external:: mpfrcmp
real (mprknd), external:: mpfrgetd2exp
integer (mpiknd) ix8
 
!  End of declaration

if (mpnw < 4 .or. mpspacer (x) < mpnw + 4 .or. mpspacer (y) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPPOLYGAMMA: Uninitialized or inadequately sized arrays')
  call mpabrt ( 579)
endif

!  End of initial input array check

mpnw1 = min (mpnw, mpnwx)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)

if (nn <= 0) then
  write (mpldb, 2)
2 format ('*** MPPOLYGAMMA: NN <= 0')
  call mpabrt ( 580)
endif

!+ call mpdmc (0.d0, 0, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrsetd (t1(1:), 0.d0, mprnd)
 
!+ call mpdmc (1.d0, 0, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfrsetd (t2(1:), 1.d0, mprnd)
 
!+ call mpcpr (x, t1, ic1, mpnw1
ic1 = mpfrcmp (x(1:), t1(1:))
 
!+ call mpcpr (x, t2, ic2, mpnw1
ic2 = mpfrcmp (x(1:), t2(1:))
 
if (ic1 <= 0 .or. ic2 >= 0) then
  write (mpldb, 3)
3 format ('*** MPPOLYGAMMA: X must be in the range (0, 1)')
  call mpabrt ( 581)
endif

!+ call mpdmc (1.d0, 0, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrsetd (t1(1:), 1.d0, mprnd)
 

do k = 1, nn
!+   call mpmuld (t1, dble(k), t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmuld (t2(1:), t1(1:), dble(k), mprnd)
 
!+   call mpeq (t2, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrset (t1(1:), t2(1:), mprnd)
 
enddo

if (mod (nn + 1, 2) == 1) then
!+   call mpneg (t1, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrneg (t2(1:), t1(1:), mprnd)
 
!+   call mpeq (t2, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrset (t1(1:), t2(1:), mprnd)
 
endif
call mphurwitzzetan (nn + 1, x, t2, mpnw1)
!+ call mpmul (t1, t2, t3, mpnw1
call mpfrsetprec (t3(1:), mpnw1)
call mpfrmul (t3(1:), t1(1:), t2(1:), mprnd)
 
!+ call mproun (t3, mpnw
call mpfrprecround (t3(1:), mpnw, mprnd)
 
!+ call mpeq (t3, y, mpnw
call mpfrsetprec (y(1:), mpnw)
call mpfrset (y(1:), t3(1:), mprnd)
 

return
end subroutine mppolygamma

subroutine mppolygammabe (nb1, nb2, berne, nn, x, y, mpnw)

!  This returns polygamma (nn, x) for nn >= 0, by calling mphurwitzzetanbe.
!  The array berne contains precomputed even Bernoulli numbers (see MPBERNER
!  above). Its dimensions must be as shown below. NB2 must be greater than
!  1.4 x precision in decimal digits.

implicit none
integer, intent(in):: nb1, nb2, nn, mpnw
integer (mpiknd), intent(in):: berne(0:nb1+5,nb2), x(0:)
integer (mpiknd), intent(out):: y(0:)
real (mprknd), parameter:: dber = 1.4d0
integer i1, i2, k, mpnw1, n1
real (mprknd) d1
integer (mpiknd) t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6)

integer, external:: mpfrsgn
integer, external:: mpfrcmp
real (mprknd), external:: mpfrgetd2exp
integer (mpiknd) ix8
 
!  End of declaration

if (mpnw < 4 .or. mpspacer (x) < mpnw + 4 .or. mpspacer (y) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPPOLYGAMMABE: Uninitialized or inadequately sized arrays')
  call mpabrt ( 582)
endif

!  End of initial input array check

!   Check if berne array has sufficient entries.

!+ call mpmdc (berne(0:nb1+5,1), d1, n1, mpnw
d1 = 2.d0 * mpfrgetd2exp (ix8, berne(1:nb1+5,1), mprnd)
if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
d1 = d1 * 2.d0 ** n1
if (mpwprecr (berne(0:nb1+5,1)) < mpnw .or. &
  abs (d1 - 1.d0 / 6.d0) > mprdfz .or. nb2 < int (dber * mpdpw * mpnw)) then
  write (mpldb, 3) int (dber * mpdpw * mpnw)
3 format ('*** MPPOLYGAMMABE: Array of even Bernoulli coefficients must be initialized'/ &
   'with at least',i8,' entries by calling MPBERNE or MPBERER.')
  call mpabrt ( 583)
endif

mpnw1 = min (mpnw, mpnwx)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)

if (nn <= 0) then
  write (mpldb, 4)
4 format ('*** MPPOLYGAMMABE: NN <= 0')
  call mpabrt ( 584)
endif

!+ if (mpsgn (x) < 0) the
if (mpfrsgn (x(1:)) < 0) then
 
  write (mpldb, 5)
5 format ('*** MPPOLYGAMMABE: X < 0')
  call mpabrt ( 585)
endif

!+ call mpdmc (1.d0, 0, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrsetd (t1(1:), 1.d0, mprnd)
 

do k = 1, nn
!+   call mpmuld (t1, dble(k), t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmuld (t2(1:), t1(1:), dble(k), mprnd)
 
!+   call mpeq (t2, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrset (t1(1:), t2(1:), mprnd)
 
enddo

if (mod (nn + 1, 2) == 1) then
!+   call mpneg (t1, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrneg (t2(1:), t1(1:), mprnd)
 
!+   call mpeq (t2, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrset (t1(1:), t2(1:), mprnd)
 
endif
call mphurwitzzetanbe (nb1, nb2, berne, nn + 1, x, t2, mpnw1)
!+ call mpmul (t1, t2, t3, mpnw1
call mpfrsetprec (t3(1:), mpnw1)
call mpfrmul (t3(1:), t1(1:), t2(1:), mprnd)
 
!+ call mproun (t3, mpnw
call mpfrprecround (t3(1:), mpnw, mprnd)
 
!+ call mpeq (t3, y, mpnw
call mpfrsetprec (y(1:), mpnw)
call mpfrset (y(1:), t3(1:), mprnd)
 

return
end subroutine mppolygammabe

subroutine mppolylogini (na, nn, arr, mpnw)

!   Initializes the MP array arr with data for mppolylogneg.
!   NN must be in the range (-nmax, -1).

implicit none
integer, intent(in):: na, nn, mpnw
integer (mpiknd), intent(out):: arr(0:na+5,1:abs(nn))
integer, parameter:: nmax = 1000
integer i1, i2, k, n, nna, mpnw1
integer (mpiknd) aa(0:mpnw+6,2,abs(nn)), t1(0:mpnw+6), t2(0:mpnw+6)

integer, external:: mpfrsgn
integer, external:: mpfrcmp
real (mprknd), external:: mpfrgetd2exp
integer (mpiknd) ix8
 
!  End of declaration

if (nn == 0 .or. abs (nn) > nmax) then
  write (mpldb, 1) -nmax
1 format ('*** MPPOLYLOGINI: N = 0 or N > ',i6)
  call mpabrt ( 586)
endif

i1 = 1000000000

do k = 1, abs (nn)
  i1 = min (i1, mpspacer (arr(0:na+5,k)))
enddo

if (mpnw < 4 .or. i1 < mpnw + 6) then
  write (mpldb, 2)
2 format ('*** MPPOLYLOGINI: Uninitialized or inadequately sized arrays')
  call mpabrt ( 587)
endif

!  End of initial input array check

nna = abs (nn)
mpnw1 = min (mpnw + 1, mpnwx)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
i1 = 2
i2 = 1
call mpinitwds (aa(0:mpnw+6,1,1), mpnw1)
call mpinitwds (aa(0:mpnw+6,2,1), mpnw1)
!+ call mpdmc (1.d0, 0, aa(0:mpnw+6,1,1), mpnw1
call mpfrsetprec (aa(1:mpnw+6,1,1), mpnw1)
call mpfrsetd (aa(1:mpnw+6,1,1), 1.d0, mprnd)
 
!+ call mpdmc (1.d0, 0, aa(0:mpnw+6,2,1), mpnw1
call mpfrsetprec (aa(1:mpnw+6,2,1), mpnw1)
call mpfrsetd (aa(1:mpnw+6,2,1), 1.d0, mprnd)
 

do k = 2, nna
  call mpinitwds (aa(0:mpnw+6,1,k), mpnw1)
  call mpinitwds (aa(0:mpnw+6,2,k), mpnw1)
!+   call mpdmc (0.d0, 0, aa(0:mpnw+6,1,k), mpnw1
  call mpfrsetprec (aa(1:mpnw+6,1,k), mpnw1)
  call mpfrsetd (aa(1:mpnw+6,1,k), 0.d0, mprnd)
 
!+   call mpdmc (0.d0, 0, aa(0:mpnw+6,2,k), mpnw1
  call mpfrsetprec (aa(1:mpnw+6,2,k), mpnw1)
  call mpfrsetd (aa(1:mpnw+6,2,k), 0.d0, mprnd)
 
enddo

do n = 2, nna
  i1 = 3 - i1
  i2 = 3 - i1

  do k = 2, n
!+     call mpmuld (aa(0:mpnw+6,i1,k-1), dble (n + 1 - k), t1, mpnw1
    call mpfrsetprec (t1(1:), mpnw1)
    call mpfrmuld (t1(1:), aa(1:mpnw+6,i1,k-1), dble (n + 1 - k), mprnd)
 
!+     call mpmuld (aa(0:mpnw+6,i1,k), dble (k), t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrmuld (t2(1:), aa(1:mpnw+6,i1,k), dble (k), mprnd)
 
!+     call mpadd (t1, t2, aa(0:mpnw+6,i2,k), mpnw1
    call mpfrsetprec (aa(1:mpnw+6,i2,k), mpnw1)
    call mpfradd (aa(1:mpnw+6,i2,k), t1(1:), t2(1:), mprnd)
 
  enddo
enddo

do k = 1, nna
!+   call mpeq (aa(0:mpnw+6,i2,k), arr(0:na+5,k), mpnw
  call mpfrsetprec (arr(1:na+5,k), mpnw)
  call mpfrset (arr(1:na+5,k), aa(1:mpnw+6,i2,k), mprnd)
 
enddo

return
end subroutine mppolylogini

subroutine mppolylogneg (na, nn, arr, x, y, mpnw)

!   This returns polylog (nn, x) for the case nn < 0. Before calling this,
!   one must call mppolylognini to initialize the array arr for this NN.
!   The dimensions of arr must be as shown below.
!   NN must be in the range (-nmax, -1).
!   The parameter nmxa is the maximum number of additional words of
!   precision needed to overcome cancelation errors when x is negative,
!   for nmax = 1000.

implicit none
integer, intent(in):: na, nn, mpnw
integer (mpiknd), intent(in):: arr(0:na+5,1:abs(nn)), x(0:)
integer (mpiknd), intent(out):: y(0:)
integer, parameter:: nmax = 1000, nmxa = 8525 / mpnbt + 1
integer i1, i2, k, mpnw1, n1, n2, nna
real (mprknd) d1, d2
integer (mpiknd) t1(0:mpnw+6+nmxa), t2(0:mpnw+6+nmxa), t3(0:mpnw+6+nmxa), &
  t4(0:mpnw+6+nmxa)

integer, external:: mpfrsgn
integer, external:: mpfrcmp
real (mprknd), external:: mpfrgetd2exp
integer (mpiknd) ix8
 
!  End of declaration

if (nn < -nmax .or. nn >= 0) then
  write (mpldb, 1) -nmax
1 format ('*** MPPOLYLOGNEG: N is <',i6,' or n >= 0.'/ &
  'For n >= 0, call mppolylogpos or polylog_pos.')
  call mpabrt ( 588)
endif

!  End of initial input array check

nna = abs (nn)
!+ call mpmdc (arr(0:na+5,1), d1, n1, mpnw
d1 = 2.d0 * mpfrgetd2exp (ix8, arr(1:na+5,1), mprnd)
if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
d1 = d1 * 2.d0 ** n1
!+ call mpmdc (arr(0:na+5,nna), d2, n2, mpnw
d2 = 2.d0 * mpfrgetd2exp (ix8, arr(1:na+5,nna), mprnd)
if (d2 == 0) then; n2 = 0; else; n2 = ix8 - 1; endif
 
d2 = d2 * 2.d0 ** n2

if (d1 /= 1.d0 .or. d2 /= 1.d0) then
  write (mpldb, 2)
2 format ('*** MPPOLYLOGNEG: Uninitialized or inadequately sized arrays'/ &
  'Call mppolylogini or polylog_ini to initialize array. See documentation.')
  call mpabrt ( 589)
endif

mpnw1 = min (mpnw + 1, mpnwx)
call mpinitwds (t1, mpnw + nmxa + 1)
call mpinitwds (t2, mpnw + nmxa + 1)
call mpinitwds (t3, mpnw + nmxa + 1)
call mpinitwds (t4, mpnw + nmxa + 1)

!+ if (mpsgn (x) < 0) the
if (mpfrsgn (x(1:)) < 0) then
 
  i1 = (nna + 1) / 2
!+   call mpmdc (arr(0:mpnw+5,i1), d1, n1, mpnw1
  d1 = 2.d0 * mpfrgetd2exp (ix8, arr(1:mpnw+5,i1), mprnd)
  if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
  mpnw1 = min (mpnw1 + (n1 + 1) / mpnbt + 1, mpnwx)
endif

!+ call mpeq (x, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrset (t1(1:), x(1:), mprnd)
 
!+ call mpeq (t1, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfrset (t2(1:), t1(1:), mprnd)
 

do k = 2, nna
!+   call mpmul (x, t1, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrmul (t3(1:), x(1:), t1(1:), mprnd)
 
!+   call mpeq (t3, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrset (t1(1:), t3(1:), mprnd)
 
!+   call mpmul (arr(0:mpnw+5,k), t1, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrmul (t3(1:), arr(1:mpnw+5,k), t1(1:), mprnd)
 
!+   call mpadd (t2, t3, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfradd (t4(1:), t2(1:), t3(1:), mprnd)
 
!+   call mpeq (t4, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrset (t2(1:), t4(1:), mprnd)
 
enddo

!+ call mpdmc (1.d0, 0, t3, mpnw1
call mpfrsetprec (t3(1:), mpnw1)
call mpfrsetd (t3(1:), 1.d0, mprnd)
 
!+ call mpsub (t3, x, t4, mpnw1
call mpfrsetprec (t4(1:), mpnw1)
call mpfrsub (t4(1:), t3(1:), x(1:), mprnd)
 
!+ call mpnpwr (t4, nna + 1, t3, mpnw1
call mpfrsetprec (t3(1:), mpnw1)
call mpfrpowsi (t3(1:), t4(1:), nna + 1, mprnd)
 
!+ call mpdiv (t2, t3, t4, mpnw1
call mpfrsetprec (t4(1:), mpnw1)
call mpfrdiv (t4(1:), t2(1:), t3(1:), mprnd)
 
!+ call mproun (t4, mpnw
call mpfrprecround (t4(1:), mpnw, mprnd)
 
!+ call mpeq (t4, y, mpnw
call mpfrsetprec (y(1:), mpnw)
call mpfrset (y(1:), t4(1:), mprnd)
 

return
end subroutine mppolylogneg

subroutine mppolylogpos (nn, x, y, mpnw)

!   This returns polylog (nn, x) for the case nn >= 0.

implicit none
integer, intent(in):: nn, mpnw
integer (mpiknd), intent(in):: x(0:)
integer (mpiknd), intent(out):: y(0:)
integer, parameter:: itrmax = 1000000
integer ic1, k, mpnw1
integer (mpiknd) t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), &
  t5(0:mpnw+6), tc1(0:9), tc2(0:9), tc3(0:9), eps (0:9)

integer, external:: mpfrsgn
integer, external:: mpfrcmp
real (mprknd), external:: mpfrgetd2exp
integer (mpiknd) ix8
 
!  End of declaration

if (mpnw < 4 .or. mpspacer (x) < mpnw + 4 .or. mpspacer (y) < mpnw + 6) then
  write (mpldb, 2)
2 format ('*** MPPOLYLOGPOS: Uninitialized or inadequately sized arrays')
  call mpabrt ( 591)
endif

!  End of initial input array check

if (nn < 0) then
  write (mpldb, 1)
1 format ('*** MPPOLYLOGPOS: N is less than zero.'/ &
  'For negative n, call mppolylogneg or polylog_neg. See documentation.')
  call mpabrt ( 590)
endif

mpnw1 = min (mpnw + 1, mpnwx)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
!+ call mpdmc (2.d0, 0, tc1, 4
call mpfrsetprec (tc1(1:), 4)
call mpfrsetd (tc1(1:), 2.d0, mprnd)
 
!+ call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4
call mpfrsetprec (eps(1:), 4)
call mpfrpowsi (eps(1:), tc1(1:), -mpnw1*mpnbt, mprnd)
 

!+ call mpabs (x, t1, mpnw
call mpfrsetprec (t1(1:), mpnw)
call mpfrabs (t1(1:), x(1:), mprnd)
 
!+ call mpdmc (1.d0, 0, t2, mpnw
call mpfrsetprec (t2(1:), mpnw)
call mpfrsetd (t2(1:), 1.d0, mprnd)
 
!+ call mpcpr (t1, t2, ic1, mpnw
ic1 = mpfrcmp (t1(1:), t2(1:))
 
if (ic1 >= 0) then
  write (mpldb, 3)
3 format ('*** MPPOLYLOGPOS: |X| must be less than one.')
  call mpabrt ( 592)
endif

if (nn == 0) then
!+   call mpdmc (1.d0, 0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrsetd (t1(1:), 1.d0, mprnd)
 
!+   call mpsub (t1, x, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrsub (t2(1:), t1(1:), x(1:), mprnd)
 
!+   call mpdiv (x, t2, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrdiv (t3(1:), x(1:), t2(1:), mprnd)
 
!+   call mproun (t3, mpnw
  call mpfrprecround (t3(1:), mpnw, mprnd)
 
!+   call mpeq (t3, y, mpnw
  call mpfrsetprec (y(1:), mpnw)
  call mpfrset (y(1:), t3(1:), mprnd)
 
else
!+   call mpeq (x, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrset (t1(1:), x(1:), mprnd)
 
!+   call mpeq (x, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrset (t2(1:), x(1:), mprnd)
 

  do k = 2, itrmax
!+     call mpmul (x, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrmul (t3(1:), x(1:), t2(1:), mprnd)
 
!+     call mpeq (t3, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrset (t2(1:), t3(1:), mprnd)
 
!+     call mpdmc (dble (k), 0, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrsetd (t3(1:), dble (k), mprnd)
 
!+     call mpnpwr (t3, nn, t4, mpnw1
    call mpfrsetprec (t4(1:), mpnw1)
    call mpfrpowsi (t4(1:), t3(1:), nn, mprnd)
 
!+     call mpdiv (t2, t4, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrdiv (t3(1:), t2(1:), t4(1:), mprnd)
 
!+     call mpadd (t1, t3, t4, mpnw1
    call mpfrsetprec (t4(1:), mpnw1)
    call mpfradd (t4(1:), t1(1:), t3(1:), mprnd)
 
!+     call mpeq (t4, t1, mpnw1
    call mpfrsetprec (t1(1:), mpnw1)
    call mpfrset (t1(1:), t4(1:), mprnd)
 

!+     call mpabs (t3, tc1, 4
    call mpfrsetprec (tc1(1:), 4)
    call mpfrabs (tc1(1:), t3(1:), mprnd)
 
!+     call mpmul (eps, t1, tc3, 4
    call mpfrsetprec (tc3(1:), 4)
    call mpfrmul (tc3(1:), eps(1:), t1(1:), mprnd)
 
!+     call mpabs (tc3, tc2, 4
    call mpfrsetprec (tc2(1:), 4)
    call mpfrabs (tc2(1:), tc3(1:), mprnd)
 
!+     call mpcpr (tc1, tc2, ic1, 4
    ic1 = mpfrcmp (tc1(1:), tc2(1:))
 
    if (ic1 <= 0) goto 100
  enddo

  write (mpldb, 4)
4 format ('*** MPPOLYLOGPOS: Loop end error')
  call mpabrt ( 593)

100 continue

!+   call mproun (t1, mpnw
  call mpfrprecround (t1(1:), mpnw, mprnd)
 
!+   call mpeq (t1, y, mpnw
  call mpfrsetprec (y(1:), mpnw)
  call mpfrset (y(1:), t1(1:), mprnd)
 
endif

return
end subroutine mppolylogpos

subroutine mpstruvehn (nu, ss, zz, mpnw)

!   This returns the StruveH function with integer arg NU and MPFR argument SS.

implicit none
integer, intent(in):: nu, mpnw
integer (mpiknd), intent(in):: ss(0:)
integer (mpiknd), intent(out):: zz(0:)
integer, parameter:: itrmax = 1000000
real (mprknd), parameter:: dmax = 1000.d0, pi = 3.1415926535897932385d0
integer ic1, mpnw1, k, n1
real (mprknd) d1
integer (mpiknd) sum(0:2*mpnw+6), td1(0:2*mpnw+6), td2(0:2*mpnw+6), tn1(0:2*mpnw+6), &
  tnm1(0:2*mpnw+6), t1(0:2*mpnw+6), t2(0:2*mpnw+6), tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)

integer, external:: mpfrsgn
integer, external:: mpfrcmp
real (mprknd), external:: mpfrgetd2exp
integer (mpiknd) ix8
 
!  End of declaration

if (mpnw < 4 .or. mpspacer (ss) < mpnw + 4 .or. mpspacer (zz) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPSTRUVEHN: Uninitialized or inadequately sized arrays')
  call mpabrt ( 594)
endif

!  Check if PI has been precomputed.

!+ call mpmdc (mppicon, d1, n1, mpnw+1
call mpfixlocr (mppicon)
d1 = 2.d0 * mpfrgetd2exp (ix8, mppicon(1:), mprnd)
if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
if (n1 /= 1 .or. abs (d1 * 2.d0**n1 - pi) > mprdfz &
  .or. mpwprecr (mppicon) < mpnw+1) then
  write (mpldb, 4) mpnw+1
4 format ('*** MPSTRUVEHN: Pi must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 597)
endif

!  End of initial input array check

if (nu < 0) then
  write (mpldb, 2)
2 format ('*** MPSTRUVEHN: NU < 0')
  call mpabrt ( 595)
endif

!+ call mpmdc (ss, d1, n1, mpnw
d1 = 2.d0 * mpfrgetd2exp (ix8, ss(1:), mprnd)
if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
d1 = abs (d1) * 2.d0**n1
if (d1 > dmax) then
  write (mpldb, 3)
3 format ('*** MPSTRUVEHN: ABS(SS) >',f8.2)
  call mpabrt ( 596)
endif

mpnw1 = min (mpnw * (1.d0 + d1 / dmax), dble (mpnwx))

call mpinitwds (sum, 2*mpnw+1)
call mpinitwds (td1, 2*mpnw+1)
call mpinitwds (td2, 2*mpnw+1)
call mpinitwds (tn1, 2*mpnw+1)
call mpinitwds (tnm1, 2*mpnw+1)
call mpinitwds (t1, 2*mpnw+1)
call mpinitwds (t2, 2*mpnw+1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
!+ call mpdmc (2.d0, 0, tc1, 4
call mpfrsetprec (tc1(1:), 4)
call mpfrsetd (tc1(1:), 2.d0, mprnd)
 
!+ call mpnpwr (tc1, -mpnw*mpnbt, eps, 4
call mpfrsetprec (eps(1:), 4)
call mpfrpowsi (eps(1:), tc1(1:), -mpnw*mpnbt, mprnd)
 

! tn1 = mpreal (1.d0, nwds1)
! tnm1 = -0.25d0 * mpreal (ss, nwds1)**2
! td1 = 0.5d0 * sqrt (mppi (nwds1))
! td2 = td1

!+ call mpdmc (1.d0, 0, tn1, mpnw1
call mpfrsetprec (tn1(1:), mpnw1)
call mpfrsetd (tn1(1:), 1.d0, mprnd)
 
!+ call mpmul (ss, ss, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrmul (t1(1:), ss(1:), ss(1:), mprnd)
 
!+ call mpmuld (t1, -0.25d0, tnm1, mpnw1
call mpfrsetprec (tnm1(1:), mpnw1)
call mpfrmuld (tnm1(1:), t1(1:), -0.25d0, mprnd)
 
!+ call mpsqrt (mppicon, t1, mpnw1
call mpfixlocr (mppicon)
call mpfrsetprec (t1(1:), mpnw1)
call mpfrsqrt (t1(1:), mppicon(1:), mprnd)
 
!+ call mpmuld (t1, 0.5d0, td1, mpnw1
call mpfrsetprec (td1(1:), mpnw1)
call mpfrmuld (td1(1:), t1(1:), 0.5d0, mprnd)
 
!+ call mpeq (td1, td2, mpnw1
call mpfrsetprec (td2(1:), mpnw1)
call mpfrset (td2(1:), td1(1:), mprnd)
 

! do k = 1, nu
!  td2 = (k + 0.5d0) * td2
! enddo

do k = 1, nu
  d1 = k + 0.5d0
!+   call mpmuld (td2, d1, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmuld (t1(1:), td2(1:), d1, mprnd)
 
!+   call mpeq (t1, td2, mpnw1
  call mpfrsetprec (td2(1:), mpnw1)
  call mpfrset (td2(1:), t1(1:), mprnd)
 
enddo

! sum = tn1 / (td1 * td2)

!+ call mpmul (td1, td2, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfrmul (t2(1:), td1(1:), td2(1:), mprnd)
 
!+ call mpdiv (tn1, t2, sum, mpnw1
call mpfrsetprec (sum(1:), mpnw1)
call mpfrdiv (sum(1:), tn1(1:), t2(1:), mprnd)
 

do k = 1, itrmax

!  tn1 = tnm1 * tn1
!  td1 = (k + 0.5d0) * td1
!  td2 = (nu + k + 0.5d0) * td2
!  t1 = tn1 / (td1 * td2)
!  sum = sum + t1

!+   call mpmul (tnm1, tn1, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmul (t1(1:), tnm1(1:), tn1(1:), mprnd)
 
!+   call mpeq (t1, tn1, mpnw1
  call mpfrsetprec (tn1(1:), mpnw1)
  call mpfrset (tn1(1:), t1(1:), mprnd)
 
  d1 = k + 0.5d0
!+   call mpmuld (td1, d1, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmuld (t1(1:), td1(1:), d1, mprnd)
 
!+   call mpeq (t1, td1, mpnw1
  call mpfrsetprec (td1(1:), mpnw1)
  call mpfrset (td1(1:), t1(1:), mprnd)
 
  d1 = nu + k + 0.5d0
!+   call mpmuld (td2, d1, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmuld (t1(1:), td2(1:), d1, mprnd)
 
!+   call mpeq (t1, td2, mpnw1
  call mpfrsetprec (td2(1:), mpnw1)
  call mpfrset (td2(1:), t1(1:), mprnd)
 
!+   call mpmul (td1, td2, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmul (t2(1:), td1(1:), td2(1:), mprnd)
 
!+   call mpdiv (tn1, t2, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrdiv (t1(1:), tn1(1:), t2(1:), mprnd)
 
!+   call mpadd (sum, t1, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfradd (t2(1:), sum(1:), t1(1:), mprnd)
 
!+   call mpeq (t2, sum, mpnw1
  call mpfrsetprec (sum(1:), mpnw1)
  call mpfrset (sum(1:), t2(1:), mprnd)
 

!  if (abs (t1) < eps) goto 100

!+   call mpabs (t1, tc1, 4
  call mpfrsetprec (tc1(1:), 4)
  call mpfrabs (tc1(1:), t1(1:), mprnd)
 
!+   call mpmul (eps, sum, tc3, 4
  call mpfrsetprec (tc3(1:), 4)
  call mpfrmul (tc3(1:), eps(1:), sum(1:), mprnd)
 
!+   call mpabs (tc3, tc2, 4
  call mpfrsetprec (tc2(1:), 4)
  call mpfrabs (tc2(1:), tc3(1:), mprnd)
 
!+   call mpcpr (tc1, tc2, ic1, 4
  ic1 = mpfrcmp (tc1(1:), tc2(1:))
 
  if (ic1 <= 0) goto 100
enddo

write (mpldb, 5)
5 format ('*** MPSTRUVEHN: Loop end error')
call mpabrt ( 598)

100 continue

! struvehn = (0.5d0 * ss)**(nu + 1) * sum

!+ call mpmuld (ss, 0.5d0, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrmuld (t1(1:), ss(1:), 0.5d0, mprnd)
 
n1 = nu + 1
!+ call mpnpwr (t1, n1, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfrpowsi (t2(1:), t1(1:), n1, mprnd)
 
!+ call mpmul (t2, sum, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrmul (t1(1:), t2(1:), sum(1:), mprnd)
 
!+ call mproun (t1, mpnw
call mpfrprecround (t1(1:), mpnw, mprnd)
 
!+ call mpeq (t1, zz, mpnw
call mpfrsetprec (zz(1:), mpnw)
call mpfrset (zz(1:), t1(1:), mprnd)
 
return
end subroutine mpstruvehn

subroutine mpzetar (ss, zz, mpnw)

!   This returns the zeta function of an MPR argument SS using an algorithm
!   due to Peter Borwein.

implicit none
integer (mpiknd), intent(in):: ss(0:)
integer (mpiknd), intent(out):: zz(0:)
integer, intent(in):: mpnw
integer, parameter:: itrmax = 1000000
real (mprknd), parameter:: dfrac = 1.d0+ceiling(mpdpw), pi = 3.1415926535897932385d0
integer i, ic1, iss, j, mpnw1, n, n1, n2
real (mprknd) d1, d2
integer (mpiknd) f1(0:mpnw+6), s(0:mpnw+6), t1(0:mpnw+6), t2(0:mpnw+6), &
  t3(0:mpnw+6), t4(0:mpnw+6), t5(0:mpnw+6), tn(0:mpnw+6), tt(0:mpnw+6), &
  tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)
real (mprknd) sgn

integer, external:: mpfrsgn
integer, external:: mpfrcmp
real (mprknd), external:: mpfrgetd2exp
integer (mpiknd) ix8
 
!  End of declaration

if (mpnw < 4 .or. mpspacer (ss) < mpnw + 4 .or. mpspacer (zz) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPZETAR: Uninitialized or inadequately sized arrays')
  call mpabrt ( 599)
endif

!  Check if PI has been precomputed.

!+ call mpmdc (mppicon, d1, n1, mpnw+1
call mpfixlocr (mppicon)
d1 = 2.d0 * mpfrgetd2exp (ix8, mppicon(1:), mprnd)
if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
if (n1 /= 1 .or. abs (d1 * 2.d0**n1 - pi) > mprdfz &
  .or. mpwprecr (mppicon) < mpnw+1) then
  write (mpldb, 4) mpnw+1
4 format ('*** MPZETAR: Pi must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 600)
endif

!  End of initial input array check

mpnw1 = min (mpnw + 1, mpnwx)
call mpinitwds (f1, mpnw1)
call mpinitwds (s, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (tn, mpnw1)
call mpinitwds (tt, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
!+ call mpdmc (2.d0, 0, tc1, 4
call mpfrsetprec (tc1(1:), 4)
call mpfrsetd (tc1(1:), 2.d0, mprnd)
 
!+ call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4
call mpfrsetprec (eps(1:), 4)
call mpfrpowsi (eps(1:), tc1(1:), -mpnw1*mpnbt, mprnd)
 

!+ call mpdmc (1.d0, 0, f1, mpnw1
call mpfrsetprec (f1(1:), mpnw1)
call mpfrsetd (f1(1:), 1.d0, mprnd)
 
!+ call mpcpr (ss, f1, ic1, mpnw1
ic1 = mpfrcmp (ss(1:), f1(1:))
 
!+ call mpinfr (ss, t1, t2, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrsetprec (t2(1:), mpnw1)
call mpfrtrunc (t1(1:), ss(1:))
call mpfrsub (t2(1:), ss(1:), t1(1:), mprnd)
 

if (ic1 == 0) then
  write (mpldb, 2)
2 format ('*** MPZETAR: argument is 1')
  call mpabrt ( 601)
!+ elseif (mpsgn (t2) == 0) the
elseif (mpfrsgn (t2(1:)) == 0) then
 

!   The argument is an integer value. Call mpzetaintr instead.

!+   call mpmdc (ss, d1, n1, mpnw
  d1 = 2.d0 * mpfrgetd2exp (ix8, ss(1:), mprnd)
  if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
  iss = d1 * 2.d0**n1
  call mpzetaintr (iss, t1, mpnw)
  goto 200
!+ elseif (mpsgn (ss) < 0) the
elseif (mpfrsgn (ss(1:)) < 0) then
 

!   If arg < 0, compute zeta(1-ss), and later apply Riemann's formula.

!+   call mpsub (f1, ss, tt, mpnw1
  call mpfrsetprec (tt(1:), mpnw1)
  call mpfrsub (tt(1:), f1(1:), ss(1:), mprnd)
 
else
!+   call mpeq (ss, tt, mpnw1
  call mpfrsetprec (tt(1:), mpnw1)
  call mpfrset (tt(1:), ss(1:), mprnd)
 
endif

!  Check if argument is large enough that computing with definition is faster.

d1 = mpnbt * mpnw * log (2.d0) / log (2.d0 * mpnbt * mpnw / 3.d0)
!+ call mpmdc (tt, d2, n2, mpnw1
d2 = 2.d0 * mpfrgetd2exp (ix8, tt(1:), mprnd)
if (d2 == 0) then; n2 = 0; else; n2 = ix8 - 1; endif
 
d2 = d2 * 2.d0 ** n2

if (d2 > d1) then

!   Evaluate the infinite series.

!+   call mpdmc (1.d0, 0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrsetd (t1(1:), 1.d0, mprnd)
 

  do i = 2, itrmax
!+     call mpdmc (dble (i), 0, t4, mpnw1
    call mpfrsetprec (t4(1:), mpnw1)
    call mpfrsetd (t4(1:), dble (i), mprnd)
 
!+     call mppower (t4, tt, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrpow (t2(1:), t4(1:), tt(1:), mprnd)
 
!+     call mpdiv (f1, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrdiv (t3(1:), f1(1:), t2(1:), mprnd)
 
!+     call mpadd (t1, t3, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfradd (t2(1:), t1(1:), t3(1:), mprnd)
 
!+     call mpeq (t2, t1, mpnw1
    call mpfrsetprec (t1(1:), mpnw1)
    call mpfrset (t1(1:), t2(1:), mprnd)
 

!+     call mpabs (t3, tc1, 4
    call mpfrsetprec (tc1(1:), 4)
    call mpfrabs (tc1(1:), t3(1:), mprnd)
 
!+     call mpmul (eps, t1, tc3, 4
    call mpfrsetprec (tc3(1:), 4)
    call mpfrmul (tc3(1:), eps(1:), t1(1:), mprnd)
 
!+     call mpabs (tc3, tc2, 4
    call mpfrsetprec (tc2(1:), 4)
    call mpfrabs (tc2(1:), tc3(1:), mprnd)
 
!+     call mpcpr (tc1, tc2, ic1, 4
    ic1 = mpfrcmp (tc1(1:), tc2(1:))
 
    if (ic1 <= 0) goto 200
  enddo

  write (mpldb, 3) 1, itrmax
3 format ('*** MPZETAR: iteration limit exceeded',2i10)
  call mpabrt ( 602)
endif

n = dfrac * mpnw1
!+ call mpdmc (2.d0, 0, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrsetd (t1(1:), 2.d0, mprnd)
 
!+ call mpnpwr (t1, n, tn, mpnw1
call mpfrsetprec (tn(1:), mpnw1)
call mpfrpowsi (tn(1:), t1(1:), n, mprnd)
 
!+ call mpneg (tn, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrneg (t1(1:), tn(1:), mprnd)
 
!+ call mpdmc (0.d0, 0, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfrsetd (t2(1:), 0.d0, mprnd)
 
!+ call mpdmc (0.d0, 0, s, mpnw1
call mpfrsetprec (s(1:), mpnw1)
call mpfrsetd (s(1:), 0.d0, mprnd)
 

sgn = 1.d0

do j = 0, 2 * n - 1
!+   call mpdmc (dble (j + 1), 0, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfrsetd (t4(1:), dble (j + 1), mprnd)
 
!+   call mppower (t4, tt, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrpow (t3(1:), t4(1:), tt(1:), mprnd)
 
!+   call mpdiv (t1, t3, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfrdiv (t4(1:), t1(1:), t3(1:), mprnd)
 
!+   call mpmuld (t4, sgn, t5, mpnw1
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfrmuld (t5(1:), t4(1:), sgn, mprnd)
 
!+   call mpadd (s, t5, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfradd (t4(1:), s(1:), t5(1:), mprnd)
 
!+   call mpeq (t4, s, mpnw1
  call mpfrsetprec (s(1:), mpnw1)
  call mpfrset (s(1:), t4(1:), mprnd)
 
  sgn = - sgn

  if (j < n - 1) then
!+     call mpdmc (0.d0, 0, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrsetd (t2(1:), 0.d0, mprnd)
 
  elseif (j == n - 1) then
!+     call mpdmc (1.d0, 0, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrsetd (t2(1:), 1.d0, mprnd)
 
  else
!+     call mpmuld (t2, dble (2 * n - j), t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrmuld (t3(1:), t2(1:), dble (2 * n - j), mprnd)
 
!+     call mpdivd (t3, dble (j + 1 - n), t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrdivd (t2(1:), t3(1:), dble (j + 1 - n), mprnd)
 
  endif
!+   call mpadd (t1, t2, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfradd (t3(1:), t1(1:), t2(1:), mprnd)
 
!+   call mpeq (t3, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrset (t1(1:), t3(1:), mprnd)
 
enddo

!+ call mpsub (f1, tt, t3, mpnw1
call mpfrsetprec (t3(1:), mpnw1)
call mpfrsub (t3(1:), f1(1:), tt(1:), mprnd)
 
!+ call mpdmc (2.d0, 0, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfrsetd (t2(1:), 2.d0, mprnd)
 
!+ call mppower (t2, t3, t4, mpnw1
call mpfrsetprec (t4(1:), mpnw1)
call mpfrpow (t4(1:), t2(1:), t3(1:), mprnd)
 
!+ call mpsub (f1, t4, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfrsub (t2(1:), f1(1:), t4(1:), mprnd)
 
!+ call mpmul (tn, t2, t3, mpnw1
call mpfrsetprec (t3(1:), mpnw1)
call mpfrmul (t3(1:), tn(1:), t2(1:), mprnd)
 
!+ call mpdiv (s, t3, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrdiv (t1(1:), s(1:), t3(1:), mprnd)
 
!+ call mpneg (t1, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfrneg (t2(1:), t1(1:), mprnd)
 
!+ call mpeq (t2, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrset (t1(1:), t2(1:), mprnd)
 

!   If original argument was negative, apply Riemann's formula.

!+ if (mpsgn (ss) < 0) the
if (mpfrsgn (ss(1:)) < 0) then
 
!+   call mpgammar (tt, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrgamma (t3(1:), tt(1:), mprnd)
 
!+   call mpmul (t1, t3, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmul (t2(1:), t1(1:), t3(1:), mprnd)
 
!+   call mpmul (mppicon, tt, t1, mpnw1
call mpfixlocr (mppicon)
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmul (t1(1:), mppicon(1:), tt(1:), mprnd)
 
!+   call mpmuld (t1, 0.5d0, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrmuld (t3(1:), t1(1:), 0.5d0, mprnd)
 
!+   call mpcssnr (t3, t4, t5, mpnw1
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfrsincos (t5(1:), t4(1:), t3(1:), mprnd)
 
!+   call mpmul (t2, t4, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmul (t1(1:), t2(1:), t4(1:), mprnd)
 
!+   call mpmuld (mppicon, 2.d0, t2, mpnw1
call mpfixlocr (mppicon)
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmuld (t2(1:), mppicon(1:), 2.d0, mprnd)
 
!+   call mppower (t2, tt, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrpow (t3(1:), t2(1:), tt(1:), mprnd)
 
!+   call mpdiv (t1, t3, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrdiv (t2(1:), t1(1:), t3(1:), mprnd)
 
!+   call mpmuld (t2, 2.d0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmuld (t1(1:), t2(1:), 2.d0, mprnd)
 
endif

200 continue

!+ call mproun (t1, mpnw
call mpfrprecround (t1(1:), mpnw, mprnd)
 
!+ call mpeq (t1, zz, mpnw
call mpfrsetprec (zz(1:), mpnw)
call mpfrset (zz(1:), t1(1:), mprnd)
 
return
end subroutine mpzetar

subroutine mpzetaintr (iss, zz, mpnw)

!   This returns the zeta function of the integer argument ISS using an algorithm
!   due to Peter Borwein.

implicit none
integer, intent(in):: iss, mpnw
integer (mpiknd), intent(out):: zz(0:)
integer, parameter:: itrmax = 1000000
real (mprknd), parameter:: dfrac = 1.d0+ceiling(mpdpw), pi = 3.1415926535897932385d0
integer i, ic1, j, mpnw1, n, n1, itt
real (mprknd) d1, sgn
integer (mpiknd) f1(0:mpnw+6), s(0:mpnw+6), t1(0:mpnw+6), t2(0:mpnw+6), &
  t3(0:mpnw+6), t4(0:mpnw+6), t5(0:mpnw+6), tn(0:mpnw+6), &
  tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)

integer, external:: mpfrsgn
integer, external:: mpfrcmp
real (mprknd), external:: mpfrgetd2exp
integer (mpiknd) ix8
 
!  End of declaration

if (mpnw < 4 .or. mpspacer (zz) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPZETAINTR: Uninitialized or inadequately sized arrays')
  call mpabrt ( 603)
endif

!  Check if PI has been precomputed.

!+ call mpmdc (mppicon, d1, n1, mpnw+1
call mpfixlocr (mppicon)
d1 = 2.d0 * mpfrgetd2exp (ix8, mppicon(1:), mprnd)
if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
if (n1 /= 1 .or. abs (d1 * 2.d0**n1 - pi) > mprdfz &
  .or. mpwprecr (mppicon) < mpnw+1) then
  write (mpldb, 4) mpnw+1
4 format ('*** MPZETAINTR: Pi must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 604)
endif

!  End of initial input array check

mpnw1 = min (mpnw + 1, mpnwx)
call mpinitwds (f1, mpnw1)
call mpinitwds (s, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (tn, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
!+ call mpdmc (2.d0, 0, tc1, 4
call mpfrsetprec (tc1(1:), 4)
call mpfrsetd (tc1(1:), 2.d0, mprnd)
 
!+ call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4
call mpfrsetprec (eps(1:), 4)
call mpfrpowsi (eps(1:), tc1(1:), -mpnw1*mpnbt, mprnd)
 

!+ call mpdmc (1.d0, 0, f1, mpnw1
call mpfrsetprec (f1(1:), mpnw1)
call mpfrsetd (f1(1:), 1.d0, mprnd)
 

if (iss == 1) then
  write (mpldb, 2)
2 format ('*** MPZETAINTR: argument is 1')
  call mpabrt ( 605)
elseif (iss == 0) then

!   Argument is zero -- result is -1/2.

!+   call mpdmc (-0.5d0, 0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrsetd (t1(1:), -0.5d0, mprnd)
 
  goto 200
elseif (iss < 0) then

!   If argument is a negative even integer, the result is zero.

  if (mod (iss, 2) == 0) then
!+     call mpdmc (0.d0, 0, t1, mpnw1
    call mpfrsetprec (t1(1:), mpnw1)
    call mpfrsetd (t1(1:), 0.d0, mprnd)
 
    goto 200
  endif

!   Otherwise if arg < 0, compute zeta(1-is), and later apply Riemann's formula.

  itt = 1 - iss
else
  itt = iss
endif

!  Check if argument is large enough that computing with definition is faster.

d1 = mpnbt * mpnw * log (2.d0) / log (2.d0 * mpnbt * mpnw / 3.d0)

if (itt > d1) then

!   Evaluate the infinite series.

!+   call mpdmc (1.d0, 0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrsetd (t1(1:), 1.d0, mprnd)
 

  do i = 2, itrmax
!+     call mpdmc (dble (i), 0, t4, mpnw1
    call mpfrsetprec (t4(1:), mpnw1)
    call mpfrsetd (t4(1:), dble (i), mprnd)
 
!+     call mpnpwr (t4, itt, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrpowsi (t2(1:), t4(1:), itt, mprnd)
 
!+     call mpdiv (f1, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrdiv (t3(1:), f1(1:), t2(1:), mprnd)
 
!+     call mpadd (t1, t3, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfradd (t2(1:), t1(1:), t3(1:), mprnd)
 
!+     call mpeq (t2, t1, mpnw1
    call mpfrsetprec (t1(1:), mpnw1)
    call mpfrset (t1(1:), t2(1:), mprnd)
 

!+     call mpabs (t3, tc1, 4
    call mpfrsetprec (tc1(1:), 4)
    call mpfrabs (tc1(1:), t3(1:), mprnd)
 
!+     call mpmul (eps, t1, tc3, 4
    call mpfrsetprec (tc3(1:), 4)
    call mpfrmul (tc3(1:), eps(1:), t1(1:), mprnd)
 
!+     call mpabs (tc3, tc2, 4
    call mpfrsetprec (tc2(1:), 4)
    call mpfrabs (tc2(1:), tc3(1:), mprnd)
 
!+     call mpcpr (tc1, tc2, ic1, 4
    ic1 = mpfrcmp (tc1(1:), tc2(1:))
 
    if (ic1 <= 0) goto 200
  enddo

  write (mpldb, 3) 1, itrmax
3 format ('*** MPZETAINTR: iteration limit exceeded',2i10)
  call mpabrt ( 606)
endif

n = dfrac * mpnw1
!+ call mpdmc (2.d0, 0, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrsetd (t1(1:), 2.d0, mprnd)
 
!+ call mpnpwr (t1, n, tn, mpnw1
call mpfrsetprec (tn(1:), mpnw1)
call mpfrpowsi (tn(1:), t1(1:), n, mprnd)
 
!+ call mpneg (tn, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrneg (t1(1:), tn(1:), mprnd)
 
!+ call mpdmc (0.d0, 0, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfrsetd (t2(1:), 0.d0, mprnd)
 
!+ call mpdmc (0.d0, 0, s, mpnw1
call mpfrsetprec (s(1:), mpnw1)
call mpfrsetd (s(1:), 0.d0, mprnd)
 

sgn = 1.d0

do j = 0, 2 * n - 1
!+   call mpdmc (dble (j + 1), 0, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfrsetd (t4(1:), dble (j + 1), mprnd)
 
!+   call mpnpwr (t4, itt, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrpowsi (t3(1:), t4(1:), itt, mprnd)
 
!+   call mpdiv (t1, t3, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfrdiv (t4(1:), t1(1:), t3(1:), mprnd)
 
!+   call mpmuld (t4, sgn, t5, mpnw1
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfrmuld (t5(1:), t4(1:), sgn, mprnd)
 
!+   call mpadd (s, t5, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfradd (t4(1:), s(1:), t5(1:), mprnd)
 
!+   call mpeq (t4, s, mpnw1
  call mpfrsetprec (s(1:), mpnw1)
  call mpfrset (s(1:), t4(1:), mprnd)
 
  sgn = - sgn

  if (j < n - 1) then
!+     call mpdmc (0.d0, 0, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrsetd (t2(1:), 0.d0, mprnd)
 
  elseif (j == n - 1) then
!+     call mpdmc (1.d0, 0, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrsetd (t2(1:), 1.d0, mprnd)
 
  else
!+     call mpmuld (t2, dble (2 * n - j), t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrmuld (t3(1:), t2(1:), dble (2 * n - j), mprnd)
 
!+     call mpdivd (t3, dble (j + 1 - n), t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrdivd (t2(1:), t3(1:), dble (j + 1 - n), mprnd)
 
  endif

!+   call mpadd (t1, t2, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfradd (t3(1:), t1(1:), t2(1:), mprnd)
 
!+   call mpeq (t3, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrset (t1(1:), t3(1:), mprnd)
 
enddo

!+ call mpdmc (2.d0, 0, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfrsetd (t2(1:), 2.d0, mprnd)
 
!+ call mpnpwr (t2, 1 - itt, t4, mpnw1
call mpfrsetprec (t4(1:), mpnw1)
call mpfrpowsi (t4(1:), t2(1:), 1 - itt, mprnd)
 
!+ call mpsub (f1, t4, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfrsub (t2(1:), f1(1:), t4(1:), mprnd)
 
!+ call mpmul (tn, t2, t3, mpnw1
call mpfrsetprec (t3(1:), mpnw1)
call mpfrmul (t3(1:), tn(1:), t2(1:), mprnd)
 
!+ call mpdiv (s, t3, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrdiv (t1(1:), s(1:), t3(1:), mprnd)
 
!+ call mpneg (t1, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfrneg (t2(1:), t1(1:), mprnd)
 
!+ call mpeq (t2, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfrset (t1(1:), t2(1:), mprnd)
 

!   If original argument was negative, apply Riemann's formula.

if (iss < 0) then
!+   call mpdmc (1.d0, 0, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrsetd (t3(1:), 1.d0, mprnd)
 
  do i = 1, itt - 1
!+     call mpmuld (t3, dble (i), t4, mpnw1
    call mpfrsetprec (t4(1:), mpnw1)
    call mpfrmuld (t4(1:), t3(1:), dble (i), mprnd)
 
!+     call mpeq (t4, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrset (t3(1:), t4(1:), mprnd)
 
  enddo

!+   call mpmul (t1, t3, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmul (t2(1:), t1(1:), t3(1:), mprnd)
 
!+   call mpmuld (mppicon, dble (itt), t1, mpnw1
call mpfixlocr (mppicon)
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmuld (t1(1:), mppicon(1:), dble (itt), mprnd)
 
!+   call mpmuld (t1, 0.5d0, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrmuld (t3(1:), t1(1:), 0.5d0, mprnd)
 
!+   call mpcssnr (t3, t4, t5, mpnw1
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfrsincos (t5(1:), t4(1:), t3(1:), mprnd)
 
!+   call mpmul (t2, t4, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmul (t1(1:), t2(1:), t4(1:), mprnd)
 
!+   call mpmuld (mppicon, 2.d0, t2, mpnw1
call mpfixlocr (mppicon)
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmuld (t2(1:), mppicon(1:), 2.d0, mprnd)
 
!+   call mpnpwr (t2, itt, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrpowsi (t3(1:), t2(1:), itt, mprnd)
 
!+   call mpdiv (t1, t3, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrdiv (t2(1:), t1(1:), t3(1:), mprnd)
 
!+   call mpmuld (t2, 2.d0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmuld (t1(1:), t2(1:), 2.d0, mprnd)
 
endif

200 continue

!+ call mproun (t1, mpnw
call mpfrprecround (t1(1:), mpnw, mprnd)
 
!+ call mpeq (t1, zz, mpnw
call mpfrsetprec (zz(1:), mpnw)
call mpfrset (zz(1:), t1(1:), mprnd)
 
return
end subroutine mpzetaintr

subroutine mpzetabe (nb1, nb2, berne, s, z, mpnw)

!  This evaluates the Riemann zeta function, using the combination of
!  the definition formula (for large s), and an Euler-Maclaurin scheme
!  (see formula 25.2.9 of the DLMF). The array berne contains precomputed
!  even Bernoulli numbers (see MPBERNER above). Its dimensions must be as
!  shown below. NB2 must be greater than 1.4 x precision in decimal digits.

implicit none
integer, intent(in):: nb1, nb2, mpnw
integer (mpiknd), intent(in):: berne(0:nb1+5,nb2), s(0:)
integer (mpiknd), intent(out):: z(0:)
integer, parameter:: itrmax = 1000000
real (mprknd), parameter:: dber = 1.5d0, dfrac = 0.6d0, &
  pi = 3.1415926535897932385d0
integer i, i1, i2, ic1, k, mpnw1, n1, n2, nn
real (mprknd) d1, d2
integer (mpiknd) t0(0:mpnw+6), t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), &
  t4(0:mpnw+6), t5(0:mpnw+6), t6(0:mpnw+6), t7(0:mpnw+6), t8(0:mpnw+6), &
  t9(0:mpnw+6), tt(0:mpnw+6), f1(0:mpnw+6), tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)

integer, external:: mpfrsgn
integer, external:: mpfrcmp
real (mprknd), external:: mpfrgetd2exp
integer (mpiknd) ix8
 
!  End of declaration

if (mpnw < 4 .or. mpspacer (s) < mpnw + 4 .or. mpspacer (z) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPZETABE: Uninitialized or inadequately sized arrays')
  call mpabrt ( 607)
endif

!  Check if PI has been precomputed.

!+ call mpmdc (mppicon, d1, n1, mpnw+1
call mpfixlocr (mppicon)
d1 = 2.d0 * mpfrgetd2exp (ix8, mppicon(1:), mprnd)
if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
if (n1 /= 1 .or. abs (d1 * 2.d0**n1 - pi) > mprdfz &
  .or. mpwprecr (mppicon) < mpnw+1) then
  write (mpldb, 5) mpnw+1
5 format ('*** MPZETABE: Pi must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 608)
endif

!  End of initial input array check

!   Check if berne array has been initialized.

!+ call mpmdc (berne(0:nb1+5,1), d1, n1, mpnw
d1 = 2.d0 * mpfrgetd2exp (ix8, berne(1:nb1+5,1), mprnd)
if (d1 == 0) then; n1 = 0; else; n1 = ix8 - 1; endif
 
d1 = d1 * 2.d0 ** n1
if (mpwprecr (berne(0:nb1+5,1)) < mpnw .or. &
  abs (d1 - 1.d0 / 6.d0) > mprdfz .or. nb2 < int (dber * mpdpw * mpnw)) then
  write (mpldb, 3) int (dber * mpdpw * mpnw)
3 format ('*** MPZETABE: Array of even Bernoulli coefficients must be initialized'/ &
   'with at least',i8,' entries.')
  call mpabrt ( 610)
endif

i = 0
k = 0
mpnw1 = min (mpnw + 1, mpnwx)
call mpinitwds (t0, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (t6, mpnw1)
call mpinitwds (t7, mpnw1)
call mpinitwds (t8, mpnw1)
call mpinitwds (t9, mpnw1)
call mpinitwds (tt, mpnw1)
call mpinitwds (f1, mpnw1)

call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
!+ call mpdmc (2.d0, 0, tc1, 4
call mpfrsetprec (tc1(1:), 4)
call mpfrsetd (tc1(1:), 2.d0, mprnd)
 
!+ call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4
call mpfrsetprec (eps(1:), 4)
call mpfrpowsi (eps(1:), tc1(1:), -mpnw1*mpnbt, mprnd)
 

!   Check if argument is 1 -- undefined.

!+ call mpdmc (1.d0, 0, t0, mpnw1
call mpfrsetprec (t0(1:), mpnw1)
call mpfrsetd (t0(1:), 1.d0, mprnd)
 
!+ call mpcpr (s, t0, ic1, mpnw1
ic1 = mpfrcmp (s(1:), t0(1:))
 
if (ic1 == 0) then
  write (mpldb, 2)
2 format ('*** MPZETABE: argument is 1')
  call mpabrt ( 609)
endif

!+ call mpdmc (1.d0, 0, f1, mpnw1
call mpfrsetprec (f1(1:), mpnw1)
call mpfrsetd (f1(1:), 1.d0, mprnd)
 

!   Check if argument is zero. If so, result is - 1/2.

!+ if (mpsgn (s) == 0) the
if (mpfrsgn (s(1:)) == 0) then
 
!+   call mpdmc (-0.5d0, 0, t1, mpnw
  call mpfrsetprec (t1(1:), mpnw)
  call mpfrsetd (t1(1:), -0.5d0, mprnd)
 
  goto 200
endif

!   Check if argument is negative.

!+ if (mpsgn (s) < 0) the
if (mpfrsgn (s(1:)) < 0) then
 

!   Check if argument is a negative even integer. If so, the result is zero.

!+   call mpmuld (s, 0.5d0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmuld (t1(1:), s(1:), 0.5d0, mprnd)
 
!+   call mpinfr (t1, t2, t3, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrtrunc (t2(1:), t1(1:))
  call mpfrsub (t3(1:), t1(1:), t2(1:), mprnd)
 
!+   if (mpsgn (t3) == 0) the
  if (mpfrsgn (t3(1:)) == 0) then
 
!+     call mpdmc (0.d0, 0, t1, mpnw1
    call mpfrsetprec (t1(1:), mpnw1)
    call mpfrsetd (t1(1:), 0.d0, mprnd)
 
    goto 200
  endif

!   Otherwise compute zeta(1-s), and later apply the reflection formula.

!+   call mpsub (f1, s, tt, mpnw1
  call mpfrsetprec (tt(1:), mpnw1)
  call mpfrsub (tt(1:), f1(1:), s(1:), mprnd)
 
else
!+   call mpeq (s, tt, mpnw1
  call mpfrsetprec (tt(1:), mpnw1)
  call mpfrset (tt(1:), s(1:), mprnd)
 
endif

!  Check if argument is large enough that computing with definition is faster.

d1 = mplogb * mpnw1 / log (32.d0 * mpnw1)
!+ call mpmdc (tt, d2, n2, mpnw1
d2 = 2.d0 * mpfrgetd2exp (ix8, tt(1:), mprnd)
if (d2 == 0) then; n2 = 0; else; n2 = ix8 - 1; endif
 
d2 = d2 * 2.d0**n2
if (d2 > d1) then
!+   call mpdmc (1.d0, 0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrsetd (t1(1:), 1.d0, mprnd)
 

  do i = 2, itrmax
!+     call mpdmc (dble (i), 0, t4, mpnw1
    call mpfrsetprec (t4(1:), mpnw1)
    call mpfrsetd (t4(1:), dble (i), mprnd)
 
!+     call mppower (t4, tt, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfrpow (t2(1:), t4(1:), tt(1:), mprnd)
 
!+     call mpdiv (f1, t2, t3, mpnw1
    call mpfrsetprec (t3(1:), mpnw1)
    call mpfrdiv (t3(1:), f1(1:), t2(1:), mprnd)
 
!+     call mpadd (t1, t3, t2, mpnw1
    call mpfrsetprec (t2(1:), mpnw1)
    call mpfradd (t2(1:), t1(1:), t3(1:), mprnd)
 
!+     call mpeq (t2, t1, mpnw1
    call mpfrsetprec (t1(1:), mpnw1)
    call mpfrset (t1(1:), t2(1:), mprnd)
 

!+     call mpabs (t3, tc1, 4
    call mpfrsetprec (tc1(1:), 4)
    call mpfrabs (tc1(1:), t3(1:), mprnd)
 
!+     call mpmul (eps, t1, tc3, 4
    call mpfrsetprec (tc3(1:), 4)
    call mpfrmul (tc3(1:), eps(1:), t1(1:), mprnd)
 
!+     call mpabs (tc3, tc2, 4
    call mpfrsetprec (tc2(1:), 4)
    call mpfrabs (tc2(1:), tc3(1:), mprnd)
 
!+     call mpcpr (tc1, tc2, ic1, 4
    ic1 = mpfrcmp (tc1(1:), tc2(1:))
 
    if (ic1 <= 0) goto 200
  enddo

  write (mpldb, 4) 1, itrmax
4 format ('*** MPZETABE: iteration limit exceeded',2i10)
  call mpabrt ( 611)
endif

!+ call mpdmc (1.d0, 0, t0, mpnw1
call mpfrsetprec (t0(1:), mpnw1)
call mpfrsetd (t0(1:), 1.d0, mprnd)
 
nn = dfrac * mpdpw * mpnw1

do k = 2, nn
!+   call mpdmc (dble (k), 0, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrsetd (t2(1:), dble (k), mprnd)
 
!+   call mppower (t2, tt, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrpow (t1(1:), t2(1:), tt(1:), mprnd)
 
!+   call mpdiv (f1, t1, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrdiv (t2(1:), f1(1:), t1(1:), mprnd)
 
!+   call mpadd (t0, t2, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfradd (t3(1:), t0(1:), t2(1:), mprnd)
 
!+   call mpeq (t3, t0, mpnw1
  call mpfrsetprec (t0(1:), mpnw1)
  call mpfrset (t0(1:), t3(1:), mprnd)
 
enddo

!+ call mpdmc (dble (nn), 0, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfrsetd (t2(1:), dble (nn), mprnd)
 
!+ call mpsub (tt, f1, t3, mpnw1
call mpfrsetprec (t3(1:), mpnw1)
call mpfrsub (t3(1:), tt(1:), f1(1:), mprnd)
 
!+ call mpmul (t1, t3, t4, mpnw1
call mpfrsetprec (t4(1:), mpnw1)
call mpfrmul (t4(1:), t1(1:), t3(1:), mprnd)
 
!+ call mpdiv (t2, t4, t3, mpnw1
call mpfrsetprec (t3(1:), mpnw1)
call mpfrdiv (t3(1:), t2(1:), t4(1:), mprnd)
 
!+ call mpadd (t0, t3, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfradd (t2(1:), t0(1:), t3(1:), mprnd)
 
!+ call mpdmc (0.5d0, 0, t3, mpnw1
call mpfrsetprec (t3(1:), mpnw1)
call mpfrsetd (t3(1:), 0.5d0, mprnd)
 
!+ call mpdiv (t3, t1, t4, mpnw1
call mpfrsetprec (t4(1:), mpnw1)
call mpfrdiv (t4(1:), t3(1:), t1(1:), mprnd)
 
!+ call mpsub (t2, t4, t0, mpnw1
call mpfrsetprec (t0(1:), mpnw1)
call mpfrsub (t0(1:), t2(1:), t4(1:), mprnd)
 

!+ call mpeq (tt, t3, mpnw1
call mpfrsetprec (t3(1:), mpnw1)
call mpfrset (t3(1:), tt(1:), mprnd)
 
d1 = 12.d0 * dble (nn)
!+ call mpmuld (t1, d1, t4, mpnw1
call mpfrsetprec (t4(1:), mpnw1)
call mpfrmuld (t4(1:), t1(1:), d1, mprnd)
 
!+ call mpdiv (t3, t4, t2, mpnw1
call mpfrsetprec (t2(1:), mpnw1)
call mpfrdiv (t2(1:), t3(1:), t4(1:), mprnd)
 
!+ call mpmuld (t1, dble (nn), t5, mpnw1
call mpfrsetprec (t5(1:), mpnw1)
call mpfrmuld (t5(1:), t1(1:), dble (nn), mprnd)
 
!+ call mpdmc (dble (nn), 0, t6, mpnw1
call mpfrsetprec (t6(1:), mpnw1)
call mpfrsetd (t6(1:), dble (nn), mprnd)
 
!+ call mpmul (t6, t6, t9, mpnw1
call mpfrsetprec (t9(1:), mpnw1)
call mpfrmul (t9(1:), t6(1:), t6(1:), mprnd)
 

do k = 2, min (nb2, itrmax)
!+   call mpdmc (dble (2 * k - 2), 0, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfrsetd (t4(1:), dble (2 * k - 2), mprnd)
 
!+   call mpadd (tt, t4, t6, mpnw1
  call mpfrsetprec (t6(1:), mpnw1)
  call mpfradd (t6(1:), tt(1:), t4(1:), mprnd)
 
!+   call mpdmc (dble (2 * k - 3), 0, t7, mpnw1
  call mpfrsetprec (t7(1:), mpnw1)
  call mpfrsetd (t7(1:), dble (2 * k - 3), mprnd)
 
!+   call mpadd (tt, t7, t8, mpnw1
  call mpfrsetprec (t8(1:), mpnw1)
  call mpfradd (t8(1:), tt(1:), t7(1:), mprnd)
 
!+   call mpmul (t6, t8, t7, mpnw1
  call mpfrsetprec (t7(1:), mpnw1)
  call mpfrmul (t7(1:), t6(1:), t8(1:), mprnd)
 
!+   call mpmul (t3, t7, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfrmul (t4(1:), t3(1:), t7(1:), mprnd)
 
!+   call mpdmc (dble (2 * k - 1), 0, t6, mpnw1
  call mpfrsetprec (t6(1:), mpnw1)
  call mpfrsetd (t6(1:), dble (2 * k - 1), mprnd)
 
!+   call mpdmc (dble (2 * k - 2), 0, t7, mpnw1
  call mpfrsetprec (t7(1:), mpnw1)
  call mpfrsetd (t7(1:), dble (2 * k - 2), mprnd)
 
!+   call mpmul (t6, t7, t8, mpnw1
  call mpfrsetprec (t8(1:), mpnw1)
  call mpfrmul (t8(1:), t6(1:), t7(1:), mprnd)
 
!+   call mpdiv (t4, t8, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrdiv (t3(1:), t4(1:), t8(1:), mprnd)
 
!+   call mpmul (t5, t9, t6, mpnw1
  call mpfrsetprec (t6(1:), mpnw1)
  call mpfrmul (t6(1:), t5(1:), t9(1:), mprnd)
 
!+   call mpeq (t6, t5, mpnw1
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfrset (t5(1:), t6(1:), mprnd)
 
!+   call mpmul (t3, berne(0:nb1+5,k), t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfrmul (t4(1:), t3(1:), berne(1:nb1+5,k), mprnd)
 
!+   call mpmuld (t5, dble (2 * k), t6, mpnw1
  call mpfrsetprec (t6(1:), mpnw1)
  call mpfrmuld (t6(1:), t5(1:), dble (2 * k), mprnd)
 
!+   call mpdiv (t4, t6, t7, mpnw1
  call mpfrsetprec (t7(1:), mpnw1)
  call mpfrdiv (t7(1:), t4(1:), t6(1:), mprnd)
 
!+   call mpadd (t2, t7, t4, mpnw1
  call mpfrsetprec (t4(1:), mpnw1)
  call mpfradd (t4(1:), t2(1:), t7(1:), mprnd)
 
!+   call mpeq (t4, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrset (t2(1:), t4(1:), mprnd)
 

!+   call mpabs (t7, tc1, 4
  call mpfrsetprec (tc1(1:), 4)
  call mpfrabs (tc1(1:), t7(1:), mprnd)
 
!+   call mpmul (eps, t2, tc3, 4
  call mpfrsetprec (tc3(1:), 4)
  call mpfrmul (tc3(1:), eps(1:), t2(1:), mprnd)
 
!+   call mpabs (tc3, tc2, 4
  call mpfrsetprec (tc2(1:), 4)
  call mpfrabs (tc2(1:), tc3(1:), mprnd)
 
!+   call mpcpr (tc1, tc2, ic1, 4
  ic1 = mpfrcmp (tc1(1:), tc2(1:))
 
  if (ic1 <= 0) goto 110
enddo

write (mpldb, 4) 2, min (nb2, itrmax)
call mpabrt ( 612)

110 continue

!+ call mpadd (t0, t2, t1, mpnw1
call mpfrsetprec (t1(1:), mpnw1)
call mpfradd (t1(1:), t0(1:), t2(1:), mprnd)
 

!   If original argument was negative, apply the reflection formula.

!+ if (mpsgn (s) < 0) the
if (mpfrsgn (s(1:)) < 0) then
 
!+   call mpgammar (tt, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrgamma (t3(1:), tt(1:), mprnd)
 
!+   call mpmul (t1, t3, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmul (t2(1:), t1(1:), t3(1:), mprnd)
 
!+   call mpmul (mppicon, tt, t1, mpnw1
call mpfixlocr (mppicon)
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmul (t1(1:), mppicon(1:), tt(1:), mprnd)
 
!+   call mpmuld (t1, 0.5d0, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrmuld (t3(1:), t1(1:), 0.5d0, mprnd)
 
!+   call mpcssnr (t3, t4, t5, mpnw1
  call mpfrsetprec (t5(1:), mpnw1)
  call mpfrsincos (t5(1:), t4(1:), t3(1:), mprnd)
 
!+   call mpmul (t2, t4, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmul (t1(1:), t2(1:), t4(1:), mprnd)
 
!+   call mpmuld (mppicon, 2.d0, t2, mpnw1
call mpfixlocr (mppicon)
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrmuld (t2(1:), mppicon(1:), 2.d0, mprnd)
 
!+   call mppower (t2, tt, t3, mpnw1
  call mpfrsetprec (t3(1:), mpnw1)
  call mpfrpow (t3(1:), t2(1:), tt(1:), mprnd)
 
!+   call mpdiv (t1, t3, t2, mpnw1
  call mpfrsetprec (t2(1:), mpnw1)
  call mpfrdiv (t2(1:), t1(1:), t3(1:), mprnd)
 
!+   call mpmuld (t2, 2.d0, t1, mpnw1
  call mpfrsetprec (t1(1:), mpnw1)
  call mpfrmuld (t1(1:), t2(1:), 2.d0, mprnd)
 
endif

200 continue

!+ call mproun (t1, mpnw
call mpfrprecround (t1(1:), mpnw, mprnd)
 
!+ call mpeq (t1, z, mpnw
call mpfrsetprec (z(1:), mpnw)
call mpfrset (z(1:), t1(1:), mprnd)
 

return
end subroutine mpzetabe

end module mpfune
