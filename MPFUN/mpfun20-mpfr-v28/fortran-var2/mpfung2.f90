!*****************************************************************************

!  MPFUN-MPFR: A thread-safe, MPFR-based arbitrary precision computation package
!  Main high-precision language interface module (module MPFUNG)
!  Variant 2: Precision level arguments ARE required.

!  Note: !> and !>> comments delimit variant differences, and are used to generate
!  the two variant files of this module using the gencodes.f90 program. Do not
!  change these comments.

!  Revision date:  7 Jan 2023

!  AUTHOR:
!    David H. Bailey
!    Lawrence Berkeley National Lab (retired) and University of California, Davis
!    Email: dhbailey@lbl.gov

!  COPYRIGHT AND DISCLAIMER:
!    All software in this package (c) 2023 David H. Bailey.
!    By downloading or using this software you agree to the copyright, disclaimer
!    and license agreement in the accompanying file DISCLAIMER.txt.

!  PURPOSE OF PACKAGE:
!    This package permits one to perform floating-point computations (real and
!    complex) to arbitrarily high numeric precision, by making only relatively
!    minor changes to existing Fortran-90 programs.  All basic arithmetic
!    operations and transcendental functions are supported, together with several
!    special functions.

!    In addition to fast execution times, one key feature of this package is a
!    100% THREAD-SAFE design, which means that user-level applications can be
!    easily converted for parallel execution, say using a threaded parallel
!    environment such as OpenMP.

!    This variant differs from the MPFUN20-Fort variant by the same author in that
!    this is based on MPFR, which presently is the fastest available low-level
!    package for high-precision floating-point computation.  At the Fortran user
!    level, most application codes written for MPFUN20-Fort may be compiled with
!    MPFUN-MPFR -- i.e., MPFUN-MPFR is plug-compatible with MPFUN20-Fort.

!  DOCUMENTATION:
!    A detailed description of this package, and instructions for compiling
!    and testing this program on various specific systems are included in the
!    README file accompanying this package, and, in more detail, in the
!    following technical paper:

!    David H. Bailey, "MPFUN2020: A new thread-safe arbitrary precision package,"
!    available at http://www.davidhbailey.com/dhbpapers/mpfun2020.pdf.

!  DESCRIPTION OF THIS MODULE (MPFUNG):
!    This module contains all high-level Fortran-90 language interfaces for
!    operations using the standard high-precision datatype.
!    There are two variants of this module:
!      mpfung1.f90	Mixed-mode operations ARE allowed; precision level
!                   arguments are NOT required in certain functions.
!      mpfung2.f90	Mixed-mode operations are NOT allowed; precision level
!                   arguments ARE required in certain functions.
!    Compile and link whichever one of these two is most appropriate for the
!    given applications. See documentation for details.

module mpfung
use mpfuna
use mpfune
use mpfunf
implicit none

!   The standard precision mp_real and mp_complex datatypes are defined here:

type mp_real
  sequence
  integer(8):: mpr(0:mpwds+5)
end type
type mp_complex
  sequence
  integer(8):: mpc(0:2*mpwds+11)
end type

!   The medium precision mp_realm and mp_complexm datatypes are defined here:

type mp_realm
  sequence
  integer (8):: mpr(0:mpwdsm+5)
end type
type mp_complexm
  sequence
  integer (8):: mpc(0:2*mpwdsm+11)
end type

!  Assignments and the five arithmetic operators:

private &
  mp_eqrr, mp_eqdr, mp_eqrd, mp_eqir, mp_eqri, mp_eqra, mp_eqrz, &
  mp_eqzr, mp_eqzz, mp_eqdz, mp_eqzd, mp_eqdcz, mp_eqzdc, &
  mp_addrr, mp_adddr, mp_addrd, mp_addir, mp_addri, mp_addzz, &
  mp_adddz, mp_addzd, mp_adddcz, mp_addzdc, mp_addrz, mp_addzr, &
  mp_subrr, mp_subdr, mp_subrd, mp_subir, mp_subri, mp_subzz, &
  mp_subdz, mp_subzd, mp_subdcz, mp_subzdc, mp_subrz, mp_subzr, &
  mp_negr, mp_negz, &
  mp_mulrr, mp_muldr, mp_mulrd, mp_mulir, mp_mulri, mp_mulzz, &
  mp_muldz, mp_mulzd, mp_muldcz, mp_mulzdc, mp_mulrz, mp_mulzr, &
  mp_divrr, mp_divdr, mp_divrd, mp_divir, mp_divri, mp_divzz, &
  mp_divdz, mp_divzd, mp_divdcz, mp_divzdc, mp_divrz, mp_divzr, &
  mp_expri, mp_exprr, mp_expzi, mp_expzz, mp_exprz, mp_expzr

!  The six comparison tests:

private &
  mp_eqtrr, mp_eqtdr, mp_eqtrd, mp_eqtir, mp_eqtri, mp_eqtzz, &
  mp_eqtdz, mp_eqtzd, mp_eqtdcz, mp_eqtzdc, mp_eqtrz, mp_eqtzr, &
  mp_netrr, mp_netdr, mp_netrd, mp_netir, mp_netri, mp_netzz, &
  mp_netdz, mp_netzd, mp_netdcz, mp_netzdc, mp_netrz, mp_netzr, &
  mp_letrr, mp_letdr, mp_letrd, mp_letir, mp_letri, &
  mp_getrr, mp_getdr, mp_getrd, mp_getir, mp_getri, &
  mp_lttrr, mp_lttdr, mp_lttrd, mp_lttir, mp_lttri, &
  mp_gttrr, mp_gttdr, mp_gttrd, mp_gttir, mp_gttri

!  Algebraic, transcendental and type conversion functions:

private &
  mp_abrt, mp_absr, mp_absz, mp_acos, mp_acosh, mp_agm, mp_aimag, &
  mp_aint, mp_airy, mp_anint, mp_asin, mp_asinh, mp_atan, &
  mp_atanh, mp_atan2, mp_ator1, mp_atorn, mp_berne, mp_bessel_i, &
  mp_bessel_in, mp_bessel_j, mp_bessel_j0, mp_bessel_j1, &
  mp_bessel_jn, mp_bessel_k, mp_bessel_kn, mp_bessel_y, &
  mp_bessel_y0, mp_bessel_y1, mp_bessel_yn, mp_binmd, mp_checkdp, &
  mp_checkqp, mp_ccos, mp_cexp, mp_clog, mp_conjg, mp_cos, &
  mp_cosh, mp_csin, mp_csqrt, mp_cssh, mp_cssn, mp_dctoz, &
  mp_dctoz2, mp_decmd, mp_digamma, mp_digamma_be, mp_dtor, &
  mp_dtor2, mp_eform, mp_egamma, mp_erf, mp_erfc, mp_exp, &
  mp_expint, mp_fixlocr, mp_fixlocrm, mp_fixlocz, mp_fixloczm, &
  mp_fform, mp_gamma, mp_hurwitz_zetan, mp_hurwitz_zetan_be, &
  mp_hypergeom_pfq, mp_hypot, mp_incgamma, mp_init, mp_initvr, &
  mp_initvz, mp_inpr, mp_inpz, mp_log, mp_log10, mp_log_gamma, &
  mp_log2, mp_max, mp_min, mp_mod, mp_mtor, mp_nrt, mp_outr, &
  mp_outz, mp_pi, mp_polygamma, mp_polygamma_be, mp_polylog_ini, &
  mp_polylog_neg, mp_polylog_pos, mp_prodd, mp_prodq, mp_qtor, &
  mp_qtor2, mp_quotd, mp_quotq, mp_rand, mp_readr1, mp_readr2, &
  mp_readr3, mp_readr4, mp_readr5, mp_readz1, mp_readz2, &
  mp_readz3, mp_readz4, mp_readz5, mp_rtod, mp_rtoq, mp_rtor, &
  mp_rtoz, mp_setwp, mp_sign, mp_sin, mp_sinh, mp_sqrt, &
  mp_struve_hn, mp_tan, mp_tanh, mp_wprec, mp_wprecz, mp_writer, &
  mp_writez, mp_zeta, mp_zeta_int, mp_zeta_be, mp_ztodc, &
  mp_ztor, mp_ztor2, mp_ztoz, mp_zmtoz

!  Operator extension interface blocks:

interface assignment (=)
  module procedure mp_eqrr
  module procedure mp_eqdr
  module procedure mp_eqir
  module procedure mp_eqrz
  module procedure mp_eqzr
  module procedure mp_eqzz
  module procedure mp_eqdz
  module procedure mp_eqdcz

!>  In variant #1, uncomment these lines:
!  module procedure mp_eqrd
!  module procedure mp_eqri
!  module procedure mp_eqra
!  module procedure mp_eqzd
!  module procedure mp_eqzdc
!>>
end interface

interface operator (+)
  module procedure mp_addrr
  module procedure mp_adddr
  module procedure mp_addrd
  module procedure mp_addir
  module procedure mp_addri
  module procedure mp_addzz
  module procedure mp_adddz
  module procedure mp_addzd
  module procedure mp_adddcz
  module procedure mp_addzdc
  module procedure mp_addrz
  module procedure mp_addzr
end interface

interface operator (-)
  module procedure mp_subrr
  module procedure mp_subdr
  module procedure mp_subrd
  module procedure mp_subir
  module procedure mp_subri
  module procedure mp_subzz
  module procedure mp_subdz
  module procedure mp_subzd
  module procedure mp_subdcz
  module procedure mp_subzdc
  module procedure mp_subrz
  module procedure mp_subzr
  module procedure mp_negr
  module procedure mp_negz
end interface

interface operator (*)
  module procedure mp_mulrr
  module procedure mp_muldr
  module procedure mp_mulrd
  module procedure mp_mulir
  module procedure mp_mulri
  module procedure mp_mulzz
  module procedure mp_muldz
  module procedure mp_mulzd
  module procedure mp_muldcz
  module procedure mp_mulzdc
  module procedure mp_mulrz
  module procedure mp_mulzr
end interface

interface operator (/)
  module procedure mp_divrr
  module procedure mp_divdr
  module procedure mp_divrd
  module procedure mp_divir
  module procedure mp_divri
  module procedure mp_divzz
  module procedure mp_divdz
  module procedure mp_divzd
  module procedure mp_divdcz
  module procedure mp_divzdc
  module procedure mp_divrz
  module procedure mp_divzr
end interface

interface operator (**)
   module procedure mp_expri
   module procedure mp_exprr
   module procedure mp_expzi
   module procedure mp_expzz
   module procedure mp_exprz
   module procedure mp_expzr
end interface

interface operator (==)
  module procedure mp_eqtrr
  module procedure mp_eqtdr
  module procedure mp_eqtrd
  module procedure mp_eqtir
  module procedure mp_eqtri
  module procedure mp_eqtzz
  module procedure mp_eqtdz
  module procedure mp_eqtzd
  module procedure mp_eqtdcz
  module procedure mp_eqtzdc
  module procedure mp_eqtrz
  module procedure mp_eqtzr
end interface

interface operator (/=)
  module procedure mp_netrr
  module procedure mp_netdr
  module procedure mp_netrd
  module procedure mp_netir
  module procedure mp_netri
  module procedure mp_netzz
  module procedure mp_netdz
  module procedure mp_netzd
  module procedure mp_netdcz
  module procedure mp_netzdc
  module procedure mp_netrz
  module procedure mp_netzr
end interface

interface operator (<=)
  module procedure mp_letrr
  module procedure mp_letdr
  module procedure mp_letrd
  module procedure mp_letir
  module procedure mp_letri
end interface

interface operator (>=)
  module procedure mp_getrr
  module procedure mp_getdr
  module procedure mp_getrd
  module procedure mp_getir
  module procedure mp_getri
end interface

interface operator (<)
  module procedure mp_lttrr
  module procedure mp_lttdr
  module procedure mp_lttrd
  module procedure mp_lttir
  module procedure mp_lttri
end interface

interface operator (>)
  module procedure mp_gttrr
  module procedure mp_gttdr
  module procedure mp_gttrd
  module procedure mp_gttir
  module procedure mp_gttri
end interface

!  MP generic function interface blocks, listed alphabetically by interface name:

interface abs
  module procedure mp_absr
  module procedure mp_absz
end interface

interface acos
  module procedure mp_acos
end interface

interface acosh
  module procedure mp_acosh
end interface

interface agm
  module procedure mp_agm
end interface

interface aimag
  module procedure mp_aimag
end interface

interface aint
  module procedure mp_aint
end interface

interface airy
  module procedure mp_airy
end interface

interface anint
  module procedure mp_anint
end interface

interface asin
  module procedure mp_asin
end interface

interface asinh
  module procedure mp_asinh
end interface

interface atan
  module procedure mp_atan
end interface

interface atanh
  module procedure mp_atanh
end interface

interface atan2
  module procedure mp_atan2
end interface

interface mpberne
  module procedure mp_berne
end interface

interface bessel_i
  module procedure mp_bessel_i
end interface

interface bessel_in
  module procedure mp_bessel_in
end interface

interface bessel_j
  module procedure mp_bessel_j
end interface

interface bessel_j0
  module procedure mp_bessel_j0
end interface

interface bessel_j1
  module procedure mp_bessel_j1
end interface

interface bessel_jn
  module procedure mp_bessel_jn
end interface

interface bessel_k
  module procedure mp_bessel_k
end interface

interface bessel_kn
  module procedure mp_bessel_kn
end interface

interface bessel_y
  module procedure mp_bessel_y
end interface

interface bessel_y0
  module procedure mp_bessel_y0
end interface

interface bessel_y1
  module procedure mp_bessel_y1
end interface

interface bessel_yn
  module procedure mp_bessel_yn
end interface

interface conjg
  module procedure mp_conjg
end interface

interface cos
  module procedure mp_cos
  module procedure mp_ccos
end interface

interface cosh
  module procedure mp_cosh
end interface

interface dble
  module procedure mp_rtod
end interface

interface dcmplx
  module procedure mp_ztodc
end interface

interface digamma
  module procedure mp_digamma
end interface

interface digamma_be
  module procedure mp_digamma_be
end interface

interface erf
  module procedure mp_erf
end interface

interface erfc
  module procedure mp_erfc
end interface

interface exp
  module procedure mp_exp
  module procedure mp_cexp
end interface

interface expint
  module procedure mp_expint
end interface

interface gamma
  module procedure mp_gamma
end interface

interface hurwitz_zetan
  module procedure mp_hurwitz_zetan
end interface

interface hurwitz_zetan_be
  module procedure mp_hurwitz_zetan_be
end interface

interface hypergeom_pfq
  module procedure mp_hypergeom_pfq
end interface

interface hypot
  module procedure mp_hypot
end interface

interface incgamma
  module procedure mp_incgamma
end interface

interface log
  module procedure mp_log
  module procedure mp_clog
end interface

interface log10
  module procedure mp_log10
end interface

interface log_gamma
  module procedure mp_log_gamma
end interface

interface max
  module procedure mp_max
end interface

interface min
  module procedure mp_min
end interface

interface mod
  module procedure mp_mod
end interface

interface mpbinmd
  module procedure mp_binmd
end interface

interface mpcmplx
  module procedure mp_dctoz
  module procedure mp_rtoz
  module procedure mp_ztoz
  module procedure mp_zmtoz
end interface

interface mpcmplxdc
  module procedure mp_dctoz2
end interface

interface mpcssh
  module procedure mp_cssh
end interface

interface mpcssn
  module procedure mp_cssn
end interface

interface mpdecmd
  module procedure mp_decmd
end interface

interface mpeform
  module procedure mp_eform
end interface

interface mpegamma
  module procedure mp_egamma
end interface

interface mpfform
  module procedure mp_fform
end interface

interface mpinit
  module procedure mp_init
end interface

interface mplog2
  module procedure mp_log2
end interface

interface mpnrt
  module procedure mp_nrt
end interface

interface mppi
  module procedure mp_pi
end interface

interface mpprod
  module procedure mp_prodd
  module procedure mp_prodq
end interface

interface mpquot
  module procedure mp_quotd
  module procedure mp_quotq
end interface

interface mprand
  module procedure mp_rand
end interface

interface mpread
  module procedure mp_readr1
  module procedure mp_readr2
  module procedure mp_readr3
  module procedure mp_readr4
  module procedure mp_readr5
  module procedure mp_readz1
  module procedure mp_readz2
  module procedure mp_readz3
  module procedure mp_readz4
  module procedure mp_readz5
end interface

interface mpreal
  module procedure mp_ator1
  module procedure mp_atorn
  module procedure mp_dtor
  module procedure mp_rtor
  module procedure mp_ztor
  module procedure mp_mtor
  module procedure mp_qtor
end interface

interface mpreald
  module procedure mp_dtor2
end interface

interface mprealq
  module procedure mp_qtor2
end interface

interface mpwprec
  module procedure mp_wprec
  module procedure mp_wprecz
end interface

interface mpwrite
  module procedure mp_writer
  module procedure mp_writez
end interface

interface polygamma
  module procedure mp_polygamma
end interface

interface polygamma_be
  module procedure mp_polygamma_be
end interface

interface polylog_ini
  module procedure mp_polylog_ini
end interface

interface polylog_neg
  module procedure mp_polylog_neg
end interface

interface polylog_pos
  module procedure mp_polylog_pos
end interface

interface qreal
  module procedure mp_rtoq
end interface

interface sign
  module procedure mp_sign
end interface

interface sin
  module procedure mp_sin
  module procedure mp_csin
end interface

interface sinh
  module procedure mp_sinh
end interface

interface sqrt
  module procedure mp_sqrt
  module procedure mp_csqrt
end interface

interface struve_hn
  module procedure mp_struve_hn
end interface

interface tan
  module procedure mp_tan
end interface

interface tanh
  module procedure mp_tanh
end interface

interface zeta
  module procedure mp_zeta
end interface

interface zeta_be
  module procedure mp_zeta_be
end interface

interface zeta_int
  module procedure mp_zeta_int
end interface

contains

!   This routine terminates execution.  Users may wish to replace the
!   default STOP with a call to a system routine that provides a traceback.

  subroutine mp_abrt (ier)
    implicit none
    integer ier
    write (mpldb, 1) ier
1   format ('*** MP_ABRT: Execution terminated, error code =',i4)
    stop
  end subroutine

!  This routine outputs an error message if iprec exceeds mpwds.

  function mp_setwp (iprec)
    integer mp_setwp
    integer, intent (in):: iprec
    if (iprec > mpwds) then
      write (mpldb, 1)
1       format ( &
        '*** MP_SETWP: requested precision level exceeds default precision.'/ &
        'Increase default precision in module MPFUNF.')
      call mp_abrt (98)
    endif
    mp_setwp = iprec
  end function

!  This routine checks if the input double precision variable has more than
!  40 significant bits; if so, an error message is output and mp_abrt is called.

  subroutine mp_checkdp (da)
    real (mprknd), intent (in):: da
    real (mprknd) d1, d2
    d1 = mpb13x * abs (da)
    d2 = abs (abs (da) + d1) - abs (d1)
    if (d2 /= abs (da)) then
      write (mpldb, 1) da
1     format ('*** MP_CHECKDP: DP value has more than 40 significant bits:', &
      1p,d25.15/'and thus very likely represents an unintended loss of accuracy.'/ &
      'Fix the issue, or else use functions mpprod, mpquot, mpreald or mpcmplxdc.'/ &
      'See documentation for details.')
      call mp_abrt (82)
    endif
  end subroutine

!  This routine checks if the input quad precision variable has more than
!  90 significant bits; if so, an error message is output and mp_abrt is called.

  subroutine mp_checkqp (qa)
    real (max (mprknd2, kind (1.0))), intent (in):: qa
    real (mprknd) d1, d2
    d1 = mpb23x * abs (qa)
    d2 = abs (abs (qa) + d1) - abs (d1)
    if (d2 /= abs (qa)) then
      write (mpldb, 1) qa
1     format ('*** MP_CHECKQP: QP value has more than 90 significant bits:'/ &
      1p,d50.35/'and thus very likely represents an unintended loss of accuracy.'/ &
      'Fix the issue, or else use functions mpprod, mpquot, mprealq or mprealqm.'/ &
      'See documentation for details.')
      call mp_abrt (82)
    endif
  end subroutine

!  These two routines are used to initialize scratch and output MP variables
!  in the routines below. The working precision level of the variable to nbt bits,
!  and the value is set to the "NAN" of MPFR.

  subroutine mp_initvr (ra, nbt)
    implicit none
    type (mp_real), intent (out):: ra
    integer, intent (in):: nbt
    ra%mpr(0) = mpwds6
    ra%mpr(1) = nbt
    ra%mpr(2) = 1
    ra%mpr(3) = mpnan
    ra%mpr(4) = loc (ra%mpr(4)) + 8
    ra%mpr(nbt/mpnbt+5) = 0
    return
  end subroutine

  subroutine mp_initvz (za, nbt)
    implicit none
    type (mp_complex), intent (out):: za
    integer, intent (in):: nbt
    integer l1
    l1 = mpwds6
    za%mpc(0) = mpwds6
    za%mpc(1) = nbt
    za%mpc(2) = 1
    za%mpc(3) = mpnan
    za%mpc(4) = loc (za%mpc(4)) + 8
    za%mpc(nbt/mpnbt+5) = 0
    za%mpc(l1) = mpwds6
    za%mpc(l1+1) = nbt
    za%mpc(l1+2) = 1
    za%mpc(l1+3) = mpnan
    za%mpc(l1+4) = loc (za%mpc(l1+4)) + 8
    za%mpc(l1+nbt/mpnbt+5) = 0
    return
  end subroutine

!  The next four subroutines are needed for most of the routines below, because
!  temporary multiprecision variables generated by Fortran compiler merely
!  copy the array, but do not correct the pointer in index 4 of the array.
!  Furthermore, calling one of these two subroutines, rather than simply
!  including the line or two of code below directly in the routine, avoids
!  error messages resulting from a conflict with the "intent (in)" attribute.

  subroutine mp_fixlocr (ra)
    implicit none
    type (mp_real):: ra
    ra%mpr(4) = loc (ra%mpr(4)) + 8
    return
  end subroutine

  subroutine mp_fixlocrm (ra)
    implicit none
    type (mp_realm):: ra
    ra%mpr(4) = loc (ra%mpr(4)) + 8
    return
  end subroutine

  subroutine mp_fixlocz (za)
    implicit none
    integer l1
    type (mp_complex):: za
    l1 = za%mpc(0)
    za%mpc(4) = loc (za%mpc(4)) + 8
    za%mpc(l1+4) = loc (za%mpc(l1+4)) + 8
    return
  end subroutine

  subroutine mp_fixloczm (za)
    implicit none
    integer l1
    type (mp_complexm):: za
    l1 = za%mpc(0)
    za%mpc(4) = loc (za%mpc(4)) + 8
    za%mpc(l1+4) = loc (za%mpc(l1+4)) + 8
    return
  end subroutine

!  Assignment routines:

  subroutine mp_eqrr (ra, rb)
    implicit none
    type (mp_real), intent (out):: ra
    type (mp_real), intent (in):: rb
    integer mpnwbt
    call mp_fixlocr (rb)
    mpnwbt = min (rb%mpr(1), mpwdsbt)
    call mp_initvr (ra, mpnwbt)
    call mpfrset (ra%mpr(1:), rb%mpr(1:), mprnd)
    return
  end subroutine

  subroutine mp_eqdr (da, rb)
    implicit none
    real (mprknd), intent (out):: da
    type (mp_real), intent (in):: rb
    real (mprknd) mpfrgetd
    external mpfrgetd
    call mp_fixlocr (rb)
    da = mpfrgetd (rb%mpr(1:), mprnd)
    return
  end subroutine

  subroutine mp_eqrd (ra, db)
    implicit none
    type (mp_real), intent (out):: ra
    real (mprknd), intent (in):: db
    integer mpnwbt
    mpnwbt = mpwdsbt
    call mp_checkdp (db)
    call mp_initvr (ra, mpnwbt)
    call mpfrsetd (ra%mpr(1:), db, mprnd)
    return
  end subroutine

  subroutine mp_eqir (ia, rb)
    implicit none
    integer, intent (out):: ia
    type (mp_real), intent (in):: rb
    real (mprknd) da
    real (mprknd) mpfrgetd
    external mpfrgetd
    call mp_fixlocr (rb)
    da = mpfrgetd (rb%mpr(1:), mprnd)
    ia = da
    return
  end subroutine

  subroutine mp_eqri (ra, ib)
    implicit none
    type (mp_real), intent (out):: ra
    integer, intent (in):: ib
    real (mprknd) db
    integer mpnwbt
    mpnwbt = mpwdsbt
    db = ib
    call mp_checkdp (db)
    call mp_initvr (ra, mpnwbt)
    call mpfrsetd (ra%mpr(1:), db, mprnd)
    return
  end subroutine

  subroutine mp_eqra (ra, ab)
    implicit none
    type (mp_real), intent (out):: ra
    character(*), intent (in):: ab
    character(1) :: chr1(len(ab)+1)
    integer i, l1, mpnw
    integer mpnwbt
    mpnwbt = mpwdsbt
    call mp_initvr (ra, mpnwbt)
    l1 = len (ab)
    do i = 1, l1
      if (ab(i:i) == 'D' .or. ab(i:i) == 'd') then
        chr1(i) = 'e'
      else
        chr1(i) = ab(i:i)
      endif
    enddo
    chr1(l1+1) = char(0)
    call mpfrsetstr (ra%mpr(1:), chr1, mprnd)
    return
  end subroutine

  subroutine mp_eqzz (za, zb)
    implicit none
    type (mp_complex), intent (out):: za
    type (mp_complex), intent (in):: zb
    integer l1, l2, mpnwbt
    l1 = zb%mpc(0)
    call mp_fixlocz (zb)
    mpnwbt = max (zb%mpc(1), zb%mpc(l1+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    l2 = mpwds6
    call mp_initvz (za, mpnwbt)
    call mpfrset (za%mpc(1:), zb%mpc(1:), mprnd)
    call mpfrset (za%mpc(l2+1:), zb%mpc(l1+1:), mprnd)
    return
  end subroutine

  subroutine mp_eqdz (da, zb)
    implicit none
    real (mprknd), intent (out):: da
    type (mp_complex), intent (in):: zb
    real (mprknd) mpfrgetd
    call mp_fixlocz (zb)
    da = mpfrgetd (zb%mpc(1:), mprnd)
    return
  end subroutine

  subroutine mp_eqzd (za, db)
    implicit none
    type (mp_complex), intent (out):: za
    real (mprknd), intent (in):: db
    real (mprknd) d1
    integer l2, mpnwbt
    mpnwbt = mpwdsbt
    l2 = mpwds6
    d1 = 0.d0
    call mp_initvz (za, mpnwbt)
    call mp_checkdp (db)
    call mpfrsetd (za%mpc(1:), db, mprnd)
    call mpfrsetd (za%mpc(l2+1:), d1, mprnd)
    return
  end subroutine

  subroutine mp_eqdcz (dca, zb)
    implicit none
    complex (kind (0.d0)), intent (out):: dca
    type (mp_complex), intent (in):: zb
    integer l1, n1, n2
    real (mprknd) d1, d2
    real (mprknd) mpfrgetd
    l1 = zb%mpc(0)
    call mp_fixlocz (zb)
    d1 = mpfrgetd (zb%mpc(1:), mprnd)
    d2 = mpfrgetd (zb%mpc(l1+1:), mprnd)
    dca = cmplx (d1, d2, mprknd)
    return
  end subroutine

  subroutine mp_eqzdc (za, dcb)
    implicit none
    type (mp_complex), intent (out):: za
    complex (kind (0.d0)), intent (in):: dcb
    integer l2, mpnwbt
    mpnwbt = mpwdsbt
    l2 = mpwds6
    call mp_initvz (za, mpnwbt)
    call mp_checkdp (dble (dcb))
    call mp_checkdp (aimag (dcb))
    call mpfrsetd (za%mpc(1:), dble (dcb), mprnd)
    call mpfrsetd (za%mpc(l2+1:), aimag (dcb), mprnd)
    return
  end subroutine

  subroutine mp_eqrz (ra, zb)
    implicit none
    type (mp_real), intent (out):: ra
    type (mp_complex), intent (in):: zb
    integer l1, mpnwbt
    l1 = zb%mpc(0)
    call mp_fixlocz (zb)
    mpnwbt = max (zb%mpc(1), zb%mpc(l1+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    call mp_initvr (ra, mpnwbt)
    call mpfrset (ra%mpr(1:), zb%mpc(1:), mprnd)
    return
  end subroutine

  subroutine mp_eqzr (za, rb)
    implicit none
    type (mp_complex), intent (out):: za
    type (mp_real), intent (in):: rb
    type (mp_real) r1
    integer l2, mpnwbt
    call mp_fixlocr (rb)
    mpnwbt = min (rb%mpr(1), mpwdsbt)
    l2 = mpwds6
    call mp_initvz (za, mpnwbt)
    call mpfrset (za%mpc(1:), rb%mpr(1:), mprnd)
    call mpfrsetd (za%mpc(l2+1:), 0.d0, mprnd)
    return
  end subroutine

!  Addition routines:

  function mp_addrr (ra, rb)
    implicit none
    type (mp_real):: mp_addrr
    type (mp_real), intent (in):: ra, rb
    integer mpnwbt
    call mp_fixlocr (ra)
    call mp_fixlocr (rb)
    mpnwbt = max (ra%mpr(1), rb%mpr(1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    call mp_initvr (mp_addrr, mpnwbt)
    call mpfradd (mp_addrr%mpr(1:), ra%mpr(1:), rb%mpr(1:), mprnd)
    return
  end function

  function mp_adddr (da, rb)
    implicit none
    type (mp_real):: mp_adddr
    real (mprknd), intent (in):: da
    type (mp_real), intent (in):: rb
    integer mpnwbt
    mpnwbt = min (rb%mpr(1), mpwdsbt)
    call mp_checkdp (da)
    call mp_fixlocr (rb)
    call mp_initvr (mp_adddr, mpnwbt)
    call mpfraddd (mp_adddr%mpr(1:), rb%mpr(1:), da, mprnd)
    return
  end function

  function mp_addrd (ra, db)
    implicit none
    type (mp_real):: mp_addrd
    type (mp_real), intent (in):: ra
    real (mprknd), intent (in):: db
    integer mpnwbt
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_checkdp (db)
    call mp_fixlocr (ra)
    call mp_initvr (mp_addrd, mpnwbt)
    call mpfraddd (mp_addrd%mpr(1:), ra%mpr(1:), db, mprnd)
    return
  end function

  function mp_addir (ia, rb)
    implicit none
    type (mp_real):: mp_addir
    integer, intent (in):: ia
    type (mp_real), intent (in):: rb
    real (mprknd) da
    integer mpnwbt
    mpnwbt = min (rb%mpr(1), mpwdsbt)
    da = ia
    call mp_checkdp (da)
    call mp_fixlocr (rb)
    call mp_initvr (mp_addir, mpnwbt)
    call mpfraddd (mp_addir%mpr(1:), rb%mpr(1:), da, mprnd)
    return
  end function

  function mp_addri (ra, ib)
    implicit none
    type (mp_real):: mp_addri
    type (mp_real), intent (in):: ra
    integer, intent (in):: ib
    real (mprknd) db
    integer mpnwbt
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    db = ib
    call mp_checkdp (db)
    call mp_fixlocr (ra)
    call mp_initvr (mp_addri, mpnwbt)
    call mpfraddd (mp_addri%mpr(1:), ra%mpr(1:), db, mprnd)
    return
  end function

  function mp_addzz (za, zb)
    implicit none
    type (mp_complex):: mp_addzz
    type (mp_complex), intent (in):: za, zb
    integer l1, l2, l3, mpnwbt
    l1 = za%mpc(0)
    call mp_fixlocz (za)
    l2 = zb%mpc(0)
    call mp_fixlocz (zb)
    mpnwbt = max (za%mpc(1), za%mpc(l1+1), zb%mpc(1), zb%mpc(l2+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    l3 = mpwds6
    call mp_initvz (mp_addzz, mpnwbt)
    call mpfradd (mp_addzz%mpc(1:), za%mpc(1:), zb%mpc(1:), mprnd)
    call mpfradd (mp_addzz%mpc(l3+1:), za%mpc(l1+1:), zb%mpc(l2+1:), mprnd)
    return
  end function

  function mp_adddz (da, zb)
    implicit none
    type (mp_complex):: mp_adddz
    real (mprknd), intent (in):: da
    type (mp_complex), intent (in):: zb
    integer l1, l2, l3, mpnwbt
    l2 = zb%mpc(0)
    call mp_checkdp (da)
    call mp_fixlocz (zb)
    mpnwbt = max (zb%mpc(1), zb%mpc(l2+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    l3 = mpwds6
    call mp_initvz (mp_adddz, mpnwbt)
    call mpfraddd (mp_adddz%mpc(1:), zb%mpc(1:), da, mprnd)
    call mpfrset (mp_adddz%mpc(l3+1:), zb%mpc(l2+1:), mprnd)
    return
  end function

  function mp_addzd (za, db)
    implicit none
    type (mp_complex):: mp_addzd
    type (mp_complex), intent (in):: za
    real (mprknd), intent (in):: db
    integer l1, l2, l3, mpnwbt
    l2 = za%mpc(0)
    call mp_checkdp (db)
    call mp_fixlocz (za)
    mpnwbt = max (za%mpc(1), za%mpc(l2+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    l3 = mpwds6
    call mp_initvz (mp_addzd, mpnwbt)
    call mpfraddd (mp_addzd%mpc(1:), za%mpc(1:), db, mprnd)
    call mpfrset (mp_addzd%mpc(l3+1:), za%mpc(l2+1:), mprnd)
    return
  end function

  function mp_adddcz (dca, zb)
    implicit none
    type (mp_complex):: mp_adddcz
    complex (kind (0.d0)), intent (in):: dca
    type (mp_complex), intent (in):: zb
    integer l1, l2, l3, mpnwbt
    l2 = zb%mpc(0)
    call mp_checkdp (dble (dca))
    call mp_checkdp (aimag (dca))
    call mp_fixlocz (zb)
    mpnwbt = max (zb%mpc(1), zb%mpc(l2+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    l3 = mpwds6
    call mp_initvz (mp_adddcz, mpnwbt)
    call mpfraddd (mp_adddcz%mpc(1:), zb%mpc(1:), dble (dca), mprnd)
    call mpfraddd (mp_adddcz%mpc(l3+1:), zb%mpc(l2+1:), aimag (dca), mprnd)
    return
  end function

  function mp_addzdc (za, dcb)
    implicit none
    type (mp_complex):: mp_addzdc
    type (mp_complex), intent (in):: za
    complex (kind (0.d0)), intent (in):: dcb
    integer l1, l2, l3, mpnwbt
    l2 = za%mpc(0)
    call mp_checkdp (dble (dcb))
    call mp_checkdp (aimag (dcb))
    call mp_fixlocz (za)
    mpnwbt = max (za%mpc(1), za%mpc(l2+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    l3 = mpwds6
    call mp_initvz (mp_addzdc, mpnwbt)
    call mpfraddd (mp_addzdc%mpc(1:), za%mpc(1:), dble (dcb), mprnd)
    call mpfraddd (mp_addzdc%mpc(l3+1:), za%mpc(l2+1:), aimag (dcb), mprnd)
    return
  end function

  function mp_addrz (ra, zb)
    implicit none
    type (mp_complex):: mp_addrz
    type (mp_real), intent (in):: ra
    type (mp_complex), intent (in):: zb
    integer l1, l2, l3, mpnwbt
    call mp_fixlocr (ra)
    l2 = zb%mpc(0)
    call mp_fixlocz (zb)
    mpnwbt = max (ra%mpr(1), zb%mpc(1), zb%mpc(l2+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    l3 = mpwds6
    call mp_initvz (mp_addrz, mpnwbt)
    call mpfradd (mp_addrz%mpc(1:), ra%mpr(1:), zb%mpc(1:), mprnd)
    call mpfrset (mp_addrz%mpc(l3+1:), zb%mpc(l2+1:), mprnd)
    return
  end function

  function mp_addzr (za, rb)
    implicit none
    type (mp_complex):: mp_addzr
    type (mp_complex), intent (in):: za
    type (mp_real), intent (in):: rb
    integer l1, l2, l3, mpnwbt
    l1 = za%mpc(0)
    call mp_fixlocz (za)
    call mp_fixlocr (rb)
    mpnwbt = max (za%mpc(1), za%mpc(l1+1), rb%mpr(1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    l3 = mpwds6
    call mp_initvz (mp_addzr, mpnwbt)
    call mpfradd (mp_addzr%mpc(1:), za%mpc(1:), rb%mpr(1:), mprnd)
    call mpfrset (mp_addzr%mpc(l3+1:), za%mpc(l1+1:), mprnd)
    return
  end function

!  Subtraction routines:

  function mp_subrr (ra, rb)
    implicit none
    type (mp_real):: mp_subrr
    type (mp_real), intent (in):: ra, rb
    integer mpnwbt
    call mp_fixlocr (ra)
    call mp_fixlocr (rb)
    mpnwbt = max (ra%mpr(1), rb%mpr(1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    call mp_initvr (mp_subrr, mpnwbt)
    call mpfrsub (mp_subrr%mpr(1:), ra%mpr(1:), rb%mpr(1:), mprnd)
    return
  end function

  function mp_subdr (da, rb)
    implicit none
    type (mp_real):: mp_subdr
    real (mprknd), intent (in):: da
    type (mp_real), intent (in):: rb
    integer mpnwbt
    call mp_checkdp (da)
    call mp_fixlocr (rb)
    mpnwbt = min (rb%mpr(1), mpwdsbt)
    call mp_initvr (mp_subdr, mpnwbt)
    call mpfrdsub (mp_subdr%mpr(1:), da, rb%mpr(1:), mprnd)
    return
  end function

  function mp_subrd (ra, db)
    implicit none
    type (mp_real):: mp_subrd
    type (mp_real), intent (in):: ra
    real (mprknd), intent (in):: db
    integer mpnwbt
    call mp_checkdp (db)
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_subrd, mpnwbt)
    call mpfrsubd (mp_subrd%mpr(1:), ra%mpr(1:), db, mprnd)
    return
  end function

  function mp_subir (ia, rb)
    implicit none
    type (mp_real):: mp_subir
    integer, intent (in):: ia
    type (mp_real), intent (in):: rb
    real (mprknd) da
    integer mpnwbt
    da = ia
    call mp_checkdp (da)
    call mp_fixlocr (rb)
    mpnwbt = min (rb%mpr(1), mpwdsbt)
    call mp_initvr (mp_subir, mpnwbt)
    call mpfrdsub (mp_subir%mpr(1:), da, rb%mpr(1:), mprnd)
    return
  end function

  function mp_subri (ra, ib)
    implicit none
    type (mp_real):: mp_subri
    type (mp_real), intent (in):: ra
    integer, intent (in):: ib
    real (mprknd) db
    integer mpnwbt
    db = ib
    call mp_checkdp (db)
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_subri, mpnwbt)
    call mpfrsubd (mp_subri%mpr(1:), ra%mpr(1:), db, mprnd)
    return
  end function

  function mp_subzz (za, zb)
    implicit none
    type (mp_complex):: mp_subzz
    type (mp_complex), intent (in):: za, zb
    integer l1, l2, l3, mpnwbt
    l1 = za%mpc(0)
    call mp_fixlocz (za)
    l2 = zb%mpc(0)
    call mp_fixlocz (zb)
    mpnwbt = max (za%mpc(1), za%mpc(l1+1), zb%mpc(1), zb%mpc(l2+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    l3 = mpwds6
    call mp_initvz (mp_subzz, mpnwbt)
    call mpfrsub (mp_subzz%mpc(1:), za%mpc(1:), zb%mpc(1:), mprnd)
    call mpfrsub (mp_subzz%mpc(l3+1:), za%mpc(l1+1:), zb%mpc(l2+1:), mprnd)
    return
  end function

  function mp_subdz (da, zb)
    implicit none
    type (mp_complex):: mp_subdz
    real (mprknd), intent (in):: da
    type (mp_complex), intent (in):: zb
    integer l1, l2, l3, mpnwbt
    l2 = zb%mpc(0)
    call mp_checkdp (da)
    call mp_fixlocz (zb)
    mpnwbt = max (zb%mpc(1), zb%mpc(l2+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    l3 = mpwds6
    call mp_initvz (mp_subdz, mpnwbt)
    call mpfrdsub (mp_subdz%mpc(1:), da, zb%mpc(1:), mprnd)
    call mpfrneg (mp_subdz%mpc(l3+1:), zb%mpc(l2+1:), mprnd)
    return
  end function

  function mp_subzd (za, db)
    implicit none
    type (mp_complex):: mp_subzd
    type (mp_complex), intent (in):: za
    real (mprknd), intent (in):: db
    integer l1, l2, l3, mpnwbt
    l2 = za%mpc(0)
    call mp_checkdp (db)
    call mp_fixlocz (za)
    mpnwbt = max (za%mpc(1), za%mpc(l2+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    l3 = mpwds6
    call mp_initvz (mp_subzd, mpnwbt)
    call mpfrsubd (mp_subzd%mpc(1:), za%mpc(1:), db, mprnd)
    call mpfrset (mp_subzd%mpc(l3+1:), za%mpc(l2+1:), mprnd)
    return
  end function

  function mp_subdcz (dca, zb)
    implicit none
    type (mp_complex):: mp_subdcz
    complex (kind (0.d0)), intent (in):: dca
    type (mp_complex), intent (in):: zb
    integer l1, l2, l3, mpnwbt
    l2 = zb%mpc(0)
    call mp_checkdp (dble (dca))
    call mp_checkdp (aimag (dca))
    call mp_fixlocz (zb)
    mpnwbt = max (zb%mpc(1), zb%mpc(l2+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    l3 = mpwds6
    call mp_initvz (mp_subdcz, mpnwbt)
    call mpfrdsub (mp_subdcz%mpc(1:), dble (dca), zb%mpc(1:), mprnd)
    call mpfrdsub (mp_subdcz%mpc(l3+1:), aimag (dca), zb%mpc(l2+1:), mprnd)
    return
  end function

  function mp_subzdc (za, dcb)
    implicit none
    type (mp_complex):: mp_subzdc
    type (mp_complex), intent (in):: za
    complex (kind (0.d0)), intent (in):: dcb
    integer l1, l2, l3, mpnwbt
    l2 = za%mpc(0)
    call mp_checkdp (dble (dcb))
    call mp_checkdp (aimag (dcb))
    call mp_fixlocz (za)
    mpnwbt = max (za%mpc(1), za%mpc(l2+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    l3 = mpwds6
    call mp_initvz (mp_subzdc, mpnwbt)
    call mpfrsubd (mp_subzdc%mpc(1:), za%mpc(1:), dble (dcb), mprnd)
    call mpfrsubd (mp_subzdc%mpc(l3+1:), za%mpc(l2+1:), aimag (dcb), mprnd)
    return
  end function

  function mp_subrz (ra, zb)
    implicit none
    type (mp_complex):: mp_subrz
    type (mp_real), intent (in):: ra
    type (mp_complex), intent (in):: zb
    integer l1, l2, l3, mpnwbt
    call mp_fixlocr (ra)
    l2 = zb%mpc(0)
    call mp_fixlocz (zb)
    mpnwbt = max (ra%mpr(1), zb%mpc(1), zb%mpc(l2+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    l3 = mpwds6
    call mp_initvz (mp_subrz, mpnwbt)
    call mpfrsub (mp_subrz%mpc(1:), ra%mpr(1:), zb%mpc(1:), mprnd)
    call mpfrset (mp_subrz%mpc(l3+1:), zb%mpc(l2+1:), mprnd)
    call mpfrneg (mp_subrz%mpc(l3+1:), mp_subrz%mpc(l3+1:), mprnd)
    return
  end function

  function mp_subzr (za, rb)
    implicit none
    type (mp_complex):: mp_subzr
    type (mp_complex), intent (in):: za
    type (mp_real), intent (in):: rb
    integer l1, l2, l3, mpnwbt
    l1 = za%mpc(0)
    call mp_fixlocz (za)
    call mp_fixlocr (rb)
    mpnwbt = max (za%mpc(1), za%mpc(l1+1), rb%mpr(1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    l3 = mpwds6
    call mp_initvz (mp_subzr, mpnwbt)
    call mpfrsub (mp_subzr%mpc(1:), za%mpc(1:), rb%mpr(1:), mprnd)
    call mpfrset (mp_subzr%mpc(l3+1:), za%mpc(l1+1:), mprnd)
    return
  end function

!  Negation routines:

  function mp_negr (ra)
    implicit none
    type (mp_real):: mp_negr
    type (mp_real), intent (in):: ra
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_negr, mpnwbt)
    call mpfrneg (mp_negr%mpr(1:), ra%mpr(1:), mprnd)
    return
  end function

  function mp_negz (za)
    implicit none
    type (mp_complex):: mp_negz
    type (mp_complex), intent (in):: za
    integer l1, l2, l3, mpnwbt
    l1 = za%mpc(0)
    call mp_fixlocz (za)
    mpnwbt = max (za%mpc(1), za%mpc(l1+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    l3 = mpwds6
    call mp_initvz (mp_negz, mpnwbt)
    call mpfrneg (mp_negz%mpc(1:), za%mpc(1:), mprnd)
    call mpfrneg (mp_negz%mpc(l3+1:), za%mpc(l1+1:), mprnd)
    return
  end function

!  Multiplication routines:

  function mp_mulrr (ra, rb)
    implicit none
    type (mp_real):: mp_mulrr
    type (mp_real), intent (in):: ra, rb
    integer mpnwbt
    call mp_fixlocr (ra)
    call mp_fixlocr (rb)
    mpnwbt = max (ra%mpr(1), rb%mpr(1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    call mp_initvr (mp_mulrr, mpnwbt)
    call mpfrmul (mp_mulrr%mpr(1:), ra%mpr(1:), rb%mpr(1:), mprnd)
    return
  end function

  function mp_muldr (da, rb)
    implicit none
    type (mp_real):: mp_muldr
    real (mprknd), intent (in):: da
    type (mp_real), intent (in):: rb
    integer mpnwbt
    call mp_checkdp (da)
    call mp_fixlocr (rb)
    mpnwbt = min (rb%mpr(1), mpwdsbt)
    call mp_initvr (mp_muldr, mpnwbt)
    call mpfrmuld (mp_muldr%mpr(1:), rb%mpr(1:), da, mprnd)
    return
  end function

  function mp_mulrd (ra, db)
    implicit none
    type (mp_real):: mp_mulrd
    type (mp_real), intent (in):: ra
    real (mprknd), intent (in):: db
    integer mpnwbt
    call mp_checkdp (db)
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_mulrd, mpnwbt)
    call mpfrmuld (mp_mulrd%mpr(1:), ra%mpr(1:), db, mprnd)
    return
  end function

  function mp_mulir (ia, rb)
    implicit none
    type (mp_real):: mp_mulir
    integer, intent (in):: ia
    type (mp_real), intent (in):: rb
    real (mprknd) da
    integer mpnwbt
    da = ia
    call mp_checkdp (da)
    call mp_fixlocr (rb)
    mpnwbt = min (rb%mpr(1), mpwdsbt)
    call mp_initvr (mp_mulir, mpnwbt)
    call mpfrmuld (mp_mulir%mpr(1:), rb%mpr(1:), da, mprnd)
    return
  end function

  function mp_mulri (ra, ib)
    implicit none
    type (mp_real):: mp_mulri
    type (mp_real), intent (in):: ra
    integer, intent (in):: ib
    real (mprknd) db
    integer mpnwbt
    db = ib
    call mp_checkdp (db)
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_mulri, mpnwbt)
    call mpfrmuld (mp_mulri%mpr(1:), ra%mpr(1:), db, mprnd)
    return
  end function

  function mp_mulzz (za, zb)
    implicit none
    type (mp_complex):: mp_mulzz
    type (mp_complex), intent (in):: za, zb
    type (mp_real) r1, r2
    integer l1, l2, l3, mpnwbt
    l1 = za%mpc(0)
    call mp_fixlocz (za)
    l2 = zb%mpc(0)
    call mp_fixlocz (zb)
    mpnwbt = max (za%mpc(1), za%mpc(l1+1), zb%mpc(1), zb%mpc(l2+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    call mp_initvr (r1, mpnwbt)
    call mp_initvr (r2, mpnwbt)
    l3 = mpwds6
    call mp_initvz (mp_mulzz, mpnwbt)
    call mpfrmul (r1%mpr(1:), za%mpc(1:), zb%mpc(1:), mprnd)
    call mpfrmul (r2%mpr(1:), za%mpc(l1+1:), zb%mpc(l2+1:), mprnd)
    call mpfrsub (mp_mulzz%mpc(1:), r1%mpr(1:), r2%mpr(1:), mprnd)
    call mpfrmul (r1%mpr(1:), za%mpc(1:), zb%mpc(l2+1:), mprnd)
    call mpfrmul (r2%mpr(1:), za%mpc(l1+1:), zb%mpc(1:), mprnd)
    call mpfradd (mp_mulzz%mpc(l3+1:), r1%mpr(1:), r2%mpr(1:), mprnd)
    return
  end function

  function mp_muldz (da, zb)
    implicit none
    type (mp_complex):: mp_muldz
    real (mprknd), intent (in):: da
    type (mp_complex), intent (in):: zb
    integer l1, l2, l3, mpnwbt
    l2 = zb%mpc(0)
    call mp_checkdp (da)
    call mp_fixlocz (zb)
    mpnwbt = max (zb%mpc(1), zb%mpc(l2+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    l3 = mpwds6
    call mp_initvz (mp_muldz, mpnwbt)
    call mpfrmuld (mp_muldz%mpc(1:), zb%mpc(1:), da, mprnd)
    call mpfrmuld (mp_muldz%mpc(l3+1:), zb%mpc(l2+1:), da, mprnd)
    return
  end function

  function mp_mulzd (za, db)
    implicit none
    type (mp_complex):: mp_mulzd
    type (mp_complex), intent (in):: za
    real (mprknd), intent (in):: db
    integer l1, l2, l3, mpnwbt
    l2 = za%mpc(0)
    call mp_checkdp (db)
    call mp_fixlocz (za)
    mpnwbt = max (za%mpc(1), za%mpc(l2+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    l3 = mpwds6
    call mp_initvz (mp_mulzd, mpnwbt)
    call mpfrmuld (mp_mulzd%mpc(1:), za%mpc(1:), db, mprnd)
    call mpfrmuld (mp_mulzd%mpc(l3+1:), za%mpc(l2+1:), db, mprnd)
    return
  end function

  function mp_muldcz (dca, zb)
    implicit none
    type (mp_complex):: mp_muldcz
    complex (kind (0.d0)), intent (in):: dca
    type (mp_complex), intent (in):: zb
    type (mp_real) r1, r2
    integer l1, l2, l3, mpnwbt
    l2 = zb%mpc(0)
    call mp_fixlocz (zb)
    mpnwbt = max (zb%mpc(1), zb%mpc(l2+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    call mp_checkdp (dble (dca))
    call mp_checkdp (aimag (dca))
    call mp_initvr (r1, mpnwbt)
    call mp_initvr (r2, mpnwbt)
    l3 = mpwds6
    call mp_initvz (mp_muldcz, mpnwbt)
    call mpfrmuld (r1%mpr(1:), zb%mpc(1:), dble (dca), mprnd)
    call mpfrmuld (r2%mpr(1:), zb%mpc(l2+1:), aimag (dca), mprnd)
    call mpfrsub (mp_muldcz%mpc(1:), r1%mpr(1:), r2%mpr(1:), mprnd)
    call mpfrmuld (r1%mpr(1:), zb%mpc(l2+1:), dble (dca), mprnd)
    call mpfrmuld (r2%mpr(1:), zb%mpc(1:), aimag (dca), mprnd)
    call mpfradd (mp_muldcz%mpc(l3+1:), r1%mpr(1:), r2%mpr(1:), mprnd)
    return
  end function

  function mp_mulzdc (za, dcb)
    implicit none
    type (mp_complex):: mp_mulzdc
    type (mp_complex), intent (in):: za
    complex (kind (0.d0)), intent (in):: dcb
    type (mp_real) r1, r2
    integer l1, l2, l3, mpnwbt
    l1 = za%mpc(0)
    call mp_fixlocz (za)
    mpnwbt = max (za%mpc(1), za%mpc(l1+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    call mp_checkdp (dble (dcb))
    call mp_checkdp (aimag (dcb))
    call mp_initvr (r1, mpnwbt)
    call mp_initvr (r2, mpnwbt)
    l3 = mpwds6
    call mp_initvz (mp_mulzdc, mpnwbt)
    call mpfrmuld (r1%mpr(1:), za%mpc(1:), dble (dcb), mprnd)
    call mpfrmuld (r2%mpr(1:), za%mpc(l1+1:), aimag (dcb), mprnd)
    call mpfrsub (mp_mulzdc%mpc(1:), r1%mpr(1:), r2%mpr(1:), mprnd)
    call mpfrmuld (r1%mpr(1:), za%mpc(l1+1:), dble (dcb), mprnd)
    call mpfrmuld (r2%mpr(1:), za%mpc(1:), aimag (dcb), mprnd)
    call mpfradd (mp_mulzdc%mpc(l3+1:), r1%mpr(1:), r2%mpr(1:), mprnd)
    return
  end function

  function mp_mulrz (ra, zb)
    implicit none
    type (mp_complex):: mp_mulrz
    type (mp_real), intent (in):: ra
    type (mp_complex), intent (in):: zb
    integer l1, l2, l3, mpnwbt
    l1 = ra%mpr(0)
    call mp_fixlocr (ra)
    l2 = zb%mpc(0)
    call mp_fixlocz (zb)
    mpnwbt = max (ra%mpr(1), zb%mpc(1), zb%mpc(l2+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    l3 = mpwds6
    call mp_initvz (mp_mulrz, mpnwbt)
    call mpfrmul (mp_mulrz%mpc(1:), ra%mpr(1:), zb%mpc(1:), mprnd)
    call mpfrmul (mp_mulrz%mpc(l3+1:), ra%mpr(1:), zb%mpc(l2+1:), mprnd)
    return
  end function

  function mp_mulzr (za, rb)
    implicit none
    type (mp_complex):: mp_mulzr
    type (mp_complex), intent (in):: za
    type (mp_real), intent (in):: rb
    integer l1, l2, l3, mpnwbt
    l1 = za%mpc(0)
    call mp_fixlocz (za)
    l2 = rb%mpr(0)
    call mp_fixlocr (rb)
    mpnwbt = max (za%mpc(1), za%mpc(l1+1), rb%mpr(1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    l3 = mpwds6
    call mp_initvz (mp_mulzr, mpnwbt)
    call mpfrmul (mp_mulzr%mpc(1:), za%mpc(1:), rb%mpr(1:), mprnd)
    call mpfrmul (mp_mulzr%mpc(l3+1:), za%mpc(l1+1:), rb%mpr(1:), mprnd)
    return
  end function

!  Division routines:

  function mp_divrr (ra, rb)
    implicit none
    type (mp_real):: mp_divrr
    type (mp_real), intent (in):: ra, rb
    integer mpnwbt
    call mp_fixlocr (ra)
    call mp_fixlocr (rb)
    mpnwbt = max (ra%mpr(1), rb%mpr(1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    call mp_initvr (mp_divrr, mpnwbt)
    call mpfrdiv (mp_divrr%mpr(1:), ra%mpr(1:), rb%mpr(1:), mprnd)
    return
  end function

  function mp_divdr (da, rb)
    implicit none
    type (mp_real):: mp_divdr
    real (mprknd), intent (in):: da
    type (mp_real), intent (in):: rb
    integer mpnwbt
    call mp_checkdp (da)
    call mp_fixlocr (rb)
    mpnwbt = min (rb%mpr(1), mpwdsbt)
    call mp_initvr (mp_divdr, mpnwbt)
    call mpfrddiv (mp_divdr%mpr(1:), da, rb%mpr(1:), mprnd)
    return
  end function

  function mp_divrd (ra, db)
    implicit none
    type (mp_real):: mp_divrd
    type (mp_real), intent (in):: ra
    real (mprknd), intent (in):: db
    integer mpnwbt
    call mp_checkdp (db)
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_divrd, mpnwbt)
    call mpfrdivd (mp_divrd%mpr(1:), ra%mpr(1:), db, mprnd)
    return
  end function

  function mp_divir (ia, rb)
    implicit none
    type (mp_real):: mp_divir
    integer, intent (in):: ia
    type (mp_real), intent (in):: rb
    real (mprknd) da
    integer mpnwbt
    da = ia
    call mp_checkdp (da)
    call mp_fixlocr (rb)
    mpnwbt = min (rb%mpr(1), mpwdsbt)
    call mp_initvr (mp_divir, mpnwbt)
    call mpfrddiv (mp_divir%mpr(1:), da, rb%mpr(1:), mprnd)
    return
  end function

  function mp_divri (ra, ib)
    implicit none
    type (mp_real):: mp_divri
    type (mp_real), intent (in):: ra
    integer, intent (in):: ib
    real (mprknd) db
    integer mpnwbt
    db = ib
    call mp_checkdp (db)
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_divri, mpnwbt)
    call mpfrdivd (mp_divri%mpr(1:), ra%mpr(1:), db, mprnd)
    return
  end function

  function mp_divzz (za, zb)
    implicit none
    type (mp_complex):: mp_divzz
    type (mp_complex), intent (in):: za, zb
    type (mp_real) r1, r2, r3
    integer l1, l2, l3, mpnwbt
    l1 = za%mpc(0)
    call mp_fixlocz (za)
    l2 = zb%mpc(0)
    call mp_fixlocz (zb)
    mpnwbt = max (za%mpc(1), za%mpc(l1+1), zb%mpc(1), zb%mpc(l2+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    call mp_initvr (r1, mpnwbt)
    call mp_initvr (r2, mpnwbt)
    call mp_initvr (r3, mpnwbt)
    l3 = mpwds6
    call mp_initvz (mp_divzz, mpnwbt)
    call mpfrmul (r1%mpr(1:), za%mpc(1:), zb%mpc(1:), mprnd)
    call mpfrmul (r2%mpr(1:), za%mpc(l1+1:), zb%mpc(l2+1:), mprnd)
    call mpfradd (mp_divzz%mpc(1:), r1%mpr(1:), r2%mpr(1:), mprnd)
    call mpfrmul (r1%mpr(1:), za%mpc(1:), zb%mpc(l2+1:), mprnd)
    call mpfrmul (r2%mpr(1:), za%mpc(l1+1:), zb%mpc(1:), mprnd)
    call mpfrsub (mp_divzz%mpc(l3+1:), r2%mpr(1:), r1%mpr(1:), mprnd)
    call mpfrmul (r1%mpr(1:), zb%mpc(1:), zb%mpc(1:), mprnd)
    call mpfrmul (r2%mpr(1:), zb%mpc(l2+1:), zb%mpc(l2+1:), mprnd)
    call mpfradd (r3%mpr(1:), r1%mpr(1:), r2%mpr(1:), mprnd)
    call mpfrdiv (mp_divzz%mpc(1:), mp_divzz%mpc(1:), r3%mpr(1:), mprnd)
    call mpfrdiv (mp_divzz%mpc(l3+1:), mp_divzz%mpc(l3+1:), r3%mpr(1:), mprnd)
    return
  end function

  function mp_divdz (da, zb)
    implicit none
    type (mp_complex):: mp_divdz
    real (mprknd), intent (in):: da
    type (mp_complex), intent (in):: zb
    real (mprknd) d1, d2
    type (mp_real) r1, r2, r3, r4, r5
    integer l1, l2, l3, mpnwbt
    l2 = zb%mpc(0)
    call mp_checkdp (da)
    call mp_fixlocz (zb)
    mpnwbt = max (zb%mpc(1), zb%mpc(l2+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    l3 = mpwds6
    call mp_initvr (r1, mpnwbt)
    call mp_initvr (r2, mpnwbt)
    call mp_initvr (r3, mpnwbt)
    call mp_initvr (r4, mpnwbt)
    call mp_initvr (r5, mpnwbt)
    call mp_initvz (mp_divdz, mpnwbt)
    call mpfrmul (r1%mpr(1:), zb%mpc(1:), zb%mpc(1:), mprnd)
    call mpfrmul (r2%mpr(1:), zb%mpc(l2+1:), zb%mpc(l2+1:), mprnd)
    call mpfradd (r3%mpr(1:), r1%mpr(1:), r2%mpr(1:), mprnd)
    call mpfrmuld (r4%mpr(1:), zb%mpc(1:), da, mprnd)
    call mpfrmuld (r5%mpr(1:), zb%mpc(l2+1:), da, mprnd)
    call mpfrdiv (mp_divdz%mpc(1:), r4%mpr(1:), r3%mpr(1:), mprnd)
    call mpfrdiv (mp_divdz%mpc(l3+1:), r5%mpr(1:), r3%mpr(1:), mprnd)
    call mpfrneg (mp_divdz%mpc(l3+1:), mp_divdz%mpc(l3+1:), mprnd)
    return
  end function

  function mp_divzd (za, db)
    implicit none
    type (mp_complex):: mp_divzd
    type (mp_complex), intent (in):: za
    real (mprknd), intent (in):: db
    integer l1, l2, l3, mpnwbt
    l2 = za%mpc(0)
    call mp_checkdp (db)
    call mp_fixlocz (za)
    mpnwbt = max (za%mpc(1), za%mpc(l2+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    l3 = mpwds6
    call mp_initvz (mp_divzd, mpnwbt)
    call mpfrdivd (mp_divzd%mpc(1:), za%mpc(1:), db, mprnd)
    call mpfrdivd (mp_divzd%mpc(l3+1:), za%mpc(l2+1:), db, mprnd)
    return
  end function

  function mp_divdcz (dca, zb)
    implicit none
    type (mp_complex):: mp_divdcz
    complex (kind (0.d0)), intent (in):: dca
    type (mp_complex), intent (in):: zb
    type (mp_real) r1, r2, r3
    integer l1, l2, l3, mpnwbt
    l2 = zb%mpc(0)
    call mp_fixlocz (zb)
    mpnwbt = max (zb%mpc(1), zb%mpc(l2+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    call mp_checkdp (dble (dca))
    call mp_checkdp (aimag (dca))
    call mp_initvr (r1, mpnwbt)
    call mp_initvr (r2, mpnwbt)
    call mp_initvr (r3, mpnwbt)
    l3 = mpwds6
    call mp_initvz (mp_divdcz, mpnwbt)
    call mpfrmuld (r1%mpr(1:), zb%mpc(1:), dble (dca), mprnd)
    call mpfrmuld (r2%mpr(1:), zb%mpc(l2+1:), aimag (dca), mprnd)
    call mpfradd (mp_divdcz%mpc(1:), r1%mpr(1:), r2%mpr(1:), mprnd)
    call mpfrmuld (r1%mpr(1:), zb%mpc(l2+1:), dble (dca), mprnd)
    call mpfrmuld (r2%mpr(1:), zb%mpc(1:), aimag (dca), mprnd)
    call mpfrsub (mp_divdcz%mpc(l3+1:), r2%mpr(1:), r1%mpr(1:), mprnd)
    call mpfrmul (r1%mpr(1:), zb%mpc(1:), zb%mpc(1:), mprnd)
    call mpfrmul (r2%mpr(1:), zb%mpc(l2+1:), zb%mpc(l2+1:), mprnd)
    call mpfradd (r3%mpr(1:), r1%mpr(1:), r2%mpr(1:), mprnd)
    call mpfrdiv (mp_divdcz%mpc(1:), mp_divdcz%mpc(1:), r3%mpr(1:), mprnd)
    call mpfrdiv (mp_divdcz%mpc(l3+1:), mp_divdcz%mpc(l3+1:), r3%mpr(1:), mprnd)
    return
  end function

  function mp_divzdc (za, dcb)
    implicit none
    type (mp_complex):: mp_divzdc
    type (mp_complex), intent (in):: za
    complex (kind (0.d0)), intent (in):: dcb
    type (mp_real) r1, r2
    real (mprknd) d3
    integer l1, l2, l3, mpnwbt
    l1 = za%mpc(0)
    call mp_fixlocz (za)
    mpnwbt = max (za%mpc(1), za%mpc(l1+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    call mp_checkdp (dble (dcb))
    call mp_checkdp (aimag (dcb))
    call mp_initvr (r1, mpnwbt)
    call mp_initvr (r2, mpnwbt)
    l3 = mpwds6
    call mp_initvz (mp_divzdc, mpnwbt)
    call mpfrmuld (r1%mpr(1:), za%mpc(1:), dble (dcb), mprnd)
    call mpfrmuld (r2%mpr(1:), za%mpc(l1+1:), aimag (dcb), mprnd)
    call mpfradd (mp_divzdc%mpc(1:), r1%mpr(1:), r2%mpr(1:), mprnd)
    call mpfrmuld (r1%mpr(1:), za%mpc(l1+1:), dble (dcb), mprnd)
    call mpfrmuld (r2%mpr(1:), za%mpc(1:), aimag (dcb), mprnd)
    call mpfrsub (mp_divzdc%mpc(l3+1:), r1%mpr(1:), r2%mpr(1:), mprnd)
    d3 = dble (dcb) ** 2 + aimag (dcb) ** 2
    call mpfrdivd (mp_divzdc%mpc(1:), mp_divzdc%mpc(1:), d3, mprnd)
    call mpfrdivd (mp_divzdc%mpc(l3+1:), mp_divzdc%mpc(l3+1:), d3, mprnd)
    return
  end function

  function mp_divrz (ra, zb)
    implicit none
    type (mp_complex):: mp_divrz
    type (mp_real), intent (in):: ra
    type (mp_complex), intent (in):: zb
    type (mp_real) r1, r2, r3
    integer l1, l2, l3, mpnwbt
    l1 = ra%mpr(0)
    call mp_fixlocr (ra)
    l2 = zb%mpc(0)
    call mp_fixlocz (zb)
    mpnwbt = max (ra%mpr(1), zb%mpc(1), zb%mpc(l2+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    call mp_initvr (r1, mpnwbt)
    call mp_initvr (r2, mpnwbt)
    call mp_initvr (r3, mpnwbt)
    l3 = mpwds6
    call mp_initvz (mp_divrz, mpnwbt)
    call mpfrmul (mp_divrz%mpc(1:), ra%mpr(1:), zb%mpc(1:), mprnd)
    call mpfrmul (r1%mpr(1:), ra%mpr(1:), zb%mpc(l2+1:), mprnd)
    call mpfrneg (mp_divrz%mpc(l3+1:), r1%mpr(1:), mprnd)
    call mpfrmul (r1%mpr(1:), zb%mpc(1:), zb%mpc(1:), mprnd)
    call mpfrmul (r2%mpr(1:), zb%mpc(l2+1:), zb%mpc(l2+1:), mprnd)
    call mpfradd (r3%mpr(1:), r1%mpr(1:), r2%mpr(1:), mprnd)
    call mpfrdiv (mp_divrz%mpc(1:), mp_divrz%mpc(1:), r3%mpr(1:), mprnd)
    call mpfrdiv (mp_divrz%mpc(l3+1:), mp_divrz%mpc(l3+1:), r3%mpr(1:), mprnd)
    return
  end function

  function mp_divzr (za, rb)
    implicit none
    type (mp_complex):: mp_divzr
    type (mp_complex), intent (in):: za
    type (mp_real), intent (in):: rb
    integer l1, l2, l3, mpnwbt
    l1 = za%mpc(0)
    call mp_fixlocz (za)
    l2 = rb%mpr(0)
    call mp_fixlocr (rb)
    mpnwbt = max (za%mpc(1), za%mpc(l1+1), rb%mpr(1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    l3 = mpwds6
    call mp_initvz (mp_divzr, mpnwbt)
    call mpfrdiv (mp_divzr%mpc(1:), za%mpc(1:), rb%mpr(1:), mprnd)
    call mpfrdiv (mp_divzr%mpc(l3+1:), za%mpc(l1+1:), rb%mpr(1:), mprnd)
    return
  end function

!  Exponentiation routines:

  function mp_expri (ra, ib)
    implicit none
    type (mp_real):: mp_expri
    type (mp_real), intent (in):: ra
    integer, intent (in):: ib
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_expri, mpnwbt)
    call mpfrpowsi (mp_expri%mpr(1:), ra%mpr(1:), ib, mprnd)
    return
  end function

  function mp_exprr (ra, rb)
    implicit none
    type (mp_real):: mp_exprr
    type (mp_real), intent (in):: ra, rb
    integer mpnwbt
    call mp_fixlocr (ra)
    call mp_fixlocr (rb)
    mpnwbt = max (ra%mpr(1), rb%mpr(1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    call mp_initvr (mp_exprr, mpnwbt)
    call mpfrpow (mp_exprr%mpr(1:), ra%mpr(1:), rb%mpr(1:), mprnd)
    return
  end function

  function mp_expzi (za, ib)
    implicit none
    type (mp_complex):: mp_expzi
    type (mp_complex), intent (in):: za
    integer, intent (in):: ib
    integer j, kk, kn, l1, l2, l3, mpnwbt, mn, nn
    type (mp_complex) z0, z1, z2
    call mp_fixlocz (za)
    l1 = za%mpc(0)
    mpnwbt = max (za%mpc(1), za%mpc(l1+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    l2 = mpwds6
    call mp_initvz (z0, mpnwbt)
    call mp_initvz (z1, mpnwbt)
    call mp_initvz (z2, mpnwbt)
    l3 = mpwds6
    call mp_initvz (mp_expzi, mpnwbt)
    nn = abs (ib)

!   Handle nn = 0, 1, 2 as special cases.

    if (nn == 0) then
      call mpfrsetd (mp_expzi%mpc(1:), 1.d0, mprnd)
      call mpfrsetd (mp_expzi%mpc(l3+1:), 0.d0, mprnd)
      goto 120
    elseif (nn == 1) then
      call mpfrset (z2%mpc(1:), za%mpc(1:), mprnd)
      call mpfrset (z2%mpc(l2+1:), za%mpc(l1+1:), mprnd)
      goto 110
    elseif (nn == 2) then
      z2 = mp_mulzz (za, za)
      goto 110
    endif

!   Determine the least integer mn such that 2^mn > nn.

    mn = log (dble (nn)) / log(2.d0) + 1.d0 + 1.d-14
    call mpfrsetd (z2%mpc(1:), 1.d0, mprnd)
    call mpfrsetd (z2%mpc(l2+1:), 0.d0, mprnd)
    call mpfrset (z0%mpc(1:), za%mpc(1:), mprnd)
    call mpfrset (z0%mpc(l2+1:), za%mpc(l1+1:), mprnd)
    kn = nn

!   Compute za^nn using the binary rule for exponentiation.

    do j = 1, mn
      kk = kn / 2
      if (kn /= 2 * kk) then
        z1 = mp_mulzz (z2, z0)
        call mpfrset (z2%mpc(1:), z1%mpc(1:), mprnd)
        call mpfrset (z2%mpc(l2+1:), z1%mpc(l2+1:), mprnd)
      endif
      kn = kk
      if (j < mn) then
        z1 = mp_mulzz (z0, z0)
        call mpfrset (z0%mpc(1:), z1%mpc(1:), mprnd)
        call mpfrset (z0%mpc(l2+1:), z1%mpc(l2+1:), mprnd)
      endif
    enddo

!   Compute reciprocal if ib is negative.

110 continue

    if (ib < 0) then
      call mpfrsetd (z1%mpc(1:), 1.d0, mprnd)
      call mpfrsetd (z1%mpc(l2+1:), 0.d0, mprnd)
      z0 = mp_divzz (z1, z2)
      call mpfrset (z2%mpc(1:), z0%mpc(1:), mprnd)
      call mpfrset (z2%mpc(l2+1:), z0%mpc(l2+1:), mprnd)
    endif

    call mpfrset (mp_expzi%mpc(1:), z2%mpc(1:), mprnd)
    call mpfrset (mp_expzi%mpc(l3+1:), z2%mpc(l2+1:), mprnd)

120 continue
    return
  end function

  function mp_expzz (za, zb)
    implicit none
    type (mp_complex):: mp_expzz
    type (mp_complex), intent (in):: za, zb
    type (mp_real):: r1, r2, r3, r4, r5, r6
    integer l1, l2, l3, mpnwbt
    call mp_fixlocz (za)
    call mp_fixlocz (zb)
    l1 = za%mpc(0)
    l2 = zb%mpc(0)
    mpnwbt = max (za%mpc(1), za%mpc(l1+1), zb%mpc(1), zb%mpc(l2+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    l3 = mpwds6
    call mp_initvr (r1, mpnwbt)
    call mp_initvr (r2, mpnwbt)
    call mp_initvr (r3, mpnwbt)
    call mp_initvr (r4, mpnwbt)
    call mp_initvr (r5, mpnwbt)
    call mp_initvr (r6, mpnwbt)
    call mp_initvz (mp_expzz, mpnwbt)
    call mpfrmul (r1%mpr(1:), za%mpc(1:), za%mpc(1:), mprnd)
    call mpfrmul (r2%mpr(1:), za%mpc(l1+1:), za%mpc(l1+1:), mprnd)
    call mpfradd (r3%mpr(1:), r1%mpr(1:), r2%mpr(1:), mprnd)
    call mpfrlog (r4%mpr(1:), r3%mpr(1:), mprnd)
    call mpfrmuld (r5%mpr(1:), r4%mpr(1:), 0.5d0, mprnd)
    call mpfrmul (r1%mpr(1:), zb%mpc(1:), r5%mpr(1:), mprnd)
    call mpfratan2 (r2%mpr(1:), za%mpc(l1+1:), za%mpc(1:), mprnd)
    call mpfrmul (r3%mpr(1:), r2%mpr(1:), zb%mpc(l2+1:), mprnd)
    call mpfrsub (r4%mpr(1:), r1%mpr(1:), r3%mpr(1:), mprnd)
    call mpfrexp (r1%mpr(1:), r4%mpr(1:), mprnd)
    call mpfrmul (r3%mpr(1:), zb%mpc(l2+1:), r5%mpr(1:), mprnd)
    call mpfrmul (r4%mpr(1:), zb%mpc(1:), r2%mpr(1:), mprnd)
    call mpfradd (r6%mpr(1:), r3%mpr(1:), r4%mpr(1:), mprnd)
    call mpfrsincos (r4%mpr(1:), r3%mpr(1:), r6%mpr(1:), mprnd)
    call mpfrmul (mp_expzz%mpc(1:), r1%mpr(1:), r3%mpr(1:), mprnd)
    call mpfrmul (mp_expzz%mpc(l3+1:), r1%mpr(1:), r4%mpr(1:), mprnd)
  end function

  function mp_exprz (ra, zb)
    implicit none
    type (mp_complex):: mp_exprz
    type (mp_real), intent (in):: ra
    type (mp_complex), intent (in):: zb
    type (mp_real):: r1, r2, r3, r4, r5
    integer l1, l2, mpnwbt
    call mp_fixlocr (ra)
    call mp_fixlocz (zb)
    l1 = zb%mpc(0)
    mpnwbt = max (ra%mpr(1), zb%mpc(1), zb%mpc(l1+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    l2 = mpwds6
    call mp_initvr (r1, mpnwbt)
    call mp_initvr (r2, mpnwbt)
    call mp_initvr (r3, mpnwbt)
    call mp_initvr (r4, mpnwbt)
    call mp_initvr (r5, mpnwbt)
    call mp_initvz (mp_exprz, mpnwbt)
    call mpfrpow (r1%mpr(1:), ra%mpr(1:), zb%mpc(1:), mprnd)
    call mpfrlog (r2%mpr(1:), ra%mpr(1:), mprnd)
    call mpfrmul (r3%mpr(1:), r2%mpr(1:), zb%mpc(l1+1:), mprnd)
    call mpfrsincos (r5%mpr(1:), r4%mpr(1:), r3%mpr(1:), mprnd)
    call mpfrmul (mp_exprz%mpc(1:), r1%mpr(1:), r4%mpr(1:), mprnd)
    call mpfrmul (mp_exprz%mpc(l2+1:), r1%mpr(1:), r5%mpr(1:), mprnd)
    return
  end function

  function mp_expzr (za, rb)
    implicit none
    type (mp_complex):: mp_expzr
    type (mp_complex), intent (in):: za
    type (mp_real), intent (in):: rb
    type (mp_real):: r1, r2, r3, r4, r5
    integer l1, l2, mpnwbt
    call mp_fixlocz (za)
    call mp_fixlocr (rb)
    l1 = za%mpc(0)
    mpnwbt = max (za%mpc(1), za%mpc(l1+1), rb%mpr(1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    l2 = mpwds6
    call mp_initvr (r1, mpnwbt)
    call mp_initvr (r2, mpnwbt)
    call mp_initvr (r3, mpnwbt)
    call mp_initvr (r4, mpnwbt)
    call mp_initvr (r5, mpnwbt)
    call mp_initvz (mp_expzr, mpnwbt)
    call mpfrmul (r1%mpr(1:), za%mpc(1:), za%mpc(1:), mprnd)
    call mpfrmul (r2%mpr(1:), za%mpc(l1+1:), za%mpc(l1+1:), mprnd)
    call mpfradd (r3%mpr(1:), r1%mpr(1:), r2%mpr(1:), mprnd)
    call mpfrlog (r4%mpr(1:), r3%mpr(1:), mprnd)
    call mpfrmuld (r5%mpr(1:), r4%mpr(1:), 0.5d0, mprnd)
    call mpfrmul (r1%mpr(1:), r5%mpr(1:), rb%mpr(1:), mprnd)
    call mpfrexp (r2%mpr(1:), r1%mpr(1:), mprnd)
    call mpfratan2 (r3%mpr(1:), za%mpc(l1+1:), za%mpc(1:), mprnd)
    call mpfrmul (r1%mpr(1:), rb%mpr(1:), r3%mpr(1:), mprnd)
    call mpfrsincos (r5%mpr(1:), r4%mpr(1:), r1%mpr(1:), mprnd)
    call mpfrmul (mp_expzr%mpc(1:), r2%mpr(1:), r4%mpr(1:), mprnd)
    call mpfrmul (mp_expzr%mpc(l2+1:), r2%mpr(1:), r5%mpr(1:), mprnd)
    return
  end function

!  Equality test routines:

  function mp_eqtrr (ra, rb)
    implicit none
    logical mp_eqtrr
    type (mp_real), intent (in):: ra, rb
    integer ic, mpfrcmp
    external mpfrcmp
    call mp_fixlocr (ra)
    call mp_fixlocr (rb)
    ic = mpfrcmp (ra%mpr(1:), rb%mpr(1:))
    if (ic == 0) then
      mp_eqtrr = .true.
    else
      mp_eqtrr = .false.
    endif
    return
  end function

  function mp_eqtdr (da, rb)
    implicit none
    logical mp_eqtdr
    real (mprknd), intent (in):: da
    type (mp_real), intent (in):: rb
    integer ic, mpfrcmpd
    external mpfrcmpd
    call mp_checkdp (da)
    call mp_fixlocr (rb)
    ic = - mpfrcmpd (rb%mpr(1:), da)
    if (ic == 0) then
      mp_eqtdr = .true.
    else
      mp_eqtdr = .false.
    endif
    return
  end function

  function mp_eqtrd (ra, db)
    implicit none
    logical mp_eqtrd
    type (mp_real), intent (in):: ra
    real (mprknd), intent (in):: db
    integer ic, mpfrcmpd
    external mpfrcmpd
    call mp_checkdp (db)
    call mp_fixlocr (ra)
    ic = mpfrcmpd (ra%mpr(1:), db)
    if (ic == 0) then
      mp_eqtrd = .true.
    else
      mp_eqtrd = .false.
    endif
    return
  end function

  function mp_eqtir (ia, rb)
    implicit none
    logical mp_eqtir
    integer, intent (in):: ia
    type (mp_real), intent (in):: rb
    real (mprknd) da
    integer ic, mpfrcmpd
    external mpfrcmpd
    da = ia
    call mp_checkdp (da)
    call mp_fixlocr (rb)
    ic = - mpfrcmpd (rb%mpr(1:), da)
    if (ic == 0) then
      mp_eqtir = .true.
    else
      mp_eqtir = .false.
    endif
    return
  end function

  function mp_eqtri (ra, ib)
    implicit none
    logical mp_eqtri
    type (mp_real), intent (in):: ra
    integer, intent (in):: ib
    real (mprknd) db
    integer ic, mpfrcmpd
    external mpfrcmpd
    db = ib
    call mp_checkdp (db)
    call mp_fixlocr (ra)
    ic = mpfrcmpd (ra%mpr(1:), db)
    if (ic == 0) then
      mp_eqtri = .true.
    else
      mp_eqtri = .false.
    endif
    return
  end function

  function mp_eqtzz (za, zb)
    implicit none
    logical mp_eqtzz
    type (mp_complex), intent (in):: za, zb
    integer ic1, ic2, l1, l2
    integer mpfrcmp
    external mpfrcmp
    call mp_fixlocz (za)
    call mp_fixlocz (zb)
    l1 = za%mpc(0)
    l2 = zb%mpc(0)
    ic1 = mpfrcmp (za%mpc(1:), zb%mpc(1:))
    ic2 = mpfrcmp (za%mpc(l1+1:), zb%mpc(l2+1:))
    if (ic1 == 0 .and. ic2 == 0) then
      mp_eqtzz = .true.
    else
      mp_eqtzz = .false.
    endif
    return
  end function

  function mp_eqtdz (da, zb)
    implicit none
    logical mp_eqtdz
    real (mprknd), intent (in):: da
    type (mp_complex), intent (in):: zb
    integer ic1, ic2, l2
    integer mpfrcmpd
    external mpfrcmpd
    call mp_checkdp (da)
    call mp_fixlocz (zb)
    l2 = zb%mpc(0)
    ic1 = - mpfrcmpd (zb%mpc(1:), da)
    if (zb%mpc(l2+3) == mpzero) then
      ic2 = 0
    else
      ic2 = 1
    endif
    if (ic1 == 0 .and. ic2 == 0) then
      mp_eqtdz = .true.
    else
      mp_eqtdz = .false.
    endif
    return
  end function

  function mp_eqtzd (za, db)
    implicit none
    logical mp_eqtzd
    type (mp_complex), intent (in):: za
    real (mprknd), intent (in):: db
    integer ic1, ic2, l1
    integer mpfrcmpd
    external mpfrcmpd
    call mp_checkdp (db)
    call mp_fixlocz (za)
    l1 = za%mpc(0)
    ic1 = mpfrcmpd (za%mpc(1:), db)
    if (za%mpc(l1+3) == mpzero) then
      ic2 = 0
    else
      ic2 = 1
    endif
    if (ic1 == 0 .and. ic2 == 0) then
      mp_eqtzd = .true.
    else
      mp_eqtzd = .false.
    endif
    return
  end function

  function mp_eqtdcz (dca, zb)
    implicit none
    logical mp_eqtdcz
    complex (kind (0.d0)), intent (in):: dca
    type (mp_complex), intent (in):: zb
    real (mprknd) da, db
    integer ic1, ic2, l2
    integer mpfrcmpd
    external mpfrcmpd
    da = dble (dca)
    db = aimag (dca)
    call mp_checkdp (da)
    call mp_checkdp (db)
    call mp_fixlocz (zb)
    l2 = zb%mpc(0)
    ic1 = - mpfrcmpd (zb%mpc(1:), da)
    ic2 = - mpfrcmpd (zb%mpc(l2+1:), db)
    if (ic1 == 0 .and. ic2 == 0) then
      mp_eqtdcz = .true.
    else
      mp_eqtdcz = .false.
    endif
    return
  end function

  function mp_eqtzdc (za, dcb)
    implicit none
    logical mp_eqtzdc
    type (mp_complex), intent (in):: za
    complex (kind (0.d0)), intent (in):: dcb
    real (mprknd) da, db
    integer ic1, ic2, l1
    integer mpfrcmpd
    external mpfrcmpd
    da = dble (dcb)
    db = aimag (dcb)
    call mp_checkdp (da)
    call mp_checkdp (db)
    call mp_fixlocz (za)
    l1 = za%mpc(0)
    ic1 = mpfrcmpd (za%mpc(1:), da)
    ic2 = mpfrcmpd (za%mpc(l1+1:), db)
    if (ic1 == 0 .and. ic2 == 0) then
      mp_eqtzdc = .true.
    else
      mp_eqtzdc = .false.
    endif
    return
  end function

  function mp_eqtrz (ra, zb)
    implicit none
    logical mp_eqtrz
    type (mp_real), intent (in):: ra
    type (mp_complex), intent (in):: zb
    integer ic1, ic2, l2, mpfrcmp
    external mpfrcmp
    call mp_fixlocr (ra)
    call mp_fixlocz (zb)
    l2 = zb%mpc(0)
    ic1 = mpfrcmp (ra%mpr(1:), zb%mpc(1:))
    if (zb%mpc(l2+3) == mpzero) then
      ic2 = 0
    else
      ic2 = 1
    endif
    if (ic1 == 0 .and. ic2 == 0) then
      mp_eqtrz = .true.
    else
      mp_eqtrz = .false.
    endif
    return
  end function

  function mp_eqtzr (za, rb)
    implicit none
    logical mp_eqtzr
    type (mp_complex), intent (in):: za
    type (mp_real), intent (in):: rb
    integer ic1, ic2, l1, mpfrcmp
    external mpfrcmp
    call mp_fixlocz (za)
    call mp_fixlocr (rb)
    l1 = za%mpc(0)
    ic1 = mpfrcmp (za%mpc(1:), rb%mpr(1:))
    if (za%mpc(l1+3) == mpzero) then
      ic2 = 0
    else
      ic2 = 1
    endif
    if (ic1 == 0 .and. ic2 == 0) then
      mp_eqtzr = .true.
    else
      mp_eqtzr = .false.
    endif
    return
  end function

!  Non-equality test routines:

  function mp_netrr (ra, rb)
    implicit none
    logical mp_netrr
    type (mp_real), intent (in):: ra, rb
    integer ic, mpfrcmp
    external mpfrcmp
    call mp_fixlocr (ra)
    call mp_fixlocr (rb)
    ic = mpfrcmp (ra%mpr(1:), rb%mpr(1:))
    if (ic /= 0) then
      mp_netrr = .true.
    else
      mp_netrr = .false.
    endif
    return
  end function

  function mp_netdr (da, rb)
    implicit none
    logical mp_netdr
    real (mprknd), intent (in):: da
    type (mp_real), intent (in):: rb
    integer ic, mpfrcmpd
    external mpfrcmpd
    call mp_checkdp (da)
    call mp_fixlocr (rb)
    ic = - mpfrcmpd (rb%mpr(1:), da)
    if (ic /= 0) then
      mp_netdr = .true.
    else
      mp_netdr = .false.
    endif
    return
  end function

  function mp_netrd (ra, db)
    implicit none
    logical mp_netrd
    type (mp_real), intent (in):: ra
    real (mprknd), intent (in):: db
    integer ic, mpfrcmpd
    external mpfrcmpd
    call mp_checkdp (db)
    call mp_fixlocr (ra)
    ic = mpfrcmpd (ra%mpr(1:), db)
    if (ic /= 0) then
      mp_netrd = .true.
    else
      mp_netrd = .false.
    endif
    return
  end function

  function mp_netir (ia, rb)
    implicit none
    logical mp_netir
    integer, intent (in):: ia
    type (mp_real), intent (in):: rb
    real (mprknd) da
    integer ic, mpfrcmpd
    external mpfrcmpd
    da = ia
    call mp_checkdp (da)
    call mp_fixlocr (rb)
    ic = - mpfrcmpd (rb%mpr(1:), da)
    if (ic /= 0) then
      mp_netir = .true.
    else
      mp_netir = .false.
    endif
    return
  end function

  function mp_netri (ra, ib)
    implicit none
    logical mp_netri
    type (mp_real), intent (in):: ra
    integer, intent (in):: ib
    real (mprknd) db
    integer ic, mpfrcmpd
    external mpfrcmpd
    db = ib
    call mp_checkdp (db)
    call mp_fixlocr (ra)
    ic = mpfrcmpd (ra%mpr(1:), db)
    if (ic /= 0) then
      mp_netri = .true.
    else
      mp_netri = .false.
    endif
    return
  end function

  function mp_netzz (za, zb)
    implicit none
    logical mp_netzz
    type (mp_complex), intent (in):: za, zb
    integer ic1, ic2, l1, l2
    integer mpfrcmp
    external mpfrcmp
    call mp_fixlocz (za)
    call mp_fixlocz (zb)
    l1 = za%mpc(0)
    l2 = zb%mpc(0)
    ic1 = mpfrcmp (za%mpc(1:), zb%mpc(1:))
    ic2 = mpfrcmp (za%mpc(l1+1:), zb%mpc(l2+1:))
    if (ic1 /= 0 .or. ic2 /= 0) then
      mp_netzz = .true.
    else
      mp_netzz = .false.
    endif
    return
  end function

  function mp_netdz (da, zb)
    implicit none
    logical mp_netdz
    real (mprknd), intent (in):: da
    type (mp_complex), intent (in):: zb
    integer ic1, ic2, l2
    integer mpfrcmpd
    external mpfrcmpd
    call mp_checkdp (da)
    call mp_fixlocz (zb)
    l2 = zb%mpc(0)
    ic1 = - mpfrcmpd (zb%mpc(1:), da)
    if (zb%mpc(l2+3) == mpzero) then
      ic2 = 0
    else
      ic2 = 1
    endif
    if (ic1 /= 0 .or. ic2 /= 0) then
      mp_netdz = .true.
    else
      mp_netdz = .false.
    endif
    return
  end function

  function mp_netzd (za, db)
    implicit none
    logical mp_netzd
    type (mp_complex), intent (in):: za
    real (mprknd), intent (in):: db
    integer ic1, ic2, l1
    integer mpfrcmpd
    external mpfrcmpd
    call mp_checkdp (db)
    call mp_fixlocz (za)
    l1 = za%mpc(0)
    ic1 = mpfrcmpd (za%mpc(1:), db)
    if (za%mpc(l1+3) == mpzero) then
      ic2 = 0
    else
      ic2 = 1
    endif
    if (ic1 /= 0 .or. ic2 /= 0) then
      mp_netzd = .true.
    else
      mp_netzd = .false.
    endif
    return
  end function

  function mp_netdcz (dca, zb)
    implicit none
    logical mp_netdcz
    complex (kind (0.d0)), intent (in):: dca
    type (mp_complex), intent (in):: zb
    real (mprknd) da, db
    integer ic1, ic2, l2
    integer mpfrcmpd
    external mpfrcmpd
    da = dble(dca)
    db = aimag(dca)
    call mp_checkdp (da)
    call mp_checkdp (db)
    call mp_fixlocz (zb)
    l2 = zb%mpc(0)
    ic1 = - mpfrcmpd (zb%mpc(1:), da)
    ic2 = - mpfrcmpd (zb%mpc(l2+1:), db)
    if (ic1 /= 0 .or. ic2 /= 0) then
      mp_netdcz = .true.
    else
      mp_netdcz = .false.
    endif
    return
  end function

  function mp_netzdc (za, dcb)
    implicit none
    logical mp_netzdc
    type (mp_complex), intent (in):: za
    complex (kind (0.d0)), intent (in):: dcb
    real (mprknd) da, db
    integer ic1, ic2, l1
    integer mpfrcmpd
    external mpfrcmpd
    da = dble (dcb)
    db = aimag (dcb)
    call mp_checkdp (da)
    call mp_checkdp (db)
    call mp_fixlocz (za)
    l1 = za%mpc(0)
    ic1 = mpfrcmpd (za%mpc(1:), da)
    ic2 = mpfrcmpd (za%mpc(l1+1:), db)
    if (ic1 /= 0 .or. ic2 /= 0) then
      mp_netzdc = .true.
    else
      mp_netzdc = .false.
    endif
    return
  end function

  function mp_netrz (ra, zb)
    implicit none
    logical mp_netrz
    type (mp_real), intent (in):: ra
    type (mp_complex), intent (in):: zb
    integer ic1, ic2, l2, mpfrcmp
    external mpfrcmp
    call mp_fixlocr (ra)
    call mp_fixlocz (zb)
    l2 = zb%mpc(0)
    ic1 = mpfrcmp (ra%mpr(1:), zb%mpc(1:))
    if (zb%mpc(l2+3) == mpzero) then
      ic2 = 0
    else
      ic2 = 1
    endif
    if (ic1 /= 0 .or. ic2 /= 0) then
      mp_netrz = .true.
    else
      mp_netrz = .false.
    endif
    return
  end function

  function mp_netzr (za, rb)
    implicit none
    logical mp_netzr
    type (mp_complex), intent (in):: za
    type (mp_real), intent (in):: rb
    integer ic1, ic2, l1, mpfrcmp
    external mpfrcmp
    call mp_fixlocz (za)
    call mp_fixlocr (rb)
    l1 = za%mpc(0)
    ic1 = mpfrcmp (za%mpc(1:), rb%mpr(1:))
    if (za%mpc(l1+3) == mpzero) then
      ic2 = 0
    else
      ic2 = 1
    endif
    if (ic1 /= 0 .or. ic2 /= 0) then
      mp_netzr = .true.
    else
      mp_netzr = .false.
    endif
    return
  end function

!  Less-than-or-equal test routines:

  function mp_letrr (ra, rb)
    implicit none
    logical mp_letrr
    type (mp_real), intent (in):: ra, rb
    integer ic, mpfrcmp
    external mpfrcmp
    call mp_fixlocr (ra)
    call mp_fixlocr (rb)
    ic = mpfrcmp (ra%mpr(1:), rb%mpr(1:))
    if (ic <= 0) then
      mp_letrr = .true.
    else
      mp_letrr = .false.
    endif
    return
  end function

  function mp_letdr (da, rb)
    implicit none
    logical mp_letdr
    real (mprknd), intent (in):: da
    type (mp_real), intent (in):: rb
    integer ic, mpfrcmpd
    external mpfrcmpd
    call mp_checkdp (da)
    call mp_fixlocr (rb)
    ic = - mpfrcmpd (rb%mpr(1:), da)
    if (ic <= 0) then
      mp_letdr = .true.
    else
      mp_letdr = .false.
    endif
    return
  end function

  function mp_letrd (ra, db)
    implicit none
    logical mp_letrd
    type (mp_real), intent (in):: ra
    real (mprknd), intent (in):: db
    integer ic, mpfrcmpd
    external mpfrcmpd
    call mp_checkdp (db)
    call mp_fixlocr (ra)
    ic = mpfrcmpd (ra%mpr(1:), db)
    if (ic <= 0) then
      mp_letrd = .true.
    else
      mp_letrd = .false.
    endif
    return
  end function

  function mp_letir (ia, rb)
    implicit none
    logical mp_letir
    integer, intent (in):: ia
    type (mp_real), intent (in):: rb
    real (mprknd) da
    integer ic, mpfrcmpd
    external mpfrcmpd
    da = ia
    call mp_checkdp (da)
    call mp_fixlocr (rb)
    ic = - mpfrcmpd (rb%mpr(1:), da)
    if (ic <= 0) then
      mp_letir = .true.
    else
      mp_letir = .false.
    endif
    return
  end function

  function mp_letri (ra, ib)
    implicit none
    logical mp_letri
    type (mp_real), intent (in):: ra
    integer, intent (in):: ib
    real (mprknd) db
    integer ic, mpfrcmpd
    external mpfrcmpd
    db = ib
    call mp_checkdp (db)
    call mp_fixlocr (ra)
    ic = mpfrcmpd (ra%mpr(1:), db)
    if (ic <= 0) then
      mp_letri = .true.
    else
      mp_letri = .false.
    endif
    return
  end function

!  Greater-than-or-equal test routines:

  function mp_getrr (ra, rb)
    implicit none
    logical mp_getrr
    type (mp_real), intent (in):: ra, rb
    integer ic, mpfrcmp
    external mpfrcmp
    call mp_fixlocr (ra)
    call mp_fixlocr (rb)
    ic = mpfrcmp (ra%mpr(1:), rb%mpr(1:))
    if (ic >= 0) then
      mp_getrr = .true.
    else
      mp_getrr = .false.
    endif
    return
  end function

  function mp_getdr (da, rb)
    implicit none
    logical mp_getdr
    real (mprknd), intent (in):: da
    type (mp_real), intent (in):: rb
    integer ic, mpfrcmpd
    external mpfrcmpd
    call mp_checkdp (da)
    call mp_fixlocr (rb)
    ic = - mpfrcmpd (rb%mpr(1:), da)
    if (ic >= 0) then
      mp_getdr = .true.
    else
      mp_getdr = .false.
    endif
    return
  end function

  function mp_getrd (ra, db)
    implicit none
    logical mp_getrd
    type (mp_real), intent (in):: ra
    real (mprknd), intent (in):: db
    integer ic, mpfrcmpd
    external mpfrcmpd
    call mp_checkdp (db)
    call mp_fixlocr (ra)
    ic = mpfrcmpd (ra%mpr(1:), db)
    if (ic >= 0) then
      mp_getrd = .true.
    else
      mp_getrd = .false.
    endif
    return
  end function

  function mp_getir (ia, rb)
    implicit none
    logical mp_getir
    integer, intent (in):: ia
    type (mp_real), intent (in):: rb
    real (mprknd) da
    integer ic, mpfrcmpd
    external mpfrcmpd
    da = ia
    call mp_checkdp (da)
    call mp_fixlocr (rb)
    ic = - mpfrcmpd (rb%mpr(1:), da)
    if (ic >= 0) then
      mp_getir = .true.
    else
      mp_getir = .false.
    endif
    return
  end function

  function mp_getri (ra, ib)
    implicit none
    logical mp_getri
    type (mp_real), intent (in):: ra
    integer, intent (in):: ib
    real (mprknd) db
    integer ic, mpfrcmpd
    external mpfrcmpd
    db = ib
    call mp_checkdp (db)
    call mp_fixlocr (ra)
    ic = mpfrcmpd (ra%mpr(1:), db)
    if (ic >= 0) then
      mp_getri = .true.
    else
      mp_getri = .false.
    endif
    return
  end function

!  Less-than test routines:

  function mp_lttrr (ra, rb)
    implicit none
    logical mp_lttrr
    type (mp_real), intent (in):: ra, rb
    integer ic, mpfrcmp
    external mpfrcmp
    call mp_fixlocr (ra)
    call mp_fixlocr (rb)
    ic = mpfrcmp (ra%mpr(1:), rb%mpr(1:))
    if (ic < 0) then
      mp_lttrr = .true.
    else
      mp_lttrr = .false.
    endif
    return
  end function

  function mp_lttdr (da, rb)
    implicit none
    logical mp_lttdr
    real (mprknd), intent (in):: da
    type (mp_real), intent (in):: rb
    integer ic, mpfrcmpd
    external mpfrcmpd
    call mp_checkdp (da)
    call mp_fixlocr (rb)
    ic = - mpfrcmpd (rb%mpr(1:), da)
    if (ic < 0) then
      mp_lttdr = .true.
    else
      mp_lttdr = .false.
    endif
    return
  end function

  function mp_lttrd (ra, db)
    implicit none
    logical mp_lttrd
    type (mp_real), intent (in):: ra
    real (mprknd), intent (in):: db
    integer ic, mpfrcmpd
    external mpfrcmpd
    call mp_checkdp (db)
    call mp_fixlocr (ra)
    ic = mpfrcmpd (ra%mpr(1:), db)
    if (ic < 0) then
      mp_lttrd = .true.
    else
      mp_lttrd = .false.
    endif
    return
  end function

  function mp_lttir (ia, rb)
    implicit none
    logical mp_lttir
    integer, intent (in):: ia
    type (mp_real), intent (in):: rb
    real (mprknd) da
    integer ic, mpfrcmpd
    external mpfrcmpd
    da = ia
    call mp_checkdp (da)
    call mp_fixlocr (rb)
    ic = - mpfrcmpd (rb%mpr(1:), da)
    if (ic < 0) then
      mp_lttir = .true.
    else
      mp_lttir = .false.
    endif
    return
  end function

  function mp_lttri (ra, ib)
    implicit none
    logical mp_lttri
    type (mp_real), intent (in):: ra
    integer, intent (in):: ib
    real (mprknd) db
    integer ic, mpfrcmpd
    external mpfrcmpd
    db = ib
    call mp_checkdp (db)
    call mp_fixlocr (ra)
    ic = mpfrcmpd (ra%mpr(1:), db)
    if (ic < 0) then
      mp_lttri = .true.
    else
      mp_lttri = .false.
    endif
    return
  end function

!  Greater-than test routines:

  function mp_gttrr (ra, rb)
    implicit none
    logical mp_gttrr
    type (mp_real), intent (in):: ra, rb
    integer ic, mpfrcmp
    external mpfrcmp
    call mp_fixlocr (ra)
    call mp_fixlocr (rb)
    ic = mpfrcmp (ra%mpr(1:), rb%mpr(1:))
    if (ic > 0) then
      mp_gttrr = .true.
    else
      mp_gttrr = .false.
    endif
    return
  end function

  function mp_gttdr (da, rb)
    implicit none
    logical mp_gttdr
    real (mprknd), intent (in):: da
    type (mp_real), intent (in):: rb
    integer ic, mpfrcmpd
    external mpfrcmpd
    call mp_checkdp (da)
    call mp_fixlocr (rb)
    ic = - mpfrcmpd (rb%mpr(1:), da)
    if (ic > 0) then
      mp_gttdr = .true.
    else
      mp_gttdr = .false.
    endif
    return
  end function

  function mp_gttrd (ra, db)
    implicit none
    logical mp_gttrd
    type (mp_real), intent (in):: ra
    real (mprknd), intent (in):: db
    integer ic, mpfrcmpd
    external mpfrcmpd
    call mp_checkdp (db)
    call mp_fixlocr (ra)
    ic = mpfrcmpd (ra%mpr(1:), db)
    if (ic > 0) then
      mp_gttrd = .true.
    else
      mp_gttrd = .false.
    endif
    return
  end function

  function mp_gttir (ia, rb)
    implicit none
    logical mp_gttir
    integer, intent (in):: ia
    type (mp_real), intent (in):: rb
    real (mprknd) da
    integer ic, mpfrcmpd
    external mpfrcmpd
    da = ia
    call mp_checkdp (da)
    call mp_fixlocr (rb)
    ic = - mpfrcmpd (rb%mpr(1:), da)
    if (ic > 0) then
      mp_gttir = .true.
    else
      mp_gttir = .false.
    endif
    return
  end function

  function mp_gttri (ra, ib)
    implicit none
    logical mp_gttri
    type (mp_real), intent (in):: ra
    integer, intent (in):: ib
    real (mprknd) db
    integer ic, mpfrcmpd
    external mpfrcmpd
    db = ib
    call mp_checkdp (db)
    call mp_fixlocr (ra)
    ic = mpfrcmpd (ra%mpr(1:), db)
    if (ic > 0) then
      mp_gttri = .true.
    else
      mp_gttri = .false.
    endif
    return
  end function

!  Algebraic and transcendental function definitions, listed alphabetically:

  function mp_absr (ra)
    implicit none
    type (mp_real):: mp_absr
    type (mp_real), intent (in):: ra
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_absr, mpnwbt)
    call mpfrabs (mp_absr%mpr(1:), ra%mpr(1:), mprnd)
    return
  end function

  function mp_absz (za)
    implicit none
    type (mp_real):: mp_absz
    type (mp_complex), intent (in):: za
    integer l1, mpnwbt
    type (mp_real) r1, r2, r3
    call mp_fixlocz (za)
    l1 = za%mpc(0)
    mpnwbt = max (za%mpc(1), za%mpc(l1+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    call mp_initvr (r1, mpnwbt)
    call mp_initvr (r2, mpnwbt)
    call mp_initvr (r3, mpnwbt)
    call mp_initvr (mp_absz, mpnwbt)
    call mpfrmul (r1%mpr(1:), za%mpc(1:), za%mpc(1:), mprnd)
    call mpfrmul (r2%mpr(1:), za%mpc(l1+1:), za%mpc(l1+1:), mprnd)
    call mpfradd (r3%mpr(1:), r1%mpr(1:), r2%mpr(1:), mprnd)
    call mpfrsqrt (mp_absz%mpr(1:), r3%mpr(1:), mprnd)
    return
  end function

  function mp_acos (ra)
    implicit none
    type (mp_real):: mp_acos
    type (mp_real), intent (in):: ra
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_acos, mpnwbt)
    call mpfracos (mp_acos%mpr(1:), ra%mpr(1:), mprnd)
    return
  end function

  function mp_acosh (ra)
    implicit none
    type (mp_real):: mp_acosh
    type (mp_real), intent (in):: ra
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_acosh, mpnwbt)
    call mpfracosh (mp_acosh%mpr(1:), ra%mpr(1:), mprnd)
    return
  end function

   function mp_agm (ra, rb)
    implicit none
    type (mp_real):: mp_agm
    type (mp_real), intent (in):: ra, rb
    integer mpnwbt
    call mp_fixlocr (ra)
    call mp_fixlocr (rb)
    mpnwbt = max (ra%mpr(1), rb%mpr(1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    call mp_initvr (mp_agm, mpnwbt)
    call mpfragm (mp_agm%mpr(1:), ra%mpr(1:), rb%mpr(1:), mprnd)
    return
  end function

  function mp_aimag (za)
    implicit none
    type (mp_real):: mp_aimag
    type (mp_complex), intent (in):: za
    integer l1, mpnwbt
    call mp_fixlocz (za)
    l1 = za%mpc(0)
    mpnwbt = max (za%mpc(1), za%mpc(l1+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    call mp_initvr (mp_aimag, mpnwbt)
    call mpfrset (mp_aimag%mpr(1:), za%mpc(l1+1:), mprnd)
    return
  end function

  function mp_aint (ra)
    implicit none
    type (mp_real):: mp_aint
    type (mp_real), intent (in):: ra
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_aint, mpnwbt)
    call mpfrtrunc (mp_aint%mpr(1:), ra%mpr(1:))
    return
  end function

  function mp_airy (ra)
    implicit none
    type (mp_real):: mp_airy
    type (mp_real), intent (in):: ra
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_airy, mpnwbt)
    call mpfrai (mp_airy%mpr(1:), ra%mpr(1:), mprnd)
    return
  end function

   function mp_anint (ra)
    implicit none
    type (mp_real):: mp_anint
    type (mp_real), intent (in):: ra
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_anint, mpnwbt)
    call mpfrround (mp_anint%mpr(1:), ra%mpr(1:))
    return
  end function

   function mp_asin (ra)
    implicit none
    type (mp_real):: mp_asin
    type (mp_real), intent (in):: ra
    type (mp_real) r1, r2, r3
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_asin, mpnwbt)
    call mpfrasin (mp_asin%mpr(1:), ra%mpr(1:), mprnd)
    return
  end function

   function mp_asinh (ra)
    implicit none
    type (mp_real):: mp_asinh
    type (mp_real), intent (in):: ra
    type (mp_real) r1, r2, r3
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_asinh, mpnwbt)
    call mpfrasinh (mp_asinh%mpr(1:), ra%mpr(1:), mprnd)
    return
  end function

   function mp_atan (ra)
    implicit none
    type (mp_real):: mp_atan
    type (mp_real), intent (in):: ra
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_atan, mpnwbt)
    call mpfratan (mp_atan%mpr(1:), ra%mpr(1:), mprnd)
    return
  end function

   function mp_atanh (ra)
    implicit none
    type (mp_real):: mp_atanh
    type (mp_real), intent (in):: ra
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_atanh, mpnwbt)
    call mpfratanh (mp_atanh%mpr(1:), ra%mpr(1:), mprnd)
    return
  end function

   function mp_atan2 (ra, rb)
    implicit none
    type (mp_real):: mp_atan2
    type (mp_real), intent (in):: ra, rb
    integer mpnwbt
    call mp_fixlocr (ra)
    call mp_fixlocr (rb)
    mpnwbt = max (ra%mpr(1), rb%mpr(1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    call mp_initvr (mp_atan2, mpnwbt)
    call mpfratan2 (mp_atan2%mpr(1:), ra%mpr(1:), rb%mpr(1:), mprnd)
    return
  end function

  function mp_ator1 (a, ib, iprec)
    implicit none
    type (mp_real):: mp_ator1
    integer, intent (in):: ib
    character(1), intent (in):: a(ib)
    character(1) a1(ib+1)
    integer i, mpnwbt

!>  In variant #1, uncomment these lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnwbt = mpnbt * mp_setwp (iprec)
!    else
!      mpnwbt = mpnbt * mpwds
!    endif
!  Otherwise in variant #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnwbt = mpnbt * mp_setwp (iprec)
!>>
    call mp_initvr (mp_ator1, mpnwbt)
    do i = 1, ib
      if (a(i) == 'D' .or. a(i) == 'd') then
        a1(i) = 'e'
      else
        a1(i) = a(i)
      endif
    enddo
    a1(ib+1) = char(0)
    call mpfrsetstr (mp_ator1%mpr(1:), a1, mprnd)
    return
  end function

  function mp_atorn (aa, iprec)
    implicit none
    character(  *), intent (in):: aa
    type (mp_real):: mp_atorn
    character(1) :: chr1(len(aa)+1)
    integer i, l1, mpnwbt

!>  In variant #1, uncomment these lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnwbt = mpnbt * mp_setwp (iprec)
!    else
!      mpnwbt = mpnbt * mpwds
!    endif
!  Otherwise in variant #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnwbt = mpnbt * mp_setwp (iprec)
!>>
    l1 = len (aa)
    call mp_initvr (mp_atorn, mpnwbt)
    do i = 1, l1
      if (aa(i:i) == 'D' .or. aa(i:i) == 'd') then
        chr1(i) = 'e'
      else
        chr1(i) = aa(i:i)
      endif
    enddo
    chr1(l1+1) = char(0)
    call mpfrsetstr (mp_atorn%mpr(1:), chr1, mprnd)
    return
  end function

  subroutine mp_berne (nb, rb, iprec)
    implicit none
    integer, intent (in):: nb
    type (mp_real), intent (out):: rb(nb)
    integer i, n1, mpnwbt, mpnw

!>  In variant #1, uncomment these lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnwbt = mpnbt * mp_setwp (iprec)
!    else
!      mpnwbt = mpnbt * mpwds
!    endif
!  Otherwise in variant #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnwbt = mpnbt * mp_setwp (iprec)
!>>
    n1 = mpwds
    mpnw = mpnwbt / mpnbt
    do i = 1, nb
      call mp_initvr (rb(i), mpnwbt)
    enddo
    call mpberner (n1, nb, rb(1)%mpr(0), mpnw)
    return
  end subroutine

  function mp_bessel_i (ra, rb)
    implicit none
    type (mp_real):: mp_bessel_i
    type (mp_real), intent (in):: ra, rb
    integer mpnwbt, mpnw
    call mp_fixlocr (ra)
    call mp_fixlocr (rb)
    mpnwbt = min (ra%mpr(1), rb%mpr(1), mpwdsbt)
    mpnw = mpnwbt / mpnbt
    call mp_initvr (mp_bessel_i, mpnwbt)
    call mpbesselir (ra%mpr(0:), rb%mpr(0:), mp_bessel_i%mpr(0:), mpnw)
    return
  end function

  function mp_bessel_in (ia, rb)
    implicit none
    type (mp_real):: mp_bessel_in
    integer, intent (in):: ia
    type (mp_real), intent (in):: rb
    integer mpnwbt, mpnw
    call mp_fixlocr (rb)
    mpnwbt = min (rb%mpr(1), mpwdsbt)
    mpnw = mpnwbt / mpnbt
    call mp_initvr (mp_bessel_in, mpnwbt)
    call mpbesselinr (ia, rb%mpr(0:), mp_bessel_in%mpr(0:), mpnw)
    return
  end function

  function mp_bessel_j (ra, rb)
    implicit none
    type (mp_real):: mp_bessel_j
    type (mp_real), intent (in):: ra, rb
    integer mpnwbt, mpnw
    call mp_fixlocr (ra)
    call mp_fixlocr (rb)
    mpnwbt = min (ra%mpr(1), rb%mpr(1), mpwdsbt)
    mpnw = mpnwbt / mpnbt
    call mp_initvr (mp_bessel_j, mpnwbt)
    call mpbesseljr (ra%mpr(0:), rb%mpr(0:), mp_bessel_j%mpr(0:), mpnw)
    return
  end function

  function mp_bessel_j0 (ra)
    implicit none
    type (mp_real):: mp_bessel_j0
    type (mp_real), intent (in):: ra
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_bessel_j0, mpnwbt)
    call mpfrj0 (mp_bessel_j0%mpr(1:), ra%mpr(1:), mprnd)
    return
  end function

  function mp_bessel_j1 (ra)
    implicit none
    type (mp_real):: mp_bessel_j1
    type (mp_real), intent (in):: ra
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_bessel_j1, mpnwbt)
    call mpfrj1 (mp_bessel_j1%mpr(1:), ra%mpr(1:), mprnd)
    return
  end function

  function mp_bessel_jn (ia, rb)
    implicit none
    type (mp_real):: mp_bessel_jn
    integer, intent (in):: ia
    type (mp_real), intent (in):: rb
    integer mpnwbt
    call mp_fixlocr (rb)
    mpnwbt = min (rb%mpr(1), mpwdsbt)
    call mp_initvr (mp_bessel_jn, mpnwbt)
    call mpfrbesseljn (mp_bessel_jn%mpr(1:), ia, rb%mpr(1:), mprnd)
    return
  end function

  function mp_bessel_k (ra, rb)
    implicit none
    type (mp_real):: mp_bessel_k
    type (mp_real), intent (in):: ra, rb
    integer mpnwbt, mpnw
    call mp_fixlocr (ra)
    call mp_fixlocr (rb)
    mpnwbt = min (ra%mpr(1), rb%mpr(1), mpwdsbt)
    mpnw = mpnwbt / mpnbt
    call mp_initvr (mp_bessel_k, mpnwbt)
    call mpbesselkr (ra%mpr(0:), rb%mpr(0:), mp_bessel_k%mpr(0:), mpnw)
    return
  end function

  function mp_bessel_kn (ia, rb)
    implicit none
    type (mp_real):: mp_bessel_kn
    integer, intent (in):: ia
    type (mp_real), intent (in):: rb
    integer mpnwbt, mpnw
    call mp_fixlocr (rb)
    mpnwbt = min (rb%mpr(1), mpwdsbt)
    mpnw = mpnwbt / mpnbt
    call mp_initvr (mp_bessel_kn, mpnwbt)
    call mpbesselknr (ia, rb%mpr(0:), mp_bessel_kn%mpr(0:), mpnw)
    return
  end function

  function mp_bessel_y (ra, rb)
    implicit none
    type (mp_real):: mp_bessel_y
    type (mp_real), intent (in):: ra, rb
    integer mpnwbt, mpnw
    call mp_fixlocr (ra)
    call mp_fixlocr (rb)
    mpnwbt = min (ra%mpr(1), rb%mpr(1), mpwdsbt)
    mpnw = mpnwbt / mpnbt
    call mp_initvr (mp_bessel_y, mpnwbt)
    call mpbesselyr (ra%mpr(0:), rb%mpr(0:), mp_bessel_y%mpr(0:), mpnw)
    return
  end function

  function mp_bessel_y0 (ra)
    implicit none
    type (mp_real):: mp_bessel_y0
    type (mp_real), intent (in):: ra
    type (mp_real) r1
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_bessel_y0, mpnwbt)
    call mpfry0 (mp_bessel_y0%mpr(1:), ra%mpr(1:), mprnd)
    return
  end function

  function mp_bessel_y1 (ra)
    implicit none
    type (mp_real):: mp_bessel_y1
    type (mp_real), intent (in):: ra
    type (mp_real) r1
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_bessel_y1, mpnwbt)
    call mpfry1 (mp_bessel_y1%mpr(1:), ra%mpr(1:), mprnd)
    return
  end function

  function mp_bessel_yn (ia, rb)
    implicit none
    type (mp_real):: mp_bessel_yn
    integer, intent (in):: ia
    type (mp_real), intent (in):: rb
    integer mpnwbt
    call mp_fixlocr (rb)
    mpnwbt = min (rb%mpr(1), mpwdsbt)
    call mp_initvr (mp_bessel_yn, mpnwbt)
    call mpfrbesselyn (mp_bessel_yn%mpr(1:), ia, rb%mpr(1:), mprnd)
    return
  end function

  subroutine mp_binmd (ra, db, ib)
    implicit none
    type (mp_real), intent (in):: ra
    real (mprknd), intent (out):: db
    integer, intent (out):: ib
    integer(8) ic8
    real (mprknd), external:: mpfrgetd2exp
    call mp_fixlocr (ra)
    db = 2.d0 * mpfrgetd2exp (ic8, ra%mpr(1:), mprnd)
    if (db == 0.d0) then
      ib = 0
    else
      ib = ic8 - 1
    endif
    return
  end subroutine

  function mp_ccos (za)
    implicit none
    type (mp_complex):: mp_ccos
    type (mp_complex), intent (in):: za
    integer l1, l2, l3, mpnwbt
    type (mp_real) r1
    type (mp_complex) z1, z2, z3
    call mp_fixlocz (za)
    l1 = za%mpc(0)
    mpnwbt = max (za%mpc(1), za%mpc(l1+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    l2 = mpwds6
    call mp_initvz (z1, mpnwbt)
    call mp_initvz (z2, mpnwbt)
    call mp_initvz (z3, mpnwbt)
    l3 = mpwds6
    call mp_initvz (mp_ccos, mpnwbt)
    call mpfrsetd (z1%mpc(1:), 1.d0, mprnd)
    call mpfrsetd (z1%mpc(l2+1:), 0.d0, mprnd)
    call mpfrset (z3%mpc(l2+1:), za%mpc(1:), mprnd)
    call mpfrset (z3%mpc(1:), za%mpc(l1+1:), mprnd)
    call mpfrneg (z3%mpc(1:), z3%mpc(1:), mprnd)
    z2 = mp_cexp (z3)
    z3 = mp_divzz (z1, z2)
    z1 = mp_addzz (z2, z3)
    call mpfrmuld (mp_ccos%mpc(1:), z1%mpc(1:), 0.5d0, mprnd)
    call mpfrmuld (mp_ccos%mpc(l3+1:), z1%mpc(l2+1:), 0.5d0, mprnd)
    return
  end function

  function mp_cexp (za)
    implicit none
    type (mp_complex):: mp_cexp
    type (mp_complex), intent (in):: za
    integer l1, l2, mpnwbt
    type (mp_real) r1
    call mp_fixlocz (za)
    l1 = za%mpc(0)
    mpnwbt = max (za%mpc(1), za%mpc(l1+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    call mp_initvr (r1, mpnwbt)
    l2 = mpwds6
    call mp_initvz (mp_cexp, mpnwbt)
    call mpfrsincos (mp_cexp%mpc(l2+1:), mp_cexp%mpc(1:), za%mpc(l1+1:), mprnd)
    call mpfrexp (r1%mpr(1:), za%mpc(1:), mprnd)
    call mpfrmul (mp_cexp%mpc(1:), mp_cexp%mpc(1:), r1%mpr(1:), mprnd)
    call mpfrmul (mp_cexp%mpc(l2+1:), mp_cexp%mpc(l2+1:), r1%mpr(1:), mprnd)
    return
  end function

  function mp_clog (za)
    implicit none
    type (mp_complex):: mp_clog
    type (mp_complex), intent (in):: za
    integer l1, l2, mpnwbt
    type (mp_real) r1, r2, r3
    call mp_fixlocz (za)
    l1 = za%mpc(0)
    mpnwbt = max (za%mpc(1), za%mpc(l1+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    call mp_initvr (r1, mpnwbt)
    call mp_initvr (r2, mpnwbt)
    call mp_initvr (r3, mpnwbt)
    l2 = mpwds6
    call mp_initvz (mp_clog, mpnwbt)
    call mpfrmul (r1%mpr(1:), za%mpc(1:), za%mpc(1:), mprnd)
    call mpfrmul (r2%mpr(1:), za%mpc(l1+1:), za%mpc(l1+1:), mprnd)
    call mpfradd (r3%mpr(1:), r1%mpr(1:), r2%mpr(1:), mprnd)
    call mpfrlog (r1%mpr(1:), r3%mpr(1:), mprnd)
    call mpfrmuld (mp_clog%mpc(1:), r1%mpr(1:), 0.5d0, mprnd)
    call mpfratan2 (mp_clog%mpc(l2+1:), za%mpc(l1+1:), za%mpc(1:), mprnd)
    return
  end function

  function mp_conjg (za)
    implicit none
    type (mp_complex):: mp_conjg
    type (mp_complex), intent (in):: za
    integer l1, l2, mpnwbt
    call mp_fixlocz (za)
    l1 = za%mpc(0)
    mpnwbt = max (za%mpc(1), za%mpc(l1+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    l2 = mpwds6
    call mp_initvz (mp_conjg, mpnwbt)
    call mpfrset (mp_conjg%mpc(1:), za%mpc(1:), mprnd)
    call mpfrset (mp_conjg%mpc(l2+1:), za%mpc(l1+1:), mprnd)
    call mpfrneg (mp_conjg%mpc(l2+1:), mp_conjg%mpc(l2+1:), mprnd)
    return
  end function

  function mp_cos (ra)
    implicit none
    type (mp_real):: mp_cos
    type (mp_real), intent (in):: ra
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_cos, mpnwbt)
    call mpfrcos (mp_cos%mpr(1:), ra%mpr(1:), mprnd)
    return
  end function

  function mp_cosh (ra)
    implicit none
    type (mp_real):: mp_cosh
    type (mp_real), intent (in):: ra
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_cosh, mpnwbt)
    call mpfrcosh (mp_cosh%mpr(1:), ra%mpr(1:), mprnd)
    return
  end function

  function mp_csin (za)
    implicit none
    type (mp_complex):: mp_csin
    type (mp_complex), intent (in):: za
    integer l1, l2, l3, mpnwbt
    type (mp_real) r1
    type (mp_complex) z1, z2, z3
    call mp_fixlocz (za)
    l1 = za%mpc(0)
    mpnwbt = max (za%mpc(1), za%mpc(l1+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    l2 = mpwds6
    call mp_initvz (z1, mpnwbt)
    call mp_initvz (z2, mpnwbt)
    call mp_initvz (z3, mpnwbt)
    l3 = mpwds6
    call mp_initvz (mp_csin, mpnwbt)
    call mpfrsetd (z1%mpc(1:), 1.d0, mprnd)
    call mpfrsetd (z1%mpc(l2+1:), 0.d0, mprnd)
    call mpfrset (z3%mpc(l2+1:), za%mpc(1:), mprnd)
    call mpfrset (z3%mpc(1:), za%mpc(l1+1:), mprnd)
    call mpfrneg (z3%mpc(1:), z3%mpc(1:), mprnd)
    z2 = mp_cexp (z3)
    z3 = mp_divzz (z1, z2)
    z1 = mp_subzz (z2, z3)
    call mpfrmuld (mp_csin%mpc(l3+1:), z1%mpc(1:), 0.5d0, mprnd)
    call mpfrmuld (mp_csin%mpc(1:), z1%mpc(l2+1:), 0.5d0, mprnd)
    call mpfrneg (mp_csin%mpc(l3+1:), mp_csin%mpc(l3+1:), mprnd)
    return
  end function

  function mp_csqrt (za)
    implicit none
    type (mp_complex):: mp_csqrt
    type (mp_complex), intent (in):: za
    integer l1, l2, mpnwbt
    type (mp_real) r1, r2, r3
    call mp_fixlocz (za)
    l1 = za%mpc(0)
    mpnwbt = max (za%mpc(1), za%mpc(l1+1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    call mp_initvr (r1, mpnwbt)
    call mp_initvr (r2, mpnwbt)
    call mp_initvr (r3, mpnwbt)
    l2 = mpwds6
    call mp_initvz (mp_csqrt, mpnwbt)
    call mpfrmul (r1%mpr(1:), za%mpc(1:), za%mpc(1:), mprnd)
    call mpfrmul (r2%mpr(1:), za%mpc(l1+1:), za%mpc(l1+1:), mprnd)
    call mpfradd (r3%mpr(1:), r1%mpr(1:), r2%mpr(1:), mprnd)
    call mpfrsqrt (r1%mpr(1:), r3%mpr(1:), mprnd)

    if (za%mpc(2) == 1) then
      call mpfradd (r2%mpr(1:), r1%mpr(1:), za%mpc(1:), mprnd)
      call mpfrsqrt (r3%mpr(1:), r2%mpr(1:), mprnd)
      call mpfrset (mp_csqrt%mpc(1:), r3%mpr(1:), mprnd)
      call mpfrdiv (mp_csqrt%mpc(l2+1:), za%mpc(l1+1:), r3%mpr(1:), mprnd)
    else
      call mpfrsub (r2%mpr(1:), r1%mpr(1:), za%mpc(1:), mprnd)
      call mpfrsqrt (r3%mpr(1:), r2%mpr(1:), mprnd)
      call mpfrabs (r2%mpr(1:), za%mpc(l1+1:), mprnd)
      call mpfrdiv (mp_csqrt%mpc(1:), r2%mpr(1:), r3%mpr(1:), mprnd)
      call mpfrset (mp_csqrt%mpc(l2+1:), r3%mpr(1:), mprnd)
      if (za%mpc(l1+2) /= 1) &
         call mpfrneg (mp_csqrt%mpc(l2+1:), mp_csqrt%mpc(l2+1:), mprnd)
    endif

    call mpfrsetd (r1%mpr(1:), 2.d0, mprnd)
    call mpfrsqrt (r2%mpr(1:), r1%mpr(1:), mprnd)
    call mpfrdiv (mp_csqrt%mpc(1:), mp_csqrt%mpc(1:), r2%mpr(1:), mprnd)
    call mpfrdiv (mp_csqrt%mpc(l2+1:), mp_csqrt%mpc(l2+1:), r2%mpr(1:), mprnd)
    return
  end function

  subroutine mp_cssh (ra, rb, rc)
    implicit none
    type (mp_real), intent (in):: ra
    type (mp_real), intent (out):: rb, rc
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (rb, mpnwbt)
    call mp_initvr (rc, mpnwbt)
    call mpfrsinhcosh (rc%mpr(1:), rb%mpr(1:), ra%mpr(1:), mprnd)
    return
  end subroutine

  subroutine mp_cssn (ra, rb, rc)
    implicit none
    type (mp_real), intent (in):: ra
    type (mp_real), intent (out):: rb, rc
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (rb, mpnwbt)
    call mp_initvr (rc, mpnwbt)
    call mpfrsincos (rc%mpr(1:), rb%mpr(1:), ra%mpr(1:), mprnd)
    return
  end subroutine

  function mp_dctoz (dca, iprec)
    implicit none
    type (mp_complex):: mp_dctoz
    complex (kind(0.d0)), intent (in):: dca
    integer l1, l2, mpnwbt

!>  In variant #1, uncomment these lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnwbt = mpnbt * mp_setwp (iprec)
!    else
!      mpnwbt = mpnbt * mpwds
!    endif
!  Otherwise in variant #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnwbt = mpnbt * mp_setwp (iprec)
!>>
    call mp_checkdp (dble (dca))
    call mp_checkdp (aimag (dca))
    l2 = mpwds6
    call mp_initvz (mp_dctoz, mpnwbt)
    call mpfrsetd (mp_dctoz%mpc(1:), dble (dca), mprnd)
    call mpfrsetd (mp_dctoz%mpc(l2+1:), aimag (dca), mprnd)
    return
  end function

  function mp_dctoz2 (dca, iprec)
    implicit none
    type (mp_complex):: mp_dctoz2
    complex (kind(0.d0)), intent (in):: dca
    integer l1, l2, mpnwbt

!>  In variant #1, uncomment these lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnwbt = mpnbt * mp_setwp (iprec)
!    else
!      mpnwbt = mpnbt * mpwds
!    endif
!  Otherwise in variant #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnwbt = mpnbt * mp_setwp (iprec)
!>>
    l2 = mpwds6
    call mp_initvz (mp_dctoz2, mpnwbt)
    call mpfrsetd (mp_dctoz2%mpc(1:), dble (dca), mprnd)
    call mpfrsetd (mp_dctoz2%mpc(l2+1:), aimag (dca), mprnd)
    return
  end function

  subroutine mp_decmd (ra, db, ib)
    implicit none
    type (mp_real), intent (in):: ra
    real (mprknd), intent (out):: db
    integer, intent (out):: ib
    integer(8) ic8
    real (mprknd), external:: mpfrgetd2exp
    real (mprknd), parameter:: alg102 = 0.301029995663981195d0
    real (mprknd) dt1, dt2
    call mp_fixlocr (ra)
    dt1 = 2.d0 * mpfrgetd2exp (ic8, ra%mpr(1:), mprnd)
    if (dt1 /= 0.d0) then
      dt2 = alg102 * (ic8 - 1) + log10 (abs (dt1))
      ib = dt2
      if (dt2 < 0.d0) ib = ib - 1
      db = sign (10.d0 ** (dt2 - ib), dt1)
    else
      db = 0.d0
      ib = 0
    endif
    return
  end subroutine

  function mp_digamma (ra)
    implicit none
    type (mp_real):: mp_digamma
    type (mp_real), intent (in):: ra
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_digamma, mpnwbt)
    call mpfrdigamma (mp_digamma%mpr(1:), ra%mpr(1:), mprnd)
    return
  end function

  function mp_digamma_be (nb, rb, rc)
    implicit none
    integer, intent (in):: nb
    type (mp_real):: mp_digamma_be
    type (mp_real), intent (in):: rb(nb), rc
    integer n1, mpnw, mpnwbt
    mpnwbt = max (int (rb(1)%mpr(1)), int (rc%mpr(1)))
    mpnwbt = min (mpnwbt, mpwdsbt)
    mpnw = mpnwbt / mpnbt
    n1 = mpwds
    call mp_fixlocr (rc)
    call mp_initvr (mp_digamma_be, mpnwbt)
    call mpdigammabe (n1, nb, rb(1)%mpr(0:), rc%mpr(0:), mp_digamma_be%mpr, mpnw)
    return
  end function

  function mp_dtor (da, iprec)
    implicit none
    type (mp_real):: mp_dtor
    real (mprknd), intent (in):: da
    integer mpnwbt
    real (mprknd), external:: mpmask13

!>  In variant #1, uncomment these lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnwbt = mpnbt * mp_setwp (iprec)
!    else
!      mpnwbt = mpnbt * mpwds
!    endif
!  Otherwise in variant #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnwbt = mpnbt * mp_setwp (iprec)
!>>
    call mp_checkdp (da)
    call mp_initvr (mp_dtor, mpnwbt)
    call mpfrsetd (mp_dtor%mpr(1:), da, mprnd)
    return
  end function

  function mp_dtor2 (da, iprec)
    implicit none
    type (mp_real):: mp_dtor2
    real (mprknd), intent (in):: da
    integer mpnwbt

!>  In variant #1, uncomment these lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnwbt = mpnbt * mp_setwp (iprec)
!    else
!      mpnwbt = mpnbt * mpwds
!    endif
!  Otherwise in variant #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnwbt = mpnbt * mp_setwp (iprec)
!>>
    call mp_initvr (mp_dtor2, mpnwbt)
    call mpfrsetd (mp_dtor2%mpr(1:), da, mprnd)
    return
  end function

  subroutine mp_eform (ra, nb, nd, b)
    implicit none
    type (mp_real), intent (in):: ra
    integer, intent (in):: nb, nd
    character(1), intent (out):: b(nb)
    character(1) b1(nd+8)
    integer ix, i, j
    integer(8) iexp
    character(16) str16

    call mp_fixlocr (ra)

!  Check for overflow of field.
    if (nb < nd + 20) then
      do i = 1, nb
        b(i) = '*'
      enddo
      return
    endif

!  Call mpfrgetstr to convert number.

    call mpfrgetstr (ra%mpr(1:), b1, nd, iexp, mprnd)

    if (b1(1) == '-') then
      b(1) = '-'
      b(2) = b1(2)
      b(3) = '.'
      do i = 1, nd - 2
        b(i+3) = b1(i+2)
      enddo
      ix = nd + 1
    else
      b(1) = b1(1)
      b(2) = '.'
      do i = 1, nd - 2
        b(i+2) = b1(i+1)
      enddo
      ix = nd
    endif
!  Insert exponent.
    b(ix+1) = 'e'
    ix = ix + 1
    write (str16, '(i16)') iexp - 1
    do i = 1, 16
      if (str16(i:i) /= ' ') goto 100
    enddo
100 continue
    do j = 1, 16 - i + 1
      b(ix+j) = str16(i+j-1:i+j-1)
    enddo
    ix = ix + 16 - i + 1
    do j = ix + 1, nb
      b(j) = ' '
    enddo

    return
  end subroutine

  function mp_egamma (iprec)
    implicit none
    type (mp_real):: mp_egamma
    integer mpnwbt

!>  In variant #1, uncomment these lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnwbt = mpnbt * mp_setwp (iprec)
!    else
!      mpnwbt = mpnbt * mpwds
!    endif
!  Otherwise in variant #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnwbt = mpnbt * mp_setwp (iprec)
!>>
    call mp_initvr (mp_egamma, mpnwbt)
    if (mpwprecr (mpegammacon) < mpnwbt / mpnbt) then
      write (mpldb, 1) mpnwbt / mpnbt
1     format ('*** MP_EGAMMA: Egamma must be precomputed to precision',i9,' words'/ &
      'by calling mpinit. See documentation for details.')
      call mp_abrt (53)
    endif
    call mpfixlocr (mpegammacon)
    call mpfrset (mp_egamma%mpr(1:), mpegammacon(1), mprnd)
  end function

  function mp_expint (ra)
    implicit none
    type (mp_real):: mp_expint
    type (mp_real), intent (in):: ra
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_expint, mpnwbt)
    call mpfreint (mp_expint%mpr(1:), ra%mpr(1:), mprnd)
    return
  end function

  function mp_erf (ra)
    implicit none
    type (mp_real):: mp_erf
    type (mp_real), intent (in):: ra
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_erf, mpnwbt)
    call mpfrerf (mp_erf%mpr(1:), ra%mpr(1:), mprnd)
    return
  end function

  function mp_erfc (ra)
    implicit none
    type (mp_real):: mp_erfc
    type (mp_real), intent (in):: ra
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_erfc, mpnwbt)
    call mpfrerfc (mp_erfc%mpr(1:), ra%mpr(1:), mprnd)
    return
  end function

  function mp_exp (ra)
    implicit none
    type (mp_real):: mp_exp
    type (mp_real), intent (in):: ra
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_exp, mpnwbt)
    call mpfrexp (mp_exp%mpr(1:), ra%mpr(1:), mprnd)
    return
  end function

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

    call mpfrgetstr (ra%mpr(1:), b1, nb, iexp, mprnd)

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

  function mp_gamma (ra)
    implicit none
    type (mp_real):: mp_gamma
    type (mp_real), intent (in):: ra
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_gamma, mpnwbt)
    call mpfrgamma (mp_gamma%mpr(1:), ra%mpr(1:), mprnd)
    return
  end function

  function mp_hurwitz_zetan (ia, rb)
    implicit none
    type (mp_real):: mp_hurwitz_zetan
    integer, intent (in):: ia
    type (mp_real), intent (in):: rb
    integer mpnwbt, mpnw
    call mp_fixlocr (rb)
    mpnwbt = min (rb%mpr(1), mpwdsbt)
    mpnw = mpnwbt / mpnbt
    call mp_initvr (mp_hurwitz_zetan, mpnwbt)
    call mphurwitzzetan (ia, rb%mpr(0:), mp_hurwitz_zetan%mpr(0:), mpnw)
    return
  end function

  function mp_hurwitz_zetan_be (nb, rb, ia, rc)
    implicit none
    type (mp_real):: mp_hurwitz_zetan_be
    integer, intent (in):: nb, ia
    type (mp_real), intent (in):: rb(nb), rc
    integer mpnwbt, mpnw, n1
    call mp_fixlocr (rc)
    mpnwbt = min (rc%mpr(1), mpwdsbt)
    mpnw = mpnwbt / mpnbt
    n1 = mpwds
    call mp_initvr (mp_hurwitz_zetan_be, mpnwbt)
    call mphurwitzzetanbe (n1, nb, rb(1)%mpr(0:), ia, rc%mpr(0:), &
      mp_hurwitz_zetan_be%mpr(0:), mpnw)
    return
  end function

  function mp_hypergeom_pfq (np, nq, aa, bb, xx)
    implicit none
    type (mp_real):: mp_hypergeom_pfq
    integer, intent (in):: np, nq
    type (mp_real), intent (in):: aa(np), bb(nq), xx
    integer k, mpnwbt, mpnw
    call mp_fixlocr (xx)
    mpnwbt = min (xx%mpr(1), mpwdsbt)
    mpnw = mpnwbt / mpnbt
    call mp_initvr (mp_hypergeom_pfq, mpnwbt)
    call mphypergeompfq (np, nq, mpwds, aa(1)%mpr(0:), bb(1)%mpr(0:), &
      xx%mpr(0:), mp_hypergeom_pfq%mpr(0:), mpnw)
    return
  end function

  function mp_hypot (ra, rb)
    implicit none
    type (mp_real):: mp_hypot
    type (mp_real), intent (in):: ra, rb
    integer mpnwbt
    call mp_fixlocr (ra)
    call mp_fixlocr (rb)
    mpnwbt = max (ra%mpr(1), rb%mpr(1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    call mp_initvr (mp_hypot, mpnwbt)
    call mpfrhypot (mp_hypot%mpr(1:), ra%mpr(1:), rb%mpr(1:), mprnd)
    return
  end function

  function mp_incgamma (ra, rb)
    implicit none
    type (mp_real):: mp_incgamma
    type (mp_real), intent (in):: ra, rb
    integer mpnwbt
    call mp_fixlocr (ra)
    call mp_fixlocr (rb)
    mpnwbt = max (ra%mpr(1), rb%mpr(1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    call mp_initvr (mp_incgamma, mpnwbt)
    call mpfrgammainc (mp_incgamma%mpr(1:), ra%mpr(1:), rb%mpr(1:), mprnd)
    return
  end function

  subroutine mp_init (iprec)
    implicit none
    integer mpnwbt

!>  In variant #1, uncomment these lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnwbt = mpnbt * mp_setwp (iprec)
!    else
!      mpnwbt = mpnbt * mpwds
!    endif
!  Otherwise in variant #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnwbt = mpnbt * mp_setwp (iprec)
!>>
    if (mpnwbt / mpnbt > mpl2pi) then
      call mpinitwds (mplog2con, mpnwbt / mpnbt)
      call mpinitwds (mppicon, mpnwbt / mpnbt)
      call mpinitwds (mpegammacon, mpnwbt / mpnbt)
      call mpfrconstlog2 (mplog2con(1), mprnd)
      call mpfrconstpi (mppicon(1), mprnd)
      call mpfrconsteuler (mpegammacon(1), mprnd)
    endif
    return
  end subroutine

  function mp_log (ra)
    implicit none
    type (mp_real):: mp_log
    type (mp_real), intent (in):: ra
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_log, mpnwbt)
    call mpfrlog (mp_log%mpr(1:), ra%mpr(1:), mprnd)
    return
  end function

  function mp_log10 (ra)
    implicit none
    type (mp_real):: mp_log10
    type (mp_real), intent (in):: ra
    type (mp_real) t1, t2
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (t1, mpnwbt)
    call mp_initvr (t2, mpnwbt)
    call mp_initvr (mp_log10, mpnwbt)
    call mpfrsetd (t1%mpr(1:), 10.d0, mprnd)
    call mpfrlog (t2%mpr(1:), t1%mpr(1:), mprnd)
    call mpfrlog (t1%mpr(1:), ra%mpr(1:), mprnd)
    call mpfrdiv (mp_log10%mpr(1:), t1%mpr(1:), t2%mpr(1:), mprnd)
    return
  end function

  function mp_log_gamma (ra)
    implicit none
    type (mp_real):: mp_log_gamma
    type (mp_real), intent (in):: ra
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_log_gamma, mpnwbt)
    call mpfrlngamma (mp_log_gamma%mpr(1:), ra%mpr(1:), mprnd)
    return
  end function

  function mp_log2 (iprec)
    implicit none
    type (mp_real):: mp_log2
    type (mp_real) qpi
    integer mpnwbt

!>  In variant #1, uncomment these lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnwbt = mpnbt * mp_setwp (iprec)
!    else
!      mpnwbt = mpnbt * mpwds
!    endif
!  Otherwise in variant #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnwbt = mpnbt * mp_setwp (iprec)
!>>
    call mp_initvr (mp_log2, mpnwbt)
    if (mpwprecr (mplog2con) < mpnwbt / mpnbt) then
      write (mpldb, 1) mpnwbt / mpnbt
1     format ('*** MP_LOG2: Log(2) must be precomputed to precision',i9,' words'/ &
      'by calling mpinit. See documentation for details.')
      call mp_abrt (53)
    endif
    call mpfixlocr (mplog2con)
    call mpfrset (mp_log2%mpr(1:), mplog2con(1), mprnd)
  end function

  function mp_max (ra, rb, rc)
    implicit none
    type (mp_real):: mp_max
    type (mp_real), intent (in):: ra, rb
    type (mp_real), optional, intent (in):: rc
    integer mpnwbt
    call mp_fixlocr (ra)
    call mp_fixlocr (rb)
    mpnwbt = max (ra%mpr(1), rb%mpr(1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    call mp_initvr (mp_max, mpnwbt)
    call mpfrmax (mp_max%mpr(1:), ra%mpr(1:), rb%mpr(1:), mprnd)
    if (present (rc)) then
    call mp_fixlocr (rc)
      call mpfrmax (mp_max%mpr(1:), mp_max%mpr(1:), rc%mpr(1:), mprnd)
    endif
    return
  end function

  function mp_min (ra, rb, rc)
    implicit none
    type (mp_real):: mp_min
    type (mp_real), intent (in):: ra, rb
    type (mp_real), optional, intent (in):: rc
    integer mpnwbt
    call mp_fixlocr (ra)
    call mp_fixlocr (rb)
    mpnwbt = max (ra%mpr(1), rb%mpr(1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    call mp_initvr (mp_min, mpnwbt)
    call mpfrmin (mp_min%mpr(1:), ra%mpr(1:), rb%mpr(1:), mprnd)
    if (present (rc)) then
      call mp_fixlocr (rc)
      call mpfrmin (mp_min%mpr(1:), mp_min%mpr(1:), rc%mpr(1:), mprnd)
    endif
  end function

  function mp_mod (ra, rb)
    implicit none
    type (mp_real):: mp_mod
    type (mp_real), intent (in):: ra, rb
    integer mpnwbt
    call mp_fixlocr (ra)
    call mp_fixlocr (rb)
    mpnwbt = max (ra%mpr(1), rb%mpr(1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    call mp_initvr (mp_mod, mpnwbt)
    call mpfrfmod (mp_mod%mpr(1:), ra%mpr(1:), rb%mpr(1:), mprnd)
    return
  end function

  function mp_mtor (ra, iprec)
    implicit none
    type (mp_real):: mp_mtor
    type (mp_realm), intent (in):: ra
    integer mpnwbt

!>  In variant #1, uncomment these lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnwbt = mpnbt * mp_setwp (iprec)
!    else
!      mpnwbt = mpnbt * mpwds
!    endif
!  Otherwise in variant #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnwbt = mpnbt * mp_setwp (iprec)
!>>
    call mp_fixlocrm (ra)
    call mp_initvr (mp_mtor, mpnwbt)
    call mpfrset (mp_mtor%mpr(1:), ra%mpr(1:), mprnd)
    return
  end function

  function mp_nrt (ra, ib)
    implicit none
    type (mp_real):: mp_nrt
    type (mp_real), intent (in):: ra
    integer, intent (in):: ib
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_nrt, mpnwbt)
    call mpfrroot (mp_nrt%mpr(1:), ra%mpr(1:), ib, mprnd)
    return
  end function

  function mp_pi (iprec)
    implicit none
    type (mp_real):: mp_pi
    integer mpnwbt

!>  In variant #1, uncomment these lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnwbt = mpnbt * mp_setwp (iprec)
!    else
!      mpnwbt = mpnbt * mpwds
!    endif
!  Otherwise in variant #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnwbt = mpnbt * mp_setwp (iprec)
!>>
    call mp_initvr (mp_pi, mpnwbt)
    if (mpwprecr (mppicon) < mpnwbt / mpnbt) then
      write (mpldb, 1) mpnwbt / mpnbt
1     format ('*** MP_PI: Pi must be precomputed to precision',i9,' words'/ &
      'by calling mpinit. See documentation for details.')
      call mp_abrt (53)
    endif
    call mpfixlocr (mppicon)
    call mpfrset (mp_pi%mpr(1:), mppicon(1), mprnd)
  end function

  function mp_polygamma (ia, rb)
    implicit none
    type (mp_real):: mp_polygamma
    integer, intent (in):: ia
    type (mp_real), intent (in):: rb
    integer mpnwbt, mpnw
    call mp_fixlocr (rb)
    mpnwbt = min (rb%mpr(1), mpwdsbt)
    mpnw = mpnwbt / mpnbt
    call mp_initvr (mp_polygamma, mpnwbt)
    call mppolygamma (ia, rb%mpr(0:), mp_polygamma%mpr(0:), mpnw)
    return
  end function

  function mp_polygamma_be (nb, rb, ia, rc)
    implicit none
    type (mp_real):: mp_polygamma_be
    integer, intent (in):: nb, ia
    type (mp_real), intent (in):: rb(nb), rc
    integer mpnwbt, mpnw, n1
    call mp_fixlocr (rc)
    mpnwbt = min (rc%mpr(1), mpwdsbt)
    mpnw = mpnwbt / mpnbt
    n1 = mpwds
    call mp_initvr (mp_polygamma_be, mpnwbt)
    call mppolygammabe (n1, nb, rb(1)%mpr(0:), ia, rc%mpr(0:), &
      mp_polygamma_be%mpr(0:), mpnw)
    return
  end function

  subroutine mp_polylog_ini (nb, rb, iprec)
    implicit none
    integer, intent (in):: nb
    type (mp_real), intent (out):: rb(abs(nb))
    integer i, n1, n2, mpnwbt, mpnw

!>  In variant #1, uncomment these lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnwbt = mpnbt * mp_setwp (iprec)
!    else
!      mpnwbt = mpnbt * mpwds
!    endif
!  Otherwise in variant #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnwbt = mpnbt * mp_setwp (iprec)
!>>
    n1 = mpwds
    n2 = abs (nb)
    mpnw = mpnwbt / mpnbt
    do i = 1, n2
      call mp_initvr (rb(i), mpnwbt)
    enddo
    call mppolylogini (n1, n2, rb(1)%mpr(0), mpnw)
    return
  end subroutine

  function mp_polylog_neg (nb, rb, rc)
    implicit none
    integer, intent (in):: nb
    type (mp_real):: mp_polylog_neg
    type (mp_real), intent (in):: rb(abs(nb)), rc
    integer n1, mpnw, mpnwbt
    mpnwbt = max (int (rb(1)%mpr(1)), int (rc%mpr(1)))
    mpnwbt = min (mpnwbt, mpwdsbt)
    mpnw = mpnwbt / mpnbt
    n1 = mpwds
    call mp_fixlocr (rc)
    call mp_initvr (mp_polylog_neg, mpnwbt)
    call mppolylogneg (n1, nb, rb(1)%mpr(0:), rc%mpr(0:), mp_polylog_neg%mpr, mpnw)
    return
  end function

  function mp_polylog_pos (ia, rb)
    implicit none
    type (mp_real):: mp_polylog_pos
    integer, intent (in):: ia
    type (mp_real), intent (in):: rb
    integer mpnwbt, mpnw
    call mp_fixlocr (rb)
    mpnwbt = min (rb%mpr(1), mpwdsbt)
    mpnw = mpnwbt / mpnbt
    call mp_initvr (mp_polylog_pos, mpnwbt)
    call mppolylogpos (ia, rb%mpr(0:), mp_polylog_pos%mpr(0:), mpnw)
    return
  end function

  function mp_prodd (ra, db)
    implicit none
    type (mp_real):: mp_prodd
    type (mp_real), intent (in):: ra
    real (mprknd), intent (in):: db
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_prodd, mpnwbt)
    call mpfrmuld (mp_prodd%mpr(1:), ra%mpr(1:), db, mprnd)
  end function

  function mp_prodq (ra, qb)
    implicit none
    type (mp_real):: mp_prodq
    type (mp_real), intent (in):: ra
    real (max (mprknd2, kind (1.))), intent (in):: qb
    real (max (mprknd2, kind (1.))):: q1
    real (mprknd) d1, d2, d3
    type (mp_real) r1, r2, r3
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_prodq, mpnwbt)
    call mp_initvr (r1, mpnwbt)
    call mp_initvr (r2, mpnwbt)
    call mp_initvr (r3, mpnwbt)
    d1 = qb
    q1 = qb - d1
    d2 = q1
    d3 = q1 - d2
    call mpfrsetd (r1%mpr(1:), d1, mprnd)
    call mpfraddd (r2%mpr(1:), r1%mpr(1:), d2, mprnd)
    call mpfraddd (r3%mpr(1:), r2%mpr(1:), d3, mprnd)
    call mpfrmul (mp_prodq%mpr(1:), ra%mpr(1:), r3%mpr(1:), mprnd)
  end function

function mp_qtor (qa, iprec)
    implicit none
    type (mp_real):: mp_qtor
    integer mpnwbt
    real (max (mprknd2, kind (1.))), intent (in):: qa
    real (max (mprknd2, kind (1.))):: q1
    real (mprknd) d1, d2, d3
    type (mp_real) r1, r2

!>  In variant #1, uncomment these lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnwbt = mpnbt * mp_setwp (iprec)
!    else
!      mpnwbt = mpnbt * mpwds
!    endif
!  Otherwise in variant #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnwbt = mpnbt * mp_setwp (iprec)
!>>
    call mp_checkqp (qa)
    call mp_initvr (mp_qtor, mpnwbt)
    call mp_initvr (r1, mpnwbt)
    call mp_initvr (r2, mpnwbt)
    d1 = qa
    q1 = qa - d1
    d2 = q1
    d3 = q1 - d2
    call mpfrsetd (r1%mpr(1:), d1, mprnd)
    call mpfraddd (r2%mpr(1:), r1%mpr(1:), d2, mprnd)
    call mpfraddd (mp_qtor%mpr(1:), r2%mpr(1:), d3, mprnd)
    return
  end function

  function mp_qtor2 (qa, iprec)
    implicit none
    type (mp_real):: mp_qtor2
    integer mpnwbt
    real (max (mprknd2, kind (1.))), intent (in):: qa
    real (max (mprknd2, kind (1.))):: q1
    real (mprknd) d1, d2, d3
    type (mp_real) r1, r2

!>  In variant #1, uncomment these lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnwbt = mpnbt * mp_setwp (iprec)
!    else
!      mpnwbt = mpnbt * mpwds
!    endif
!  Otherwise in variant #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnwbt = mpnbt * mp_setwp (iprec)
!>>
    call mp_initvr (mp_qtor2, mpnwbt)
    call mp_initvr (r1, mpnwbt)
    call mp_initvr (r2, mpnwbt)
    d1 = qa
    q1 = qa - d1
    d2 = q1
    d3 = q1 - d2
    call mpfrsetd (r1%mpr(1:), d1, mprnd)
    call mpfraddd (r2%mpr(1:), r1%mpr(1:), d2, mprnd)
    call mpfraddd (mp_qtor2%mpr(1:), r2%mpr(1:), d3, mprnd)
    return
  end function

  function mp_quotd (ra, db)
    implicit none
    type (mp_real):: mp_quotd
    type (mp_real), intent (in):: ra
    real (mprknd), intent (in):: db
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_quotd, mpnwbt)
    call mpfrdivd (mp_quotd%mpr(1:), ra%mpr(1:), db, mprnd)
  end function

  function mp_quotq (ra, qb)
    implicit none
    type (mp_real):: mp_quotq
    type (mp_real), intent (in):: ra
    real (max (mprknd2, kind (1.))), intent (in):: qb
    real (max (mprknd2, kind (1.))):: q1
    real (mprknd) d1, d2, d3
    type (mp_real) r1, r2, r3
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_quotq, mpnwbt)
    call mp_initvr (r1, mpnwbt)
    call mp_initvr (r2, mpnwbt)
    call mp_initvr (r3, mpnwbt)
    d1 = qb
    q1 = qb - d1
    d2 = q1
    d3 = q1 - d2
    call mpfrsetd (r1%mpr(1:), d1, mprnd)
    call mpfraddd (r2%mpr(1:), r1%mpr(1:), d2, mprnd)
    call mpfraddd (r3%mpr(1:), r2%mpr(1:), d3, mprnd)
    call mpfrdiv (mp_quotq%mpr(1:), ra%mpr(1:), r3%mpr(1:), mprnd)
  end function

  function mp_rand (ra)
    implicit none
    type (mp_real):: mp_rand
    type (mp_real), intent (in):: ra
    type (mp_real) r1, r2
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_rand, mpnwbt)
    call mp_initvr (r1, mpnwbt)
    call mp_initvr (r2, mpnwbt)
    call mpfrmuld (r1%mpr(1:), ra%mpr(1:), mprandx, mprnd)
    call mpfrtrunc (r2%mpr(1:), r1%mpr(1:))
    call mpfrsetd (r1%mpr(1:), mprandx, mprnd)
    call mpfrfms (mp_rand%mpr(1:), ra%mpr(1:), r1%mpr(1:), r2%mpr(1:), mprnd)
    return
  end function

!   Five variations of read are necessary due to Fortran rules about optional arguments.

  subroutine mp_readr1 (iu, r1, iprec)
    implicit none
    integer, intent (in):: iu
    type (mp_real), intent (out):: r1
    integer mpnw, mpnwbt

!>  In variant #1, uncomment these lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnwbt = mpnbt * mp_setwp (iprec)
!    else
!      mpnwbt = mpnbt * mpwds
!    endif
!  Otherwise in variant #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnwbt = mpnbt * mp_setwp (iprec)
!>>
    mpnw = mpnwbt / mpnbt
    call mp_inpr (iu, r1, mpnw)
    return
  end subroutine

  subroutine mp_readr2 (iu, r1, r2, iprec)
    implicit none
    integer, intent (in):: iu
    type (mp_real), intent (out):: r1, r2
    integer mpnw, mpnwbt

!>  In variant #1, uncomment these lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnwbt = mpnbt * mp_setwp (iprec)
!    else
!      mpnwbt = mpnbt * mpwds
!    endif
!  Otherwise in variant #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnwbt = mpnbt * mp_setwp (iprec)
!>>
    mpnw = mpnwbt / mpnbt
    call mp_inpr (iu, r1, mpnw)
    call mp_inpr (iu, r2, mpnw)
    return
  end subroutine

  subroutine mp_readr3 (iu, r1, r2, r3, iprec)
    implicit none
    integer, intent (in):: iu
    type (mp_real), intent (out):: r1, r2, r3
    integer mpnw, mpnwbt

!>  In variant #1, uncomment these lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnwbt = mpnbt * mp_setwp (iprec)
!    else
!      mpnwbt = mpnbt * mpwds
!    endif
!  Otherwise in variant #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnwbt = mpnbt * mp_setwp (iprec)
!>>
    mpnw = mpnwbt / mpnbt
    call mp_inpr (iu, r1, mpnw)
    call mp_inpr (iu, r2, mpnw)
    call mp_inpr (iu, r3, mpnw)
    return
  end subroutine

  subroutine mp_readr4 (iu, r1, r2, r3, r4, iprec)
    implicit none
    integer, intent (in):: iu
    type (mp_real), intent (out):: r1, r2, r3, r4
    integer mpnw, mpnwbt

!>  In variant #1, uncomment these lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnwbt = mpnbt * mp_setwp (iprec)
!    else
!      mpnwbt = mpnbt * mpwds
!    endif
!  Otherwise in variant #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnwbt = mpnbt * mp_setwp (iprec)
!>>
    mpnw = mpnwbt / mpnbt
    call mp_inpr (iu, r1, mpnw)
    call mp_inpr (iu, r2, mpnw)
    call mp_inpr (iu, r3, mpnw)
    call mp_inpr (iu, r4, mpnw)
    return
  end subroutine

  subroutine mp_readr5 (iu, r1, r2, r3, r4, r5, iprec)
    implicit none
    integer, intent (in):: iu
    type (mp_real), intent (out):: r1, r2, r3, r4, r5
    integer mpnw, mpnwbt

!>  In variant #1, uncomment these lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnwbt = mpnbt * mp_setwp (iprec)
!    else
!      mpnwbt = mpnbt * mpwds
!    endif
!  Otherwise in variant #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnwbt = mpnbt * mp_setwp (iprec)
!>>
    mpnw = mpnwbt / mpnbt
    call mp_inpr (iu, r1, mpnw)
    call mp_inpr (iu, r2, mpnw)
    call mp_inpr (iu, r3, mpnw)
    call mp_inpr (iu, r4, mpnw)
    call mp_inpr (iu, r5, mpnw)
    return
  end subroutine

  subroutine mp_readz1 (iu, z1, iprec)
    implicit none
    integer, intent (in):: iu
    type (mp_complex), intent (out):: z1
    integer mpnw, mpnwbt

!>  In variant #1, uncomment these lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnwbt = mpnbt * mp_setwp (iprec)
!    else
!      mpnwbt = mpnbt * mpwds
!    endif
!  Otherwise in variant #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnwbt = mpnbt * mp_setwp (iprec)
!>>
    mpnw = mpnwbt / mpnbt
    call mp_inpz (iu, z1, mpnw)
    return
  end subroutine

  subroutine mp_readz2 (iu, z1, z2, iprec)
    implicit none
    integer, intent (in):: iu
    type (mp_complex), intent (out):: z1, z2
    integer mpnw, mpnwbt

!>  In variant #1, uncomment these lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnwbt = mpnbt * mp_setwp (iprec)
!    else
!      mpnwbt = mpnbt * mpwds
!    endif
!  Otherwise in variant #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnwbt = mpnbt * mp_setwp (iprec)
!>>
    mpnw = mpnwbt / mpnbt
    call mp_inpz (iu, z1, mpnw)
    call mp_inpz (iu, z2, mpnw)
    return
  end subroutine

  subroutine mp_readz3 (iu, z1, z2, z3, iprec)
    implicit none
    integer, intent (in):: iu
    type (mp_complex), intent (out):: z1, z2, z3
    integer mpnw, mpnwbt

!>  In variant #1, uncomment these lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnwbt = mpnbt * mp_setwp (iprec)
!    else
!      mpnwbt = mpnbt * mpwds
!    endif
!  Otherwise in variant #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnwbt = mpnbt * mp_setwp (iprec)
!>>
    mpnw = mpnwbt / mpnbt
    call mp_inpz (iu, z1, mpnw)
    call mp_inpz (iu, z2, mpnw)
    call mp_inpz (iu, z3, mpnw)
    return
  end subroutine

  subroutine mp_readz4 (iu, z1, z2, z3, z4, iprec)
    implicit none
    integer, intent (in):: iu
    type (mp_complex), intent (out):: z1, z2, z3, z4
    integer mpnw, mpnwbt

!>  In variant #1, uncomment these lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnwbt = mpnbt * mp_setwp (iprec)
!    else
!      mpnwbt = mpnbt * mpwds
!    endif
!  Otherwise in variant #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnwbt = mpnbt * mp_setwp (iprec)
!>>
    mpnw = mpnwbt / mpnbt
    call mp_inpz (iu, z1, mpnw)
    call mp_inpz (iu, z2, mpnw)
    call mp_inpz (iu, z3, mpnw)
    call mp_inpz (iu, z4, mpnw)
    return
  end subroutine

  subroutine mp_readz5 (iu, z1, z2, z3, z4, z5, iprec)
    implicit none
    integer, intent (in):: iu
    type (mp_complex), intent (out):: z1, z2, z3, z4, z5
    integer mpnw, mpnwbt

!>  In variant #1, uncomment these lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnwbt = mpnbt * mp_setwp (iprec)
!    else
!      mpnwbt = mpnbt * mpwds
!    endif
!  Otherwise in variant #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnwbt = mpnbt * mp_setwp (iprec)
!>>
    mpnw = mpnwbt / mpnbt
    call mp_inpz (iu, z1, mpnw)
    call mp_inpz (iu, z2, mpnw)
    call mp_inpz (iu, z3, mpnw)
    call mp_inpz (iu, z4, mpnw)
    call mp_inpz (iu, z5, mpnw)
    return
  end subroutine

  function mp_rtod (ra)
    implicit none
    real (mprknd):: mp_rtod
    type (mp_real), intent (in):: ra
    real (mprknd) mpfrgetd
    external mpfrgetd
    call mp_fixlocr (ra)
    mp_rtod = mpfrgetd (ra%mpr(1:), mprnd)
    return
  end function

  function mp_rtoq (ra)
    implicit none
    real (max (mprknd2, kind (1.))):: mp_rtoq
    type (mp_real), intent (in):: ra
    type (mp_real) r1, r2
    real (mprknd) mpfrgetd, d1, d2, d3
    external mpfrgetd
    integer mpnwbt
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_fixlocr (ra)
    call mp_initvr (r1, mpnwbt)
    call mp_initvr (r2, mpnwbt)
    d1 = mpfrgetd (ra%mpr(1:), mprnd)
    mp_rtoq = d1
    call mpfrsubd (r1%mpr(1:), ra%mpr(1:), d1, mprnd)
    d2 = mpfrgetd (r1%mpr(1:), mprnd)
    mp_rtoq = mp_rtoq + d2
    call mpfrsubd (r2%mpr(1:), r1%mpr(1:), d2, mprnd)
    d3 = mpfrgetd (r2%mpr(1:), mprnd)
    mp_rtoq = mp_rtoq + d3
    return
  end function

  function mp_rtor (ra, iprec)
    implicit none
    type (mp_real):: mp_rtor
    type (mp_real), intent (in):: ra
    integer mpnwbt

!>  In variant #1, uncomment these lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnwbt = mpnbt * mp_setwp (iprec)
!    else
!      mpnwbt = mpnbt * mpwds
!    endif
!  Otherwise in variant #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnwbt = mpnbt * mp_setwp (iprec)
!>>
    call mp_fixlocr (ra)
    call mp_initvr (mp_rtor, mpnwbt)
    call mpfrset (mp_rtor%mpr(1:), ra%mpr(1:), mprnd)
    return
  end function

  function mp_rtoz (ra, rb, iprec)
    implicit none
    type (mp_complex):: mp_rtoz
    type (mp_real), intent (in):: ra, rb
    integer l1, mpnwbt

!>  In variant #1, uncomment these lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnwbt = mpnbt * mp_setwp (iprec)
!    else
!      mpnwbt = mpnbt * mpwds
!    endif
!  Otherwise in variant #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnwbt = mpnbt * mp_setwp (iprec)
!>>
    call mp_fixlocr (ra)
    call mp_fixlocr (rb)
    l1 = mpwds6
    call mp_initvz (mp_rtoz, mpnwbt)
    call mpfrset (mp_rtoz%mpc(1:), ra%mpr(1:), mprnd)
    call mpfrset (mp_rtoz%mpc(l1+1:), rb%mpr(1:), mprnd)
    return
  end function

  function mp_sign (ra, rb)
    implicit none
    type (mp_real):: mp_sign
    type (mp_real), intent (in):: ra, rb
    integer mpnwbt
    call mp_fixlocr (ra)
    call mp_fixlocr (rb)
    mpnwbt = max (ra%mpr(1), rb%mpr(1))
    mpnwbt = min (mpnwbt, mpwdsbt)
    call mp_initvr (mp_sign, mpnwbt)
    if (rb%mpr(2) == 1) then
      call mpfrset (mp_sign%mpr(1:), ra%mpr(1:), mprnd)
    else
      call mpfrneg (mp_sign%mpr(1:), ra%mpr(1:), mprnd)
    endif
    return
  end function

  function mp_sin (ra)
    implicit none
    type (mp_real):: mp_sin
    type (mp_real), intent (in):: ra
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_sin, mpnwbt)
    call mpfrsin (mp_sin%mpr(1:), ra%mpr(1:), mprnd)
    return
  end function

  function mp_sinh (ra)
    implicit none
    type (mp_real):: mp_sinh
    type (mp_real), intent (in):: ra
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_sinh, mpnwbt)
    call mpfrsinh (mp_sinh%mpr(1:), ra%mpr(1:), mprnd)
    return
  end function

  function mp_sqrt (ra)
    implicit none
    type (mp_real):: mp_sqrt
    type (mp_real), intent (in):: ra
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_sqrt, mpnwbt)
    call mpfrsqrt (mp_sqrt%mpr(1:), ra%mpr(1:), mprnd)
    return
  end function

  function mp_struve_hn (nu, ra)
    implicit none
    integer, intent(in):: nu
    type (mp_real):: mp_struve_hn
    type (mp_real), intent (in):: ra
    integer mpnwbt, mpnw
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    mpnw = mpnwbt / mpnbt
    call mp_initvr (mp_struve_hn, mpnwbt)
    call mpstruvehn (nu, ra%mpr(0:), mp_struve_hn%mpr(0:), mpnw)
    return
  end function

  function mp_tan (ra)
    implicit none
    type (mp_real):: mp_tan
    type (mp_real), intent (in):: ra
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_tan, mpnwbt)
    call mpfrtan (mp_tan%mpr(1:), ra%mpr(1:), mprnd)
    return
  end function

  function mp_tanh (ra)
    implicit none
    type (mp_real):: mp_tanh
    type (mp_real), intent (in):: ra
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_tanh, mpnwbt)
    call mpfrtanh (mp_tanh%mpr(1:), ra%mpr(1:), mprnd)
    return
  end function

!   Return working precision level (in words) of input MP value.

  function mp_wprec (ra)
    implicit none
    integer mp_wprec
    type (mp_real), intent (in):: ra
    integer mpnwbt
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    mp_wprec = mpnwbt / mpnbt
    return
  end function

  function mp_wprecz (za)
    implicit none
    integer mp_wprecz
    type (mp_complex), intent (in):: za
    integer l1, mpnwbt
    l1 = za%mpc(0)
    mpnwbt = max (int (za%mpc(1)), int (za%mpc(l1+1)))
    mp_wprecz = min (mpnwbt, mpwdsbt) / mpnbt
    return
  end function

!   Output routines.

  subroutine mp_writer (iu, ln, ld, r1, r2, r3, r4, r5)
    implicit none
    integer, intent (in):: iu, ln, ld
    type (mp_real), intent (in):: r1, r2, r3, r4, r5
    optional:: r2, r3, r4, r5

    call mp_outr (iu, ln, ld, r1)
    if (present (r2)) then
      call mp_outr (iu, ln, ld, r2)
    endif
    if (present (r3)) then
      call mp_outr (iu, ln, ld, r3)
    endif
    if (present (r4)) then
      call mp_outr (iu, ln, ld, r4)
    endif
    if (present (r5)) then
      call mp_outr (iu, ln, ld, r5)
    endif

    return
  end subroutine

  subroutine mp_writez (iu, ln, ld, z1, z2, z3, z4, z5)
    implicit none
    integer, intent (in):: iu, ln, ld
    type (mp_complex), intent (in):: z1, z2, z3, z4, z5
    optional:: z2, z3, z4, z5

    call mp_outz (iu, ln, ld, z1)
    if (present (z2)) then
      call mp_outz (iu, ln, ld, z2)
    endif
    if (present (z3)) then
      call mp_outz (iu, ln, ld, z3)
    endif
    if (present (z4)) then
      call mp_outz (iu, ln, ld, z4)
    endif
    if (present (z5)) then
      call mp_outz (iu, ln, ld, z5)
    endif

    return
  end subroutine

  function mp_zeta (ra)
    implicit none
    type (mp_real):: mp_zeta
    type (mp_real), intent (in):: ra
    integer mpnwbt
    call mp_fixlocr (ra)
    mpnwbt = min (ra%mpr(1), mpwdsbt)
    call mp_initvr (mp_zeta, mpnwbt)
    call mpfrzeta (mp_zeta%mpr(1:), ra%mpr(1:), mprnd)
    return
  end function

  function mp_zeta_int (ia, iprec)
    implicit none
    type (mp_real):: mp_zeta_int
    integer, intent (in):: ia
    integer mpnwbt, mpnw

!>  In variant #1, uncomment these lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnwbt = mpnbt * mp_setwp (iprec)
!    else
!      mpnwbt = mpnbt * mpwds
!    endif
!  Otherwise in variant #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnwbt = mpnbt * mp_setwp (iprec)
!>>
    mpnw = mpnwbt / mpnbt
    call mp_initvr (mp_zeta_int, mpnwbt)
    call mpzetaintr (ia, mp_zeta_int%mpr(0:), mpnw)
    return
  end function

  function mp_zeta_be (nb, rb, rc)
    implicit none
    integer, intent (in):: nb
    type (mp_real):: mp_zeta_be
    type (mp_real), intent (in):: rb(nb), rc
    integer n1, mpnw, mpnwbt
    mpnwbt = max (int (rb(1)%mpr(1)), int (rc%mpr(1)))
    mpnwbt = min (mpnwbt, mpwdsbt)
    mpnw = mpnwbt / mpnbt
    n1 = mpwds
    call mp_fixlocr (rc)
    call mp_initvr (mp_zeta_be, mpnwbt)
    call mpzetabe (n1, nb, rb(1)%mpr(0:), rc%mpr(0:), mp_zeta_be%mpr, mpnw)
    return
  end function

  function mp_zmtoz (za, iprec)
    implicit none
    type (mp_complex):: mp_zmtoz
    type (mp_complexm), intent (in):: za
    integer l1, l2, mpnwbt

!>  In variant #1, uncomment these lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnwbt = mpnbt * mp_setwp (iprec)
!    else
!      mpnwbt = mpnbt * mpwds
!    endif
!  Otherwise in variant #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnwbt = mpnbt * mp_setwp (iprec)
!>>
    l1 = za%mpc(0)
    call mp_fixloczm (za)
    l2 = mpwds6
    call mp_initvz (mp_zmtoz, mpnwbt)
    call mpfrset (mp_zmtoz%mpc(1:), za%mpc(1:), mprnd)
    call mpfrset (mp_zmtoz%mpc(l2+1:), za%mpc(l1+1:), mprnd)
    return
  end function

  function mp_ztodc (za)
    implicit none
    complex (kind(0.d0)):: mp_ztodc
    type (mp_complex), intent (in):: za
    integer l1
    real (mprknd) d1, d2, mpfrgetd
    external mpfrgetd
    call mp_fixlocz (za)
    l1 = za%mpc(0)
    d1 = mpfrgetd (za%mpc(1:), mprnd)
    d2 = mpfrgetd (za%mpc(l1+1:), mprnd)
    mp_ztodc = cmplx (d1, d2, mprknd)
    return
  end function

  function mp_ztor (za, iprec)
    implicit none
    type (mp_real):: mp_ztor
    type (mp_complex), intent (in):: za
    integer l1, mpnwbt

!>  In variant #1, uncomment these lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnwbt = mpnbt * mp_setwp (iprec)
!    else
!      mpnwbt = mpnbt * mpwds
!    endif
!  Otherwise in variant #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnwbt = mpnbt * mp_setwp (iprec)
!>>
    l1 = za%mpc(0)
    call mp_fixlocz (za)
    call mp_initvr (mp_ztor, mpnwbt)
    call mpfrset (mp_ztor%mpr(1:), za%mpc(1:), mprnd)
    return
  end function

  function mp_ztor2 (za)
    implicit none
    type (mp_real):: mp_ztor2
    type (mp_complex), intent (in):: za
    integer l1, mpnwbt
    mpnwbt = mpnbt * mpwds
    l1 = za%mpc(0)
    call mp_fixlocz (za)
    call mp_initvr (mp_ztor2, mpnwbt)
    call mpfrset (mp_ztor2%mpr(1:), za%mpc(1:), mprnd)
    return
  end function

  function mp_ztoz (za, iprec)
    implicit none
    type (mp_complex):: mp_ztoz
    type (mp_complex), intent (in):: za
    integer l1, l2, mpnwbt

!>  In variant #1, uncomment these lines:
!    integer, optional, intent (in):: iprec
!    if (present (iprec)) then
!      mpnwbt = mpnbt * mp_setwp (iprec)
!    else
!      mpnwbt = mpnbt * mpwds
!    endif
!  Otherwise in variant #2, uncomment these two lines:
    integer, intent (in):: iprec
    mpnwbt = mpnbt * mp_setwp (iprec)
!>>
    l1 = za%mpc(0)
    call mp_fixlocz (za)
    l2 = mpwds6
    call mp_initvz (mp_ztoz, mpnwbt)
    call mpfrset (mp_ztoz%mpc(1:), za%mpc(1:), mprnd)
    call mpfrset (mp_ztoz%mpc(l2+1:), za%mpc(l1+1:), mprnd)
    return
  end function

  subroutine mp_inpr (iu, a, mpnw)

!   This routine reads the MPR number A from logical unit IU.  The digits of A
!   may span more than one line, provided that a "\" appears at the end of
!   a line to be continued (any characters after the "\" on the same line
!   are ignored).  Individual input lines may not exceed 2048 characters in
!   length, although this limit can be changed in the system parameters
!   (parameter mpnstr) in module MPFUNA.  Embedded blanks are allowed anywhere.
!   An exponent with "e" or "d" may optionally follow the numeric value.

!   A scratch array below (CHR1) holds character data for input to MPFR_SET_STR.
!   It is dimensioned MPNWBT * (MPNDPW + 1) + 1000 (see below).
!   If more nonblank input characters than this are input, they are ignored.

  implicit none
  integer i, i1, iu, lnc1, lncx, ln1, mpnw
  character(mpnstr) line1
  character(18), parameter:: validc = ' 0123456789+-.dDeE'
  character(1) chr1(mpnw*(mpndpw+1)+1001)
  type (mp_real) a
  integer mpnwbt

  mpnwbt = mpnw * mpnbt
  call mp_initvr (a, mpnwbt)
  lnc1 = 0
  lncx = mpnw * (mpndpw + 1) + 1000

100 continue

  read (iu, '(a)', end = 200) line1

!   Find the last nonblank character.

  do i = mpnstr, 1, -1
    if (line1(i:i) /= ' ') goto 110
  enddo

!   Input line is blank -- ignore.

  goto 100

110 continue

  ln1 = i

!   Scan input line, looking for valid characters.

  do i = 1, ln1
    if (line1(i:i) == '\') goto 100
    i1 = index (validc, line1(i:i))
    if (i1 == 0 .and. line1(i:i) /= ' ') then
      write (6, 2) line1(i:i)
2     format ('*** mp_inpr: Invalid input character = ',a)
      call mp_abrt (87)
    elseif (line1(i:i) /= ' ') then
      if (lnc1 < lncx) then
        lnc1 = lnc1 + 1
        if (line1(i:i) == 'D' .or. line1(i:i) == 'd') then
          chr1(lnc1) = 'e'
        else
          chr1(lnc1) = line1(i:i)
        endif
      endif
    endif
  enddo

  chr1(lnc1+1) = char(0)
  call mpfrsetstr (a%mpr(1:), chr1, mprnd)

  goto 300

200  continue

  write (mpldb, 4)
4 format ('*** mp_inpr: End-of-file encountered.')
  call mp_abrt (72)

300 return
  end subroutine

  subroutine mp_inpz (iu, a, mpnw)
    implicit none
    integer, intent (in):: iu, mpnw
    type (mp_complex), intent (out):: a
    type (mp_real) r1, r2
    call mp_inpr (iu, r1, mpnw)
    call mp_inpr (iu, r2, mpnw)
    a = mp_rtoz (r1, r2, mpnw)
    return
  end subroutine

  subroutine mp_outr (iu, ln, nd, a)

!   This routine writes MPR number A to logical unit IU in E(LN,ND) format.
!   This is output on mpoutln characters per line.  The value of mpoutln is set
!   in the system parameters at the start of module MPFUNA.

  implicit none
  integer i, iu, ln, ln1, nd
  character(1) chr1(ln)
  character(32) cform1, cform2
  type (mp_real) a

  call mp_fixlocr (a)
  call mp_eform (a, ln, nd, chr1)

  write (cform1, 1) mpoutln
  1 format ('(',i8,'a1)')
  write (cform2, 2) mpoutln
  2 format ('(',i8,'a1,"\")')

  if (ln <= mpoutln) then
    write (iu, fmt = cform1) (chr1(i), i = 1, ln)
  elseif (mod (ln, mpoutln) == 0) then
    ln1 = mpoutln * (ln / mpoutln) - mpoutln
    write (iu, fmt = cform2) (chr1(i), i = 1, ln1)
    write (iu, fmt = cform1) (chr1(i), i = ln1 + 1, ln)
  else
    ln1 = mpoutln * (ln / mpoutln)
    write (iu, fmt = cform2) (chr1(i), i = 1, ln1)
    write (iu, fmt = cform1) (chr1(i), i = ln1 + 1, ln)
  endif

  return
  end subroutine

  subroutine mp_outz (iu, ln, nd, a)
    implicit none
    integer, intent (in):: iu, ln, nd
    type (mp_complex), intent (in):: a
    call mp_outr (iu, ln, nd, mp_ztor2 (a))
    call mp_outr (iu, ln, nd, mp_aimag (a))
    return
  end subroutine

end module mpfung

