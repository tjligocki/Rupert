
MPFUN2015: A thread-safe arbitrary precision package
MPFUN-Fort (all-Fortran) version

Revision date:  24 Apr 2016

AUTHOR:
   David H. Bailey
   Lawrence Berkeley National Lab (retired) and University of California, Davis
   Email: dhbailey@lbl.gov
   
COPYRIGHT AND DISCLAIMER:
  All software in this package (c) 2016 David H. Bailey.
  By downloading or using this software you agree to the copyright, disclaimer
  and license agreement in the accompanying file DISCLAIMER.txt.

1. PURPOSE OF PACKAGE:
  This system permits one to perform floating-point computations (real and
  complex) to arbitrarily high numeric precision, by making only relatively
  minor changes to existing Fortran-90 programs.  All basic arithmetic
  operations and transcendental functions are supported, together with several
  special functions.
  
  In addition to fast execution times, one key feature of this package is a
  100% THREAD-SAFE design, which means that user-level applications can be
  easily converted for parallel execution, say using a threaded parallel
  environment such as OpenMP.  There are NO global shared variables (except
  static compile-time data), and NO initialization is necessary unless
  extremely high precision (> 19,500 digits) is required.
  
  A related software system by the same author is MPFUN-MPFR, which is
  similar to MPFUN-Fort in its user interface, but which calls the MPFR
  package for all low-level functions and operations.  The MPFUN-MPFR version
  is significantly faster on most operations than MPFUN-Fort, but installation
  is more complicated.  It is available from the same site as this package:
  http://www.davidhbailey.com/dhbsoftware/

2. DOCUMENTATION:
  A detailed description of this software, with instructions for writing Fortran
  code to use the package, is available in this technical paper:
   
  David H. Bailey, "MPFUN2015: A thread-safe arbitrary precision package," 
  available at http://www.davidhbailey.com/dhbpapers/mpfun2015.pdf.
   
3. INSTALLATION, COMPILATION AND LINKING:
  Installation, compilation and linking is relatively straightforward,
  provided that you have a Unix-based system, such as Linux or Apple OSX, with
  a command-line interface (such as the Terminal application of Apple OSX).  
  
  For Apple OSX systems, you first must install the "Command Line Tools"
  package, which is available (for free) from the Apple Developer website.
  Login (or register, if this your first access), using your Apple ID at:
    https://developer.apple.com/devcenter/mac/index.action
  Then click on "View all downloads" and select "Command Line Tools" for
  your particular version of the MAC OSX operating system.  Install the
  downloaded package on your system.

  The gfortran compiler, which is highly recommended for this package, is
  available for a variety of systems at the website
    https://gcc.gnu.org/wiki/GFortranBinaries
  When you attempt to install gfortran, you may get the message "gfortran.pkg
  can't be opened because it is from an unidentified developer".  If so, go to
  the Security page of your Mac's System Preferences, and click "Open anyway".

  The MPFUN-Fort package also works with IBM's xlf_r compiler, Intel's ifort
  compiler and Portland Group's pgf90.

  Once this has been done, it is easy to install the MPFUN-Fort package.  To do
  this, download the file "MPFUN-Fort-v07.tar.gz" (or whatever is the latest
  version) into your Documents folder.  If it is not decompressed by your browser,
  type "gunzip MPFUN-Fort-v07.tar.gz" (or the equivalent name for the latest
  version).  Then type "tar xfv MPFUN-Fort-v07.tar (or the equivalent name).  
  This should create the directory MPFUN-Fort-v07 (or the equivalent name).
  
  For both MPFUN-Fort and MPFUN-MPFR, there are actually two variants of the
  software, both of which are included in the distribution file:
    
  Variant 1: This is recommended for basic applications that do not dynamically
    change the precision level (or do so only rarely).
  Variant 2: This is recommended for more sophisticated applications that
    dynamically change the precision level.

  The two variants of the packages correspond to two variants of module
  MPFUNG, the high-level language interface module.  Compile/link scripts are
  available in the fortran directory of the MPFUN2015 software for the gfortran
  compiler and other supported compilers.  For the gfortran and ifort compilers,
  which support the real*16 datatype, the respective compile/link scripts
  include the proper modules.

  For example, to compile variant 1 of the library with the GNU gfortran
  compiler, go to the fortran directory and type
    ./gnu-complib1.scr
  To compile and link the application program prog.f90 for variant 1, producing
  the executable file prog, type
    ./gnu-complink1.scr prog

  NOTE: For both compilers, the very first time you compile the library (using
  either the complib1.scr or complib2.scr scripts), you may see "fatal" errors,
  such as various modules not found.  This is normal -- just repeat the library
  compile scripts. The library compile scripts invoke the compiler twice for
  this reason.

  Seven test programs are included in the fortran directory of the package.  The
  script mpfun-tests.scr, which is included in the distribution package for each
  version, compiles variant 2 of the library, then compiles, links and runs all
  seven of the test programs. 

4. CODING INSTRUCTIONS AND USAGE:

  Here is a brief summary of Fortran coding instructions.  For full details,
  see the documentation paper mentioned above.
  
  To use either version from a Fortran program, first set the parameter mpipl,
  the "default" precision level in digits, which is the maximum precision
  level to be used for subsequent computation, and is used to specify the amount
  of storage required for multiprecision data.  mpipl is set in a parameter
  statement at the start of module MPFUNF, which is in file mpfunf.f90.  In the
  code as distributed, mpipl is set to 1200 digits (sufficient to run the seven
  test programs), but it can be set to any level greater than or equal to 30 digits.
  mpipl is automatically converted to mantissa words by the formula 
    mpwds = int (mpipl / mpdpw + 2),
  where mpdpw is a system parameter, and where int () means truncate to integer.
  For MPFUN-Fort, mpdpw is log_{10} (2^{48}) = 14.44943979187..., whereas for
  MPFUN-MPFR it is  log_{10}(2^{64}) =  19.26591972249... (both values are double
  precision approximations).  The resulting parameter mpwds is the internal default
  precision level, in words.  All computations are performed to mpwds precision
  unless the user, within an application code, specifies a lower precision level.
  
  After setting the value of mpipl in module MPFUNF, compile the appropriate
  version of the library, using one of the scripts mentioned above.
  
  Next, place the following line in every subprogram of the user's application
  code that contains a multiprecision variable or array, at the beginning of the
  declaration section, before any implicit or type statements:
    use mpmodule

  To designate a variable or array as multiprecision real (MPR) in your
  application code, use the Fortran-90 type statement with the type "mp_real",
  as in this example:
    type (mp_real) a, b(m), c(m,n)
  Similarly, to designate a variable or array as multiprecision complex
  (MPC), use a type statement with "mp_complex".

  Thereafter when one of these variables or arrays appears in code, e.g.,
     d = a + b(i) * sqrt(3.d0 - c(i,j))
  the proper multiprecision routines are automatically called by the
  Fortran compiler.
  
  Most common mixed-mode combinations (arithmetic operations, comparisons and
  assignments) involving MPR, MPC, double precision (DP) and integer arguments
  are supported.  A complete list of supported mixed-mode operations is given
  in the documentation paper.  However, there are some hazards.

  For example, the code r1 = 3.14159d0, where r1 is MPR, does NOT produce
  the true multiprecision equivalent of 3.14159, unless the numerical
  value is a modest-sized whole number or exact binary fraction.  Similarly,
  the code r2 = r1 + 3.d0 * sqrt (2.d0), where r1 and r2 are MPR, does NOT
  produce the true multiprecision value, since the expression
  3.d0 * sqrt (2.d0) will be performed in double precision (according to
  standard Fortran-90 precedence rules).  See documentation for instructions
  on handling these situations.  

  Input/output of MP variables or array elements is done using the
  subroutines mpread and mpwrite.  See documentation for details.

  Standard Fortran intrinsics are supported with MPR and MPC arguments,
  and they operate similarly to the standard double precision (DP) and
  double complex (DC) equivalents.  A complete list of supported functions and
  subroutines is given in the documentation paper.  

  Some applications do not need to change the working precision from the
  initially-defined default level, whereas others need to change frequently.
  Accordingly, for both MPFUN-Fort and MPFUN-MPFR, there are two variants of
  the language interface module MPFUNG:

  Variant 1: This is recommended for basic applications that do not dynamically
    change the precision level (or do so only rarely).
  Variant 2: This is recommended for more sophisticated applications that
    dynamically change the precision level.

  See documentation for full details on the differences between these two
  variants.
  
  Several application programs (tpphix3.f90, tpslq1.f90, tpslqm1.f90,
  tpslqm2.f90, tpslqm3.f90, tpphix3.f90, tquadts.f90 and tquadtsp.f90) are
  included, together with corresponding output files for comparison with user
  results. If, after compiling the library and running each of these programs,
  the results in these reference output files can be reproduced (except for
  timings, iteration counts, etc.), then one can be fairly confident that the
  software is working properly.  Full descriptions of these programs are
  included in the documentation paper.



