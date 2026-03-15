*****************************************************************************

MPFUN2020: A thread-safe arbitrary precision package
MPFUN20-MPFR version

Revision date: 3 Feb 2023

AUTHOR:
David H. Bailey
Lawrence Berkeley National Lab (retired) and University of California, Davis
Email: dhbailey@lbl.gov
 
COPYRIGHT AND DISCLAIMER:
All software in this package (c) 2023 David H. Bailey. By downloading or using this software you agree to the copyright, disclaimer and license agreement in the accompanying file DISCLAIMER.txt.

FULL DOCUMENTATION:
 David H. Bailey, "MPFUN2020: A new thread-safe arbitrary precision package," 
https://www.davidhbailey.com/dhbpapers/mpfun2020.pdf

INDEX OF THIS README FILE:

I. PURPOSE OF PACKAGE
II. INSTALLING CODING ENVIRONMENT (FOR MAC OS X SYSTEMS)
III. INSTALLING FORTRAN COMPILER (IF NEEDED)
IV. DOWNLOADING AND INSTALLING GMP AND MPFR
V. DOWNLOADING AND COMPILING MPFUN20-MPFR
VI. BRIEF SUMMARY OF CODING INSTRUCTIONS AND USAGE
VII. SAMPLE APPLICATION PROGRAMS AND TESTS
VIII. RECENTLY TESTED PLATFORMS
IX. RECENT UPDATES
X. APPENDIX: TRANSCENDENTAL, SPECIAL AND MISCELLANEOUS FUNCTIONS

+++++


I. PURPOSE OF PACKAGE

This package permits one to perform floating-point computations (real and complex) to arbitrarily high numeric precision, by making only relatively minor changes to existing Fortran-90 programs (mostly changes to type statements). All basic arithmetic operations and transcendental functions are supported, together with numerous special functions.

The package comes in two versions: one completely self-contained, all-Fortran version that is simple to install; and one version based on the MPFR package that is more complicated to install but runs somewhat faster on most applications. Both versions are completely thread-safe, which means that user-level applications can be easily converted for parallel execution, say by using a threaded parallel environment such as OpenMP. Both versions also detect, and provide means to overcome, accuracy problems rooted in the usage of inexact double-precision constants and expressions. A high-level Fortran-90 interface, supporting both multiprecision real and complex datatypes, is provided for each, so that most users need only to make minor changes to existing double-precision code. The two versions are "plug-compatible" in the sense that applications written for one also run with the other (provided a simple guideline is followed). 

The two versions of this package are:

MPFUN20-Fort: This is an all-Fortran version based on 8-byte integer arithmetic. It includes support for a medium precision datatype, which results in faster execution on very large problems, and features FFT-based multiplication to accelerate very high precision computations. It compiles in just a few seconds on any system with a Fortran-2008 compliant compiler (examples include the GNU gfortran compiler, the Intel ifort compiler and the NAG Fortran compiler).

MPFUN20-MPFR: This is virtually identical to MPFUN20-Fort in its user interface, but it calls the MPFR package for all low-level functions and operations. The MPFUN20-MPFR version is faster than MPFUN20-Fort on most applications, particularly those that involve transcendental functions. However, installation of MPFUN20-MPFR is significantly more complicated (because the GMP and MPFR packages must first be installed, usually requiring administrator privilege).

What follows are the instructions for MPFUN20-MPFR.


II. INSTALLING CODING ENVIRONMENT (FOR MAC OS X SYSTEMS)

For Apple Mac OS X systems (highly recommended for MPFUN20-Fort), first install the latest supported version of Xcode, which is available for free from the Apple Developer website:

  https://developer.apple.com/

Here click on Account, then enter your Apple ID and password, then go to 

  https://developer.apple.com/download/more/?=xcode

From this list, select Xcode. As of the above date, the latest version of Xcode is 13.3. Download this package on your system, double click to decompress, and place the resulting Xcode app somewhere on your system (typically in the Applications folder). Double-click on the app to run Xcode, allowing it install additional components, then quit Xcode. Open a terminal window, using the Terminal application in the Utilities folder, and type

  xcode-select --install

which installs various command-line tools. The entire process of downloading Xcode and installing command-line tools takes roughly 25 minutes. When this is completed, you should be ready to continue with the installation.


III. INSTALLING FORTRAN COMPILER (IF NEEDED)

Running MPFUN20-MPFR is relatively straightforward, provided that one has a Unix-based system, such as Linux or Apple OS X, and a Fortran-2008 compliant compiler. These requirements are met by the GNU gfortran compiler, the Intel ifort compiler, the NAG nagfor compiler, IBM's xlf, PGI's pgf90 and others.

The gfortran compiler (highly recommended for MPFUN20-Fort) is available free for a variety of systems at this website:

  https://gcc.gnu.org/wiki/GFortranBinaries

For Apple Mac OS X systems, download the installer file here:

  https://github.com/fxcoudert/gfortran-for-macOS/releases

The gfortran compiler is normally placed in /usr/local/lib and /usr/local/bin. Thus before one uses gfortran, one must insert a line in one's shell initialization file (if the Z shell is used, as on most Apple OS X systems, the shell initialization file is ~/.zshrc). The line to be included is:

  PATH=/usr/local/lib:/usr/local/bin:$PATH

The following line is also recommended for gfortran compiler users:

  GFORTRAN_UNBUFFERED_ALL=yes; export GFORTRAN_UNBUFFERED_ALL

The following line is recommended for inclusion in the shell initialization file, no matter what compiler is used (it prevents stack overflow system errors):

  ulimit -s unlimited

On most Unix systems (including Apple Mac OS X systems), the shell initialization file must be manually executed upon initiating a terminal shell, typically by typing "source .zshrc".


IV. DOWNLOADING AND INSTALLING GMP AND MPFR

To use MPFUN20-MPFR, you must first install the GMP and MPFR packages, as follows (this presumes you have a command-line interface and gcc installed):
  
1. Using a browser (e.g., Apple Safari or Firefox), download the file "gmp-6.2.1.tar.xz" (or whatever is latest version) from https://gmplib.org. The .xz file is more likely to be decompressed by your browser than the .lz file. Move this file to a suitable spot on your system, typically to the Documents folder.

2. In the Documents folder (or wherever the tar file was moved), type 

  tar -xf gmp-6.2.1.tar.xz

(or whatever version was downloaded). This should create the directory "gmp-6.2.1" (or a similar name).  

3. Change directory to gmp-6.2.1 (or similar name) and type "./configure", followed by "make", then "make check". All tests should pass. Then type "make install". On Apple OS X systems and some others, you may need to type instead "sudo make install", which will request your computer system's admin password. This should place several files, including libgmp.10.dylib, libgma.a, libgmp.dylib and libgmp.la, in /usr/local/lib. This process typically takes 3-10 minutes.

4. Using a browser (e.g., Apple Safari or Firefox), download "mpfr-4.1.0.tar.xz" (or whatever is the latest version) from https://www.mpfr.org/mpfr-current/. Move this file to a suitable spot on your system, typically in your Documents folder.

5. In the Documents folder (or wherever the tar file was moved), type

  tar -xf mpfr-4.1.0.tar.xz

(or whatever version was downloaded). This should create the directory "mpfr-4.1.0" (or a similar name).  

6. Change directory to mpfr-4.1.0 (or similar name) and type ./configure, followed by "make", then "make check". All tests should pass. Then type "make install". On Apple OS X systems and some others, you may need to type instead "sudo make install", which will request your computer system's admin password. This should place several files, including libmpfr.4.dylib, libmpfr.6.dydlib, libmpfr.a, libmpfr.dylib and libmpfr.la, in /usr/local/lib. This process typically takes 3-10 minutes.

7. To test the installation, place the C program "sample" (beginning with the line "#include <stdio.h>") from the URL https://www.mpfr.org/sample.html into a file "sample.c" (located anywhere within the Documents folder), then compile by typing

  gcc -o sample sample.c -lmpfr -lgmp -L/usr/local/lib

Then when you type "./sample", you should see the single line of output given at the bottom of the URL https://www.mpfr.org/sample.html.


V. DOWNLOADING AND COMPILING MPFUN20-MPFR

From the website https://www.davidhbailey.com/dhbsoftware, download the file "mpfun20-mpfr-vnn.tar.gz" (replace "vnn" by whatever is the current version on the website, such as "v22"). If the file is not decompressed by your browser, use gunzip at the shell level to do this. Some browsers (such as the Apple Safari browser) do not drop the ".gz" suffix after decompression; if so, remove this suffix manually at the shell level. Then type

  tar xfv mpfun20-mpfr-vnn.tar  

(where again "vnn" is replaced by the downloaded version). This should create the directory and unpack all files.

The MPFUN20-MPFR software comes in two variants, which are in directories fortran-var1 and fortran-var2, respectively:

Variant 1: This is recommended for beginning users and for basic applications that do not dynamically change the working precision level (or do so only rarely).

Variant 2: This is recommended for more sophisticated applications that dynamically change the working precision level. It does not allow some mixed-mode combinations, and requires one to explicitly specify a working precision parameter for some functions. However, in the present author's experience, these restrictions result in less overall effort to produce a debugged, efficient application code.

See documentation paper for additional details on the differences between these two variants. The Fortran source files and scripts required for each of these variants are in the respective directories fortran-var1 and fortran-var2.

Compile/link scripts are available for the GNU gfortran and the Intel ifort compilers. These scripts automatically select the proper source files from the package for compilation and employ the appropriate compiler flags. For example, to compile Variant 1 of the library using the GNU gfortran compiler, go to the fortran-var1 directory and type

  ./gnu-complib1.scr

Then to compile and link the application program tpslq1.f90 for variant 1, using the GNU gfortran compiler, producing the executable file tpslq1, type

  ./gnu-complink1.scr tpslq1

To execute the program, with output to tpslq1.txt, type

  ./tpslq1 > tpslq1.txt

These scripts assume that the user program is in the same directory as the library files; this can easily be changed by editing the script files.

Several sample test programs, together with reference output files, are included in the fortran-var1 and fortran-var2 directories -- see Section VIII below.


VI. BRIEF SUMMARY OF CODING INSTRUCTIONS AND USAGE

What follows is a brief summary of Fortran coding instructions. For full details, see the documentation paper:

David H. Bailey, "MPFUN2020: A new thread-safe arbitrary precision package," 
https://www.davidhbailey.com/dhbpapers/mpfun2020.pdf

First set the parameter mpipl, the default standard precision level in digits, which is the maximum precision level to be used for subsequent computation, and is used to specify the amount of storage required for multiprecision data. mpipl is set in a parameter statement in file mpfunf.f90 in the fortran-var1 or fortran-var2 directory of the software. In the code as distributed, mpipl is set to 2500 digits (sufficient to run each of the test programs), but it can be set to any level greater than 50 digits. mpipl is automatically converted to mantissa words by the formula 

  mpwds = int (mpipl / mpdpw + 2)

where mpdpw (digits per word) is a system parameter (approx. 18.0617997398 for MPFUN20-Fort and 19.2659197224 in MPFUN-MPFR) set in file mpfuna.f90. The resulting parameter mpwds is the internal default precision level, in words. All subsequent computations are performed to mpwds words precision unless the user, within an application code, specifies a lower precision.

After setting the value of mpipl, compile the library, using one of the scripts mentioned above (e.g., gnu-complib1.scr if using the GNU gfortran compiler or intel-complib1.scr if using the Intel compiler).

Next, place the following line in every subprogram of the user's application code that contains a multiprecision variable or array, at the beginning of the declaration section, before any implicit or type statements:

  use mpmodule

To designate a variable or array as multiprecision real (MPR) in an application code, use the Fortran-90 type statement with the type "mp_real", as in this example:

  type (mp_real) a, b(m), c(m,n)

Similarly, to designate a variable or array as multiprecision complex (MPC), use a type statement with "mp_complex".

Thereafter when one of these variables or arrays appears in code, e.g.,

  d = a + b(i) * sqrt(3.d0 - c(i,j))

the proper multiprecision routines are automatically called by the Fortran compiler.

Most common mixed-mode combinations (arithmetic operations, comparisons and assignments) involving MPR, MPC, double precision (DP) and integer arguments are supported, although restrictions apply if one uses Variant 2 of the MPFUN20-Fort software. A complete list of supported mixed-mode operations is given in the documentation paper.

Users should be aware, however, that there are some hazards in this type of programming, inherent in conventions adopted by all Fortran compilers. For example, the code r1 = 3.14159d0, where r1 is MPR, does NOT produce the true multiprecision equivalent of 3.14159. In fact, the software will flag such usage with a run-time error. To obtain the full MPR converted value, write this as r1 = '3.14159', or, if using variant 2, as r1 = mpreal ('3.14159', nwds), where nwds is the level of working precision to be assigned to r1. Similarly, the code r2 = r1 + 3.d0 * sqrt (2.d0), where r1 and r2 are MPR, does NOT produce the true multiprecision value one might expect, since the expression 3.d0 * sqrt (2.d0) will be performed in double precision, according to Fortran-90 precedence rules. In fact, the above line of code will result in a run-time error. To obtain the fully accurate result, write this as r2 = r1 + 3.d0 * sqrt (mpreal (2.d0)), or, if using variant 2, as r2 = r1 + 3.d0 * sqrt (mpreal (2.d0, nwds)), where nwds is the level of working precision. See documentation paper for details.

Input and output of MPR and MPC data are performed using the subroutines mpread and mpwrite. For example, to output the variable r1 in E format to Fortran unit 6 (standard output), to 100-digit accuracy, in a field of width 120 characters, use the line of code

  call mpwrite (6, 120, 100, r1)

The second argument (120 in the above example) must be at least 20 larger than the third argument (100 in the above example). To read the variable r1 from Fortran unit 5 (standard input), use the line of code

  call mpread (5, r1)

or, if using variant 2, as

  call mpread (5, r1, nwds)

where nwds is the level of working precision to be assigned to r1. See documentation paper for details such as formatting.

Most Fortran-2008 intrinsic functions are supported with MPR and MPC arguments, as appropriate, and numerous special functions are also supported. A complete list of supported functions and subroutines is summarized in the Appendix below (section X), and also in the documentation paper. 


VII. SAMPLE APPLICATION PROGRAMS AND TESTS

The current release of the software includes a set of sample application programs in the fortran-var1 and fortran-var2 directories (the files are identical between directories):

testmpfun.f90  Tests most arithmetic, transcendental and special functions.
tpslq1.f90   Performs the standard 1-level PSLQ integer relation algorithm.
tpslqm1.f90  Performs the 1-level multipair PSLQ integer relation algorithm.
tpslqm2.f90  Performs the 2-level multipair PSLQ integer relation algorithm.
tpslqm3.f90  Performs the 3-level multipair PSLQ integer relation algorithm.
tpphix3.f90  Performs a Poisson polynomial application, using 3-level multipair PSLQ.
tquad.f90   Evaluates a set of definite integrals, using tanh-sinh, exp-sinh and sinh-sinh algorithms.
tquadgs.f90  Evaluates a set of definite integrals, using Gaussian quadrature.

Corresponding reference output files (e.g., tpphix3.ref.txt) are also included for each of the above programs.

In addition, the fortran-var1 and fortran-var2 directories include test scripts that compile the library and run each of the above sample programs above (except tquadgs.f90, which takes considerably more run time). In directory fortran-var1, these scripts are:

gnu-mpfun-tests1.scr
intel-mpfun-tests1.scr
nag-mpfun-tests1.scr

and the same scripts in directory fortran-var2, except for 2 instead of 1 in the filenames. For each test program, the script outputs either TEST PASSED or TEST FAILED. If all tests pass, then one can be fairly confident that the MPFUN2020 software and underlying compilers are working properly. Full descriptions of these application programs are included in the documentation paper:

David H. Bailey, "MPFUN2020: A new thread-safe arbitrary precision package," 
https://www.davidhbailey.com/dhbpapers/mpfun2020.pdf


VIII. RECENTLY TESTED PLATFORMS AND NOTES:

1. gfortran compiler, version 12.0.0, on an Apple MacBook Pro, OS X version 12.3.1, with an Apple Silicon (ARM) M1 Pro processor.
2. gfortran compiler, version 12.0.0, on an Apple Mac Studio, OS X version 12.3.1, with an Apple Silicon (ARM) M1 Max processor.
3. gfortran compiler, version 10.2.1, on a Debian Linux system, version 4.19.152-1, with an Intel processor.
4. NAG nagfor compiler, version 7.0, on a Debian Linux system, version 4.19.152-1, with an Intel processor.
5. Intel ifort compiler, version 2021.4.0, on a Debian Linux system, version 4.19.152-1, with an Intel processor. 

NOTE for Intel ifort: Due to what is evidently a compiler bug, one must make these changes to mpfunh1.f90 and mpfunh2.f90 before compiling: replace "mp_polylog_ini" (3 instances) with "mp_polylog_inim"; and replace "mp_polylog_neg" (6 instances) with "mp_polylog_negm".


IX. RECENT UPDATES:

3 Mar 2022: Revised mpfune.f90, mpfung1.f90, mpfung2.f90, mpfunh1.f90 and mpfunh2.f90 to implement the hypergeom_pfq function; updated testmpfun.f90.

18 Apr 2022: Fixed a bug in mpfung1.f90, mpfung2.f90, mpfunh1.f90 and mpfunh2.f90.

14 May 2022: Implemented the special functions hurwitz_zetan_be and polygamma_be in mpfune.f90, with corresponding changes in mpfung1.f90, mpfung2.f90, mpfunh1.f90, mpfunh2.f90 and testmpfun.f90.

19 May 2022: Implemented the special functions bessel_i, bessel_j, bessel_k, bessel_y in mpfune.f90 with corresponding changes in mpfung1.f90, mpfung2.f90, mpfunh1.f90, mpfunh2.f90, mpmodule.f90 and testmpfun.f90.

7 Jan 2023: Fixed a problem with mpeformat and mpfformat; inserted intent statements in all routines that did not already have them; changed all parameter statements to object oriented style.

13 Jan 2023: Made significant improvements to tpslq1.f90, tpslqm1.f90, tpslqm2.f90, tpslqm3.f90 and tpphix3.f90.

2 Feb 2023: Modified the mpfuna.f90 and mpfune.f90 codes to facilitate translation to DQFUN; fixed three minor bugs.


X. APPENDIX: TRANSCENDENTAL, SPECIAL AND MISCELLANEOUS FUNCTIONS:

As mentioned above, the MPFUN20-Fort and MPFUN20-MPFR packages support most Fortran-2008 intrinsics, including all the well-known transcendentals (e.g., sin, exp, log, etc.), and, in addition, a set of 30 special functions (e.g., Bessel functions, gamma function, zeta function, etc.). Further, the package includes a set of I/O and conversion functions, such as functions to convert between double precision and multiprecision real or between double complex and multiprecision complex. These functions and subroutines are listed below (a few are listed in more than one table). Some additional functions to work with the medium precision datatype are listed in the documentation paper.

In these tables, "F" denotes function, "S" denotes subroutine, "MPR" denotes multiprecision real, "MPC" denotes multiprecision complex, "DP" denotes double precision, "DC" denotes double complex, "Int" denotes integer and "QP" denotes IEEE quad precision (if supported by the compiler). The variable names r1,r2,r3 are MPR; z1 is MPC; d1 is DP; dc1 is DC; q1 is QP; i1,i2,i3,n,nb,np,nq are integers; s1 is character(1); sn is character(n); rb is an MPR vector of length nb; ss is an MPR vector of length n; aa is an MPR vector of length np; and bb is an MPR vector of length nq.

1. Standard Fortran-2008 transcendental functions:

Type   Name              Description
MPR    abs(r1)           Absolute value
MPR    abs(z1)           Absolute value of complex arg
MPR    acos(r1)          Inverse cosine
MPR    acosh(r1)         Inverse hyperbolic cosine
MPR    aimag(z1)         Imaginary part of complex arg
MPR    aint(r1)          Truncates to integer
MPR    anint(r1)         Rounds to closest integer
MPR    asin(r1)          Inverse sine
MPR    asinh(r1)         Inverse hyperbolic sine
MPR    atan(r1)          Inverse tangent
MPR    atan2(r1,r2)      Arctangent with two args
MPR    atanh(r1)         Inverse hyperbolic tangent
MPR    bessel_j0(r1)     Bessel function of the first kind, order 0
MPR    bessel_j1(r1)     Bessel function of the first kind, order 1
MPR    bessel_jn(n,r1)   Besel function of the first kind, order n
MPR    bessel_y0(r1)     Bessel function of the second kind, order 0
MPR    bessel_y1(r1)     Bessel function of the second kind, order 1
MPR    bessel_yn(n,r1)   Besel function of the second kind, order n
MPC    conjg(z1)         Complex conjugate
MPR    cos(r1)           Cosine of real arg
MPC    cos(z1)           Cosine of complex arg
MPR    cosh(r1)          Hyperbolic cosine
DP     dble(r1)          Converts MPR argument to DP
DC     dcmplx(z1)        Converts MPC argument to DC
MPR    erf(r1)           Error function
MPR    erfc(r1)          Complementary error function
MPR    exp(r1)           Exponential function of real arg
MPC    exp(z1)           Exponential function of complex arg
MPR    gamma(r1)         Gamma function
MPR    hypot(r1,r2)      Hypotenuse of two args
MPR    log(r1)           Natural logarithm of real arg
MPC    log(z1)           Natural logarithm of complex arg
MPR    log10(r1)         Base-10 logarithm of real arg
MPR    log_gamma(r1)     Log gamma function
MPR    max(r1,r2)        Maximum of two (or three) args
MPR    min(r1,r2)        Minimum of two (or three) args
MPR    mod(r1,r2)        Mod function = r1 - r2*aint(r1/r2)
MPR    sign(r1,r2)       Transfers sign from r2 to r1
MPR    sin(r1)           Sine function of real arg
MPC    sin(z1)           Sine function of complex arg
MPR    sinh(r1)          Hyperbolic sine
MPR    sqrt(r1)          Square root of real arg
MPC    sqrt(z1)          Square root of complex arg
MPR    tan(r1)           Tangent function
MPR    tanh(r1)          Hyperbolic tangent function

2. Special functions:

Type   Name                 Description 
F:MPR  agm(r1,r2)           Arithmetic-geometric mean 
F:MPR  airy(r1)             Airy function [1] 
S      mpberne(nb,rb)       Initialize array rb of length nb with first nb even 
                              Bernoulli numbers [2][3] 
F:MPR  bessel_i(r1,r2)      BesselI function, order r1, of r2  
F:MPR  bessel_in(n,r1)      BesselI function, order n, of r1 
F:MPR  bessel_j(r1,r2)      BesselJ function, order r1, of r2  
F:MPR  bessel_jn(n,r1)      BesselJ function, order n, of r1 
F:MPR  bessel_i(r1,r2)      BesselK function, order r1, of r2  
F:MPR  bessel_kn(n,r1)      BesselK function, order n, of r1 
F:MPR  bessel_y(r1,r2)      BesselY function, order r1, of r2  
F:MPR  bessel_yn(n,r1)      BesselY function, order n, of r1 
F:MPR  digamma(r1)          Digamma function of r1 [1] 
F:MPR  digamma_be(nb,rb,r1) Digamma function of r1, using nb even Bernoulli 
                              numbers in rb [3] 
F:MPR  erf(r1)              Error function 
F:MPR  erfc(r1)             Complementary error function 
F:MPR  expint(r1)           Exponential integral function 
F:MPR  gamma(r1)            Gamma function 
F:MPR  hurwitz_zetan(k,r1)  Hurwitz zeta function, order n >= 2, of r1 [4] 
F:MPR  hurwitz_zetan_be (nb,rb,n,r1)  Hurwitz zeta function, order n >= 2, of r1, 
                              using nb even Bernoulli numbers in rb [3] 
F:MPR  hypergeom_pfq(np,nq,aa,bb,r1)  Hypergeometric pFq function of aa, bb and r1; 
                              dimensions are aa(np) and bb(nq) [5] 
F:MPR  incgamma(r1,r2)      Incomplete gamma function [6] 
F:MPR  polygamma(k,r1)      Polygamma function, order k >= 1, of r1 [4] 
F:MPR  polygamma(nb,rb,k,r1)  Polygamma function, order k >= 1, of r1, using 
                              nb even Bernoulli numbers in rb [3] 
S      polylog_ini(n,ss)    Initialize array ss, of size |n|, for computing
                              polylogarithms of order n when n < 0 [2][6]
F:MPR  polylog_neg(n,ss,r1) Polylogarithm function of r1, for n < 0, using 
                              precomputed data in ss [7]
F:MPR  polylog_pos(n,r1)    Polylogarithm function, order n >= 0, of r1 [8]
F:MPR  struve_hn(n,r1)      StruveH function, order n >= 0, of r1 [9]
F:MPR  zeta(r1)             Zeta function of r1
F:MPR  zeta_be(nb,rb,r1)    Zeta function of r1, using nb even Bernoulli
                              numbers in rb [3]
F:MPR  zeta_int(n)          Zeta function of integer argument n [2]

Notes:
[1]: Only available with MPFUN20-MPFR. 
[2]: In variant 1, an integer precision level argument (mantissa words) may optionally be added as the final argument; this argument is required in variant 2.
[3]: For most applications, set nb > 2X precision in decimal digits; see mpberne above.
[4]: For hurwitz_zetan and polygamma, the argument r1 is limited to the range (0, 1).
[5]: For hypergeom_pfq, the integers np and nq must not exceed 10.
[6]: For incgamma, r1 must not be zero, and must not be negative unless r2 = 0.
[7]: For polylog_ini and polylog_neg, the integer n is limited to the range [-1000, -1]. 
[8]: For polylog_pos, the argument r1 is limited to the range (-1,1). 
[9]: For struve_hn, the argument r1 is limited to the range [-1000, 1000].

3. Miscellaneous I/O, conversion and transcendental functions:

Type   Name                 Description
F:MPC  mpcmplx(r1,r2)       Converts (r1,r2) to MPC [1]   
F:MPC  mpcmplx(dc1)         Converts DC arg to MPC [1]
F:MPC  mpcmplx(z1)          Converts MPC arg to MPC [1]
F:MPC  mpcmplxdc(dc1)       Converts DC to MPC, without checking [1][2]
S      mpcssh(r1,r2,r3)     Returns both cosh and sinh of r1, in the same 
                              time as calling just cosh or just sinh
S      mpcssn(r1,r2,r3)     Returns both cos and sin of r1, in the same
                              time as calling just cos or just sin
S      mpdecmd(r1,d1,i1)    Converts r1 to the form d1*10^i1
S      mpeform(r1,i1,i2,s1) Converts r1 to char(1) string in Ei1.i2
                              suitable for output
S      mpfform(r1,i1,i2,s1) Converts r1 to char(1) string in Fi1.i2
                              suitable for output
F:MPR  mpegamma()           Returns Euler's gamma constant [1]
S      mpinit()             Initializes for extra-high precision [1]
F:MPR  mplog2()             Returns log(2) [1]
F:MPR  mpnrt(r1,i1)         Returns the i1-th root of r1
F:MPR  mppi()               Returns pi [1]
F:MPR  mpprod(r1,d1)        Returns r1*d1, without checking [2]
F:MPR  mpquot(r1,d1)        Returns r1/d1, without checking [2]
F:MPR  mprand(r1)           Returns pseudorandom number, based on r1
                              Start with an irrational, say r1 = mplog2()
                              Typical iterated usage: r1 = mprand(r1)
S      mpread(i1,r1)        Inputs r1 from Fortran unit i1; up to five
                              MPR args may be listed [1]
S      mpread(i1,z1)        Inputs z1 from Fortran unit i1; up to five
                              MPC args may be listed [1]
F:MPR  mpreal(r1)           Converts MPR arg to MPR [1]
F:MPR  mpreal(z1)           Converts MPC arg to MPR [1]
F:MPR  mpreal(d1)           Converts DP arg to MPR [1][2]
F:MPR  mpreal(q1)           Converts QP arg to MPR [1]
F:MPR  mpreal(s1,i1)        Converts char(1) string of length i1 to MPR [1]
F:MPR  mpreal(sn)           Converts char(n) string to MPR [1]
F:MPR  mpreald(d1)          Converts DP arg to MPR, without checking [1][2]
F:MPR  mprealq(q1)          Converts QP arg to MPR, without checking [1][2]
F:Int  mpwprec(r1)          Returns precision in words assigned to r1
F:Int  mpwprec(z1)          Returns precision in words assigned to z1
S      mpwrite(i1,i2,i3,r1) Outputs r1 in Ei2.i3 format to unit i1; up to
                              five MPR args may be listed
S      mpwrite(i1,i2,i3,z1) Outputs z1 in Ei2.i3 format to unit i1; up to
                              five MPC args may be listed
F:QP   qreal(r1)            Converts MPR arg to quad precision (if supported)

Notes:
[1]: In variant 1, an integer precision level argument (mantissa words) may optionally be added as the final argument; this argument is required in variant 2.
[2]: These do not check DP, DC or QP values.
