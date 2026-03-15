#include "rupert.h"

main(int argc, char** argv) {
  int m,n;
  int i1,i2;

  m = atoi(argv[1]);
  n = atoi(argv[2]);

  if (m > n) {
    fprintf(stderr,"m (%d) must be <= n (%d)\n",m,n);
    exit(1);
  }

  Rupert original(m,n);
  Rupert cur(m,n);
  Rupert next(m,n);

  original.Read();

  cur = original;

  i1 = n-2;
  i2 = n-1;
  for (int k = 0; k < 1000000; k++) {
    int iMin,iMax;

    cur.Compute();

    if (mpfr_cmp(cur.GetSum()[i1].Get(),cur.GetSum()[i2].Get()) < 0) {
      iMin = i1;
      iMax = i2;
    } else {
      iMin = i2;
      iMax = i1;
    }

    Number minSum,maxSum,diffSum;

    minSum = cur.GetSum()[iMin];
    maxSum = cur.GetSum()[iMax];
    mpfr_sub(diffSum.Get(),maxSum.Get(),minSum.Get(),GMP_RNDN);

    Number sMin,sMinSum;
    Number sCen,sCenSum;
    Number sMax,sMaxSum;

    mpfr_set_d(sCen.Get(),0.0,GMP_RNDN);
    sCenSum = cur.GetMaxSum();

    mpfr_set_d(sMin.Get(),-0.0001,GMP_RNDN);
    cur.Transform(next,iMin,iMax,sMin);

    next.Compute();
    sMinSum = next.GetMaxSum();

    mpfr_set_d(sMax.Get(),0.0001,GMP_RNDN);
    cur.Transform(next,iMin,iMax,sMax);

    next.Compute();
    sMaxSum = next.GetMaxSum();

    Number oldDelta1,oldDelta2;
    mpfr_set_d(oldDelta1.Get(),0.0,GMP_RNDN);
    mpfr_set_d(oldDelta2.Get(),0.0,GMP_RNDN);

    for (int i = 0; i < 1000; i++) {
      Number delta1,delta2;

      mpfr_sub(delta1.Get(),sCen.Get(),sMin.Get(),GMP_RNDN);
      mpfr_sub(delta2.Get(),sMax.Get(),sCen.Get(),GMP_RNDN);

      if (mpfr_cmp(delta1.Get(),oldDelta1.Get()) == 0 &&
          mpfr_cmp(delta2.Get(),oldDelta2.Get()) == 0) {
        break;
      }

#if 0
mpfr_out_str(stdout,10,0,sMin.Get(),GMP_RNDN);
printf("\n");
mpfr_out_str(stdout,10,0,sMinSum.Get(),GMP_RNDN);
printf("\n");
mpfr_out_str(stdout,10,0,sCen.Get(),GMP_RNDN);
printf("\n");
mpfr_out_str(stdout,10,0,sCenSum.Get(),GMP_RNDN);
printf("\n");
mpfr_out_str(stdout,10,0,sMax.Get(),GMP_RNDN);
printf("\n");
mpfr_out_str(stdout,10,0,sMaxSum.Get(),GMP_RNDN);
printf("\n");
printf("\n");
mpfr_out_str(stdout,10,0,delta1.Get(),GMP_RNDN);
printf("\n");
mpfr_out_str(stdout,10,0,delta2.Get(),GMP_RNDN);
printf("\n");
printf("\n");
fflush(stdout);
#endif

      oldDelta1 = delta1;
      oldDelta2 = delta2;

      if (mpfr_cmp(delta1.Get(),delta2.Get()) > 0) {
        Number s,sSum;

        mpfr_add(s.Get(),sMin.Get(),sCen.Get(),GMP_RNDN);
        mpfr_div_ui(s.Get(),s.Get(),2,GMP_RNDN);

        cur.Transform(next,iMin,iMax,s);
        next.Compute();

        sSum = next.GetMaxSum();

        if (mpfr_cmp(sSum.Get(),sCenSum.Get()) < 0) {
          sMax = sCen;
          sMaxSum = sCenSum;

          sCen = s;
          sCenSum = sSum;
        } else if (mpfr_cmp(sSum.Get(),sCenSum.Get()) == 0) {
          sMax = sCen;
          sMaxSum = sCenSum;

          sCen = s;
          sCenSum = sSum;
      
          break;
        } else {
          sMin = s;
          sMinSum = sSum;
        }
      } else {
        Number s,sSum;

        mpfr_add(s.Get(),sCen.Get(),sMax.Get(),GMP_RNDN);
        mpfr_div_ui(s.Get(),s.Get(),2,GMP_RNDN);

        cur.Transform(next,iMin,iMax,s);
        next.Compute();

        sSum = next.GetMaxSum();

        if (mpfr_cmp(sSum.Get(),sCenSum.Get()) < 0) {
          sMin = sCen;
          sMinSum = sCenSum;

          sCen = s;
          sCenSum = sSum;
        } else if (mpfr_cmp(sSum.Get(),sCenSum.Get()) == 0) {
          sMin = sCen;
          sMinSum = sCenSum;

          sCen = s;
          sCenSum = sSum;
      
          break;
        } else {
          sMax = s;
          sMaxSum = sSum;
        }
      }
    }

    printf("Iter: %d %d ",iMin,iMax);
    mpfr_out_str(stdout,10,0,sCen.Get(),GMP_RNDN);
    printf("\n");
    printf("          ");
    mpfr_out_str(stdout,10,0,sCenSum.Get(),GMP_RNDN);
    printf("\n");
    fflush(stdout);

    cur.Transform(next,iMin,iMax,sCen);
    next.Compute();

    cur = next;

    i2 = (i2 + 1) % n;

    if (i2 == i1) {
      i1 = (i1 + 1) % n;
      i2 = (i1 + 1) % n;
    }
  }

  cur.Write();

#if 0
  printData(o,m,n);
#endif
#if 0
  printData(a,m,n);
#endif

#if 0
  minSum = 1.0/sqrt(9.0/8.0);
  printf("%30.27Lf\n",minSum);
#endif

  exit(0);
}
