#include "rupert.h"

void computeAdvance(Rupert& next, Rupert& cur,
                    vector<Number>& grad, Number& scale,
                    int m, int n);

main(int argc, char** argv) {
  int m,n;

  m = atoi(argv[1]);
  n = atoi(argv[2]);

  if (m > n) {
    fprintf(stderr,"m (%d) must be <= n (%d)\n",m,n);
    exit(1);
  }

  Rupert cur(m,n),next(m,n);
  cur.Read();
  cur.Write();

  vector<Number> grad((n*(n-1))/2);
  Number tmp,mag;
  double dt = 1e-45;

  for (int k = 0; k < 200; k++) {
    int index = 0;
    mpfr_set_d(mag.Get(),0.0,GMP_RNDN);
    for (int i = 0; i < n; i++) {
      for (int j = i+1; j < n; j++) {
        cur.Gradient(grad[index],i,j,dt);
        printf("%3d %3d   ",i,j);
        mpfr_out_str(stdout,10,0,grad[index].Get(),GMP_RNDN);
        printf("\n");
        fflush(stdout);
        mpfr_mul(tmp.Get(),grad[index].Get(),grad[index].Get(),GMP_RNDN);
        mpfr_add(mag.Get(),mag.Get(),tmp.Get(),GMP_RNDN);

        index++;
      }
    }
    mpfr_sqrt(mag.Get(),mag.Get(),GMP_RNDN);
    printf("Length:   ");
    mpfr_out_str(stdout,10,0,mag.Get(),GMP_RNDN);
    printf("\n");
    printf("\n");

    Number sMin,sMinSum;
    Number sCen,sCenSum;
    Number sMax,sMaxSum;

    mpfr_set_d(tmp.Get(),0.0,GMP_RNDN);
    if (mpfr_cmp(mag.Get(),tmp.Get()) == 0) {
#if 0
      if (k == 0) {
        mpfr_set_d(tmp.Get(),1e-60,GMP_RNDN);

        for (int i = 0; i < n; i++) {
          for (int j = i+1; j < n; j++) {
            cur.Transform(next,i,j,tmp);
            cur = next;
          }
        }

        cur.Compute();

        continue;
      } else {
#endif
        break;
#if 0
      }
#endif
    } else {
      mpfr_set_d(sCen.Get(),0.0,GMP_RNDN);
      sCenSum = cur.GetMaxSum();

      mpfr_set_d(sMin.Get(),-0.0001,GMP_RNDN);
      computeAdvance(next,cur,grad,sMin,m,n);

      next.Compute();
      sMinSum = next.GetMaxSum();

      mpfr_set_d(sMax.Get(),0.0001,GMP_RNDN);
      computeAdvance(next,cur,grad,sMax,m,n);

      next.Compute();
      sMaxSum = next.GetMaxSum();

      Number oldDelta1,oldDelta2;
      mpfr_set_d(oldDelta1.Get(),0.0,GMP_RNDN);
      mpfr_set_d(oldDelta2.Get(),0.0,GMP_RNDN);

      for (int i = 0; i < 1000; i++) {
        Number delta1,delta2;

#if 0
        mpfr_out_str(stdout,10,0,sMin.Get(),GMP_RNDN);
        printf("\n");
        printf("   ");
        mpfr_out_str(stdout,10,0,sMinSum.Get(),GMP_RNDN);
        printf("\n");

        mpfr_out_str(stdout,10,0,sCen.Get(),GMP_RNDN);
        printf("\n");
        printf("   ");
        mpfr_out_str(stdout,10,0,sCenSum.Get(),GMP_RNDN);
        printf("\n");

        mpfr_out_str(stdout,10,0,sMax.Get(),GMP_RNDN);
        printf("\n");
        printf("   ");
        mpfr_out_str(stdout,10,0,sMaxSum.Get(),GMP_RNDN);
        printf("\n");
        printf("\n");
#endif

        mpfr_sub(delta1.Get(),sCen.Get(),sMin.Get(),GMP_RNDN);
        mpfr_sub(delta2.Get(),sMax.Get(),sCen.Get(),GMP_RNDN);

        if (mpfr_cmp(delta1.Get(),oldDelta1.Get()) == 0 &&
            mpfr_cmp(delta2.Get(),oldDelta2.Get()) == 0) {
          printf("Hi 1...\n");
          break;
        }

        oldDelta1 = delta1;
        oldDelta2 = delta2;

        if (mpfr_cmp(delta1.Get(),delta2.Get()) > 0) {
          Number s,sSum;

          mpfr_add(s.Get(),sMin.Get(),sCen.Get(),GMP_RNDN);
          mpfr_div_ui(s.Get(),s.Get(),2,GMP_RNDN);

          computeAdvance(next,cur,grad,s,m,n);
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

            printf("Hi 2...\n");
            break;
          } else {
            sMin = s;
            sMinSum = sSum;
          }
        } else {
          Number s,sSum;

          mpfr_add(s.Get(),sCen.Get(),sMax.Get(),GMP_RNDN);
          mpfr_div_ui(s.Get(),s.Get(),2,GMP_RNDN);

          computeAdvance(next,cur,grad,s,m,n);
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
          
            printf("Hi 3...\n");
            break;
          } else {
            sMax = s;
            sMaxSum = sSum;
          }
        }
      }

      mpfr_out_str(stdout,10,0,sMin.Get(),GMP_RNDN);
      printf("\n");
      printf("   ");
      mpfr_out_str(stdout,10,0,sMinSum.Get(),GMP_RNDN);
      printf("\n");

      mpfr_out_str(stdout,10,0,sCen.Get(),GMP_RNDN);
      printf("\n");
      printf("   ");
      mpfr_out_str(stdout,10,0,sCenSum.Get(),GMP_RNDN);
      printf("\n");

      mpfr_out_str(stdout,10,0,sMax.Get(),GMP_RNDN);
      printf("\n");
      printf("   ");
      mpfr_out_str(stdout,10,0,sMaxSum.Get(),GMP_RNDN);
      printf("\n");
      printf("\n");
      fflush(stdout);

      computeAdvance(next,cur,grad,sCen,m,n);

      cur = next;
      cur.Compute();
    }
  }

  cur.Compute();
  cur.Write();

  exit(0);
}

void computeAdvance(Rupert& next, Rupert& cur,
                    vector<Number>& grad, Number& scale,
                    int m, int n)
{
  vector<Number> gradScale((n*(n-1))/2);
  gradScale = grad;
  int index;

  index = 0;
  for (int i = 0; i < n; i++) {
    for (int j = i+1; j < n; j++) {
      mpfr_mul(gradScale[index].Get(),grad[index].Get(),scale.Get(),GMP_RNDN);
      index++;
    }
  }

  next = cur;

  index = 0;
  for (int i = 0; i < n; i++) {
    for (int j = i+1; j < n; j++) {
      Rupert temp(m,n);

      next.Transform(temp,i,j,gradScale[index]);
      next = temp;

      index++;
    }
  }
}
