#include "rupert.h"

void matrixMultiply(Rupert& m1, Rupert& m2, Rupert& m3, int n)
{
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      mpfr_set_d(m1.GetMatrix()[i][j].Get(),0.0,GMP_RNDN);

      for (int k = 0; k < n; k++) {
        Number tmp;

        mpfr_mul(tmp.Get(),
                 m2.GetMatrix()[i][k].Get(),
                 m3.GetMatrix()[k][j].Get(),GMP_RNDN);
        mpfr_add(m1.GetMatrix()[i][j].Get(),
                 m1.GetMatrix()[i][j].Get(),
                 tmp.Get(),GMP_RNDN);
      }
    }
  }
}

main(int argc, char** argv) { int m,n;

  m = 2;
  n = 3;

  if (m > n) {
    fprintf(stderr,"m (%d) must be <= n (%d)\n",m,n);
    exit(1);
  }

  Rupert best(m,n);
  best.Read();
  best.Compute();
  best.Write();

  Rupert exact(m,n);

  Number tmp1,tmp2,tmp3,tmp4;

  mpfr_set_ui(tmp1.Get(),8,GMP_RNDN);
  mpfr_div_ui(tmp1.Get(),tmp1.Get(),9,GMP_RNDN);
  mpfr_sqrt(tmp1.Get(),tmp1.Get(),GMP_RNDN);
  mpfr_set(exact.GetMatrix()[0][0].Get(),tmp1.Get(),GMP_RNDN);

  mpfr_set_d(exact.GetMatrix()[0][1].Get(),0.0,GMP_RNDN);

  mpfr_set_ui(tmp1.Get(),1,GMP_RNDN);
  mpfr_div_ui(tmp1.Get(),tmp1.Get(),3,GMP_RNDN);
  mpfr_set(exact.GetMatrix()[0][2].Get(),tmp1.Get(),GMP_RNDN);

  mpfr_set_ui(tmp1.Get(),5,GMP_RNDN);
  mpfr_div_ui(tmp1.Get(),tmp1.Get(),90,GMP_RNDN);
  mpfr_sqrt(tmp1.Get(),tmp1.Get(),GMP_RNDN);
  mpfr_set(exact.GetMatrix()[1][0].Get(),tmp1.Get(),GMP_RNDN);

  mpfr_set_ui(tmp1.Get(),1,GMP_RNDN);
  mpfr_div_ui(tmp1.Get(),tmp1.Get(),2,GMP_RNDN);
  mpfr_sqrt(tmp1.Get(),tmp1.Get(),GMP_RNDN);
  mpfr_set(exact.GetMatrix()[1][1].Get(),tmp1.Get(),GMP_RNDN);

  mpfr_set_ui(tmp1.Get(),2,GMP_RNDN);
  mpfr_div_ui(tmp1.Get(),tmp1.Get(),3,GMP_RNDN);
  mpfr_neg(tmp1.Get(),tmp1.Get(),GMP_RNDN);
  mpfr_set(exact.GetMatrix()[1][2].Get(),tmp1.Get(),GMP_RNDN);

  mpfr_set_ui(tmp1.Get(),5,GMP_RNDN);
  mpfr_div_ui(tmp1.Get(),tmp1.Get(),90,GMP_RNDN);
  mpfr_sqrt(tmp1.Get(),tmp1.Get(),GMP_RNDN);
  mpfr_neg(tmp1.Get(),tmp1.Get(),GMP_RNDN);
  mpfr_set(exact.GetMatrix()[2][0].Get(),tmp1.Get(),GMP_RNDN);

  mpfr_set_ui(tmp1.Get(),1,GMP_RNDN);
  mpfr_div_ui(tmp1.Get(),tmp1.Get(),2,GMP_RNDN);
  mpfr_sqrt(tmp1.Get(),tmp1.Get(),GMP_RNDN);
  mpfr_set(exact.GetMatrix()[2][1].Get(),tmp1.Get(),GMP_RNDN);

  mpfr_set_ui(tmp1.Get(),2,GMP_RNDN);
  mpfr_div_ui(tmp1.Get(),tmp1.Get(),3,GMP_RNDN);
  mpfr_set(exact.GetMatrix()[2][2].Get(),tmp1.Get(),GMP_RNDN);

  exact.Compute();
  // exact.Write();

  printf("\n");
  printf("Best:   ");
  mpfr_out_str(stdout,10,0,best.GetMaxSum().Get(),GMP_RNDN);
  printf("\n");

  // exact.Write();

  printf("Exact:  ");
  mpfr_out_str(stdout,10,0,exact.GetMaxSum().Get(),GMP_RNDN);
  printf("\n");

  Number diff;
  mpfr_sub(diff.Get(),best.GetMaxSum().Get(),exact.GetMaxSum().Get(),GMP_RNDN);

  printf("Diff:   ");
  mpfr_out_str(stdout,10,0,diff.Get(),GMP_RNDN);
  printf("\n");

  Rupert invBest(m,n);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      mpfr_set(invBest.GetMatrix()[j][i].Get(),
               best   .GetMatrix()[i][j].Get(),GMP_RNDN);
    }
  }

  // invBest.Write();

  Rupert delta(m,n);

  matrixMultiply(delta,exact,invBest,n);
  delta.Write();

  Rupert direct(m,n);

  matrixMultiply(direct,delta,best,n);
  direct.Compute();

  printf("\n");
  printf("Direct: ");
  mpfr_out_str(stdout,10,0,direct.GetMaxSum().Get(),GMP_RNDN);
  printf("\n");

  // exact.Write();

  printf("Exact:  ");
  mpfr_out_str(stdout,10,0,exact.GetMaxSum().Get(),GMP_RNDN);
  printf("\n");

  mpfr_sub(diff.Get(),direct.GetMaxSum().Get(),exact.GetMaxSum().Get(),GMP_RNDN);

  printf("Diff:   ");
  mpfr_out_str(stdout,10,0,diff.Get(),GMP_RNDN);
  printf("\n");
  printf("\n");

  Number angle12;

  mpfr_set(tmp1.Get(),delta.GetMatrix()[1][1].Get(),GMP_RNDN);
  mpfr_neg(tmp2.Get(),delta.GetMatrix()[1][2].Get(),GMP_RNDN);

  mpfr_atan2(angle12.Get(),tmp2.Get(),tmp1.Get(),GMP_RNDN);
  mpfr_neg(angle12.Get(),angle12.Get(),GMP_RNDN);

  Rupert rotate12(m,n);

  delta.Transform(rotate12,1,2,angle12);
  // rotate12.Write();

  Number angle02;

  mpfr_set(tmp1.Get(),rotate12.GetMatrix()[0][0].Get(),GMP_RNDN);
  mpfr_neg(tmp2.Get(),rotate12.GetMatrix()[0][2].Get(),GMP_RNDN);

  mpfr_atan2(angle02.Get(),tmp2.Get(),tmp1.Get(),GMP_RNDN);
  mpfr_neg(angle02.Get(),angle02.Get(),GMP_RNDN);

  Rupert rotate02(m,n);

  rotate12.Transform(rotate02,0,2,angle02);
  // rotate02.Write();

  Number angle01;

  mpfr_set(tmp1.Get(),rotate02.GetMatrix()[0][0].Get(),GMP_RNDN);
  mpfr_neg(tmp2.Get(),rotate02.GetMatrix()[0][1].Get(),GMP_RNDN);

  mpfr_atan2(angle01.Get(),tmp2.Get(),tmp1.Get(),GMP_RNDN);
  mpfr_neg(angle01.Get(),angle01.Get(),GMP_RNDN);

  Rupert rotate01(m,n);

  rotate02.Transform(rotate01,0,1,angle01);
  // rotate01.Write();

  Rupert temp(best);
  Rupert bester(best);
  // bester.Write();

  printf("Start:  ");
  mpfr_out_str(stdout,10,0,bester.GetMaxSum().Get(),GMP_RNDN);
  printf("\n");

  mpfr_neg(angle01.Get(),angle01.Get(),GMP_RNDN);
  bester.Transform(temp,0,1,angle01);
  bester = temp;
  bester.Compute();
  // bester.Write();

  printf("Rot01:  ");
  mpfr_out_str(stdout,10,0,bester.GetMaxSum().Get(),GMP_RNDN);
  printf("\n");

  mpfr_neg(angle02.Get(),angle02.Get(),GMP_RNDN);
  bester.Transform(temp,0,2,angle02);
  bester = temp;
  bester.Compute();
  // bester.Write();

  printf("Rot02:  ");
  mpfr_out_str(stdout,10,0,bester.GetMaxSum().Get(),GMP_RNDN);
  printf("\n");

  mpfr_neg(angle12.Get(),angle12.Get(),GMP_RNDN);

  bester.Transform(temp,1,2,angle12);
  bester = temp;
  bester.Compute();
  // bester.Write();

  printf("Rot12:  ");
  mpfr_out_str(stdout,10,0,bester.GetMaxSum().Get(),GMP_RNDN);
  printf("\n");
  printf("\n");

  // exact.Write();
  printf("Exact:  ");
  mpfr_out_str(stdout,10,0,exact.GetMaxSum().Get(),GMP_RNDN);
  printf("\n");
  printf("\n");

  mpfr_sub(diff.Get(),bester.GetMaxSum().Get(),exact.GetMaxSum().Get(),GMP_RNDN);

  printf("Diff:   ");
  mpfr_out_str(stdout,10,0,diff.Get(),GMP_RNDN);
  printf("\n");
  printf("\n");

  printf("Angle 0,1: ");
  mpfr_out_str(stdout,10,0,angle01.Get(),GMP_RNDN);
  printf("\n");

  printf("Angle 0,2: ");
  mpfr_out_str(stdout,10,0,angle02.Get(),GMP_RNDN);
  printf("\n");

  printf("Angle 1,2: ");
  mpfr_out_str(stdout,10,0,angle12.Get(),GMP_RNDN);
  printf("\n");
  printf("\n");

  vector<Number> grad((n*(n-1))/2);

  printf("Best gradient:\n");
  int ng=0;
  for (int i = 0; i < n; i++) {
    for (int j = i+1; j < n; j++) {
      Number deriv;

      double dt = 1e-45;
      best.Gradient(deriv,i,j,dt);

      grad[ng] = deriv;
      ng++;
    
      printf("  %3d %3d %5.0le  ",i,j,dt);
      mpfr_out_str(stdout,10,0,deriv.Get(),GMP_RNDN);
      printf("\n");
    }
  }
  printf("\n");

  printf("Bester gradient:\n");
  for (int i = 0; i < n; i++) {
    for (int j = i+1; j < n; j++) {
      Number deriv;

      double dt = 1e-45;
      bester.Gradient(deriv,i,j,dt);
    
      printf("  %3d %3d %5.0le  ",i,j,dt);
      mpfr_out_str(stdout,10,0,deriv.Get(),GMP_RNDN);
      printf("\n");
    }
  }
  printf("\n");

  printf("Exact gradient:\n");
  for (int i = 0; i < n; i++) {
    for (int j = i+1; j < n; j++) {
      Number deriv;

      double dt = 1e-45;
      exact.Gradient(deriv,i,j,dt);
    
      printf("  %3d %3d %5.0le  ",i,j,dt);
      mpfr_out_str(stdout,10,0,deriv.Get(),GMP_RNDN);
      printf("\n");
    }
  }
  printf("\n");

  printf("Best graph:\n");
  double range = M_PI;
  for (double dt = -range; dt <= 1.001*range; dt += 2.0*range/1000.0) {
    Number scale;
    vector<Number> gradScale((n*(n-1))/2);

    mpfr_set_d(scale.Get(),dt,GMP_RNDN);

    for (int i = 0; i < ng; i++) {
      mpfr_mul(gradScale[i].Get(),grad[i].Get(),scale.Get(),GMP_RNDN);
    }

    temp = best;
    ng = 0;
    for (int i = 0; i < n; i++) {
      for (int j = i+1; j < n; j++) {
        Rupert temp2(m,n);

        temp.Transform(temp2,i,j,gradScale[ng]);

        temp = temp2;
        ng++;
      }
    }

    temp.Compute();

    printf("%13.6le ",dt);
    Number maxDiff;
    mpfr_sub(maxDiff.Get(),temp.GetMaxSum().Get(),best.GetMaxSum().Get(),GMP_RNDN);
    mpfr_out_str(stdout,10,0,maxDiff.Get(),GMP_RNDN);
    printf("\n");
  }
  printf("\n");

  printf("Best surface:\n");
  for (double dt1 = -range; dt1 <= 1.001*range; dt1 += 2.0*range/100.0) {
    for (double dt2 = -range; dt2 <= 1.001*range; dt2 += 2.0*range/100.0) {
      Number angle01;
      Number angle02;
      Number angle12;

      mpfr_set_d(angle01.Get(),dt1,GMP_RNDN);
      mpfr_set_d(angle02.Get(),0.0,GMP_RNDN);
      mpfr_set_d(angle12.Get(),dt2,GMP_RNDN);

      temp = best;

      Rupert temp2(m,n);

      temp.Transform(temp2,0,1,angle01);
      temp = temp2;

      temp.Transform(temp2,0,2,angle02);
      temp = temp2;

      temp.Transform(temp2,1,2,angle12);
      temp = temp2;

      temp.Compute();

      printf("%13.6le %13.6le ",dt1,dt2);
      Number maxDiff;
      mpfr_sub(maxDiff.Get(),temp.GetMaxSum().Get(),best.GetMaxSum().Get(),GMP_RNDN);
      mpfr_out_str(stdout,10,0,maxDiff.Get(),GMP_RNDN);
      printf("\n");
    }
    printf("\n");
  }
  printf("\n");

  printf("Best volume:\n");
  for (double dt1 = -range; dt1 <= 1.001*range; dt1 += 2.0*range/100.0) {
    for (double dt2 = -range; dt2 <= 1.001*range; dt2 += 2.0*range/100.0) {
      for (double dt3 = -range; dt3 <= 1.001*range; dt3 += 2.0*range/100.0) {
        Number angle01;
        Number angle02;
        Number angle12;

        mpfr_set_d(angle01.Get(),dt1,GMP_RNDN);
        mpfr_set_d(angle02.Get(),dt2,GMP_RNDN);
        mpfr_set_d(angle12.Get(),dt3,GMP_RNDN);

        temp = best;

        Rupert temp2(m,n);

        temp.Transform(temp2,0,1,angle01);
        temp = temp2;

        temp.Transform(temp2,0,2,angle02);
        temp = temp2;

        temp.Transform(temp2,1,2,angle12);
        temp = temp2;

        temp.Compute();

        // printf("%13.6le %13.6le %13.6le ",dt1,dt2,dt3);
        Number maxDiff;
        mpfr_sub(maxDiff.Get(),temp.GetMaxSum().Get(),best.GetMaxSum().Get(),GMP_RNDN);
        mpfr_out_str(stdout,10,0,maxDiff.Get(),GMP_RNDN);
        printf("\n");
      }
      // printf("\n");
    }
    // printf("\n");
  }
  printf("\n");

  exit(0);
}
