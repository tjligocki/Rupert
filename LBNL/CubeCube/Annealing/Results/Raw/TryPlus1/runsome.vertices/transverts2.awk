BEGIN {
  i = 0;
}
i < NF {
  i = NR-1
  n = NF;
  for (j = 0; j < n; j++) {
    matrix[i,j] = $(j+1);
  }
  next;
}
NF == 1 {
  max = $1;
}
END {
  for (i = 0; i < n; i++) {
    printf("    %2d   ",i+1);
  }
  printf("\n\n");

  for (i = 0; i < n; i++) {
    v[i] = 0.0;
  }

  for (i = 0; i < m; i++) {
    v[i] = -1.0;
  }

  transVerts(matrix,v,w0);
  printVerts(w0,n);
  print "";

  for (i = 0; i < m; i++) {
    v[i] = 1.0;

    transVerts(matrix,v,w);
    printVerts(w,n);

    v[i] = -1.0;
  }
  print "";

  for (i = 0; i < m; i++) {
    v[i] = 1.0;

    transVerts(matrix,v,w);
    subVerts(w,w0);
    printVerts(w,n);

    v[i] = -1.0;
  }
}
function genVerts(v,n,m,matrix) {
  v[m] = -1.0;

  if (m == 0) {
    transVerts(matrix,v,w);
    printVerts(w,n);
  } else {
    genVerts(v,n,m-1,matrix);
  }

  v[m] = 1.0;

  if (m == 0) {
    transVerts(matrix,v,w);
    printVerts(w,n);
  } else {
    genVerts(v,n,m-1,matrix);
  }
}
function transVerts(matrix,v,w, i) {
  for (i = 0; i < n; i++) {
    w[i] = 0.0;
    for (j = 0; j < n; j++) {
      w[i] = w[i] + matrix[i,j]*v[j];
    }
  }
}
function subVerts(v,w, i) {
  for (i = 0; i < n; i++) {
    v[i] = v[i] - w[i];
  }
}
function printVerts(v,n, i) {
  for (i = 0; i < n; i++) {
    printf("%8.5lf ",v[i]/max);
  }
  printf("\n");
}
