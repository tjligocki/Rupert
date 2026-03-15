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
    v[i] = 0.0;
  }

  for (i = 0; i < m; i++) {
    v[i] = -1.0;
  }

  transVerts(matrix,v,w);
  printVerts(w,n);

  for (i = 0; i < m; i++) {
    v[i] = 1.0;

    transVerts(matrix,v,w);
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
function printVerts(v,n, i) {
  for (i = 0; i < n; i++) {
    printf("%10.7lf ",v[i]/max);
  }
  printf("\n");
}
