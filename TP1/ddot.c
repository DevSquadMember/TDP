#include <stdio.h>
#include <stdlib.h>

double cblas_ddot(const int n, const double* dx, const int incx, const double* dy, const int incy) {
  int i;
  double res = 0;
  if (incx == 1 && incy == 1) {
    for (i = 0 ; i < n ; ++i) {
      res += dx[i] * dy[i];
    } 
  } else {
    int x = 0, y = 0;
    for (i = 0 ; i < n ; ++i) {
      res += dx[x] * dy[y];
      x += incx;
      y += incy;
    }
  }
  return res;
}

