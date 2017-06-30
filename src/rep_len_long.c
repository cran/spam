#include <stdint.h>

void rep_len_long(double *x, int64_t *lx, int64_t *length, double *out) {
  int64_t i;
  for(i=0; i < *length; i++){
    out[i] = x[i%*lx];
  }
}
