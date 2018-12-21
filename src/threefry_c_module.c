#include <threefry.h>
#include <math.h>

#define DOUBLE_MULT 5.421010862427522e-20

int threefry_int32(threefry2x64_ctr_t *ctr, threefry2x64_key_t *key) {
  threefry2x64_ctr_t rand = threefry2x64(*ctr, *key);
  ctr->v[0]++;
  return rand.v[0];
}

double threefry_c_double(threefry2x64_ctr_t *ctr, threefry2x64_key_t *key) {
  threefry2x64_ctr_t rand = threefry2x64(*ctr, *key);
  ctr->v[0]++;
  return rand.v[0]*DOUBLE_MULT;
}

double threefry_c_normal(threefry2x64_ctr_t *ctr, threefry2x64_key_t *key) {
  threefry2x64_ctr_t rand;
  double u1, u2, radius;
  int flag;

  flag = 1;
  while (flag) {
    rand = threefry2x64(*ctr, *key);
    ctr->v[0]++;
    u1 = 2*rand.v[0]*DOUBLE_MULT - 1;
    u2 = 2*rand.v[1]*DOUBLE_MULT - 1;
    radius = u1*u1 + u2*u2;
    if (radius > 0 && radius < 1) flag = 0;
  }

  return u1 * sqrt(-2*log(radius)/radius) ;
}
