#include <threefry.h>

#define DOUBLE_MULT 5.421010862427522e-20

int threefry_int32(threefry2x64_ctr_t *ctr, threefry2x64_key_t *key) {
  threefry2x64_ctr_t rand = threefry2x64(*ctr, *key);
  return rand.v[0];
}

double threefry_c_double(threefry2x64_ctr_t *ctr, threefry2x64_key_t *key) {
  threefry2x64_ctr_t rand = threefry2x64(*ctr, *key);
  return rand.v[0]*DOUBLE_MULT;
}
