#ifndef _SCLOUDPLUS_MATRIX_H_
#define _SCLOUDPLUS_MATRIX_H_
#include <stdint.h>

void scloudplus_add(uint16_t *in0, uint16_t *in1, int len, uint16_t *out);

void scloudplus_sub(uint16_t *in0, uint16_t *in1, int len, uint16_t *out);

void scloudplus_mul_cs(uint16_t *C, uint16_t *S, uint16_t *out);

void scloudplus_mul_add_sb_e(const uint16_t *S, const uint16_t *B,
							 const uint16_t *E, uint16_t *out);
#endif
