#ifndef _SCLOUDPLUS_SAMPLE_H_
#define _SCLOUDPLUS_SAMPLE_H_
#include "aes.h"
#include "param.h"
void scloudplus_mul_add_as_e(const uint8_t *seedA, const uint16_t *S,
							 const uint16_t *E, uint16_t *B);
void scloudplus_mul_add_sa_e(const uint8_t *seedA, const uint16_t *S,
							 uint16_t *E, uint16_t *C);
void scloudplus_sampleeta1(uint8_t *seed, uint16_t *matrixE);
void scloudplus_sampleeta2(uint8_t *seed, uint16_t *matrixe1,
						   uint16_t *matrixe2);
void scloudplus_samplepsi(uint8_t *seed, uint16_t *matrixs);
void scloudplus_samplephi(uint8_t *seed, uint16_t *matrixs);
#endif