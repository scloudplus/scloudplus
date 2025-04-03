#ifndef _SCLOUDPLUS_PKE_H_
#define _SCLOUDPLUS_PKE_H_
#include <stdint.h>
void scloudplus_pkekeygen(uint8_t *pk, uint8_t *sk);
void scloudplus_pkeenc(uint8_t *pk, uint8_t *m, uint8_t *r, uint8_t *ctx);
void scloudplus_pkedec(uint8_t *sk, uint8_t *ctx, uint8_t *m);
#endif