#ifndef _SCLOUDPLUS_KEM_H_
#define _SCLOUDPLUS_KEM_H_
#include <stdint.h>
void scloud_kemkeygen(uint8_t *pk, uint8_t *sk);
void scloud_kemencaps(uint8_t *pk, uint8_t *ctx, uint8_t *ss);
void scloud_kemdecaps(uint8_t *sk, uint8_t *ctx, uint8_t *ss);
#endif