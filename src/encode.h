#ifndef _SCLOUDPLUS_ENCODE_H_
#define _SCLOUDPLUS_ENCODE_H_
#include <stdint.h>
void scloudplus_msgencode(const uint8_t *msg, uint16_t *matrixM);
void scloudplus_msgdecode(const uint16_t *matrixM, uint8_t *msg);
void scloudplus_packsk(uint16_t *S, uint8_t *sk);
void scloudplus_unpacksk(uint8_t *sk, uint16_t *S);
void scloudplus_packpk(uint16_t *B, uint8_t *pk);
void scloudplus_unpackpk(uint8_t *pk, uint16_t *B);
void scloudplus_compressc1(uint16_t *C, uint16_t *out);
void scloudplus_decompressc1(uint16_t *in, uint16_t *C);
void scloudplus_compressc2(uint16_t *C, uint16_t *out);
void scloudplus_decompressc2(uint16_t *in, uint16_t *C);
void scloudplus_packc1(uint16_t *C, uint8_t *out);
void scloudplus_unpackc1(uint8_t *in, uint16_t *C);
void scloudplus_packc2(uint16_t *C, uint8_t *out);
void scloudplus_unpackc2(uint8_t *in, uint16_t *C);
#endif