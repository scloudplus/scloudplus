#ifndef _SCLOUDPLUS_UTIL_H_
#define _SCLOUDPLUS_UTIL_H_
#include <stddef.h>
#include <stdint.h>
int8_t scloudplus_verify(const uint8_t *a, const uint8_t *b, size_t len);
void scloudplus_cmov(uint8_t *r, const uint8_t *a, const uint8_t *b, size_t len,
					 int8_t bl);
#endif