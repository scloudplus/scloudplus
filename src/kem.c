#include "../include/kem.h"
#include "../include/pke.h"
#include "../include/param.h"
#include "../include/hash.h"
#include "../include/util.h"
#include "../include/random.h"
#include <string.h>
#include <stdlib.h>
void scloud_kemkeygen(uint8_t *pk, uint8_t *sk)
{
	uint8_t z[32];
	randombytes(z, 32);
	scloudplus_pkekeygen(pk, sk);
	memcpy(sk + scloudplus_pke_sk, pk, scloudplus_pk);
	scloudplus_H(sk + scloudplus_pke_sk + scloudplus_pk, pk, scloudplus_pk);
	memcpy(sk + scloudplus_kem_sk - 32, z, 32);
}
void scloud_kemencaps(uint8_t *pk, uint8_t *ctx, uint8_t *ss)
{
	uint8_t *kc = (uint8_t *)malloc(sizeof(uint8_t) * (scloudplus_ctx + 32));
	uint8_t m[scloudplus_ss + 32], rk[64];
	randombytes(m, scloudplus_ss);
	scloudplus_H(m + scloudplus_ss, pk, scloudplus_pk);
	scloudplus_G(rk, m, scloudplus_ss + 32);
	scloudplus_pkeenc(pk, m, rk, ctx);
	memcpy(kc, rk + 32, 32);
	memcpy(kc + 32, ctx, scloudplus_ctx);
	scloudplus_K(ss, scloudplus_ss, kc, scloudplus_ctx + 32);
}
void scloud_kemdecaps(uint8_t *sk, uint8_t *ctx, uint8_t *ss)
{
	uint8_t m[scloudplus_ss + 32], rk[64];
	uint8_t *ctx1 = (uint8_t *)malloc(sizeof(uint8_t) * (scloudplus_ctx + 32));
	scloudplus_pkedec(sk, ctx, m);
	memcpy(m + scloudplus_ss, sk + scloudplus_pke_sk + scloudplus_pk, 32);
	scloudplus_G(rk, m, scloudplus_ss + 32);
	scloudplus_pkeenc(sk + scloudplus_pke_sk, m, rk, ctx1 + 32);
	int8_t bl = scloudplus_verify(ctx, ctx1 + 32, scloudplus_ctx);
	memcpy(ctx1 + 32, ctx, scloudplus_ctx);
	scloudplus_cmov(ctx1, rk + 32, sk + scloudplus_kem_sk - 32, 32, bl);
	scloudplus_K(ss, scloudplus_ss, ctx1, scloudplus_ctx + 32);
}