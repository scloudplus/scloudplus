#include "pke.h"
#include "param.h"
#include "encode.h"
#include "matrix.h"
#include "sample.h"
#include "fips202.h"
#include "random.h"
#include <stdlib.h>
#include <string.h>
void scloudplus_pkekeygen(uint8_t *pk, uint8_t *sk)
{
	uint16_t *S =
		(uint16_t *)malloc(sizeof(uint16_t) * scloudplus_n * scloudplus_nbar);
	uint16_t *E =
		(uint16_t *)malloc(sizeof(uint16_t) * scloudplus_m * scloudplus_nbar);
	uint16_t *B =
		(uint16_t *)malloc(sizeof(uint16_t) * scloudplus_m * scloudplus_nbar);
	uint8_t alpha[32], seed[80];
	uint8_t *seedA = seed;
	uint8_t *r1 = seed + 16;
	uint8_t *r2 = seed + 48;
	randombytes(alpha, 32);
	scloudplus_F(seed, 80, alpha, 32);
	scloudplus_samplepsi(r1, S);
	scloudplus_sampleeta1(r2, E);
	scloudplus_mul_add_as_e(seedA, S, E, B);
	scloudplus_packpk(B, pk);
	memcpy(pk + scloudplus_pk - 16, seedA, 16);
	scloudplus_packsk(S, sk);
	free(S);
	free(E);
	free(B);
}

void scloudplus_pkeenc(uint8_t *pk, uint8_t *m, uint8_t *r, uint8_t *ctx)
{
	uint16_t *S1 =
		(uint16_t *)malloc(sizeof(uint16_t) * scloudplus_mbar * scloudplus_m);
	uint16_t *E1 =
		(uint16_t *)malloc(sizeof(uint16_t) * scloudplus_mbar * scloudplus_n);
	uint16_t *E2 =
		(uint16_t *)malloc(sizeof(uint16_t) * scloudplus_mbar * scloudplus_nbar);
	uint16_t *mu0 =
		(uint16_t *)malloc(sizeof(uint16_t) * scloudplus_mbar * scloudplus_nbar);
	uint16_t *C1 =
		(uint16_t *)malloc(sizeof(uint16_t) * scloudplus_mbar * scloudplus_n);
	uint16_t *C2 =
		(uint16_t *)malloc(sizeof(uint16_t) * scloudplus_mbar * scloudplus_nbar);
	uint16_t *B =
		(uint16_t *)malloc(sizeof(uint16_t) * scloudplus_m * scloudplus_nbar);
	uint8_t seed[64];
	uint8_t *seedA = pk + scloudplus_pk - 16;
	uint8_t *r1 = seed;
	uint8_t *r2 = seed + 32;
	scloudplus_F(seed, 64, r, 32);
	scloudplus_samplephi(r1, S1);
	scloudplus_sampleeta2(r2, E1, E2);
	scloudplus_msgencode(m, mu0);
	scloudplus_unpackpk(pk, B);
	scloudplus_mul_add_sa_e(seedA, S1, E1, C1);
	scloudplus_mul_add_sb_e(S1, B, E2, C2);
	scloudplus_add(C2, mu0, scloudplus_mbar * scloudplus_nbar, C2);
	scloudplus_compressc1(C1, C1);
	scloudplus_compressc2(C2, C2);
	scloudplus_packc1(C1, ctx);
	scloudplus_packc2(C2, ctx + scloudplus_c1);
	free(S1);
	free(E1);
	free(E2);
	free(mu0);
	free(C1);
	free(C2);
	free(B);
}

void scloudplus_pkedec(uint8_t *sk, uint8_t *ctx, uint8_t *m)
{
	uint16_t *S =
		(uint16_t *)malloc(sizeof(uint16_t) * scloudplus_n * scloudplus_nbar);
	uint16_t *C1 =
		(uint16_t *)malloc(sizeof(uint16_t) * scloudplus_mbar * scloudplus_n);
	uint16_t *C2 =
		(uint16_t *)malloc(sizeof(uint16_t) * scloudplus_mbar * scloudplus_nbar);
	uint16_t *D =
		(uint16_t *)malloc(sizeof(uint16_t) * scloudplus_mbar * scloudplus_nbar);
	scloudplus_unpacksk(sk, S);
	scloudplus_unpackc1(ctx, C1);
	scloudplus_unpackc2(ctx + scloudplus_c1, C2);
	scloudplus_decompressc1(C1, C1);
	scloudplus_decompressc2(C2, C2);
	scloudplus_mul_cs(C1, S, D);
	scloudplus_sub(C2, D, scloudplus_mbar * scloudplus_nbar, D);
	scloudplus_msgdecode(D, m);
	free(S);
	free(C1);
	free(C2);
	free(D);
}
