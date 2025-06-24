#include "../include/ds_benchmark.h"
#include "../include/encode.h"
#include "../include/kem.h"
#include "../include/random.h"
#include "../include/param.h"
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define KEM_TEST_ITERATIONS 100
#define KEM_BENCH_SECONDS 1
#if (scloudplus_l == 128)
#define SYSTEM_NAME "scloudplus 128"
#elif (scloudplus_l == 192)
#define SYSTEM_NAME "scloudplus 192"
#elif (scloudplus_l == 256)
#define SYSTEM_NAME "scloudplus 256"
#endif
static int kem_test(const char *named_parameters, int iterations)
{
	uint8_t pk[scloudplus_pk];
	uint8_t sk[scloudplus_kem_sk];
	uint8_t ctx[scloudplus_ctx];
	uint8_t ssa[scloudplus_ss];
	uint8_t ssb[scloudplus_ss];
	scloud_kemkeygen(pk, sk);
	scloud_kemencaps(pk, ctx, ssa);
	scloud_kemdecaps(sk, ctx, ssb);

	printf("====================================================================="
		   "========================================================\n");
	printf("Testing correctness of key encapsulation mechanism (KEM),system %s,"
		   "tests for %d iterations\n",
		   named_parameters, iterations);
	printf("====================================================================="
		   "========================================================\n");

	for (int i = 0; i < KEM_TEST_ITERATIONS; i++)
	{

		scloud_kemkeygen(pk, sk);
		scloud_kemencaps(pk, ctx, ssa);
		scloud_kemdecaps(sk, ctx, ssb);
		if (memcmp(ssa, ssb, scloudplus_ss) != 0)
		{
			printf("\n");
			for (int i = 0; i < scloudplus_ss; i++)
			{
				printf("%d ", ssa[i]);
			}
			printf("\n");
			for (int i = 0; i < scloudplus_ss; i++)
			{
				printf("%d ", ssb[i]);
			}
			printf("wrong id is %d", i);
			printf("\n");
			return false;
		}
	}
	printf("Tests PASSED. All session keys matched.\n");

	return true;
}
static void kem_bench(const int seconds)
{
	uint8_t pk[scloudplus_pk];
	uint8_t sk[scloudplus_kem_sk];
	uint8_t ctx[scloudplus_ctx];
	uint8_t ssa[scloudplus_ss];
	uint8_t ssb[scloudplus_ss];

	TIME_OPERATION_SECONDS({ scloud_kemkeygen(pk, sk); }, "Key generation", seconds);

	scloud_kemkeygen(pk, sk);
	TIME_OPERATION_SECONDS({ scloud_kemencaps(pk, ctx, ssa); }, "KEM encapsulate", seconds);

	scloud_kemencaps(pk, ctx, ssa);
	TIME_OPERATION_SECONDS({ scloud_kemdecaps(sk, ctx, ssb); }, "KEM decapsulate", seconds);

	TIME_OPERATION_SECONDS(
		{
			scloud_kemencaps(pk, ctx, ssa);
			scloud_kemdecaps(sk, ctx, ssb);
		},
		"KEM enc and decapsulate", seconds);
}

int main()
{
	int OK = true;

	OK = kem_test(SYSTEM_NAME, KEM_TEST_ITERATIONS);
	if (OK != true)
	{
		goto exit;
	}

	PRINT_TIMER_HEADER
	kem_bench(KEM_BENCH_SECONDS);
exit:
	return (OK == true) ? EXIT_SUCCESS : EXIT_FAILURE;
}