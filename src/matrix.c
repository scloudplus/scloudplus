#include "matrix.h"
#include "param.h"
#include <string.h>
void scloudplus_add(uint16_t *in0, uint16_t *in1, int len, uint16_t *out)
{
	for (int i = 0; i < len; i++)
	{
		out[i] = (in0[i] + in1[i]) & 0xFFF;
	}
}
void scloudplus_sub(uint16_t *in0, uint16_t *in1, int len, uint16_t *out)
{
	for (int i = 0; i < len; i++)
	{
		out[i] = (in0[i] - in1[i]) & 0xFFF;
	}
}

void scloudplus_mul_cs(uint16_t *C, uint16_t *S, uint16_t *out)
{
	memset(out, 0, scloudplus_mbar * scloudplus_nbar * 2);
	for (int i = 0; i < scloudplus_mbar; i++)
	{
		for (int j = 0; j < scloudplus_nbar; j++)
		{
			for (int k = 0; k < scloudplus_n; k++)
			{
				out[i * scloudplus_nbar + j] +=
					C[i * scloudplus_n + k] * (uint16_t)S[j * scloudplus_n + k];
			}
		}
	}
}
void scloudplus_mul_add_sb_e(const uint16_t *S, const uint16_t *B,
							 const uint16_t *E, uint16_t *out)
{
	memcpy(out, E, scloudplus_mbar * scloudplus_nbar * 2);
	for (int i = 0; i < scloudplus_mbar; i++)
	{
		for (int j = 0; j < scloudplus_nbar; j++)
		{
			for (int k = 0; k < scloudplus_m; k++)
			{
				out[i * scloudplus_nbar + j] +=
					(uint16_t)S[i * scloudplus_m + k] * B[k * scloudplus_nbar + j];
			}
		}
	}
}
