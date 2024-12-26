#include "encode.h"
#include "param.h"
#include <string.h>
#include <complex.h>
#include <math.h>
/**
 * @brief Compute the vector v from the message vector m.

 * @param m The input message vector of length µ.
 * @param v The output complex vector of length 16.
 *
 * Algorithm 2: The Labeling Method
 * Input: A message vector m ∈ {0, 1}µ
 * Input: Positive integers n, τ such that n = 2^k ≥ 4, ⌊ k/2 ⌋ ≤ τ , and µ ≤ τn
 − n/4 (k − 1)
 * Output: A lattice vector x ∈ C, where C is the lattice code defined based on
 the
 * nested lattices 2^τ · Z[i]^n/2 ⊆ BWn
 * Write m′ = (u0 , u1 , . . . , u_n/2 −1 ) such that uj ∈ {0, 1}^2τ −wH (j)
 * Compute v = (v0 , . . . , v_n/2 −1 ), where vj = f_2τ −wH (j) (uj ) for 0 ≤ j
 < n/2
 */

static inline void compute_v(const uint8_t *m, double complex v[16])
{
	uint8_t A[6] = {0};
	uint8_t B[20] = {0};
	uint8_t C[6] = {0};
#if (scloudplus_tau == 3)
	A[0] = (m[0] >> 0) & 0x07;
	A[1] = (m[0] >> 3) & 0x07;
	A[2] = ((m[0] >> 6) & 0x03) | ((m[1] << 2) & 0x04);
	A[3] = (m[1] >> 1) & 0x07;
	A[4] = (m[1] >> 4) & 0x07;
	A[5] = ((m[1] >> 7) & 0x01) | ((m[2] << 1) & 0x06);

	for (int i = 0; i < 3; ++i)
	{
		B[i] = (m[2] >> (2 + 2 * i)) & 0x03;
	}

	for (int i = 0; i < 4; ++i)
	{
		B[3 + i] = (m[3] >> (2 * i)) & 0x03;
		B[7 + i] = (m[4] >> (2 * i)) & 0x03;
		B[11 + i] = (m[5] >> (2 * i)) & 0x03;
		B[15 + i] = (m[6] >> (2 * i)) & 0x03;
	}
	B[19] = m[7] & 0x03;
	C[0] = (m[7] >> 2) & 0x01;
	C[1] = (m[7] >> 3) & 0x01;
	C[2] = (m[7] >> 4) & 0x01;
	C[3] = (m[7] >> 5) & 0x01;
	C[4] = (m[7] >> 6) & 0x01;
	C[5] = (m[7] >> 7) & 0x01;
#elif (scloudplus_tau == 4)
	A[0] = m[0] & 0x0F;
	A[1] = (m[0] >> 4) & 0x0F;
	A[2] = m[1] & 0x0F;
	A[3] = (m[1] >> 4) & 0x0F;
	A[4] = m[2] & 0x0F;
	A[5] = (m[2] >> 4) & 0x0F;

	B[0] = m[3] & 0x07;
	B[1] = (m[3] >> 3) & 0x07;
	B[2] = ((m[3] >> 6) & 0x03) | ((m[4] << 2) & 0x04);
	B[3] = (m[4] >> 1) & 0x07;
	B[4] = (m[4] >> 4) & 0x07;
	B[5] = ((m[4] >> 7) & 0x01) | ((m[5] << 1) & 0x06);
	B[6] = (m[5] >> 2) & 0x07;
	B[7] = (m[5] >> 5) & 0x07;

	B[8] = m[6] & 0x07;
	B[9] = (m[6] >> 3) & 0x07;
	B[10] = ((m[6] >> 6) & 0x03) | ((m[7] << 2) & 0x04);
	B[11] = (m[7] >> 1) & 0x07;
	B[12] = (m[7] >> 4) & 0x07;
	B[13] = ((m[7] >> 7) & 0x01) | ((m[8] << 1) & 0x06);
	B[14] = (m[8] >> 2) & 0x07;
	B[15] = (m[8] >> 5) & 0x07;

	B[16] = m[9] & 0x07;
	B[17] = (m[9] >> 3) & 0x07;
	B[18] = ((m[9] >> 6) & 0x03) | ((m[10] << 2) & 0x04);
	B[19] = (m[10] >> 1) & 0x07;

	C[0] = (m[10] >> 4) & 0x03;
	C[1] = (m[10] >> 6) & 0x03;
	C[2] = m[11] & 0x03;
	C[3] = (m[11] >> 2) & 0x03;
	C[4] = (m[11] >> 4) & 0x03;
	C[5] = (m[11] >> 6) & 0x03;

#endif
	uint8_t D[32] = {A[0], A[1], A[2], B[0], A[3], B[1], B[2], B[3],
					 A[4], B[4], B[5], B[6], B[7], B[8], B[9], C[0],
					 A[5], B[10], B[11], B[12], B[13], B[14], B[15], C[1],
					 B[16], B[17], B[18], C[2], B[19], C[3], C[4], C[5]};

	for (int i = 0; i < 16; ++i)
	{
		v[i] = D[2 * i] + D[2 * i + 1] * I;
	}
}

/**
 * @brief Compute the lattice vector w from the complex vector v.
 * @param v The input complex vector of length 16.
 * @param w The output lattice vector of length 32.
 *
 * for l from 1 to k − 1 do
 *     Write v = (w1 , w2 , . . . , w_n/2^l ), where wj ∈ Z[i]^2
 *     Update v ← (w1 , w1 + ϕw2 , w3 , w3 + ϕw4 , . . . , w_n/2^l −1 , w_n/2^l
 * −1 + ϕw_n/2^l ) end for Compute w = [v]_2^τ return w
 */
static inline void compute_w(const double complex v[16], uint16_t w[32])
{
	double complex tmp[16];
	for (int i = 0; i < 16; i++)
	{
		tmp[i] = v[i];
	}

	for (int i = 0; i < 8; i++)
	{
		tmp[2 * i + 1] = tmp[2 * i] + tmp[2 * i + 1] * (1 + I);
	}
	for (int i = 0; i < 4; i++)
	{
		tmp[4 * i + 2] = tmp[4 * i] + tmp[4 * i + 2] * (1 + I);
		tmp[4 * i + 3] = tmp[4 * i + 1] + tmp[4 * i + 3] * (1 + I);
	}
	for (int i = 0; i < 2; i++)
	{
		tmp[8 * i + 4] = tmp[8 * i] + tmp[8 * i + 4] * (1 + I);
		tmp[8 * i + 5] = tmp[8 * i + 1] + tmp[8 * i + 5] * (1 + I);
		tmp[8 * i + 6] = tmp[8 * i + 2] + tmp[8 * i + 6] * (1 + I);
		tmp[8 * i + 7] = tmp[8 * i + 3] + tmp[8 * i + 7] * (1 + I);
	}
	for (int i = 0; i < 8; i++)
	{
		tmp[8 + i] = tmp[i] + tmp[8 + i] * (1 + I);
	}

#if (scloudplus_tau == 3)
	for (int i = 0; i < 16; i++)
	{
		w[2 * i] =
			((uint16_t)((int)creal(tmp[i]) & 0x7) * (1 << (scloudplus_logq - 3))) &
			0xFFF;
		w[2 * i + 1] =
			((uint16_t)((int)cimag(tmp[i]) & 0x7) * (1 << (scloudplus_logq - 3))) &
			0xFFF;
	}
#elif (scloudplus_tau == 4)
	for (int i = 0; i < 16; i++)
	{
		w[2 * i] =
			((uint16_t)((int)creal(tmp[i]) & 0xF) * (1 << (scloudplus_logq - 4))) &
			0xFFF;
		w[2 * i + 1] =
			((uint16_t)((int)cimag(tmp[i]) & 0xF) * (1 << (scloudplus_logq - 4))) &
			0xFFF;
	}
#endif
}
/**
 * @brief Reduce the lattice vector w to the complex vector v.
 * @param inout The input/output complex vector of length 16.
 *
 * This function reduces the lattice vector w to the complex vector v by
 * applying the modulus operation based on the value of scloudplus_tau.
 */
static inline void reduce_w(double complex inout[16])
{
	double mod, sub, real, imag;
#if (scloudplus_tau == 3)
	inout[0] = (double)((int)creal(inout[0]) & 0x7) +
			   (double)((int)cimag(inout[0]) & 0x7) * I;
	inout[3] = (double)((int)creal(inout[3]) & 0x3) +
			   (double)((int)cimag(inout[3]) & 0x3) * I;
	inout[5] = (double)((int)creal(inout[5]) & 0x3) +
			   (double)((int)cimag(inout[5]) & 0x3) * I;
	inout[6] = (double)((int)creal(inout[6]) & 0x3) +
			   (double)((int)cimag(inout[6]) & 0x3) * I;
	inout[9] = (double)((int)creal(inout[9]) & 0x3) +
			   (double)((int)cimag(inout[9]) & 0x3) * I;
	inout[10] = (double)((int)creal(inout[10]) & 0x3) +
				(double)((int)cimag(inout[10]) & 0x3) * I;
	inout[12] = (double)((int)creal(inout[12]) & 0x3) +
				(double)((int)cimag(inout[12]) & 0x3) * I;
	inout[15] = (double)((int)creal(inout[15]) & 0x1) +
				(double)((int)cimag(inout[15]) & 0x1) * I;

	mod = ((int)cimag(inout[1]) & 0x3);
	sub = mod - (int)cimag(inout[1]);
	inout[1] = (double)(((int)creal(inout[1]) + (int)sub) & 0x7) + mod * I;

	mod = ((int)cimag(inout[2]) & 0x3);
	sub = mod - (int)cimag(inout[2]);
	inout[2] = (double)(((int)creal(inout[2]) + (int)sub) & 0x7) + mod * I;

	mod = ((int)cimag(inout[4]) & 0x3);
	sub = mod - (int)cimag(inout[4]);
	inout[4] = (double)(((int)creal(inout[4]) + (int)sub) & 0x7) + mod * I;

	mod = ((int)cimag(inout[8]) & 0x3);
	sub = mod - (int)cimag(inout[8]);
	inout[8] = (double)(((int)creal(inout[8]) + (int)sub) & 0x7) + mod * I;

	mod = ((int)cimag(inout[7]) & 0x1);
	sub = mod - (int)cimag(inout[7]);
	inout[7] = (double)(((int)creal(inout[7]) + (int)sub) & 0x3) + mod * I;

	mod = ((int)cimag(inout[11]) & 0x1);
	sub = mod - (int)cimag(inout[11]);
	inout[11] = (double)(((int)creal(inout[11]) + (int)sub) & 0x3) + mod * I;

	mod = ((int)cimag(inout[13]) & 0x1);
	sub = mod - (int)cimag(inout[13]);
	inout[13] = (double)(((int)creal(inout[13]) + (int)sub) & 0x3) + mod * I;

	mod = ((int)cimag(inout[14]) & 0x1);
	sub = mod - (int)cimag(inout[14]);
	inout[14] = (double)(((int)creal(inout[14]) + (int)sub) & 0x3) + mod * I;
#elif (scloudplus_tau == 4)
	inout[0] = (double)((int)creal(inout[0]) & 0xF) +
			   (double)((int)cimag(inout[0]) & 0xF) * I;
	inout[3] = (double)((int)creal(inout[3]) & 0x7) +
			   (double)((int)cimag(inout[3]) & 0x7) * I;
	inout[5] = (double)((int)creal(inout[5]) & 0x7) +
			   (double)((int)cimag(inout[5]) & 0x7) * I;
	inout[6] = (double)((int)creal(inout[6]) & 0x7) +
			   (double)((int)cimag(inout[6]) & 0x7) * I;
	inout[9] = (double)((int)creal(inout[9]) & 0x7) +
			   (double)((int)cimag(inout[9]) & 0x7) * I;
	inout[10] = (double)((int)creal(inout[10]) & 0x7) +
				(double)((int)cimag(inout[10]) & 0x7) * I;
	inout[12] = (double)((int)creal(inout[12]) & 0x7) +
				(double)((int)cimag(inout[12]) & 0x7) * I;
	inout[15] = (double)((int)creal(inout[15]) & 0x3) +
				(double)((int)cimag(inout[15]) & 0x3) * I;

	mod = ((int)cimag(inout[1]) & 0x7);
	sub = mod - (int)cimag(inout[1]);
	inout[1] = (double)(((int)creal(inout[1]) + (int)sub) & 0xF) + mod * I;

	mod = ((int)cimag(inout[2]) & 0x7);
	sub = mod - (int)cimag(inout[2]);
	inout[2] = (double)(((int)creal(inout[2]) + (int)sub) & 0xF) + mod * I;

	mod = ((int)cimag(inout[4]) & 0x7);
	sub = mod - (int)cimag(inout[4]);
	inout[4] = (double)(((int)creal(inout[4]) + (int)sub) & 0xF) + mod * I;

	mod = ((int)cimag(inout[8]) & 0x7);
	sub = mod - (int)cimag(inout[8]);
	inout[8] = (double)(((int)creal(inout[8]) + (int)sub) & 0xF) + mod * I;

	mod = ((int)cimag(inout[7]) & 0x3);
	sub = mod - (int)cimag(inout[7]);
	inout[7] = (double)(((int)creal(inout[7]) + (int)sub) & 0x7) + mod * I;

	mod = ((int)cimag(inout[11]) & 0x3);
	sub = mod - (int)cimag(inout[11]);
	inout[11] = (double)(((int)creal(inout[11]) + (int)sub) & 0x7) + mod * I;

	mod = ((int)cimag(inout[13]) & 0x3);
	sub = mod - (int)cimag(inout[13]);
	inout[13] = (double)(((int)creal(inout[13]) + (int)sub) & 0x7) + mod * I;

	mod = ((int)cimag(inout[14]) & 0x3);
	sub = mod - (int)cimag(inout[14]);
	inout[14] = (double)(((int)creal(inout[14]) + (int)sub) & 0x7) + mod * I;
#endif
}

/**
 * @brief Recover the message vector m from the complex vector v.
 * @param v The input complex vector of length 16.
 * @param m The output message vector of length µ.
 *
 * This function recovers the message vector m from the complex vector v
 * by extracting the real and imaginary parts of v and combining them
 * based on the value of scloudplus_tau.
 */
static inline void recover_m(const double complex v[16], uint8_t *m)
{
	int A[6] = {0, 1, 2, 4, 8, 16};
	int B[20] = {3, 5, 6, 7, 9, 10, 11, 12, 13, 14,
				 17, 18, 19, 20, 21, 22, 24, 25, 26, 28};
	int C[6] = {15, 23, 27, 29, 30, 31};
	uint16_t vecv[32];
	for (int i = 0; i < 16; i++)
	{
		vecv[2 * i] = round(creal(v[i]));
		vecv[2 * i + 1] = round(cimag(v[i]));
	}
#if (scloudplus_tau == 3)
	memset(m, 0, 8);
	for (int i = 5; i >= 0; i--)
	{
		m[7] = (m[7] << 1 | vecv[C[i]]);
	}
	m[7] = (m[7] << 2) | vecv[B[19]];
	m[6] = (m[6] | vecv[B[18]]) << 2;
	m[6] = (m[6] | vecv[B[17]]) << 2;
	m[6] = (m[6] | vecv[B[16]]) << 2;
	m[6] = (m[6] | vecv[B[15]]) << 0;
	m[5] = (m[5] | vecv[B[14]]) << 2;
	m[5] = (m[5] | vecv[B[13]]) << 2;
	m[5] = (m[5] | vecv[B[12]]) << 2;
	m[5] = (m[5] | vecv[B[11]]) << 0;
	m[4] = (m[4] | vecv[B[10]]) << 2;
	m[4] = (m[4] | vecv[B[9]]) << 2;
	m[4] = (m[4] | vecv[B[8]]) << 2;
	m[4] = (m[4] | vecv[B[7]]) << 0;
	m[3] = (m[3] | vecv[B[6]]) << 2;
	m[3] = (m[3] | vecv[B[5]]) << 2;
	m[3] = (m[3] | vecv[B[4]]) << 2;
	m[3] = (m[3] | vecv[B[3]]) << 0;
	m[2] = (m[2] | vecv[B[2]]) << 2;
	m[2] = (m[2] | vecv[B[1]]) << 2;
	m[2] = (m[2] | vecv[B[0]]) << 2;
	m[2] = m[2] | (vecv[A[5]] >> 1);
	m[1] = m[1] | (vecv[A[5]] << 7);
	m[1] = m[1] | (vecv[A[4]] << 4);
	m[1] = m[1] | (vecv[A[3]] << 1);
	m[1] = m[1] | (vecv[A[2]] >> 2);
	m[0] = m[0] | (vecv[A[2]] << 6);
	m[0] = m[0] | (vecv[A[1]] << 3);
	m[0] = m[0] | (vecv[A[0]] << 0);
#elif (scloudplus_tau == 4)
	memset(m, 0, 12);
	m[11] =
		(vecv[C[5]] << 6) | (vecv[C[4]] << 4) | (vecv[C[3]] << 2) | (vecv[C[2]]);
	m[10] = (vecv[C[1]] << 6) | (vecv[C[0]] << 4) | (vecv[B[19]] << 1) |
			((vecv[B[18]]) >> 2);
	m[9] = (vecv[B[18]] << 6) | (vecv[B[17]] << 3) | vecv[B[16]];
	m[8] = (vecv[B[15]] << 5) | (vecv[B[14]] << 2) | (vecv[B[13]] >> 1);
	m[7] = (vecv[B[13]] << 7) | (vecv[B[12]] << 4) | (vecv[B[11]] << 1) |
		   (vecv[B[10]] >> 2);
	m[6] = (vecv[B[10]] << 6) | (vecv[B[9]] << 3) | vecv[B[8]];
	m[5] = (vecv[B[7]] << 5) | (vecv[B[6]] << 2) | (vecv[B[5]] >> 1);
	m[4] = (vecv[B[5]] << 7) | (vecv[B[4]] << 4) | (vecv[B[3]] << 1) |
		   (vecv[B[2]] >> 2);
	m[3] = (vecv[B[2]] << 6) | (vecv[B[1]] << 3) | vecv[B[0]];
	m[2] = (vecv[A[5]] << 4) | (vecv[A[4]]);
	m[1] = (vecv[A[3]] << 4) | (vecv[A[2]]);
	m[0] = (vecv[A[1]] << 4) | (vecv[A[0]]);
#endif
}

/**
 * @brief Recover the complex vector v from the lattice vector w.
 * @param w The input lattice vector of length 16.
 * @param v The output complex vector of length 16.
 *
 * This function recovers the complex vector v from the lattice vector w
 * by applying the inverse operations of the labeling method.
 */
static inline void recover_v(const double complex w[16], double complex v[16])
{
	double complex tmp[16] = {0};
	for (int i = 0; i < 16; i++)
	{
		tmp[i] = w[i];
	}
	for (int i = 0; i < 8; i++)
	{
		tmp[8 + i] = ((tmp[8 + i] - tmp[i]) * (1 - I)) * 0.5;
	}
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			tmp[8 * i + 4 + j] =
				((tmp[8 * i + 4 + j] - tmp[8 * i + j]) * (1 - I)) * 0.5;
		}
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			tmp[4 * i + 2 + j] =
				((tmp[4 * i + 2 + j] - tmp[4 * i + j]) * (1 - I)) * 0.5;
		}
	}
	for (int i = 0; i < 8; i++)
	{
		tmp[2 * i + 1] = ((tmp[2 * i + 1] - tmp[2 * i]) * (1 - I)) * 0.5;
	}
	reduce_w(tmp);
	for (int i = 0; i < 16; i++)
	{
		v[i] = tmp[i];
	}
}
/**
 * @brief Computes the Euclidean distance between two sets of complex numbers.
 *
 * This function calculates the Euclidean distance between two sets of complex
 * numbers by summing the squared differences of their real and imaginary parts
 * and then taking the square root of the sum.
 *
 * @param set1 Pointer to the first set of complex numbers.
 * @param set2 Pointer to the second set of complex numbers.
 * @param size The number of complex numbers in each set.
 * @return The Euclidean distance between the two sets of complex numbers.
 */
double euclidean_distance(double complex *set1, double complex *set2,
						  int size)
{
	double sum = 0.0;
	for (int i = 0; i < size; i++)
	{
		double real_diff = creal(set1[i]) - creal(set2[i]);
		double imag_diff = cimag(set1[i]) - cimag(set2[i]);

		sum += real_diff * real_diff + imag_diff * imag_diff;
	}
	return sqrt(sum);
}
/**
 * @brief Perform the Babai's nearest plane algorithm for decoding using the BDD
 * Algorithm for Barnes-Wall lattices.
 *
 * This function implements the BDD Algorithm for Barnes-Wall lattices to decode
 * a complex vector t into a lattice vector y. It recursively divides the input
 * vector into smaller parts, computes intermediate results, and selects the
 * closest lattice vector based on the Euclidean distance.
 *
 * BDD Algorithm for Barnes-Wall lattices
 * Input: A target vector t ∈ C^n/2 and the lattice BWn
 * Output: A lattice vector y ∈ BWn
 * 1: if n = 2 then
 * 2:     return ⌊t⌉
 * 3: else
 * 4:     Write t = (t1, t2) such that t1, t2 ∈ C^n/4
 * 5:     Compute y1 = BDD(t1, BWn/2), y2 = BDD(t2, BWn/2)
 * 6:     Compute z1 = BDD(ϕ^−1(t2 − y1), BWn/2), z2 = BDD(ϕ^−1(t1 − y2), BWn/2)
 * 7:     Compute x = (y1, y1 + ϕz1), x′ = (y2 + ϕz2, y2)
 * 8:     if ||x − t|| < ||x′ − t|| then
 * 9:         return x
 * 10:    else
 * 11:        return x′
 * 12:    end if
 * 13: end if
 *
 * @param t The input complex vector of length n.
 * @param y The output lattice vector of length n.
 * @param n The length of the input and output vectors.
 */

void bddbwn(double complex *t, double complex *y, int n)
{
	int tlen = n >> 1;
	int halftlen = tlen >> 1;

	if (n == 2)
	{
		y[0] = round(creal(t[0])) + round(cimag(t[0])) * I;
		return;
	}
	double complex t1[halftlen], t2[halftlen];
	for (int i = 0; i < halftlen; i++)
	{
		t1[i] = t[i];
		t2[i] = t[i + halftlen];
	}

	double complex y1[halftlen], y2[halftlen];
	bddbwn(t1, y1, tlen);
	bddbwn(t2, y2, tlen);

	double complex z1[halftlen], z2[halftlen], z1in[halftlen], z2in[halftlen];
	for (int i = 0; i < halftlen; i++)
	{
		z1in[i] = (t2[i] - y1[i]) * (1 - I) * 0.5;
		z2in[i] = (t1[i] - y2[i]) * (1 - I) * 0.5;
	}
	bddbwn(z1in, z1, tlen);
	bddbwn(z2in, z2, tlen);

	for (int i = 0; i < halftlen; i++)
	{
		z1[i] = z1[i] * (1 + I);
		z2[i] = z2[i] * (1 + I);
	}
	double complex out1[tlen], out2[tlen];
	for (int i = 0; i < halftlen; i++)
	{
		out1[i] = y1[i];
		out1[halftlen + i] = y1[i] + z1[i];
		out2[i] = y2[i] + z2[i];
		out2[halftlen + i] = y2[i];
	}
	double d1 = euclidean_distance(out1, t, tlen);
	double d2 = euclidean_distance(out2, t, tlen);

	if (d1 < d2)
	{
		for (int i = 0; i < tlen; i++)
		{
			y[i] = out1[i];
		}
	}
	else
	{
		for (int i = 0; i < tlen; i++)
		{
			y[i] = out2[i];
		}
	}
}

void scloudplus_msgencode(const uint8_t *msg, uint16_t *matrixM)
{
	uint8_t *msgPtr = (uint8_t *)msg;
	uint16_t *matrixMPtr = matrixM;
	memset(matrixM, 0, scloudplus_mbar * scloudplus_nbar * sizeof(uint16_t));
	double complex v[16];
	for (int i = 0; i < scloudplus_subm; i++)
	{
		compute_v(msgPtr, v);
		compute_w(v, matrixMPtr);
		msgPtr += scloudplus_mu >> 3;
		matrixMPtr += 32;
	}
}

void scloudplus_msgdecode(const uint16_t *matrixM, uint8_t *msg)
{
	uint8_t *msgPtr = msg;

	double complex w[16] = {0 + 0 * I};
	double complex v[16] = {0 + 0 * I};
	for (int i = 0; i < scloudplus_subm; i++)
	{
		for (int j = 0; j < 16; j++)
		{
			v[j] = ((double)matrixM[32 * i + 2 * j] * (1 << scloudplus_tau) *
					0.000244140625) +
				   ((double)matrixM[32 * i + 2 * j + 1] * (1 << scloudplus_tau) *
					0.000244140625) *
					   I;
			w[j] = 0 + 0 * I;
		}
		bddbwn(v, w, 32);
		recover_v(w, v);
		recover_m(v, msgPtr);
		msgPtr += scloudplus_mu >> 3;
	}
}

uint32_t readu16tou32(const uint16_t *ptr)
{
	return ((uint32_t)ptr[0]) | ((uint32_t)ptr[1] << 16);
}

void scloudplus_packpk(uint16_t *B, uint8_t *pk)
{
	uint16_t *ptrin = B;
	uint8_t *ptrout = pk;
	uint32_t temp;

	for (int i = 0; i < scloudplus_m * scloudplus_nbar; i = i + 2)
	{
		temp = readu16tou32(ptrin);
		temp = (temp & 0xFFF) ^ ((temp >> 4) & 0xFFF000);
		*(uint32_t *)ptrout = temp;
		ptrin = ptrin + 2;
		ptrout = ptrout + 3;
	}
}

void scloudplus_unpackpk(uint8_t *pk, uint16_t *B)
{
	uint8_t *ptrin = pk;
	uint16_t *ptrout = B;
	for (int i = 0; i < scloudplus_m * scloudplus_nbar; i = i + 2)
	{
		*ptrout = *(uint16_t *)ptrin & 0xFFF;
		*(ptrout + 1) = (*(uint16_t *)(ptrin + 1) >> 4) & 0xFFF;
		ptrin = ptrin + 3;
		ptrout = ptrout + 2;
	}
}

void scloudplus_packsk(uint16_t *S, uint8_t *sk)
{
	uint16_t *ptrin = S;
	uint8_t *ptrout = sk;
	uint8_t temp;
	for (int i = 0; i < scloudplus_n * scloudplus_nbar; i = i + 4)
	{
		temp = 0;
		temp = (*ptrin & 0x03);
		temp = ((*(ptrin + 1) << 2) & 0x0C) ^ temp;
		temp = ((*(ptrin + 2) << 4) & 0x30) ^ temp;
		temp = ((*(ptrin + 3) << 6) & 0xC0) ^ temp;
		*ptrout = temp;
		ptrin = ptrin + 4;
		ptrout = ptrout + 1;
	}
}

void scloudplus_unpacksk(uint8_t *sk, uint16_t *S)
{
	uint8_t *ptrin = sk;
	uint16_t *ptrout = S;
	int8_t temp;
	for (int i = 0; i < scloudplus_n * scloudplus_nbar; i = i + 4)
	{
		temp = *ptrin;
		*ptrout = (int16_t)((temp & 0x03) << 14) >> 14;
		*(ptrout + 1) = (int16_t)(((temp >> 2) & 0x03) << 14) >> 14;
		*(ptrout + 2) = (int16_t)(((temp >> 4) & 0x03) << 14) >> 14;
		*(ptrout + 3) = (int16_t)(((temp >> 6) & 0x03) << 14) >> 14;
		ptrin = ptrin + 1;
		ptrout = ptrout + 4;
	}
}
void scloudplus_compressc1(uint16_t *C, uint16_t *out)
{
#if (scloudplus_l == 128)
	for (int i = 0; i < scloudplus_mbar * scloudplus_n; i++)
	{
		out[i] = ((((uint32_t)(C[i] & 0xFFF) << 9) + 2048) >> 12) & 0x1FF;
	}
#elif (scloudplus_l == 192)
	memcpy(out, C, scloudplus_mbar * scloudplus_n * sizeof(uint16_t));
#elif (scloudplus_l == 256)
	for (int i = 0; i < scloudplus_mbar * scloudplus_n; i++)
	{
		out[i] = ((((uint32_t)(C[i] & 0xFFF) << 10) + 2048) >> 12) & 0x3FF;
	}
#endif
}
void scloudplus_decompressc1(uint16_t *in, uint16_t *C)
{
#if (scloudplus_l == 128)
	for (int i = 0; i < scloudplus_mbar * scloudplus_n; i++)
	{
		C[i] = ((uint32_t)((in[i] & 0x1FF) << 12) + 256) >> 9;
	}
#elif (scloudplus_l == 192)
	memcpy(C, in, scloudplus_mbar * scloudplus_n * sizeof(uint16_t));
#elif (scloudplus_l == 256)
	for (int i = 0; i < scloudplus_mbar * scloudplus_n; i++)
	{
		C[i] = ((uint32_t)((in[i] & 0x3FF) << 12) + 512) >> 10;
	}
#endif
}
void scloudplus_compressc2(uint16_t *C, uint16_t *out)
{
	uint32_t temp, remainder;
#if (scloudplus_l == 128 || scloudplus_l == 256)
	for (int i = 0; i < scloudplus_mbar * scloudplus_nbar; i++)
	{
		temp = ((((uint32_t)(C[i] & 0xFFF) << 7) + 2048) >> 12);
		remainder = (((uint32_t)(C[i] & 0xFFF) << 7) + 2048) % 6144;
		out[i] = (temp - ((!remainder) && 1)) & 0x7F;
	}
#elif (scloudplus_l == 192)
	for (int i = 0; i < scloudplus_mbar * scloudplus_nbar; i++)
	{
		temp = ((((uint32_t)(C[i] & 0xFFF) << 10) + 2048) >> 12);
		remainder = (((uint32_t)(C[i] & 0xFFF) << 10) + 2048) % 6144;
		out[i] = (temp - ((!remainder) && 1)) & 0x3FF;
	}
#endif
}
void scloudplus_decompressc2(uint16_t *in, uint16_t *C)
{
#if (scloudplus_l == 128 || scloudplus_l == 256)
	for (int i = 0; i < scloudplus_mbar * scloudplus_nbar; i++)
	{
		C[i] = ((uint32_t)((in[i] & 0x7F) << 12) + 64) >> 7;
	}
#elif (scloudplus_l == 192)
	for (int i = 0; i < scloudplus_mbar * scloudplus_nbar; i++)
	{
		C[i] = ((uint32_t)((in[i] & 0x3FF) << 12) + 512) >> 10;
	}
#endif
}

void scloudplus_packc1(uint16_t *C, uint8_t *out)
{
#if (scloudplus_l == 128)
	int inlen = scloudplus_mbar * scloudplus_n;
	uint8_t *ptrin = (uint8_t *)C;
	for (int i = 0; i < inlen; i++)
	{
		out[i] = ptrin[2 * i];
	}
	for (int i = 0; i < (inlen >> 3); i++)
	{
		for (int j = 0; j < 8; j++)
		{
			out[inlen + i] = (out[inlen + i] << 1) | ptrin[16 * i + 2 * j + 1];
		}
	}
#elif (scloudplus_l == 192)
	uint16_t *ptrin = C;
	uint8_t *ptrout = out;
	uint32_t temp;
	for (int i = 0; i < scloudplus_mbar * scloudplus_n; i = i + 2)
	{
		temp = readu16tou32(ptrin);
		temp = (temp & 0xFFF) ^ ((temp >> 4) & 0xFFF000);
		*(uint32_t *)ptrout = temp;
		ptrin = ptrin + 2;
		ptrout = ptrout + 3;
	}
#elif (scloudplus_l == 256)
	int inlen = scloudplus_mbar * scloudplus_n;
	uint8_t *ptrin = (uint8_t *)C;
	for (int i = 0; i < inlen; i++)
	{
		out[i] = ptrin[2 * i];
	}
	for (int i = 0; i < (inlen >> 2); i++)
	{
		for (int j = 0; j < 4; j++)
		{
			out[inlen + i] = (out[inlen + i] << 2) | ptrin[8 * i + 2 * j + 1];
		}
	}
#endif
}
void scloudplus_unpackc1(uint8_t *in, uint16_t *C)
{
#if (scloudplus_l == 128)
	int outlen = scloudplus_mbar * scloudplus_n;
	for (int i = 0; i < outlen; i++)
	{
		C[i] = (uint16_t)in[i];
	}
	for (int i = 0; i < (outlen >> 3); i++)
	{
		C[8 * i] = C[8 * i] | (((uint16_t)in[outlen + i] << 1) & 0x100);
		C[8 * i + 1] = C[8 * i + 1] | (((uint16_t)in[outlen + i] << 2) & 0x100);
		C[8 * i + 2] = C[8 * i + 2] | (((uint16_t)in[outlen + i] << 3) & 0x100);
		C[8 * i + 3] = C[8 * i + 3] | (((uint16_t)in[outlen + i] << 4) & 0x100);
		C[8 * i + 4] = C[8 * i + 4] | (((uint16_t)in[outlen + i] << 5) & 0x100);
		C[8 * i + 5] = C[8 * i + 5] | (((uint16_t)in[outlen + i] << 6) & 0x100);
		C[8 * i + 6] = C[8 * i + 6] | (((uint16_t)in[outlen + i] << 7) & 0x100);
		C[8 * i + 7] = C[8 * i + 7] | (((uint16_t)in[outlen + i] << 8) & 0x100);
	}
#elif (scloudplus_l == 192)
	uint8_t *ptrin = in;
	uint16_t *ptrout = C;
	for (int i = 0; i < scloudplus_mbar * scloudplus_n; i = i + 2)
	{
		*ptrout = *(uint16_t *)ptrin & 0xFFF;
		*(ptrout + 1) = (*(uint16_t *)(ptrin + 1) >> 4) & 0xFFF;
		ptrin = ptrin + 3;
		ptrout = ptrout + 2;
	}
#elif (scloudplus_l == 256)
	int outlen = scloudplus_mbar * scloudplus_n;
	for (int i = 0; i < outlen; i++)
	{
		C[i] = (uint16_t)in[i];
	}
	for (int i = 0; i < (outlen >> 2); i++)
	{
		C[4 * i] = C[4 * i] | (((uint16_t)in[outlen + i] << 2) & 0x300);
		C[4 * i + 1] = C[4 * i + 1] | (((uint16_t)in[outlen + i] << 4) & 0x300);
		C[4 * i + 2] = C[4 * i + 2] | (((uint16_t)in[outlen + i] << 6) & 0x300);
		C[4 * i + 3] = C[4 * i + 3] | (((uint16_t)in[outlen + i] << 8) & 0x300);
	}
#endif
}

void scloudplus_packc2(uint16_t *C, uint8_t *out)
{
#if (scloudplus_l == 128)
	uint16_t *ptrin = C;
	uint8_t *ptrout = out;
	int inlen = scloudplus_mbar * scloudplus_nbar;
	for (int i = 0; i < inlen; i = i + 8)
	{
		*ptrout = ((*ptrin) & 0x7F) | (*(ptrin + 1) << 7);					// 7+1
		*(ptrout + 1) = ((*(ptrin + 1) >> 1) & 0x3F) | (*(ptrin + 2) << 6); // 6+2
		*(ptrout + 2) = ((*(ptrin + 2) >> 2) & 0x1F) | (*(ptrin + 3) << 5); // 5+3
		*(ptrout + 3) = ((*(ptrin + 3) >> 3) & 0x0F) | (*(ptrin + 4) << 4); // 4+4
		*(ptrout + 4) = ((*(ptrin + 4) >> 4) & 0x07) | (*(ptrin + 5) << 3); // 3+5
		*(ptrout + 5) = ((*(ptrin + 5) >> 5) & 0x03) | (*(ptrin + 6) << 2); // 2+6
		*(ptrout + 6) = ((*(ptrin + 6) >> 6) & 0x01) | (*(ptrin + 7) << 1); // 1+7
		ptrin = ptrin + 8;
		ptrout = ptrout + 7;
	}
#elif (scloudplus_l == 192)
	int inlen = scloudplus_mbar * scloudplus_nbar;
	uint8_t *ptrin = (uint8_t *)C;
	for (int i = 0; i < inlen; i++)
	{
		out[i] = ptrin[2 * i];
	}
	for (int i = 0; i < (inlen >> 2); i++)
	{
		for (int j = 0; j < 4; j++)
		{
			out[inlen + i] = (out[inlen + i] << 2) | ptrin[8 * i + 2 * j + 1];
		}
	}
#elif (scloudplus_l == 256)
	uint16_t *ptrin = C;
	uint8_t *ptrout = out;
	int inlen = scloudplus_mbar * scloudplus_nbar -
				((scloudplus_mbar * scloudplus_nbar) & 0x7);
	for (int i = 0; i < inlen; i = i + 8)
	{
		*ptrout = ((*ptrin) & 0x7F) | (*(ptrin + 1) << 7);					// 7+1
		*(ptrout + 1) = ((*(ptrin + 1) >> 1) & 0x3F) | (*(ptrin + 2) << 6); // 6+2
		*(ptrout + 2) = ((*(ptrin + 2) >> 2) & 0x1F) | (*(ptrin + 3) << 5); // 5+3
		*(ptrout + 3) = ((*(ptrin + 3) >> 3) & 0x0F) | (*(ptrin + 4) << 4); // 4+4
		*(ptrout + 4) = ((*(ptrin + 4) >> 4) & 0x07) | (*(ptrin + 5) << 3); // 3+5
		*(ptrout + 5) = ((*(ptrin + 5) >> 5) & 0x03) | (*(ptrin + 6) << 2); // 2+6
		*(ptrout + 6) = ((*(ptrin + 6) >> 6) & 0x01) | (*(ptrin + 7) << 1); // 1+7
		ptrin = ptrin + 8;
		ptrout = ptrout + 7;
	}
	*ptrout = ((*ptrin) & 0x7F) | (*(ptrin + 1) << 7);
	*(ptrout + 1) = ((*(ptrin + 1) >> 1) & 0x3F) | (*(ptrin + 2) << 6);
	*(ptrout + 2) = ((*(ptrin + 2) >> 2) & 0x1F) | (*(ptrin + 3) << 5);
	*(ptrout + 3) = ((*(ptrin + 3) >> 3) & 0x0F);
#endif
}
void scloudplus_unpackc2(uint8_t *in, uint16_t *C)
{
#if (scloudplus_l == 128)
	uint8_t *ptrin = in;
	uint16_t *ptrout = C;
	int outlen = scloudplus_mbar * scloudplus_nbar;
	for (int i = 0; i < outlen; i = i + 8)
	{
		*ptrout = *ptrin & 0x7F;
		*(ptrout + 1) = (*(uint16_t *)ptrin >> 7) & 0x7F;
		*(ptrout + 2) = (*(uint16_t *)(ptrin + 1) >> 6) & 0x7F;
		*(ptrout + 3) = (*(uint16_t *)(ptrin + 2) >> 5) & 0x7F;
		*(ptrout + 4) = (*(uint16_t *)(ptrin + 3) >> 4) & 0x7F;
		*(ptrout + 5) = (*(uint16_t *)(ptrin + 4) >> 3) & 0x7F;
		*(ptrout + 6) = (*(uint16_t *)(ptrin + 5) >> 2) & 0x7F;
		*(ptrout + 7) = (*(ptrin + 6) >> 1) & 0x7F;
		ptrin = ptrin + 7;
		ptrout = ptrout + 8;
	}
#elif (scloudplus_l == 192)
	int outlen = scloudplus_mbar * scloudplus_nbar;
	for (int i = 0; i < outlen; i++)
	{
		C[i] = (uint16_t)in[i];
	}
	for (int i = 0; i < (outlen >> 2); i++)
	{
		C[4 * i] = C[4 * i] | (((uint16_t)in[outlen + i] << 2) & 0x300);
		C[4 * i + 1] = C[4 * i + 1] | (((uint16_t)in[outlen + i] << 4) & 0x300);
		C[4 * i + 2] = C[4 * i + 2] | (((uint16_t)in[outlen + i] << 6) & 0x300);
		C[4 * i + 3] = C[4 * i + 3] | (((uint16_t)in[outlen + i] << 8) & 0x300);
	}
#elif (scloudplus_l == 256)
	uint8_t *ptrin = in;
	uint16_t *ptrout = C;
	int outlen = scloudplus_mbar * scloudplus_nbar -
				 ((scloudplus_mbar * scloudplus_nbar) & 0x7);
	for (int i = 0; i < outlen; i = i + 8)
	{
		*ptrout = *ptrin & 0x7F;
		*(ptrout + 1) = (*(uint16_t *)ptrin >> 7) & 0x7F;
		*(ptrout + 2) = (*(uint16_t *)(ptrin + 1) >> 6) & 0x7F;
		*(ptrout + 3) = (*(uint16_t *)(ptrin + 2) >> 5) & 0x7F;
		*(ptrout + 4) = (*(uint16_t *)(ptrin + 3) >> 4) & 0x7F;
		*(ptrout + 5) = (*(uint16_t *)(ptrin + 4) >> 3) & 0x7F;
		*(ptrout + 6) = (*(uint16_t *)(ptrin + 5) >> 2) & 0x7F;
		*(ptrout + 7) = (*(ptrin + 6) >> 1) & 0x7F;
		ptrin = ptrin + 7;
		ptrout = ptrout + 8;
	}
	*ptrout = *ptrin & 0x7F;
	*(ptrout + 1) = (*(uint16_t *)ptrin >> 7) & 0x7F;
	*(ptrout + 2) = (*(uint16_t *)(ptrin + 1) >> 6) & 0x7F;
	*(ptrout + 3) = (*(uint16_t *)(ptrin + 2) >> 5) & 0x7F;
#endif
}
