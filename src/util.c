#include "util.h"

/**
 * @brief Compares two byte arrays in constant time.
 *
 * This function performs a constant-time comparison of two byte arrays
 * to prevent timing attacks. It compares each byte of the arrays and
 * accumulates the result using bitwise OR operations. The final result
 * is either 0 if the arrays are equal or 0xFF if they are not.
 *
 * @param a Pointer to the first byte array.
 * @param b Pointer to the second byte array.
 * @param len Length of the byte arrays to compare.
 * @return 0 if the arrays are equal, 0xFF if they are not.
 */
int8_t scloudplus_verify(const uint8_t *a, const uint8_t *b, size_t len)
{
	size_t i;
	uint8_t r = 0;

	for (i = 0; i < len; i++)
	{
		r |= a[i] ^ b[i];
	}

	r = (-(int8_t)(r >> 1) | -(int8_t)(r & 1)) >> (8 * sizeof(uint8_t) - 1);
	return (int8_t)r;
}
/**
 * @brief Conditionally moves data from two source arrays to a destination
 * array.
 *
 * This function performs a conditional move operation on two input arrays `a`
 * and `b` based on the value of the `bl` parameter. For each element in the
 * arrays, if `bl` is non-zero, the corresponding element from array `b` is
 * copied to the destination array `r`. Otherwise, the element from array `a` is
 * copied to `r`.
 *
 * @param r Pointer to the destination array where the result will be stored.
 * @param a Pointer to the first source array.
 * @param b Pointer to the second source array.
 * @param len The number of elements to process in the arrays.
 * @param bl The condition flag that determines which source array to use for
 * each element.
 */

void scloudplus_cmov(uint8_t *r, const uint8_t *a, const uint8_t *b, size_t len,
					 int8_t bl)
{
	for (size_t i = 0; i < len; i++)
	{
		r[i] = (~bl & a[i]) | (bl & b[i]);
	}
}