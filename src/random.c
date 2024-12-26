#include "random.h"
#include <stdio.h>
#include "config.h"
#ifdef USE_OPENSSL
#include <openssl/rand.h>
int randombytes(
	unsigned char *random_array,
	unsigned int nbytes)
{
	return RAND_bytes(random_array, (int)nbytes);
}
#else
#if defined(_WIN32) || defined(_WIN64)
#define OS_WINDOWS
#elif defined(__linux__)
#define OS_LINUX
#endif
#ifdef OS_WINDOWS
#include <Wincrypt.h>
#include <Windows.h>

int randombytes(unsigned char *buffer, unsigned int size)
{
	HCRYPTPROV prov = 0;

	if (!CryptAcquireContext(&prov, NULL, NULL, PROV_RSA_FULL,
							 CRYPT_VERIFYCONTEXT | CRYPT_SILENT))
	{
		fprintf(stderr, "Error during CryptAcquireContext!\n");
		return -1;
	}

	if (!CryptGenRandom(prov, (DWORD)size, buffer))
	{
		fprintf(stderr, "Error during CryptGenRandom!\n");
		CryptReleaseContext(prov, 0);
		return -1;
	}

	CryptReleaseContext(prov, 0);
	return 0;
}
#endif

#ifdef OS_LINUX

#include <fcntl.h>
#include <stdlib.h>
#include <unistd.h>

static int lock = -1;

static __inline void delay(unsigned int count)
{
	while (count--)
	{
	}
}

int randombytes(
	unsigned char *random_array,
	unsigned int nbytes)
{ // Generation of "nbytes" of random values

	int r, n = nbytes, count = 0;

	if (lock == -1)
	{
		do
		{
			lock = open("/dev/urandom", O_RDONLY);
			if (lock == -1)
			{
				delay(0xFFFFF);
			}
		} while (lock == -1);
	}

	while (n > 0)
	{
		do
		{
			r = read(lock, random_array + count, n);
			if (r == -1)
			{
				delay(0xFFFF);
			}
		} while (r == -1);
		count += r;
		n -= r;
	}

	return 0;
}
#endif
#endif