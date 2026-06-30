# Scloud+ Implementations

This directory contains the submitted compact Scloud+ implementation layout.

It contains the portable Reference tier and the optimized AVX2/NEON tier used
for the submitted KEM instances.

```text
Reference_Implementation/   portable Reference tier
Optimized_Implementation/   AVX2 or NEON Optimized tier
_shared/                    shared API_PKC facade and algorithm kernels
```

Each security level has a single build entry under
`{Reference_Implementation,Optimized_Implementation}/Scloudplus-*/kem/`.
Those leaf directories intentionally contain only the build files and the
level-specific `parameters.h`.  Shared code includes
`scloudplus_param_common.h`, which imports the leaf `parameters.h` through the
build entry include path.

`parameters.h` contains only the constants fixed by the security level:
dimensions, sampler parameters, message-code parameters, byte lengths, and API
sizes.  Primitive family and backend choices are supplied by CMake compile
definitions, or equivalently by a direct Makefile/embedded build that provides
the same definitions.

## Build Options

Select the primitive family with `SCLOUDPLUS_FAMILY`; it defaults to `AES` when
omitted.

- Reference tier: backend is fixed to `REF`.
- Optimized tier: `SCLOUDPLUS_BACKEND=AUTO` selects `AVX2` on x86_64/AMD64 and
  `NEON` on AArch64/ARM64.  Pass `SCLOUDPLUS_BACKEND=AVX2` or
  `SCLOUDPLUS_BACKEND=NEON` for explicit backend testing.

The public-matrix byte format is fixed to packed10 in this submission, so it is
not a build option.

Example:

```sh
cmake -S Reference_Implementation/Scloudplus-128/kem \
      -B build/ref128-aes \
      -DCMAKE_BUILD_TYPE=Release \
      -DSCLOUDPLUS_FAMILY=AES
cmake --build build/ref128-aes
./build/ref128-aes/test_scloudplus
```

Public targets exposed by each leaf CMake entry:

- `test_scloudplus`
- `kem_loop_scloudplus`
- `tamper_scloudplus`
- `kat_kem_scloudplus`
- `verify_kat_kem`
- `bench_scloudplus`
- `scloudplus_kem`

For AArch64 SM3 optimized builds, the compiler must support the ARM SM3/SM4
intrinsics used by `neon/sm3_aarch64.c`.  The optimized CMake entries pass
`-march=armv8.2-a+sm4` on non-Darwin AArch64 builds; this is the minimum target
feature set needed by GCC 10-compatible SM3/SM4 intrinsic headers.  Darwin
arm64 builds use the fallback path in `neon/sm3_aarch64.c` because Apple clang
does not provide the same intrinsic interface.
