# Scloud+ Minimal Implementation

This branch is a compact implementation package for Scloud+ KEM.  It keeps a
single public-matrix byte format and only the Reference and Optimized tiers.

## Layout

```text
Implementations/
  Reference_Implementation/Scloudplus-{128,192,256,384,512}/kem/
  Optimized_Implementation/Scloudplus-{128,192,256,384,512}/kem/
  _shared/
Self_Evaluation/
Test_Vectors/
```

All concrete builds use the same public-matrix representation; historical
algorithm names and KAT filenames keep the `packed10` tag for compatibility
with earlier submission material.

## Build

Reference example:

```sh
cmake -S Implementations/Reference_Implementation/Scloudplus-256/kem \
      -B build/ref-256-sm3 \
      -DCMAKE_BUILD_TYPE=Release \
      -DSCLOUDPLUS_FAMILY=SM3
cmake --build build/ref-256-sm3
./build/ref-256-sm3/test_scloudplus
```

Optimized x86_64 example:

```sh
cmake -S Implementations/Optimized_Implementation/Scloudplus-512/kem \
      -B build/opt-512-shake \
      -DCMAKE_BUILD_TYPE=Release \
      -DSCLOUDPLUS_FAMILY=SHAKE
cmake --build build/opt-512-shake
./build/opt-512-shake/test_scloudplus
```

Optimized AArch64 example:

```sh
cmake -S Implementations/Optimized_Implementation/Scloudplus-512/kem \
      -B build/opt-512-aes \
      -DCMAKE_BUILD_TYPE=Release \
      -DSCLOUDPLUS_FAMILY=AES
cmake --build build/opt-512-aes
./build/opt-512-aes/test_scloudplus
```

## CMake Options

- `SCLOUDPLUS_FAMILY`: `AES`, `SHAKE`, or `SM3`; defaults to `AES`.
- `SCLOUDPLUS_BACKEND`: Optimized tier only, `AUTO`, `AVX2`, or `NEON`;
  defaults to `AUTO`.

The Reference tier always uses the portable reference backend.  In the
Optimized tier, `AUTO` selects `AVX2` on x86_64/AMD64 and `NEON` on
AArch64/ARM64.  Set `SCLOUDPLUS_BACKEND` explicitly only when you need to test
one backend deliberately.  The public matrix format is fixed and is not a user
option.

## Direct Compiler Use

CMake is only a convenience wrapper.  Each `kem/` directory has a checked-in
level-specific `parameters.h`; the shared `scloudplus_param_common.h` header
derives the family/backend-dependent constants from compiler definitions.
Embedded projects or Makefiles can compile the same sources directly by passing
the family/backend macros.  For example, a portable AES-128 correctness build
can be compiled with:

```sh
cc -std=c99 -O2 \
  -DSCLOUDPLUS_FAMILY_AES -DSCLOUDPLUS_TIER_REFERENCE \
  -DSCLOUDPLUS_REF_FAMILY_AES \
  -I Implementations/Reference_Implementation/Scloudplus-128/kem \
  -I Implementations/_shared/scloudplus_core/include \
  -I Implementations/_shared/api_pkc \
  Implementations/_shared/scloudplus_core/common/aes_reference.c \
  Implementations/_shared/scloudplus_core/common/encode.c \
  Implementations/_shared/scloudplus_core/common/hash_aes_shake.c \
  Implementations/_shared/scloudplus_core/common/kem.c \
  Implementations/_shared/scloudplus_core/ref/matrix_reference.c \
  Implementations/_shared/scloudplus_core/common/pack.c \
  Implementations/_shared/scloudplus_core/common/pke.c \
  Implementations/_shared/scloudplus_core/common/random.c \
  Implementations/_shared/scloudplus_core/common/sample.c \
  Implementations/_shared/scloudplus_core/common/util.c \
  Implementations/Reference_Implementation/Scloudplus-128/kem/KEM_AlgorithmInstance.c \
  Self_Evaluation/tests/test_scloudplus.c \
  -o test_scloudplus
```

## Targets

Each `kem/` entry builds:

- `test_scloudplus`: correctness tests.
- `kem_loop_scloudplus`: repeated KEM loop smoke test.
- `tamper_scloudplus`: ciphertext tamper and implicit-rejection test.
- `kat_kem_scloudplus`: deterministic KAT generator.
- `verify_kat_kem`: deterministic KAT replay verifier.
- `bench_scloudplus`: local benchmark entry.
- `scloudplus_kem`: shared KEM library.

Canonical KAT files are stored in `Test_Vectors/`.
