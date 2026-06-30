# Scloud+ Core

Shared Scloud+ KEM implementation code.

- `common/`: platform-independent KEM/PKE flow, sampling, encoding, hash
  adapters, random input, and utility code.  It does not contain matrix
  kernels, serialization backends, AES/SM3 primitive backends, SIMD code, or
  test drivers.
- `ref/`: portable reference/scalar backends used by Reference AES/SHAKE/SM3
  builds, and by optimized builds only where no SIMD serialization backend is
  selected.
- `avx2/`: x86-specific optimized modules, including `aes_aesni.c`, AVX2
  matrix kernels, SHAKE128x4 wrappers, and AVX2 SM3 helpers.  The `keccak4x/`
  subdirectory is third-party Keccak Code Package support code used by the
  SHAKE128x4 wrapper.
- `neon/`: AArch64 NEON/ARM-Crypto optimized modules, including
  `aes_armcrypto.c`, `sm3_aarch64.c`, NEON matrix kernels, and NEON
  serialization.
- `include/`: public core interfaces used across modules. Backend-local helper
  headers live inside their backend directories.

Each leaf `kem/` entry provides a checked-in `parameters.h` containing only the
constants fixed by that security level.  Shared code includes
`include/scloudplus_param_common.h`, which imports the leaf `parameters.h` via
the build entry include path and adds family/backend selection macros.  The
public matrix byte format is fixed in this branch, so there is no separate
A-format build option.  Internal headers are included directly from `include/`;
they are not generated or copied.

The core directory intentionally does not contain test drivers.  Correctness
tests and benchmark programs live under `Self_Evaluation/`.

Reference AES and SM3 builds use explicitly named portable sources,
`ref/aes_reference.c` and `ref/sm3_reference.c`.  Optimized AES and SM3
builds use backend-named sources such as `avx2/aes_aesni.c`,
`avx2/sm3_avx2.c`, `neon/aes_armcrypto.c`, and
`neon/sm3_aarch64.c`.

## Source Responsibilities

The shared core is intentionally split by responsibility, not by concrete
parameter set.

Message label/delabel and MsgEnc/MsgDec are deliberately not split by backend:
Reference, AVX2, and NEON builds all use `common/encode.c`.  Backend
directories should optimize matrix generation/products, AES/SM3/SHAKE helpers,
serialization, or row unpacking, not maintain separate codec semantics.

| Path | Responsibility |
| --- | --- |
| `common/kem.c`, `common/pke.c` | Family-independent KEM and PKE control flow. |
| `common/encode.c` | Shared BW/RBW label, delabel, message encoding, and integer BDD decoding used by both Reference and Optimized builds. |
| `common/sample.c` | Short-noise sampling shared by all builds. |
| `common/hash_aes_shake.c` | AES/SHAKE F/G/H/K adapter and Keccak/SHAKE implementation. |
| `common/hash_sm3_portable.c` | Portable SM3 F/G/H/K adapter used by Reference SM3 and AArch64 SM3 builds. |
| `common/aes.h`, `common/sm3_backend.h`, `common/modarith.h` | Backend-internal primitive and arithmetic interfaces shared by ref/AVX2/NEON sources. |
| `ref/pack_reference.c` | Scalar public key, ciphertext, and short secret-key serialization. |
| `ref/aes_reference.c` | Portable AES-128 CTR backend for Reference AES builds. |
| `ref/sm3_reference.c` | Portable ICCS/NGCC auxfunc SM3 backend for Reference SM3 builds. |
| `ref/matrix_reference.c` | Portable packed10 public-matrix generation and matrix products for all Reference families. |
| `avx2/aes_aesni.c` | AES-NI AES-128 CTR backend for Optimized AES AVX2 builds. |
| `avx2/avx2_util.h`, `avx2/fips202x4.h` | Backend-local AVX2 helper declarations. |
| `avx2/matrix_avx2.c` | AVX2 packed10 matrix backend for AES/SHAKE/SM3 Optimized builds. |
| `avx2/fips202x4.c` | Scloud+ SHAKE128x4 wrapper used by the AVX2 SHAKE matrix backend. |
| `avx2/hash_sm3_avx2.c`, `avx2/sm3_avx2.c` | Optimized SM3 hash dispatch and AVX2 pseudo-XOF helpers. |
| `avx2/keccak4x/` | Third-party Keccak Code Package files used only by `avx2/fips202x4.c`. |
| `neon/aes_armcrypto.c` | ARM Crypto AES-128 CTR backend for Optimized AES NEON builds. |
| `neon/pack_neon.c` | NEON public key, ciphertext, and short secret-key serialization. |
| `neon/matrix_neon.c`, `neon/matrix_neon_internal.h` | NEON packed10 matrix backend and shared NEON row/unpack helpers. |
| `neon/sample_neon.c` | NEON BD6/BD12 rejection-sampling helpers used by optimized AArch64 builds. |
| `neon/sm3_aarch64.c` | AArch64 SM3 backend using SM3 crypto instructions when available, with an in-file scalar fallback for systems such as macOS where those instructions are not exposed. |
