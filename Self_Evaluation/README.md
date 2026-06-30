# Self Evaluation

This directory contains lightweight local helpers used by the CMake entries
under `Implementations`.

## Local Tests

Each `kem/` entry under `Implementations` exposes the following test targets:

- `test_scloudplus`: functional correctness checks.
- `kem_loop_scloudplus`: repeated KEM loop smoke test.
- `tamper_scloudplus`: ciphertext tamper and implicit-rejection check.
- `verify_kat_kem`: deterministic replay against `Test_Vectors/`.

Example:

```sh
cmake -S Implementations/Reference_Implementation/Scloudplus-128/kem \
      -B build/ref-128-aes \
      -DCMAKE_BUILD_TYPE=Release \
      -DSCLOUDPLUS_FAMILY=AES
cmake --build build/ref-128-aes --target test_scloudplus verify_kat_kem
./build/ref-128-aes/test_scloudplus
./build/ref-128-aes/verify_kat_kem
```

## Local Benchmark Helper

`benchmark/measure_instance.sh` expects an already configured CMake build
directory, builds `bench_scloudplus`, and records local timing output.

```sh
cmake -S Implementations/Optimized_Implementation/Scloudplus-256/kem \
      -B build/opt-256-aes-avx2 \
      -DCMAKE_BUILD_TYPE=Release \
      -DSCLOUDPLUS_FAMILY=AES \
      -DSCLOUDPLUS_BACKEND=AVX2
Self_Evaluation/benchmark/measure_instance.sh build/opt-256-aes-avx2 all 0.25 3
```

These local measurements are for developer smoke checks.  Official performance
reports should be generated with the required benchmark harness on the target
platforms.
