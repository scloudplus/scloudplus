# Scloud+ Optimized Tier

This tier contains five security-level entries: `Scloudplus-128`,
`Scloudplus-192`, `Scloudplus-256`, `Scloudplus-384`, and `Scloudplus-512`.

Each level has a single `kem/` build entry and a level-specific
`kem/parameters.h`.  The public-matrix byte format is fixed and is
intentionally not exposed as a build option.

The Optimized tier supports `SCLOUDPLUS_BACKEND=AVX2` on x86_64 and
`SCLOUDPLUS_BACKEND=NEON` on AArch64.  The default
`SCLOUDPLUS_BACKEND=AUTO` selects the backend from the target architecture.
