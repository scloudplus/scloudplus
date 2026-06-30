# Scloud+ Reference Tier

This tier contains five security-level entries: `Scloudplus-128`,
`Scloudplus-192`, `Scloudplus-256`, `Scloudplus-384`, and `Scloudplus-512`.

Each level has a single `kem/` build entry and a level-specific
`kem/parameters.h`.  The public-matrix byte format is fixed and is
intentionally not exposed as a build option.

The Reference tier uses the portable REF backend.
