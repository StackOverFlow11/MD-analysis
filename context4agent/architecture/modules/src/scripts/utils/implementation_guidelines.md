# md_analysis.scripts.utils — Implementation Guidelines

## Responsibilities

- Provide bijective index mapping between CP2K XYZ atom ordering and VASP POSCAR element-grouped ordering.
- Encode/decode mapping metadata in POSCAR comment lines for downstream traceability.
- Reorder per-atom arrays between the two orderings.

## Dependencies

- **numpy**: permutation computation, array manipulation
- **ase**: `Atoms` objects, VASP POSCAR I/O (`ase.io.write`)
- **stdlib**: `base64`, `re`, `urllib.parse` for comment line encoding

No dependency on other `md_analysis` sub-packages.

## Key Design Decisions

1. **Single-direction storage**: Only `poscar_to_xyz` is stored in the POSCAR comment line; the inverse is reconstructed in O(N) at decode time. This halves the storage cost.
2. **Stable sort**: `np.argsort(kind='stable')` preserves intra-group relative order (matching VASP convention).
3. **dtype selection**: uint16 for systems ≤ 65535 atoms, uint32 otherwise, to minimize base64 payload size.
4. **URL encoding**: Source path uses `urllib.parse.quote` to handle spaces/special characters safely.
5. **Permutation validation**: On decode, the array is verified to be a valid permutation of `0..n-1`.

## Package Structure

```
scripts/
  __init__.py             # docstring only, no re-exports
  utils/
    __init__.py           # re-exports all IndexMapper public symbols
    IndexMapper.py        # core implementation
```

## Sync Rules

Changes to this module must update:
- This file (`implementation_guidelines.md`)
- `interface_exposure.md` in the same directory
- `CLAUDE.md` project map (if public API changes)
