# md_analysis.scripts.utils — Interface Exposure

## Public API

| Symbol                      | Module           | Description                                                      |
|-----------------------------|------------------|------------------------------------------------------------------|
| `IndexMap`                  | IndexMapper.py   | Frozen dataclass: bijective XYZ ↔ POSCAR index mapping           |
| `IndexMapError`             | IndexMapper.py   | Base exception for index-mapping operations                      |
| `IndexMapParseError`        | IndexMapper.py   | POSCAR comment line decode failure                               |
| `compute_index_map`         | IndexMapper.py   | Compute permutation from XYZ atom order to POSCAR element-grouped order |
| `write_poscar_with_map`     | IndexMapper.py   | Write POSCAR with reordered atoms + encoded mapping in comment line |
| `read_index_map_from_poscar`| IndexMapper.py   | Decode IndexMap from POSCAR first line                           |
| `encode_comment_line`       | IndexMapper.py   | Encode IndexMap → POSCAR comment string                          |
| `decode_comment_line`       | IndexMapper.py   | Decode POSCAR comment string → IndexMap                          |
| `remap_array`               | IndexMapper.py   | Reorder first axis of ndarray by mapping direction               |

## POSCAR Comment Line Format

```
md_analysis::v1 frame=<int> source=<urlencoded> n=<int> order=<csv> p2x=<base64>
```

- `p2x`: `poscar_to_xyz` permutation, base64-encoded (uint16 for n≤65535, uint32 otherwise)
- Inverse permutation (`xyz_to_poscar`) is reconstructed at decode time in O(N)

## IndexMap Dataclass Fields

| Field            | Type               | Description                                        |
|------------------|--------------------|----------------------------------------------------|
| `xyz_to_poscar`  | `np.ndarray`       | shape (N,), XYZ index → POSCAR index               |
| `poscar_to_xyz`  | `np.ndarray`       | shape (N,), POSCAR index → XYZ index               |
| `frame`          | `int`              | 0-based trajectory frame number                    |
| `source`         | `str`              | Source XYZ file path (metadata)                    |
| `n_atoms`        | `int`              | Total atom count                                   |
| `element_order`  | `tuple[str, ...]`  | POSCAR element grouping order                      |

## Stability

- **Experimental** — API may evolve as downstream scripts are implemented.
- Not re-exported from `md_analysis` top-level `__init__.py`.
