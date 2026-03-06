#!/share/home/chem-wangyg/shaofl/local/src/installs/miniconda/envs/md_env/bin/python3
import argparse
import math
import os
from typing import List, Tuple


def build_lattice_vectors(a: float, b: float, c: float,
                          alpha_deg: float, beta_deg: float, gamma_deg: float) -> Tuple[Tuple[float, float, float],
                                                                                          Tuple[float, float, float],
                                                                                          Tuple[float, float, float]]:
    """Return lattice vectors (a_vec, b_vec, c_vec) for a triclinic cell in Å.

    Uses the conventional CIF definition of angles:
    - alpha = angle between b and c
    - beta  = angle between a and c
    - gamma = angle between a and b
    """
    alpha = math.radians(alpha_deg)
    beta = math.radians(beta_deg)
    gamma = math.radians(gamma_deg)

    a_vec = (a, 0.0, 0.0)
    b_x = b * math.cos(gamma)
    b_y = b * math.sin(gamma)
    b_vec = (b_x, b_y, 0.0)

    c_x = c * math.cos(beta)
    # Avoid division by zero for degenerate gamma; for valid cells sin(gamma) > 0
    sin_gamma = math.sin(gamma) if math.sin(gamma) != 0.0 else 1e-12
    c_y = c * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / sin_gamma
    # Ensure numerical stability for c_z
    c_z_sq = c * c - c_x * c_x - c_y * c_y
    c_z = math.sqrt(max(c_z_sq, 0.0))
    c_vec = (c_x, c_y, c_z)

    return a_vec, b_vec, c_vec


def frac_to_cart(frac: Tuple[float, float, float],
                 a_vec: Tuple[float, float, float],
                 b_vec: Tuple[float, float, float],
                 c_vec: Tuple[float, float, float]) -> Tuple[float, float, float]:
    u, v, w = frac
    # Map fractional coordinates to [0, 1) range
    u = u % 1.0
    v = v % 1.0
    w = w % 1.0
    
    x = u * a_vec[0] + v * b_vec[0] + w * c_vec[0]
    y = u * a_vec[1] + v * b_vec[1] + w * c_vec[1]
    z = u * a_vec[2] + v * b_vec[2] + w * c_vec[2]
    return x, y, z


def parse_cif(file_path: str):
    """Parse a minimal subset of CIF: cell lengths/angles and atom fractional coords.

    Returns:
        dict with keys: name (str), a, b, c, alpha, beta, gamma (floats),
        atoms: List[Tuple[symbol(str), (u, v, w)]]
    """
    with open(file_path, 'r', encoding='utf-8') as f:
        lines = [ln.strip() for ln in f if ln.strip()]

    name = None
    a = b = c = None
    alpha = beta = gamma = None

    # Extract block name
    for ln in lines:
        if ln.lower().startswith('data_'):
            name = ln.split('data_', 1)[1].strip() or None
            break

    # Extract cell parameters
    def get_value(prefix: str) -> float:
        for ln in lines:
            if ln.lower().startswith(prefix):
                # value may be followed by comments/quotes; take first token after key
                parts = ln.split()
                if len(parts) >= 2:
                    try:
                        return float(parts[1].strip("'\""))
                    except ValueError:
                        # Some CIFs place the value in the next token or with parentheses
                        for token in parts[1:]:
                            token_clean = token.strip("'\"")
                            try:
                                return float(token_clean)
                            except ValueError:
                                continue
        raise ValueError(f'Missing CIF key: {prefix}')

    a = get_value('_cell_length_a')
    b = get_value('_cell_length_b')
    c = get_value('_cell_length_c')
    alpha = get_value('_cell_angle_alpha')
    beta = get_value('_cell_angle_beta')
    gamma = get_value('_cell_angle_gamma')

    # Locate atom_site loop
    atoms: List[Tuple[str, Tuple[float, float, float]]] = []
    loop_start_idx = -1
    for i, ln in enumerate(lines):
        if ln.lower() == 'loop_':
            # Check if the subsequent header lines include atom_site fields
            header_start = i + 1
            header_fields: List[str] = []
            j = header_start
            while j < len(lines) and lines[j].startswith('_'):
                header_fields.append(lines[j])
                j += 1
            if any(h.lower().startswith('_atom_site_') for h in header_fields):
                loop_start_idx = i
                break

    if loop_start_idx == -1:
        raise ValueError('Could not find atom_site loop in CIF.')

    # Collect header fields for atom_site loop
    header_fields = []
    idx = loop_start_idx + 1
    while idx < len(lines) and lines[idx].startswith('_'):
        header_fields.append(lines[idx])
        idx += 1

    # Map required columns
    def find_col(name_fragment: str) -> int:
        for k, h in enumerate(header_fields):
            if h.lower().startswith(name_fragment):
                return k
        raise ValueError(f'Missing column in atom_site loop: {name_fragment}')

    sym_col = find_col('_atom_site_type_symbol')
    fx_col = find_col('_atom_site_fract_x')
    fy_col = find_col('_atom_site_fract_y')
    fz_col = find_col('_atom_site_fract_z')

    # Read data rows until next loop_ or data_/tag block
    while idx < len(lines):
        ln = lines[idx]
        if ln.lower() == 'loop_' or ln.startswith('_') or ln.lower().startswith('data_'):
            break
        parts = ln.split()
        if len(parts) <= max(sym_col, fx_col, fy_col, fz_col):
            # Not a valid atom row; stop at first non-conforming line to avoid consuming next sections
            break
        symbol = parts[sym_col].strip("'\"")
        try:
            u = float(parts[fx_col])
            v = float(parts[fy_col])
            w = float(parts[fz_col])
        except ValueError:
            # If numeric tokens are quoted
            u = float(parts[fx_col].strip("'\""))
            v = float(parts[fy_col].strip("'\""))
            w = float(parts[fz_col].strip("'\""))
        # Normalize symbol like "Cu1" -> "Cu" if type_symbol occasionally contains label
        # but prefer CIF-provided type_symbol which is usually clean (e.g., "Cu", "O")
        atoms.append((symbol, (u, v, w)))
        idx += 1

    return {
        'name': name,
        'a': a,
        'b': b,
        'c': c,
        'alpha': alpha,
        'beta': beta,
        'gamma': gamma,
        'atoms': atoms,
    }


def write_xyz(out_path: str, title: str,
              atoms_cart: List[Tuple[str, Tuple[float, float, float]]]) -> None:
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write(f"{len(atoms_cart)}\n")
        f.write(f"{title}\n")
        for sym, (x, y, z) in atoms_cart:
            # Match sample style: right-align one-letter symbols in 2 chars, 6 decimals
            f.write(f"{sym:>2}    {x:0.6f}    {y:0.6f}    {z:0.6f}\n")


def convert_cif_to_xyz(in_path: str, out_path: str = None) -> str:
    data = parse_cif(in_path)
    a_vec, b_vec, c_vec = build_lattice_vectors(
        data['a'], data['b'], data['c'], data['alpha'], data['beta'], data['gamma']
    )

    atoms_cart: List[Tuple[str, Tuple[float, float, float]]] = []
    for sym, (u, v, w) in data['atoms']:
        x, y, z = frac_to_cart((u, v, w), a_vec, b_vec, c_vec)
        atoms_cart.append((sym, (x, y, z)))

    title = data['name'] or os.path.splitext(os.path.basename(in_path))[0]
    if out_path is None:
        out_path = os.path.splitext(in_path)[0] + '.xyz'
    write_xyz(out_path, title, atoms_cart)
    return out_path


def main():
    parser = argparse.ArgumentParser(description='Convert CIF (fractional coords) to XYZ (Cartesian Å).')
    parser.add_argument('input', help='Path to input CIF file')
    parser.add_argument('-o', '--output', help='Path to output XYZ file (default: input with .xyz)')
    args = parser.parse_args()

    out = convert_cif_to_xyz(args.input, args.output)
    print(f'Wrote XYZ: {out}')


if __name__ == '__main__':
    main()


