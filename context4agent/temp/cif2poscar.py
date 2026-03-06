#!/share/home/chem-wangyg/shaofl/local/src/installs/miniconda/envs/md_env/bin/python3
import argparse
import sys
import math
import re
from collections import OrderedDict

def parse_args():
    parser = argparse.ArgumentParser(description="Convert CIF to VASP POSCAR format.")
    parser.add_argument("-i", "--input", default="out.cif", help="Input CIF file path (default: out.cif)")
    parser.add_argument("-o", "--output", default="POSCAR", help="Output POSCAR file path (default: POSCAR)")
    return parser.parse_args()

class CifToPoscar:
    def __init__(self, input_file):
        self.input_file = input_file
        self.cell_params = {
            'a': None, 'b': None, 'c': None,
            'alpha': 90.0, 'beta': 90.0, 'gamma': 90.0
        }
        self.atoms = OrderedDict()  # {Element: [[x, y, z], ...]}
        self.atom_types_order = []  # To preserve order of appearance
        
    def read_cif(self):
        try:
            with open(self.input_file, 'r') as f:
                lines = f.readlines()
        except FileNotFoundError:
            print(f"Error: Input file '{self.input_file}' not found.")
            sys.exit(1)

        in_loop = False
        loop_headers = []
        
        # Regex for cleaning data
        # Matches numbers (int or float)
        float_re = r"[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?"
        
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            # Skip comments and empty lines
            if not line or line.startswith('#'):
                i += 1
                continue

            # Parse Cell Parameters
            if line.startswith('_cell_length_a'):
                self.cell_params['a'] = float(line.split()[1])
            elif line.startswith('_cell_length_b'):
                self.cell_params['b'] = float(line.split()[1])
            elif line.startswith('_cell_length_c'):
                self.cell_params['c'] = float(line.split()[1])
            elif line.startswith('_cell_angle_alpha'):
                self.cell_params['alpha'] = float(line.split()[1])
            elif line.startswith('_cell_angle_beta'):
                self.cell_params['beta'] = float(line.split()[1])
            elif line.startswith('_cell_angle_gamma'):
                self.cell_params['gamma'] = float(line.split()[1])
            
            # Detect Loop
            elif line.startswith('loop_'):
                in_loop = True
                loop_headers = []
                i += 1
                continue
            
            # Process Loop Block
            if in_loop:
                if line.startswith('_'):
                    # Reading headers
                    header_name = line.split()[0] # Handle cases like '_atom_site_label'
                    loop_headers.append(header_name)
                else:
                    # Reading data
                    # Check if this loop contains atom data
                    if '_atom_site_fract_x' in loop_headers:
                        self._parse_atom_line(line, loop_headers)
                    else:
                        # End of loop or different loop (e.g. symmetry loop)
                        # If the line starts with loop_ or _, we let the main loop handle it
                        if line.startswith('loop_') or line.startswith('_'):
                            in_loop = False
                            continue 
            i += 1

    def _parse_atom_line(self, line, headers):
        parts = line.split()
        if len(parts) < len(headers):
            return # Skip incomplete lines

        # Map headers to indices
        idx_map = {name: i for i, name in enumerate(headers)}
        
        # Identify columns
        # Priority: _atom_site_type_symbol > _atom_site_label
        type_col = idx_map.get('_atom_site_type_symbol')
        label_col = idx_map.get('_atom_site_label')
        
        x_col = idx_map.get('_atom_site_fract_x')
        y_col = idx_map.get('_atom_site_fract_y')
        z_col = idx_map.get('_atom_site_fract_z')
        
        if x_col is None or y_col is None or z_col is None:
            return

        # Extract Element Type
        element = "X"
        if type_col is not None:
            element = parts[type_col]
        elif label_col is not None:
            element = parts[label_col]
        
        # Clean element string (e.g. "Cu1" -> "Cu")
        # Remove digits
        element = re.sub(r'\d+', '', element)
        
        # Coordinates
        try:
            x = float(parts[x_col])
            y = float(parts[y_col])
            z = float(parts[z_col])
        except ValueError:
            return # Skip if coords aren't numbers (e.g. brackets)

        # Store
        if element not in self.atoms:
            self.atoms[element] = []
            self.atom_types_order.append(element)
        
        self.atoms[element].append((x, y, z))

    def calculate_lattice_vectors(self):
        # Convert degrees to radians
        a = self.cell_params['a']
        b = self.cell_params['b']
        c = self.cell_params['c']
        alpha = math.radians(self.cell_params['alpha'])
        beta = math.radians(self.cell_params['beta'])
        gamma = math.radians(self.cell_params['gamma'])

        # VASP Convention
        # a is along x-axis
        ax = a
        ay = 0.0
        az = 0.0

        # b is in xy-plane
        bx = b * math.cos(gamma)
        by = b * math.sin(gamma)
        bz = 0.0

        # c is generic
        cx = c * math.cos(beta)
        cy = c * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / math.sin(gamma)
        
        # cz = sqrt(c^2 - cx^2 - cy^2)
        cz_sq = c**2 - cx**2 - cy**2
        cz = math.sqrt(cz_sq) if cz_sq > 0 else 0.0
        
        return [
            (ax, ay, az),
            (bx, by, bz),
            (cx, cy, cz)
        ]

    def write_poscar(self, output_path):
        vectors = self.calculate_lattice_vectors()
        
        with open(output_path, 'w') as f:
            # Line 1: Comment
            f.write(f"Generated from {self.input_file} by cif2poscar.py\n")
            
            # Line 2: Scale factor
            f.write("   1.00000000000000\n")
            
            # Lines 3-5: Lattice Vectors
            for v in vectors:
                f.write(f"    {v[0]:19.16f} {v[1]:19.16f} {v[2]:19.16f}\n")
            
            # Line 6: Element Symbols
            f.write(" " + " ".join(self.atom_types_order) + "\n")
            
            # Line 7: Atom Counts
            counts = [str(len(self.atoms[elem])) for elem in self.atom_types_order]
            f.write("    " + "    ".join(counts) + "\n")
            
            # Line 8: Direct/Cartesian
            f.write("Direct\n")
            
            # Lines 9+: Coordinates
            for elem in self.atom_types_order:
                for coord in self.atoms[elem]:
                    f.write(f"  {coord[0]:19.16f} {coord[1]:19.16f} {coord[2]:19.16f}\n")

def main():
    args = parse_args()
    
    converter = CifToPoscar(args.input)
    converter.read_cif()
    
    # Check if we parsed anything
    if not converter.cell_params['a']:
        print("Error: Could not parse cell parameters from CIF.")
        sys.exit(1)
        
    converter.write_poscar(args.output)
    print(f"Successfully converted {args.input} to {args.output}")

if __name__ == "__main__":
    main()

