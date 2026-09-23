import sys
import argparse


def get_formatted_atom_name(atom_name):
    """
    PDB atom-name format (columns 13-16).

    - 4 characters: start at column 13 and fill all 4 columns.
      Example: C1_1
    - 3 characters or fewer: start at column 14 with one leading space.
      Example: _N__
    """
    atom_name = atom_name.strip()

    if len(atom_name) >= 4:
        return f"{atom_name[:4]:<4}"
    else:
        return f" {atom_name:<3}"


def format_pdb_line(
    record,
    serial,
    name,
    alt,
    res,
    chain,
    seq,
    icode,
    x,
    y,
    z,
    occ,
    tf,
    el,
    ch,
):
    """
    Generate a fixed-width PDB-format line.
    Strict column alignment is used to avoid formatting errors.
    """

    # Atom name: columns 13-16
    fmt_name = get_formatted_atom_name(name)

    # Residue name: columns 18-20
    fmt_res = f"{res.strip():>3}"

    return (
        f"{record:<6}{serial:>5} {fmt_name}{alt}{fmt_res} "
        f"{chain}{seq:>4}{icode}   "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{occ:>6.2f}{tf:>6.2f}          "
        f"{el:>2}{ch:>2}\n"
    )


def process_pdb(args):
    print(f"[*] Input: {args.input}")
    print(f"[*] Output: {args.output}")
    print(f"[*] Target: Chain {args.chain}, Residue {args.residue}")
    print(
        f"[*] Action: Changing residue name to '{args.new_name}' "
        "and record type to 'HETATM'"
    )

    header_lines = []
    pre_target_lines = []

    # Atoms that will be merged into the modified residue
    target_residue_atoms = []

    post_target_lines = []
    other_chain_lines = []

    found_target = False
    found_ligand = False

    with open(args.input, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):

                # --------------------------------------------------
                # PDB parsing
                # --------------------------------------------------
                atom_name = line[12:16]
                alt_loc = line[16]
                res_name = line[17:21].strip()
                chain_id = line[21]

                try:
                    res_seq = int(line[22:26])
                except Exception:
                    res_seq = 0

                insert_code = line[26]

                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                try:
                    occ = float(line[54:60])
                except Exception:
                    occ = 1.0

                try:
                    tf = float(line[60:66])
                except Exception:
                    tf = 0.0

                el = line[76:78].strip() if len(line) > 76 else ""
                ch = line[78:80].strip() if len(line) > 78 else ""

                # --------------------------------------------------
                # Processing logic
                # --------------------------------------------------

                # 1. Ligand atoms
                if res_name.startswith(args.ligand):
                    found_ligand = True

                    target_residue_atoms.append(
                        {
                            "rec": "HETATM",
                            "name": atom_name,
                            "alt": alt_loc,
                            "res": args.new_name,
                            "chain": args.chain,
                            "seq": args.residue,
                            "icode": insert_code,
                            "x": x,
                            "y": y,
                            "z": z,
                            "occ": occ,
                            "tf": tf,
                            "el": el,
                            "ch": ch,
                            "is_backbone": False,
                        }
                    )

                # 2. Target protein residue
                elif chain_id == args.chain and res_seq == args.residue:
                    found_target = True

                    target_residue_atoms.append(
                        {
                            "rec": "HETATM",
                            "name": atom_name,
                            "alt": alt_loc,
                            "res": args.new_name,
                            "chain": chain_id,
                            "seq": res_seq,
                            "icode": insert_code,
                            "x": x,
                            "y": y,
                            "z": z,
                            "occ": occ,
                            "tf": tf,
                            "el": el,
                            "ch": ch,
                            "is_backbone": True,
                        }
                    )

                # 3. All other atoms
                elif chain_id == args.chain and res_seq < args.residue:
                    pre_target_lines.append(line)

                elif chain_id == args.chain and res_seq > args.residue:
                    post_target_lines.append(line)

                else:
                    other_chain_lines.append(line)

            elif line.startswith(("TER", "END", "CONECT", "MASTER")):
                continue

            else:
                header_lines.append(line)

    if not found_target:
        print(
            f"[!] Warning: Target residue "
            f"{args.chain}:{args.residue} was not found."
        )

    if not found_ligand:
        print(f"[!] Warning: Ligand '{args.ligand}' was not found.")

    # ----------------------------------------------------------
    # Writing output
    # Backbone atoms are written first, followed by ligand atoms
    # ----------------------------------------------------------

    modified_backbone = [
        atom for atom in target_residue_atoms if atom["is_backbone"]
    ]

    modified_ligand = [
        atom for atom in target_residue_atoms if not atom["is_backbone"]
    ]

    current_serial = 1

    with open(args.output, "w") as out:

        # Header
        out.writelines(header_lines)

        # 1. Atoms before the target residue
        for line in pre_target_lines:
            out.write(line[:6] + f"{current_serial:5d}" + line[11:])
            current_serial += 1

        # 2. Merged target residue
        # Backbone first
        for atom in modified_backbone:
            out.write(
                format_pdb_line(
                    atom["rec"],
                    current_serial,
                    atom["name"],
                    atom["alt"],
                    atom["res"],
                    atom["chain"],
                    atom["seq"],
                    atom["icode"],
                    atom["x"],
                    atom["y"],
                    atom["z"],
                    atom["occ"],
                    atom["tf"],
                    atom["el"],
                    atom["ch"],
                )
            )
            current_serial += 1

        # Ligand atoms second
        for atom in modified_ligand:
            out.write(
                format_pdb_line(
                    atom["rec"],
                    current_serial,
                    atom["name"],
                    atom["alt"],
                    atom["res"],
                    atom["chain"],
                    atom["seq"],
                    atom["icode"],
                    atom["x"],
                    atom["y"],
                    atom["z"],
                    atom["occ"],
                    atom["tf"],
                    atom["el"],
                    atom["ch"],
                )
            )
            current_serial += 1

        # 3. Atoms after the target residue
        for line in post_target_lines:
            out.write(line[:6] + f"{current_serial:5d}" + line[11:])
            current_serial += 1

        out.write(
            f"TER   {current_serial:5d}      "
            f"{args.new_name:>3} {args.chain}{args.residue:4d}\n"
        )
        current_serial += 1

        # 4. Other chains
        for line in other_chain_lines:
            out.write(line[:6] + f"{current_serial:5d}" + line[11:])
            current_serial += 1

        out.write("END\n")

    print("[*] Processing complete.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Merge a ligand into a target protein residue "
            "as a single HETATM residue."
        )
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input PDB file",
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output PDB file",
    )

    parser.add_argument(
        "-c",
        "--chain",
        default="A",
        help="Target chain ID (default: A)",
    )

    parser.add_argument(
        "-r",
        "--residue",
        type=int,
        required=True,
        help="Target residue number (e.g. 20)",
    )

    parser.add_argument(
        "-n",
        "--new_name",
        default="D4K",
        help="New residue name, up to 3 characters (default: D4K)",
    )

    parser.add_argument(
        "-l",
        "--ligand",
        default="LIG",
        help="Ligand residue name prefix in the input PDB (default: LIG)",
    )

    args = parser.parse_args()

    if len(args.new_name) > 3:
        print(
            "[!] Warning: New residue name is longer than 3 characters. "
            "It will be truncated to comply with the PDB format."
        )
        args.new_name = args.new_name[:3]

    process_pdb(args)
