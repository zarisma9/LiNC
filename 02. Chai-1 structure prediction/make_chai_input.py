

import argparse
import csv
import sys
from pathlib import Path


FASTA_NAME = "input_sequence.fa"
RESTRAINT_NAME = "input_covalent_bond.restraints"
RUN_NAME = "run.py"

RESTRAINT_HEADER = [
    "chainA",
    "res_idxA",
    "chainB",
    "res_idxB",
    "connection_type",
    "confidence",
    "min_distance_angstrom",
    "max_distance_angstrom",
    "comment",
    "restraint_id",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Generate Chai-1 input_sequence.fa, "
            "input_covalent_bond.restraints, and run.py."
        )
    )

    parser.add_argument(
        "--default-path",
        required=True,
        help="Directory where all generated files will be written.",
    )

    parser.add_argument(
        "--chain",
        action="append",
        nargs=3,
        metavar=("CHAIN_ID", "TYPE", "SEQUENCE"),
        required=True,
        help=(
            "Add one chain/entity. TYPE must be protein, dna, rna, or ligand. "
            "Repeat --chain as many times as needed."
        ),
    )

    parser.add_argument(
        "--bond",
        action="append",
        nargs=10,
        metavar=(
            "CHAIN_A",
            "RES_IDX_A",
            "CHAIN_B",
            "RES_IDX_B",
            "CONNECTION_TYPE",
            "CONFIDENCE",
            "MIN_DIST",
            "MAX_DIST",
            "COMMENT",
            "RESTRAINT_ID",
        ),
        required=True,
        help=(
            "Add one covalent-bond restraint. "
            "Repeat --bond as many times as needed."
        ),
    )

    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing generated files.",
    )

    return parser.parse_args()


def validate(args):
    allowed_types = {"protein", "dna", "rna", "ligand"}

    chains = []
    chain_ids = set()

    for chain_id, entity_type, sequence in args.chain:
        chain_id = chain_id.strip()
        entity_type = entity_type.strip().lower()
        sequence = "".join(sequence.split())

        if not chain_id:
            raise ValueError("Chain ID cannot be empty.")

        if chain_id in chain_ids:
            raise ValueError(f"Duplicate chain ID: {chain_id}")

        if entity_type not in allowed_types:
            raise ValueError(
                f"Unsupported entity type '{entity_type}' for chain {chain_id}. "
                f"Allowed types: {', '.join(sorted(allowed_types))}"
            )

        if not sequence:
            raise ValueError(f"Sequence/SMILES cannot be empty for chain {chain_id}.")

        chain_ids.add(chain_id)
        chains.append((chain_id, entity_type, sequence))

    bonds = []

    for values in args.bond:
        (
            chain_a,
            res_idx_a,
            chain_b,
            res_idx_b,
            connection_type,
            confidence,
            min_dist,
            max_dist,
            comment,
            restraint_id,
        ) = values

        chain_a = chain_a.strip()
        chain_b = chain_b.strip()

        if chain_a not in chain_ids:
            raise ValueError(
                f"Bond restraint references undefined chainA '{chain_a}'."
            )

        if chain_b not in chain_ids:
            raise ValueError(
                f"Bond restraint references undefined chainB '{chain_b}'."
            )

        try:
            confidence_f = float(confidence)
            min_dist_f = float(min_dist)
            max_dist_f = float(max_dist)
        except ValueError as exc:
            raise ValueError(
                "CONFIDENCE, MIN_DIST, and MAX_DIST must be numeric."
            ) from exc

        if not 0.0 <= confidence_f <= 1.0:
            raise ValueError("CONFIDENCE must be between 0.0 and 1.0.")

        if min_dist_f <= 0 or max_dist_f <= 0:
            raise ValueError("Distance values must be greater than 0.")

        if min_dist_f > max_dist_f:
            raise ValueError("MIN_DIST cannot be greater than MAX_DIST.")

        if not res_idx_a.strip() or not res_idx_b.strip():
            raise ValueError("RES_IDX_A and RES_IDX_B cannot be empty.")

        if not connection_type.strip():
            raise ValueError("CONNECTION_TYPE cannot be empty.")

        if not restraint_id.strip():
            raise ValueError("RESTRAINT_ID cannot be empty.")

        bonds.append(
            [
                chain_a,
                res_idx_a.strip(),
                chain_b,
                res_idx_b.strip(),
                connection_type.strip(),
                confidence.strip(),
                min_dist.strip(),
                max_dist.strip(),
                comment.strip(),
                restraint_id.strip(),
            ]
        )

    return chains, bonds


def build_fasta(chains):
    lines = []

    for chain_id, entity_type, sequence in chains:
        lines.append(f">{entity_type}|{chain_id}")
        lines.append(sequence)

    return "\n".join(lines) + "\n"


def build_run_py(default_path):
    default_path_repr = repr(str(default_path))

    template = '''#!/usr/bin/env python3

import argparse
import logging
import shutil
from pathlib import Path

import numpy as np
import torch
from chai_lab.chai1 import run_inference


logging.basicConfig(level=logging.INFO)


def parse_args():
    parser = argparse.ArgumentParser(description="Run Chai-1 inference.")
    parser.add_argument(
        "--default-path",
        default=__DEFAULT_PATH__,
        help="Directory containing input_sequence.fa and input_covalent_bond.restraints.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    default_path = Path(args.default_path).expanduser().resolve()
    fasta_path = default_path / "input_sequence.fa"
    constraint_path = default_path / "input_covalent_bond.restraints"
    output_dir = default_path / "output"

    if not fasta_path.is_file():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    if not constraint_path.is_file():
        raise FileNotFoundError(
            f"Covalent-bond restraint file not found: {constraint_path}"
        )

    if output_dir.exists():
        logging.warning("Removing old output directory: %s", output_dir)
        shutil.rmtree(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    candidates = run_inference(
        fasta_file=fasta_path,
        output_dir=output_dir,
        constraint_path=constraint_path,
        num_trunk_recycles=3,
        num_diffn_timesteps=200,
        seed=99,
        device="cuda:0",
        use_esm_embeddings=True,
        num_diffn_samples=1,
    )

    cif_paths = candidates.cif_paths
    agg_scores = [
        ranking_data.aggregate_score.item()
        for ranking_data in candidates.ranking_data
    ]

    torch.save(candidates.pae, output_dir / "pae.pt")
    torch.save(candidates.pde, output_dir / "pde.pt")
    torch.save(candidates.plddt, output_dir / "plddt.pt")

    scores = np.load(output_dir / "scores.model_idx_0.npz")

    logging.info("Generated CIF files: %s", cif_paths)
    logging.info("Aggregate scores: %s", agg_scores)
    logging.info("Score keys: %s", list(scores.keys()))


if __name__ == "__main__":
    main()
'''
    return template.replace("__DEFAULT_PATH__", default_path_repr)


def ensure_targets_available(output_dir, force):
    targets = [
        output_dir / FASTA_NAME,
        output_dir / RESTRAINT_NAME,
        output_dir / RUN_NAME,
    ]

    existing = [path for path in targets if path.exists()]

    if existing and not force:
        names = ", ".join(path.name for path in existing)
        raise FileExistsError(
            f"Refusing to overwrite existing file(s): {names}. "
            "Use --force to overwrite."
        )

    return targets


def main():
    args = parse_args()

    try:
        chains, bonds = validate(args)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)

    output_dir = Path(args.default_path).expanduser().resolve()

    try:
        fasta_path, restraint_path, run_path = ensure_targets_available(
            output_dir, args.force
        )
    except FileExistsError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(3)

    # No output files are created until all required inputs pass validation.
    output_dir.mkdir(parents=True, exist_ok=True)

    fasta_text = build_fasta(chains)
    run_text = build_run_py(output_dir)

    fasta_path.write_text(fasta_text, encoding="utf-8")

    with restraint_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(RESTRAINT_HEADER)
        writer.writerows(bonds)

    run_path.write_text(run_text, encoding="utf-8")
    run_path.chmod(0o755)

    print("Generated files:")
    print(f"  {fasta_path}")
    print(f"  {restraint_path}")
    print(f"  {run_path}")


if __name__ == "__main__":
    main()
