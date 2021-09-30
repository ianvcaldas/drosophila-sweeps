from pathlib import Path

import numpy as np


def vcftools_to_ms(input_file, output_file, pos_file=None):
    """Convert a .012 file (from vcftools) into an ms file.
    Also requires a path to a .012.pos file with the right positions.
    If pos_file is None, look for a .012.pos file in the same folder as input_file.
    """
    input_file = Path(input_file)
    output_file = Path(output_file)
    if pos_file is None:
        pos_file = Path(str(input_file) + ".pos")
    else:
        pos_file = Path(pos_file)
    with open(pos_file) as f:
        positions = [line.strip().split("\t")[1] for line in f]
    with open(output_file, "w") as out:
        out.write(str(input_file) + "\n\n\n//\n")
        out.write(f"segsites: {len(positions)}\n")
        out.write(f"positions: {' '.join(positions)}\n")
        with open(input_file, "r") as f:
            for line in f:
                out.write(
                    "".join(line.strip().split("\t")[1:]).replace("2", "1") + "\n"
                )


def ms_to_numpy(ms_file):
    """Reads in a ms file and outputs the positions and the genotypes.
    Genotypes are a numpy array of 0s and 1s with shape (num_segsites, num_samples).
    """
    with open(ms_file, "r") as f:
        # Read in number of segregating sites and positions
        for line in f:
            if line.startswith("segsites"):
                num_segsites = int(line.strip().split()[1])
                if num_segsites == 0:
                    # Shape of data array for 0 segregating sites should be (0, 1)
                    return np.array([]), np.array([], ndmin=2, dtype=np.uint8).T
            elif line.startswith("positions"):
                positions = np.array([float(x) for x in line.strip().split()[1:]])
                break
            else:
                continue
        # Now read in the data
        data = np.array([list(line.strip()) for line in f], dtype=np.uint8)
    return positions, data.T


def treeseq_to_ms(trees, ms_file, normalize_positions=True, header_line="", seed=""):
    header = "\n".join([header_line, str(seed), "", "//"])
    num_segsites = trees.num_sites
    segsites_line = f"segsites: {num_segsites}"
    positions = [trees.site(mut.site).position for mut in trees.mutations()]
    if normalize_positions:
        positions = [pos / trees.sequence_length for pos in positions]
    positions.sort()
    positions_line = "positions: " + " ".join(f"{pos:.6f}" for pos in positions)
    haplotypes_block = "\n".join(hap for hap in trees.haplotypes())
    with open(ms_file, "w") as f:
        f.write("\n".join([header, segsites_line, positions_line, haplotypes_block]))
