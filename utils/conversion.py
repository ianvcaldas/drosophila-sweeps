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
    """Returns positions and genotypes of a ms file as numpy arrays. Genotypes are an
    array of zeros and ones, in the format positions x samples.
    """
    try:
        with open(ms_file, "rb") as f:
            line = f.readline()
            while line.strip() != b"//":
                line = f.readline()
            data = np.genfromtxt(
                (line for line in f), skip_header=2, dtype=np.uint8, delimiter=1
            )
    except StopIteration:  # catches files wth 0 segregating sites
        data = np.array([])
    if len(data.shape) == 1:  # correct shape for 0 or 1 segregating site
        data = data.reshape((-1, 1))
    # remove all the monomorphic sites:
    uniq_cols = [len(np.unique(col)) > 1 for col in data.T]
    data = data[:, uniq_cols]
    positions = np.array([])  # another provision for 0 segregating sites
    with open(ms_file, "r") as f:
        for line in f:
            if line.startswith("positions:"):
                positions = np.array([float(x) for x in line.split()[1:]])
                # Only keep polymorphic positions:
                positions = positions[uniq_cols]
    if data.shape[1] != positions.shape[0]:
        raise FeatureError(
            "Data is shape {data.shape}, but"
            "positions are shape {positions.shape}. We"
            "expect the number of columns in data and the"
            "number of elements in positions to be equal."
        )
    assert data.shape[1] == positions.shape[0], print(data.shape, positions.shape)
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
