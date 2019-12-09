"""
Detects and trim abnormally long branches using TreeShrink v1.3.2

TreeShrink must be installed and on path
Phyx (pxrr) must be installed and on path
"""

import sys
import os
import shutil
import multiprocessing as mp
from pathlib import Path
from tempfile import mkstemp


def run_treeshrink(tree_filename, quantile, out_filename):
    """run tree_shrink"""
    # Compose file names
    tree_filename = Path(tree_filename)
    cluster_id = tree_filename.stem
    tree_extension = tree_filename.suffix
    tmp_dir = Path(f"/tmp/{cluster_id}_ts")
    tmp_file = tmp_dir / f"{cluster_id}_{quantile}{tree_extension}"

    # run treeshrink
    cmd = (
        f"run_treeshrink.py "
        f"--tree {tree_filename} "
        f"--centroid "
        f"--mode per-gene "
        f"--quantiles {quantile} "
        f"--outdir {tmp_dir}"
    )
    sys.stderr.write(cmd)
    os.system(cmd)

    # Move results to out_filename
    shutil.move(tmp_file, out_filename)
    shutil.rmtree(tmp_dir)


def remove_single_quotes_from_tree(filename_in, filename_out):
    """Clean the quotes from file"""
    with open(filename_in, "r") as ts_file_in, \
        open(filename_out, "w") as ts_file_out:
        for line in ts_file_in.readlines():
            ts_file_out.write(line.replace("'", ""))


def unroot_tree(filename_in, filename_out, executable="pxrr"):
    """Unroot the tree with pxrr"""
    cmd = f"{executable} -u -t {filename_in} -o {filename_out}"
    sys.stderr.write(cmd + "\n")
    os.system(cmd)


def process_tree(filename_in, filename_out, quantile):
    """treeshrink + clean_quotes + unroot"""

    # create temporary files
    ts_file = mkstemp()
    clean_file = mkstemp()

    # Run steps
    run_treeshrink(filename_in, quantile, ts_file[1])
    remove_single_quotes_from_tree(ts_file[1], clean_file[1])
    unroot_tree(clean_file[1], filename_out)

    # Clean
    os.remove(ts_file[1])
    os.remove(clean_file[1])



if __name__ == "__main__":
    if len(sys.argv) != 6:
        sys.stderr.write(
            """
            python tree_shrink_wrapper.py in_dir in_ext quantile outdir cores
            """
        )
        sys.exit(0)

    IN_DIR = Path(sys.argv[1])
    TREE_FILE_ENDING = sys.argv[2]
    QUANTILE = sys.argv[3]
    OUT_DIR = Path(sys.argv[4])
    CORES = int(sys.argv[5])

    TREES_IN = [
        IN_DIR / x
        for x in os.listdir(IN_DIR)
        if x.endswith(TREE_FILE_ENDING)
    ]

    os.makedirs(OUT_DIR, exist_ok=True)

    TREES_OUT = [
        (OUT_DIR / x).with_suffix(".ts.tt")
        for x in os.listdir(IN_DIR)
        if x.endswith(TREE_FILE_ENDING)
    ]

    POOL = mp.Pool(CORES)
    POOL.starmap(process_tree, zip(TREES_IN, TREES_OUT, QUANTILE))
    POOL.close()
    POOL.join()
