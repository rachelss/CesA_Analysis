# label_foreground branches.py

import re
import os

# Input tree file

input_file = "RAxML_bestTree.CesA_transx_alignment_aliview_edit.output"
output_dir = "branch_foreground_trees"
os.makedirs(output_dir, exist_ok=True)

# Read the full Newick tree

with open(input_file, "r") as f:
    tree = f.read().strip()

# Extract tip labels from the tree using regex

labels = re.findall(r'([^\(\):,]+):', tree)
labels = sorted(set(labels))

print(f"Found {len(labels)} labels.")

# For each label, insert #1 on its branch

for label in labels:

    # Find the exact label's branch (e.g., "tiplabel:0.123456")
    match = re.search(rf'({re.escape(label)}):([\d\.Ee+-]+)', tree)
    if not match:
        print(f"Could not find branch for {label}")
        continue
    original = match.group(0)         # e.g., 'TipName:0.0123'
    branch_length = match.group(2)    # e.g., '0.0123'
    modified = f"{label}:{branch_length}#1"
    new_tree = tree.replace(original, modified)

    # Write new tree to a file

    with open(os.path.join(output_dir, f"tree_{label}.nwk"), "w") as out_f:
        out_f.write(new_tree + ";\n")

print("Done creating labeled trees.")


