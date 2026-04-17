import re
import pandas as pd
import os

# ==============================
# CONFIG
# ==============================

FMAP_OUTPUT = "results/09_motif_mechanism_output/01_fmap_output/fmap_output.txt"
METADATA = "results/09_motif_mechanism_input/01_fmap_input/motif_selected_sequences.csv"
OUTPUT_FILE = "results/10_motif_mechanism_output/01_fmap_output/fmap_final_results.csv"

# ==============================
# PARSER
# ==============================

def parse_fmap_output(file):

    data = []

    peptide_id = None
    sequence = None
    membrane_energy = None
    fpen = None

    helices = {}
    current_helix_idx = None

    with open(file) as f:

        for line in f:
            line = line.strip()

            # ----------------------
            # NEW PEPTIDE
            # ----------------------
            m = re.match(r"\d+\s+\d+\.\d+\s+\w+\s+(P\d+)", line)
            if m:

                # guardar anterior
                for h in helices.values():
                    data.append(build_row(peptide_id, sequence, membrane_energy, fpen, h))

                peptide_id = m.group(1)
                sequence = None
                membrane_energy = None
                fpen = None
                helices = {}
                current_helix_idx = None
                continue

            # ----------------------
            # SEQUENCE
            # ----------------------
            if peptide_id and re.fullmatch(r"[A-Z]+", line):
                sequence = line
                continue

            # ----------------------
            # MEMBRANE ENERGY
            # ----------------------
            m = re.search(r"binding energy of peptide:\s*(-?\d+\.?\d*)", line)
            if m:
                membrane_energy = float(m.group(1))
                continue

            # ----------------------
            # FPEN
            # ----------------------
            m = re.match(r"fpen=\s*(-?\d+\.?\d*)", line)
            if m:
                fpen = float(m.group(1))
                continue

            # ----------------------
            # HELICES (solution)
            # ----------------------
            m = re.match(
                r"Solution NMR detectable helix\s+(\d+):\s+(\d+)\s*-\s*(\d+).*?(-?\d+\.?\d*)\s+and\s+(-?\d+\.?\d*)",
                line
            )
            if m:
                idx = int(m.group(1))

                helices[idx] = {
                    "start": int(m.group(2)),
                    "end": int(m.group(3)),
                    "stability_water": float(m.group(4)),
                    "stability_bound": float(m.group(5)),
                    "transfer_energy": None,
                    "tilt_angle": None,
                    "depth_thickness": None,
                    "is_false": "FALSE"
                }

                current_helix_idx = idx
                continue

            # ----------------------
            # FALSE HELIX
            # ----------------------
            if "false helix" in line.lower():
                if current_helix_idx in helices:
                    helices[current_helix_idx]["is_false"] = "TRUE"
                continue

            # ----------------------
            # TRANSFER ENERGY
            # ----------------------
            m = re.search(
                r"alpha-helix\s+(\d+).*?:\s*(-?\d+\.?\d*).*?tilt angle:\s*(-?\d+\.?\d*)\.?,\s*depth/thickness:\s*(-?\d+\.?\d*)",
                line
            )
            if m:
                idx = int(m.group(1))

                if idx in helices:
                    helices[idx]["transfer_energy"] = float(m.group(2))
                    helices[idx]["tilt_angle"] = float(m.group(3))
                    helices[idx]["depth_thickness"] = float(m.group(4))
                continue

    # último péptido
    for h in helices.values():
        data.append(build_row(peptide_id, sequence, membrane_energy, fpen, h))

    return pd.DataFrame(data)

# ==============================
# BUILD ROW
# ==============================

def build_row(peptide_id, sequence, membrane_energy, fpen, h):

    return {
        "Peptide_ID": peptide_id,
        "Sequence": sequence,
        "Membrane_Binding_Energy": membrane_energy,
        "Helix_Start": h.get("start"),
        "Helix_End": h.get("end"),
        "Stability_Water": h.get("stability_water"),
        "Stability_Bound_Coil": h.get("stability_bound"),
        "Transfer_Energy": h.get("transfer_energy"),
        "Tilt_Angle": h.get("tilt_angle"),
        "Depth_Thickness": h.get("depth_thickness"),
        "Fpen": fpen,
        "Is_False_Helix": h.get("is_false", "FALSE")
    }

# ==============================
# MAIN
# ==============================

if __name__ == "__main__":

    print("Parsing FMAP output...")

    fmap_df = parse_fmap_output(FMAP_OUTPUT)

    print("Loading metadata...")
    meta_df = pd.read_csv(METADATA)

    print("Merging...")

    final_df = fmap_df.merge(meta_df, on="Peptide_ID", how="left")

    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)

    final_df.to_csv(OUTPUT_FILE, index=False)

    print(f"✅ Saved to {OUTPUT_FILE}")
    print(f"Total helices: {len(final_df)}")