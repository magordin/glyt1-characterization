import numpy as np
from pymol import cmd

DATASET = {
    # hGlyT1
    "8WFL": (
        "124",
        "120-131",
        "resi 179-209 or resi 290-307 or resi 454-488 or resi 491-516",
        "A",
    ),
    "8WFK": (
        "124",
        "120-131",
        "resi 179-209 or resi 290-307 or resi 454-488 or resi 491-516",
        "A",
    ),
    # hGlyT2
    "9HUE": (
        "215",
        "211-224",
        "resi 261-301 or resi 393-411 or resi 560-593 or resi 602-620",
        "A",
    ),
    "9HUF": (
        "215",
        "211-224",
        "resi 261-301 or resi 393-411 or resi 560-593 or resi 602-620",
        "A",
    ),
    "9W8J": (
        "215",
        "211-224",
        "resi 261-301 or resi 393-411 or resi 560-593 or resi 602-620",
        "A",
    ),
    "9W8Q": (
        "215",
        "211-224",
        "resi 261-301 or resi 393-411 or resi 560-593 or resi 602-620",
        "A",
    ),
    # hDAT (SLC6A3)
    "8VBY": (
        "84",
        "80-91",
        "resi 139-170 or resi 238-258 or resi 404-438 or resi 441-464",
        "A",
    ),
    "8Y2G": (
        "84",
        "80-91",
        "resi 139-170 or resi 238-258 or resi 404-438 or resi 441-464",
        "A",
    ),
    "9EO4": (
        "84",
        "80-91",
        "resi 139-170 or resi 238-258 or resi 404-438 or resi 441-464",
        "A",
    ),
    "9J6R": (
        "84",
        "80-91",
        "resi 139-170 or resi 238-258 or resi 404-438 or resi 441-464",
        "A",
    ),
    "9J6S": (
        "84",
        "80-91",
        "resi 139-170 or resi 238-258 or resi 404-438 or resi 441-464",
        "A",
    ),
    "9JKM": (
        "84",
        "80-91",
        "resi 139-170 or resi 238-258 or resi 404-438 or resi 441-464",
        "A",
    ),
    # hNET
    **{
        pdb: (
            "80",
            "76-87",
            "resi 135-166 or resi 234-254 or resi 400-434 or resi 437-460",
            "A",
        )
        for pdb in [
            "8WTW",
            "8Y8Z",
            "8Y90",
            "8Y91",
            "8Y92",
            "8Y93",
            "8Y94",
            "8YR2",
            "8Z1L",
            "8ZP1",
            "8ZP2",
            "9K6X",
            "9VWR",
        ]
    },
    # hSERT
    "7TXT": (
        "103",
        "99-110",
        "resi 158-189 or resi 257-277 or resi 423-457 or resi 460-483",
        "S",
    ),
    **{
        pdb: (
            "103",
            "99-110",
            "resi 158-189 or resi 257-277 or resi 423-457 or resi 460-483",
            "A",
        )
        for pdb in [
            "5I6X",
            "5I6Z",
            "5I71",
            "5I73",
            "5I74",
            "5I75",
            "6AWN",
            "6AWO",
            "6AWP",
            "6AWQ",
            "6DZW",
            "6DZY",
            "6VRH",
            "6VRK",
            "6VRL",
            "6W2B",
            "6W2C",
            "7LIA",
            "7LWD",
            "9HCO",
            "9IY7",
            "9VWS",
        ]
    },
    # hB0AT1 (SLC6A19)
    "8I92": (
        "56",
        "52-63",
        "resi 112-142 or resi 195-213 or resi 414-444 or resi 456-474",
        "D",
    ),
    "8I93": (
        "56",
        "52-63",
        "resi 112-142 or resi 195-213 or resi 414-444 or resi 456-474",
        "D",
    ),
    **{
        pdb: (
            "56",
            "52-63",
            "resi 112-142 or resi 195-213 or resi 414-444 or resi 456-474",
            "A",
        )
        for pdb in [
            "6M17",
            "6M18",
            "6M1D",
            "7DWX",
            "7V61",
            "8WBY",
            "8WBZ",
            "9KXT",
            "9KXU",
            "9KXV",
            "9KXW",
            "9KXY",
            "9KXZ",
            "9KY0",
            "9KY1",
        ]
    },
}

def run_phase_b_trp_shift(reference_pdb="8WFL", atom_name="CA"):
  """Executes Core Alignment and extracts Trp Rotamers + Trp residue displacement relative to reference_pdb."""
  if reference_pdb not in DATASET:
    print(f"Error: Reference PDB {reference_pdb} not found in dataset mapping.")
    return

  ref_trp_resi, ref_tm1b, ref_core, ref_chain = DATASET[reference_pdb]

  print("\n" + "=" * 80)
  print(
      f"PHASE B: TRP DISPLACEMENT EXTRACTION ({atom_name} Atom Relative to"
      f" {reference_pdb})"
  )
  print("=" * 80)

  # 1. Superimpose core helices (explicitly locked to protein chain)
  print("\n[Step 1] Aligning Core Scaffold Helices (cmd.super)...")
  for obj, (_, _, core_sel, chain_id) in DATASET.items():
    if obj == reference_pdb:
      continue

    tgt_sel = f"{obj} and chain {chain_id} and polymer.protein and ({core_sel}) and name CA"
    ref_sel = f"{reference_pdb} and chain {ref_chain} and polymer.protein and ({ref_core}) and name CA"

    # Fallback if chain A does not exist in target PDB
    if cmd.count_atoms(tgt_sel) == 0:
      tgt_sel = f"{obj} and polymer.protein and ({core_sel}) and name CA"

    cmd.super(tgt_sel, ref_sel)

  print("Core Alignment Complete.\n")

  # 2. Extract Reference Trp Atom Coordinate
  ref_atom_sel = f"first ({reference_pdb} and chain {ref_chain} and resi {ref_trp_resi} and name {atom_name})"
  ref_trp_coords = cmd.get_coords(ref_atom_sel)

  if ref_trp_coords is None or len(ref_trp_coords) == 0:
    # Try without chain filter
    ref_atom_sel = f"first ({reference_pdb} and resi {ref_trp_resi} and name {atom_name})"
    ref_trp_coords = cmd.get_coords(ref_atom_sel)

  ref_trp_pos = ref_trp_coords[0]

  # Header
  print(
      f"{'PDB':<8} | {'Chi1 (°)':<10} | {'Chi2 (°)':<10} | {'Trp Shift (Å)':<15}"
  )
  print("-" * 55)

  results = {}

  for obj, (trp_resi, _, _, chain_id) in DATASET.items():
    prefix = f"{obj} and (chain {chain_id} or not chain *) and polymer.protein and resi {trp_resi}"

    # Fallback if specific chain not found
    if cmd.count_atoms(prefix) == 0:
      prefix = f"{obj} and polymer.protein and resi {trp_resi}"

    # A. Dihedrals using exact single-atom selections
    try:
      n_sel = f"first ({prefix} and name N)"
      ca_sel = f"first ({prefix} and name CA)"
      cb_sel = f"first ({prefix} and name CB)"
      cg_sel = f"first ({prefix} and name CG)"
      cd1_sel = f"first ({prefix} and name CD1)"

      chi1 = cmd.get_dihedral(n_sel, ca_sel, cb_sel, cg_sel) % 360
      chi2 = cmd.get_dihedral(ca_sel, cb_sel, cg_sel, cd1_sel) % 360
    except Exception:
      chi1, chi2 = np.nan, np.nan

    # B. Trp Atom Displacement relative to reference
    tgt_atom_sel = f"first ({prefix} and name {atom_name})"
    tgt_trp_coords = cmd.get_coords(tgt_atom_sel)

    if tgt_trp_coords is not None and len(tgt_trp_coords) > 0:
      trp_shift = np.linalg.norm(tgt_trp_coords[0] - ref_trp_pos)
    else:
      trp_shift = np.nan

    print(f"{obj:<8} | {chi1:<10.1f} | {chi2:<10.1f} | {trp_shift:<15.2f}")

    results[obj] = {
        "chi1": round(float(chi1), 1) if not np.isnan(chi1) else None,
        "chi2": round(float(chi2), 1) if not np.isnan(chi2) else None,
        "trp_shift": (
            round(float(trp_shift), 2) if not np.isnan(trp_shift) else None
        ),
    }

  print("=" * 55 + "\n")
  print(f"pdbs = {list(results.keys())}")
  print(f"chi1 = {[v['chi1'] for v in results.values()]}")
  print(f"trp_shift = {[v['trp_shift'] for v in results.values()]}")

  return results


cmd.extend("run_phase_b_trp_shift", run_phase_b_trp_shift)