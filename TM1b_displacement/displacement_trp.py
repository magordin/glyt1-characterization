import numpy as np
from pymol import cmd

# Dataset mapping: (Trp Residue, TM1b Range, Core Alignment Selection)
DATASET = {
    # hGlyT1
    "8WFL": (
        "70",
        "66-77",
        "resi 125-155 or resi 236-253 or resi 400-434 or resi 437-462",
    ),
    "8WFK": (
        "70",
        "66-77",
        "resi 125-155 or resi 236-253 or resi 400-434 or resi 437-462",
    ),
    # hGlyT2
    "9HUE": (
        "215",
        "211-224",
        "resi 261-301 or resi 393-411 or resi 560-593 or resi 602-620",
    ),
    "9HUF": (
        "215",
        "211-224",
        "resi 261-301 or resi 393-411 or resi 560-593 or resi 602-620",
    ),
    "9W8J": (
        "215",
        "211-224",
        "resi 261-301 or resi 393-411 or resi 560-593 or resi 602-620",
    ),
    # hDAT
    "8VBY": (
        "84",
        "80-91",
        "resi 139-170 or resi 238-258 or resi 404-438 or resi 441-464",
    ),
    "9J6R": (
        "84",
        "80-91",
        "resi 139-170 or resi 238-258 or resi 404-438 or resi 441-464",
    ),
    "9J6S": (
        "84",
        "80-91",
        "resi 139-170 or resi 238-258 or resi 404-438 or resi 441-464",
    ),
    "9EO4": (
        "84",
        "80-91",
        "resi 139-170 or resi 238-258 or resi 404-438 or resi 441-464",
    ),
    "8Y2G": (
        "84",
        "80-91",
        "resi 139-170 or resi 238-258 or resi 404-438 or resi 441-464",
    ),
}


def run_phase_b_trp_shift(reference_pdb="8WFL", atom_name="CA"):
  """Executes Core Alignment and extracts Trp Rotamers + Trp residue displacement relative to reference_pdb."""
  if reference_pdb not in DATASET:
    print(f"Error: Reference PDB {reference_pdb} not found in dataset mapping.")
    return

  ref_trp_resi, ref_tm1b, ref_core = DATASET[reference_pdb]

  print("\n" + "=" * 80)
  print(
      f"PHASE B: TRP DISPLACEMENT EXTRACTION ({atom_name} Atom Relative to"
      f" {reference_pdb})"
  )
  print("=" * 80)

  # 1. Superimpose core helices to reference structure
  print("\nAligning Core Scaffold Helices")
  for obj, (_, _, core_sel) in DATASET.items():
    if obj == reference_pdb:
      continue
    cmd.super(
        f"{obj} and ({core_sel}) and name CA",
        f"{reference_pdb} and ({ref_core}) and name CA",
    )
  print("Core Alignment Complete.\n")

  # 2. Extract Reference Trp Atom Coordinate
  ref_trp_coords = cmd.get_coords(
      f"{reference_pdb} and resi {ref_trp_resi} and name {atom_name}"
  )
  if ref_trp_coords is None or len(ref_trp_coords) == 0:
    print(f"Error: Reference atom {atom_name} for {reference_pdb} not found.")
    return
  ref_trp_pos = ref_trp_coords[0]

  print(
      f"{'PDB':<8} | {'Chi1 (°)':<10} | {'Chi2 (°)':<10} | {'Trp Shift (Å)':<15}"
  )
  print("-" * 55)

  results = {}

  for obj, (trp_resi, _, _) in DATASET.items():
    # A. Dihedrals
    try:
      chi1 = (
          cmd.get_dihedral(
              f"{obj} and resi {trp_resi} and name N",
              f"{obj} and resi {trp_resi} and name CA",
              f"{obj} and resi {trp_resi} and name CB",
              f"{obj} and resi {trp_resi} and name CG",
          )
          % 360
      )

      chi2 = (
          cmd.get_dihedral(
              f"{obj} and resi {trp_resi} and name CA",
              f"{obj} and resi {trp_resi} and name CB",
              f"{obj} and resi {trp_resi} and name CG",
              f"{obj} and resi {trp_resi} and name CD1",
          )
          % 360
      )
    except Exception:
      chi1, chi2 = np.nan, np.nan

    # B. Trp Atom Displacement relative to Reference
    tgt_trp_coords = cmd.get_coords(
        f"{obj} and resi {trp_resi} and name {atom_name}"
    )
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

  print("# ----------------------------------------")
  print(f"pdbs = {list(results.keys())}")
  print(f"chi1 = {[v['chi1'] for v in results.values()]}")
  print(f"trp_shift = {[v['trp_shift'] for v in results.values()]}")
  print("# ----------------------------------------\n")

  return results


cmd.extend("run_phase_b_trp_shift", run_phase_b_trp_shift)