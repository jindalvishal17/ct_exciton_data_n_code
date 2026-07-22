# P3HT / O-IDTBR Charge-Transfer Exciton Model

Companion code for:

> Vishal Jindal, Michael J. Janik, and Scott T. Milner,
> *"First-Principles Modeling of Interfacial Charge Transfer Exciton States in Organic Solar Cells"*
> (Journal of Chemical Theory and Computation).

This repository implements the tight-binding model used in the paper to compute
charge-transfer (CT) exciton energies at the donor–acceptor interface between
poly(3-hexylthiophene) (P3HT) and the non-fullerene acceptor O-IDTBR.

## What this code does

Starting from a single molecular-dynamics frame of an amorphous P3HT:IDTBR
blend, the code:

1. **Imports** the atomistic configuration (`npt.gro`).
2. **Locates** aromatic ring centres and the 237 BT–thiophene close contacts
   (centre-to-centre distance < 5 Å, chain ends excluded).
3. **Builds** the smeared (Gaussian) 14×14 Coulomb interaction kernel
   (Eq. 6–7).
4. **Computes** the interchain hopping `t_ab` via a machine-learning model
   trained on ring-pair geometry (Fig. 6).
5. **Computes** the intrachain hopping matrix elements from DFT dimer energies
   using the energy-splitting method (Figs. 4–5, Table 2).
6. **Assembles** the onsite, hopping, and Coulomb parameters (Tables 1–2).
7. **Solves** for polaron states on P3HT and IDTBR (initial guesses).
8. **Builds** the 14×14 CT-exciton tight-binding Hamiltonian (Eq. 4).
9. **Minimizes** the total CT exciton energy
   `E_CT = E_K + E_C + E_P` (Eq. 12) for all 237 contacts, then decomposes it
   into one-body (`E_K`) and two-body (`E_C + E_P`) contributions for the
   correlation/outlier analysis (Figs. 8–13, Table 3).

## Repository layout

```
.
├── CT_Exciton_Master.nb        # entry-point notebook (documents + runs all modules)
├── src/
│   ├── 00_config.wl                    # paths & physical constants
│   ├── 01_import_configuration.wl      # MD frame import, PBC
│   ├── 02_ring_centers_close_contacts.wl
│   ├── 03_coulomb_interaction_matrix.wl
│   ├── 04_interchain_hopping_ML.wl
│   ├── 05_intrachain_hopping_DFT.wl
│   ├── 06_gather_parameters.wl
│   ├── 07_polaron_initial_condition.wl
│   ├── 08_hamiltonian_matrix.wl
│   └── 09_ct_exciton_energy_analysis.wl
├── data/                        # place input data here (see DATA_MANIFEST.md)
└── results/                     # generated outputs (ML input file, plots, etc.)
```

## Requirements

- **Wolfram Mathematica** (developed on 13.1; any recent version should work).
  The `.wl` files are plain-text Wolfram Language scripts and can also be run
  head-less with `wolframscript`.

## How to run

1. Place the required input data files in `data/` (see
   [DATA_MANIFEST.md](DATA_MANIFEST.md)).
2. Open `CT_Exciton_Master.nb` in Mathematica.
3. Evaluate the **Setup** cell, then evaluate the module cells **00 → 09** in
   order. Module 09 minimizes the exciton energy for all 237 contacts and may
   take several minutes.

Alternatively, run head-less:

```bash
wolframscript -file src/00_config.wl
wolframscript -file src/01_import_configuration.wl
# ... through 09
```
(When running head-less, set `rootDir` to the repository root at the top of
`00_config.wl` instead of `NotebookDirectory[]`.)

## Key outputs

| Symbol | Description | Paper figure/table |
|---|---|---|
| `ExcitonEnergies2` | CT exciton energy per contact | Fig. 9 |
| `ExcitonEnergies3` | IDTBR-localized exciton energy per contact | Fig. 8 |
| `results` | per-contact `{E_CT, E_C+E_P, E_K, E_C, E_P, E_Ke, E_Kh, ψ_e, ψ_h}` | Figs. 11–13, Table 3 |
| `keepIdx` / `RemoveIds` | main-population vs. outlier indices | Table 3 |
| `table3Main` / `table3Outliers` | mean Coulomb / polarization energies | Table 3 |
| `lm2`, `residuals2` | refined `E_CT` vs `E_1b` linear fit | Fig. 11 |

## Notes on this reorganization

- The original single large notebook (`CT-Exciton-Final-Code-v8.nb`) has been
  split into ten commented modules that map one-to-one onto the Methods and
  Results sections of the paper. Each module header cites the relevant equation,
  figure, or table.
- All file paths were converted from absolute machine-specific paths to
  **relative** paths under `data/` so the repository is portable.
- **Reconstructed code:** In the original notebook the total-energy function
  `CTexctionE` and the initial wavefunctions `psiEfull` / `psiHfull` were
  referenced but not defined. `CTexctionE` has been reconstructed exactly from
  Eq. (12) and the component decomposition in `CTSolve` (which was present);
  `psiEfull` / `psiHfull` are given localized polaron-based initial guesses.
  These reconstructions are flagged in `src/09_ct_exciton_energy_analysis.wl`.

## Data availability

The input data files (MD frame, ML predictions, DFT energies) are **not**
included in this repository. See [DATA_MANIFEST.md](DATA_MANIFEST.md) for the
complete list, expected locations, and formats. Add the data files before
attempting to reproduce the results.

## License

See [LICENSE](LICENSE).
