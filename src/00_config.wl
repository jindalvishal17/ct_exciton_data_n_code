(* ::Package:: *)

(*============================================================================*)
(*  00_config.wl                                                              *)
(*  Global configuration: paths and physical constants for the                *)
(*  P3HT / O-IDTBR charge-transfer (CT) exciton tight-binding model.          *)
(*                                                                            *)
(*  Companion to:                                                             *)
(*    Jindal, Janik & Milner,                                                  *)
(*    "First-Principles Modeling of Interfacial Charge Transfer Exciton       *)
(*     States in Organic Solar Cells" (JCTC).                                  *)
(*                                                                            *)
(*  This file defines directory locations (as RELATIVE paths so the repo is   *)
(*  portable) and the physical constants used throughout the calculation.     *)
(*  It should be evaluated first; every other module assumes these symbols.   *)
(*============================================================================*)

(*----------------------- Directory layout ----------------------------------*)
(* All paths are relative to the repository root, which is set as the        *)
(* working directory from the master notebook. Adjust rootDir only if you     *)
(* run a module standalone from a different working directory.               *)

rootDir    = NotebookDirectory[];          (* set by the master notebook      *)
dataDir    = FileNameJoin[{rootDir, "data"}];
resultsDir = FileNameJoin[{rootDir, "results"}];
SetDirectory[rootDir];

(*----------------------- Physical constants --------------------------------*)
(* Unit conversions.                                                          *)
(*   fac   = Hartree -> eV            (atomic-orbital energy units -> eV)    *)
(*   eFac  = (Hartree * Bohr) -> eV*A  (Coulomb 1/R with R in Angstroms)     *)
(*           27.2114 eV/Hartree  *  0.529177 A/Bohr                          *)
fac   = 27.2114;          (* Hartree to eV                                  *)
eFac  = 27.2 * 0.529;     (* eV*A ; used in the smeared-Coulomb integrals    *)

(* Dielectric medium.                                                        *)
(*   eps = relative permittivity eps_r of the amorphous blend (Eq. 8,10,12). *)
eps = 3.0;

(*----------------------- System definition ----------------------------------*)
(* The CT exciton lives on a 14-site chain:                                  *)
(*   sites 1..7  = seven-thiophene P3HT segment (central Th in close contact) *)
(*   sites 8..14 = the 7 aromatic moieties of an O-IDTBR molecule             *)
(*                (Rh - BT - Th - Ph - Th - BT - Rh ; see Fig. 1 of the paper)*)
nSites = 14;
n1 = 7;   (* number of P3HT sites                                            *)
n2 = 7;   (* number of IDTBR sites                                           *)

(* Number of donor-acceptor close contacts identified from the MD frame     *)
(* (BT-thiophene center-to-center distance < 5 A, chain ends excluded).       *)
nContacts = 237;

(*----------------------- Convenience ---------------------------------------*)
(* Quiet heavy messages from FindMinimum when optimizing many configurations *)
(* (uncomment the next line if you want a silent batch run).                 *)
(* Off[FindMinimum::cvmit]; *)
