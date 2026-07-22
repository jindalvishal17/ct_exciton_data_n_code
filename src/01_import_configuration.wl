(* ::Package:: *)

(*============================================================================*)
(*  01_import_configuration.wl                                                *)
(*  Import a single MD frame of the amorphous P3HT:O-IDTBR blend.            *)
(*                                                                            *)
(*  Paper ref:  Section 2.3 "Donor-acceptor interface".                       *)
(*  The simulated melt contains 128 P3HT chains (24 monomers each) and 300    *)
(*  O-IDTBR molecules in a periodic box (7 x 7 x 27.5 nm^3), equilibrated     *)
(*  by virtual-site coarse graining over 2 microseconds (Agarwala et al.).    *)
(*                                                                            *)
(*  Input : data/npt.gro  (GROMACS .gro frame, atom coordinates in nm)        *)
(*  Output: {pbc3D, atomP3HT, atomIDTBR}                                     *)
(*    pbc3D       = periodic box dimensions (Angstroms)                       *)
(*    atomP3HT    = {name, x, y, z} for every P3HT atom   (Angstroms)         *)
(*    atomIDTBR   = {name, x, y, z} for every IDTBR atom  (Angstroms)         *)
(*                                                                            *)
(*  NOTE: coordinates are converted nm -> Angstrom (x10) so that subsequent   *)
(*  Coulomb / distance calculations use eV and Angstrom consistently.         *)
(*============================================================================*)

Clear[importFrame];

(* importFrame[fName] reads a .gro frame and returns the box, P3HT atoms and  *)
(* IDTBR atoms. The atom-index slicing below is hard-wired to the specific    *)
(* npt.gro frame layout (128 chains x 24 monomers of P3HT = 11 atoms each,    *)
(* followed by 300 IDTBR molecules of 88 atoms each).                         *)
importFrame[fName_] := Module[{frame, nAtoms, lBox3D, atoms1, atoms2, atoms,
                              newNames, idtbrAtoms, idtbrNames},
  frame = Import[fName, "Table"];

  nAtoms = frame[[2, 1]];                       (* atom count, header line 2 *)
  lBox3D = 10 * frame[[-1]];                   (* box vector (nm -> Angstrom) *)

  (* ---- P3HT : 128 chains x 24 monomers, 11 atoms per monomer ---- *)
  (* rows 3..10001 hold the first block of P3HT atoms (cols 2,4,5,6) *)
  atoms1 = frame[[3 ;; 10001, {2, 4, 5, 6}]];
  atoms2 = frame[[10002 ;; 33794, {2, 3, 4, 5}]];
  atoms  = Join[atoms1, atoms2];
  atoms[[All, 2 ;; -1]] *= 10;                  (* nm -> Angstrom           *)
  (* assign coarse-grained atom names: 11 atoms per thiophene monomer *)
  newNames = Flatten[Table[
      {"CA","CA","CA","CA","S","C2","C2","C2","C2","C2","C2"},
      128 * 24]];
  atoms[[All, 1]] = newNames;

  (* ---- IDTBR : 300 molecules, 88 atoms each ---- *)
  idtbrAtoms = frame[[33795 ;; -2, {2, 3, 4, 5}]];
  idtbrAtoms[[All, 2 ;; -1]] *= 10;             (* nm -> Angstrom           *)
  (* full 88-atom-per-molecule atom-name template for O-IDTBR *)
  idtbrNames = Flatten[Table[{
      "C","C","N","S","C","C","S","O","C","C","C","C","C","C","C","N","N","S",
      "C","C","C","C","S","C","C","C","C","C","C","C","C","C","C","C","C","C",
      "C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","S","C",
      "C","C","C","C","C","C","C","C","C","N","N","S","C","S","C","O","N","C",
      "C","C","S","C","C","C","C","C","C","C","C","C","C","C","C","C"}, 300]];
  idtbrAtoms[[All, 1]] = idtbrNames;

  {lBox3D, atoms, idtbrAtoms}
];

(*---- Load the equilibrated frame ----*)
(* Place npt.gro in the data/ directory (see DATA_MANIFEST.md).              *)
{pbc3D, atomP3HT, atomIDTBR} = importFrame[FileNameJoin[{dataDir, "npt.gro"}]];
