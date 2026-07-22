(* ::Package:: *)

(*============================================================================*)
(*  03_coulomb_interaction_matrix.wl                                          *)
(*  Smeared (Gaussian) Coulomb interaction matrix for the 14 CT-exciton      *)
(*  sites of a single donor-acceptor contact.                                 *)
(*                                                                            *)
(*  Paper ref:  Section 2.2 "Coulomb energy term", Eqs. (5)-(7), Table 1.     *)
(*  Each monomer orbital is approximated as a 3D Gaussian charge cloud of     *)
(*  width sigma. The direct Coulomb integral between sites i and j is         *)
(*                                                                            *)
(*        D_ij = (1/R_ij) Erf[ R_ij / (2 sigma_ij) ],    sigma_ij = (sigma_i + sigma_j)/2   *)
(*                                                                            *)
(*  (Eq. 6). The onsite (i==j) limit is 1/(Sqrt[Pi] sigma_i). The full        *)
(*  Coulomb energy EC (Eq. 7) combines an onsite term Es,i with the offsite   *)
(*  Gaussian integrals; here we pre-build the 14x14 kernel matCoulIntCT and  *)
(*  the onsite Coulomb vector VcCT, which are contracted with the            *)
(*  wavefunction amplitudes in module 09 (Eq. 7 and Eq. 12).                  *)
(*                                                                            *)
(*  Monomer size parameters deltaCT (proportional to the Table 1 smearing     *)
(*  lengths sigma_k) for the 14 sites in order                                *)
(*      [Th Th Th Th Th Th Th | Rh BT Th Ph Th BT Rh].                       *)
(*============================================================================ *)

Clear[deltaP3HT, deltaIDTBR, deltaCT, matCoulIntCT, IDTBRringsCenter,
      getCloseContactPosn];

(* Monomer size parameters (Angstrom). deltaP3HT = thiophene value for all 7  *)
(* P3HT sites; deltaIDTBR lists the 7 IDTBR moieties (Rh BT Th Ph Th BT Rh).  *)
deltaP3HT  = Table[4., {7}];
deltaIDTBR = {6., 4.4, 4., 4.15, 4., 4.4, 6.};
deltaCT    = Join[deltaP3HT, deltaIDTBR];

(*----------------------- Contact geometry ----------------------------------*)
(* IDTBRringsCenter[n] -> 7 centres (Rh BT Th Ph Th BT Rh) of molecule n.     *)
IDTBRringsCenter[n_] := {
   ringCenter[rh1Atoms[n]], ringCenter[btd1Atoms[n]], ringCenter[idt1Atoms[n]],
   ringCenter[idbAtoms[n]],  ringCenter[idt2Atoms[n]], ringCenter[btd2Atoms[n]],
   ringCenter[rh2Atoms[n]]};

(* getCloseContactPosn[n] -> 14 site centres for contact n:                    *)
(*   sites 1..7 = seven thiophene centres of the P3HT window (contact Th +/- 3) *)
(*   sites 8..14 = the 7 IDTBR aromatic centres.                              *)
getCloseContactPosn[nThpair_] := Module[{a, b, seg1, seg2},
  {a, b} = GetCloseContactRingPair[nThpair];
  seg1 = Table[ringCenter[p3htAtoms[i]], {i, a - 3, a + 3}];
  seg2 = IDTBRringsCenter[b];
  Join[seg1, seg2]
];

(*----------------------- Smeared Coulomb matrix (14x14) --------------------*)
(* matCoulIntCT[posns] builds the full D_ij kernel (Eq. 6) for the 14 site    *)
(* positions of one contact, including the onsite (i==j) Gaussian limit.     *)
(* Distances use the minimum image (dpbc). Units: eV*A from eFac.            *)
matCoulIntCT[posns_] := Module[{mat, i, j, R},
  mat = eFac * Outer[0 &, posns, posns, 1];
  Do[
    If[i != j,
       R = Norm[dpbc[posns[[i]] - posns[[j]]]];
       mat[[i, j]] = (eFac / R) *
         Erf[R / ((deltaCT[[i]] + deltaCT[[j]]) / 2)]],
    {i, 14}, {j, 14}];
  Do[mat[[i, i]] = eFac / (Sqrt[Pi] (deltaCT[[i]] / 2)), {i, 14}];
  mat
];

(*----------------------- Onsite Coulomb parameter vector --------------------*)
(* VcCT = Es,k (Table 1, "Coulomb interaction parameter") for the 14 sites,   *)
(* in order [Th x7 | Rh BT Th Ph Th BT Rh]. These multiply the onsite         *)
(* |a_h|^2 |a_e|^2 products in Eq. 7.                                         *)
VcP3HT  = Table[4.717, {7}];                       (* thiophene Es = 4.72 eV  *)
Vcidtbr = {3.854, 4.309, 4.717, 4.914, 4.717, 4.309, 3.854};  (* Rh BT Th Ph Th BT Rh *)
VcCT    = Join[VcP3HT, Vcidtbr];

(* Sanity check: typical spacing of adjacent IDTBR rings (A). *)
(* Table[Norm[IDTBRringsCenter[130][[i]] - IDTBRringsCenter[130][[i+1]]], {i, 6}] *)
