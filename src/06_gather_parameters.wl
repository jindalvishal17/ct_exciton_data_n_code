(* ::Package:: *)

(*============================================================================*)
(*  06_gather_parameters.wl                                                   *)
(*  Assemble the tight-binding parameters (onsite, hopping, Coulomb) for the  *)
(*  14-site CT-exciton Hamiltonian.                                          *)
(*                                                                            *)
(*  Paper ref:  Tables 1 & 2, and the Hamiltonian of Eq. (4).                *)
(*                                                                            *)
(*  Site ordering for the 14-site model (Fig. 1):                            *)
(*      P3HT :  1..7  = Th Th Th Th Th Th Th          (7-thiophene window)    *)
(*      IDTBR:  8..14 = Rh BT Th Ph Th BT Rh                                 *)
(*  The interchain contact t12 couples the central P3HT thiophene (site 4)   *)
(*  to the left BT (site 9, "Left" contact) or the right BT (site 13,        *)
(*  "Right" contact) of the IDTBR.                                            *)
(*============================================================================ *)

(*======================= P3HT (donor) parameters ========================*)
(* Onsite energies (eV): HOMO level = -cation energy, LUMO = anion energy.    *)
eH = Table[-eCationThio, {7}];          (* HOMO onsite (hole) ; Th = -8.88 eV *)
eL = Table[ eAnionThio,  {7}];          (* LUMO onsite (electron); Th = 1.71 eV*)

(* Intrachain hopping along the 7-thiophene window. tH[i] returns the six      *)
(* hopping terms flanking thiophene i (i-3 .. i+2), drawn from the per-bond    *)
(* P3HT Th-Th distribution (module 05).                                       *)
tH[i_] := Table[-tHoleP3HT[[k]], {k, i - 3, i + 2}];   (* hole hopping      *)
tL[i_] := Table[-tElecP3HT[[k]], {k, i - 3, i + 2}];   (* electron hopping  *)

(* Onsite Coulomb parameter Es,k (Table 1) for the 14 sites.                 *)
VcP3HT  = Table[4.717, {7}];
Vcidtbr = {3.854, 4.309, 4.717, 4.914, 4.717, 4.309, 3.854}; (* Rh BT Th Ph Th BT Rh *)
VcCT    = Join[VcP3HT, Vcidtbr];

(*======================= IDTBR (acceptor) parameters ====================*)
(* Split the per-contact IDTBR dimer hoppings into left-BT and right-BT sets. *)
tHoleRhBTleft  = tHoleRhBT[[1 ;; 300]];   tHoleRhBTright  = tHoleRhBT[[301 ;; 600]];
tHoleThBTleft  = tHoleThBT[[1 ;; 300]];   tHoleThBTright  = tHoleThBT[[301 ;; 600]];
tElecRhBTleft  = tElecRhBT[[1 ;; 300]];   tElecRhBTright  = tElecRhBT[[301 ;; 600]];
tElecThBTleft  = tElecThBT[[1 ;; 300]];   tElecThBTright  = tElecThBT[[301 ;; 600]];

(* Onsite energies for the 7 IDTBR sites (Table 1).                            *)
eHidtbr = -{eCationRh, eCationBTD, eCationThio, eCationBen,
            eCationThio, eCationBTD, eCationRh};     (* HOMO onsite (hole)   *)
eLidtbr =  {eAnionRh,  eAnionBTD,  eAnionThio,  eAnionBen,
            eAnionThio, eAnionBTD,  eAnionRh};        (* LUMO onsite (electron)*)

(* Intrachain hopping along IDTBR (6 bonds). The central Th-Ph couplings are  *)
(* the rigid constants tHoleThBen / tElecThBen; the BT couplings vary per     *)
(* contact (index i).                                                         *)
tHidtbr[i_] := -{tHoleRhBTleft[[i]], tHoleThBTleft[[i]], tHoleThBen, tHoleThBen,
                 tHoleThBTright[[i]], tHoleRhBTright[[i]]};
tLidtbr[i_] :=  {tElecRhBTleft[[i]], tElecThBTleft[[i]], tElecThBen, tElecThBen,
                 tElecThBTright[[i]], tElecRhBTright[[i]]};

(*======================= Interchain contact hopping =====================*)
(* t12[k] = t_ab for the k-th close contact (ML-predicted, module 04).         *)
t12[kthPair_] := tInterchain[[kthPair]];

(* Re-import the ML t_ab here as well so this module is self-contained if      *)
(* run out of order.                                                          *)
If[! ValueQ[tInterchain],
   tInterchain = Import[FileNameJoin[{dataDir, "predictedt600.csv"}]] // Flatten];
