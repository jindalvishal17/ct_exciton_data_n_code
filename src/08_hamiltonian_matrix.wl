(* ::Package:: *)

(*============================================================================*)
(*  08_hamiltonian_matrix.wl                                                  *)
(*  Build the 14x14 tight-binding Hamiltonian for a CT exciton at one         *)
(*  donor-acceptor contact (Eq. 4 of the paper).                             *)
(*                                                                            *)
(*  The Hamiltonian is tridiagonal for the two chains (7 P3HT sites + 7       *)
(*  IDTBR sites) plus a single off-diagonal interchain coupling t12 between    *)
(*  the central P3HT thiophene (site 4 = Ceiling[7/2]) and the contacting BT  *)
(*  ring of IDTBR. The BT can be the LEFT BT (IDTBR site 9, "Left" contact)  *)
(*  or the RIGHT BT (IDTBR site 13, "Right" contact).                        *)
(*                                                                            *)
(*  H[i,j] =   eps_i                         if i == j                        *)
(*          -  t_i                           if j == i+1  (intrachain)        *)
(*          -  t_j                           if j == i-1  (intrachain)        *)
(*          -  t12                           at the interchain contact pair    *)
(*============================================================================ *)

Clear[elemCTLeft, elemCTRight, elemCT, matCTp3ht, matCTexciton,
      GetRightorLeft];

(*----------------------- Element functions --------------------------------*)
(* elemCTLeft  : interchain link between site 4 (P3HT) and site n1+2=9 (left BT). *)
(* elemCTRight : interchain link between site 4 (P3HT) and site n1+6=13 (right BT).*)
elemCTLeft[i_, j_, e1_, e2_, t1_, t2_, t12_, n1_, n2_] := Which[
  i == j && i <= n1,                     e1[[i]],
  (j == i + 1) && i < n1,               -t1[[i]],
  (j == i - 1) && i <= n1,              -t1[[j]],
  i == j && n1 < i <= n1 + n2,           e2[[i - n1]],
  (j == i + 1) && n1 < i < n1 + n2,     -t2[[i - n1]],
  (j == i - 1) && n1 + 1 < i <= n1 + n2,-t2[[j - n1]],
  (i == Ceiling[n1/2]) && (j == n1 + 2), -t12,            (* interchain link *)
  (j == Ceiling[n1/2]) && (i == n1 + 2), -t12,
  True, 0];

elemCTRight[i_, j_, e1_, e2_, t1_, t2_, t12_, n1_, n2_] := Which[
  i == j && i <= n1,                     e1[[i]],
  (j == i + 1) && i < n1,               -t1[[i]],
  (j == i - 1) && i <= n1,              -t1[[j]],
  i == j && n1 < i <= n1 + n2,           e2[[i - n1]],
  (j == i + 1) && n1 < i < n1 + n2,     -t2[[i - n1]],
  (j == i - 1) && n1 + 1 < i <= n1 + n2,-t2[[j - n1]],
  (i == Ceiling[n1/2]) && (j == n1 + 6), -t12,           (* interchain link *)
  (j == Ceiling[n1/2]) && (i == n1 + 6), -t12,
  True, 0];

(* Standalone P3HT-only 7x7 Hamiltonian (no interchain term).                 *)
elemCT[i_, j_, e1_, t1_] := Which[
  i == j && i <= n1,        e1[[i]],
  (j == i + 1) && i < n1,  -t1[[i]],
  (j == i - 1) && i <= n1, -t1[[j]],
  True, 0];
matCTp3ht[e1_, t1_] := Array[elemCT[#1, #2, e1, t1] &, {7, 7}];

(*----------------------- Full 14x14 CT Hamiltonian -----------------------*)
(* GetRightorLeft[i] -> 1 if the contacting BT is the LEFT BT, 2 if RIGHT.     *)
GetRightorLeft[i_] := Module[{ThioNo, BTDno},
  {ThioNo, BTDno} = ThioBTDpairList[[i]];
  If[BTDno > 300, 2, 1]
];

(* matCTexciton[RorL, eP3HT, eIDTBR, tP3HT, tIDTBR, t12] -> 14x14 matrix.    *)
matCTexciton[RorL_, eP3HT_, eIDTBR_, tP3HT_, tIDTBR_, tCloseContact_] :=
  If[RorL == 2,
     Array[elemCTRight[#1, #2, eP3HT, eIDTBR, tP3HT, tIDTBR, tCloseContact, 7, 7] &, {14, 14}],
  If[RorL == 1,
     Array[elemCTLeft[#1, #2,  eP3HT, eIDTBR, tP3HT, tIDTBR, tCloseContact, 7, 7] &, {14, 14}]]];

(* Example: inspect one Hamiltonian (contact 4).                              *)
(* matCTexciton[1, eH, eHidtbr, tH[8], tHidtbr[8], -tInterchain[[4]]] // MatrixForm *)
