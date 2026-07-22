(* ::Package:: *)

(*============================================================================*)
(*  02_ring_centers_close_contacts.wl                                         *)
(*  Locate aromatic ring centers and find P3HT-IDTBR close contacts.         *)
(*                                                                            *)
(*  Paper ref:  Section 2.3 "Donor-acceptor interface" and Fig. 1.           *)
(*  The most probable contact is between the benzothiadiazole (BT) ring of   *)
(*  O-IDTBR and a thiophene (Th) ring of P3HT (steric shielding by the octyl  *)
(*  side chains on the IDT core). We keep BT-Th pairs with center-to-center  *)
(*  distance < 5 A, then exclude thiophenes within 3 monomers of a chain     *)
(*  end, leaving 237 donor-acceptor contacts (ThioBTDpairList).               *)
(*                                                                            *)
(*  Each P3HT segment used later is a 7-thiophene window centred on the       *)
(*  contact thiophene; each IDTBR is reduced to its 7 aromatic sites:        *)
(*     Rh - BT - Th - Ph - Th - BT - Rh   (Fig. 1).                           *)
(*============================================================================ *)

Clear[p3htAtoms, rh1Atoms, btd1Atoms, idt1Atoms, idbAtoms, rh2Atoms,
      btd2Atoms, idt2Atoms, pbc, dpbc, nearestImage, ringCenter,
      nChain, isEnd, findCloseContacts, GetCloseContactRingPair,
      btdAtoms, thioAtoms];

(*----------------------- Ring atom selectors -------------------------------*)
(* Each function returns the Cartesian coords of the atoms in the k-th ring.  *)
(* P3HT: 11 atoms per monomer; the first 5 are the thiophene ring.            *)
p3htAtoms[k_] := atomP3HT[[11 (k - 1) + 1 ;; 11 (k - 1) + 5]][[All, 2 ;; -1]];

(* IDTBR: 88 atoms per molecule. The 7 aromatic sites are indexed relative to *)
(* the start of molecule k (see atom-name template in 01_import_configuration).*)
rh1Atoms[k_]  := atomIDTBR[[88 (k - 1) + 1  ;; 88 (k - 1) + 5 ]][[All, 2 ;; -1]]; (* left rhodanine *)
btd1Atoms[k_] := atomIDTBR[[88 (k - 1) + 10 ;; 88 (k - 1) + 15]][[All, 2 ;; -1]]; (* left BT       *)
idt1Atoms[k_] := atomIDTBR[[88 (k - 1) + 19 ;; 88 (k - 1) + 23]][[All, 2 ;; -1]]; (* left Th of IDT*)
idbAtoms[k_]  := atomIDTBR[[{88 (k-1)+24, 88 (k-1)+25, 88 (k-1)+27,
                            88 (k-1)+44, 88 (k-1)+45, 88 (k-1)+47}]][[All, 2;;-1]]; (* Ph center *)
rh2Atoms[k_]  := atomIDTBR[[{88 (k-1)+66, 88 (k-1)+67, 88 (k-1)+68,
                            88 (k-1)+70, 88 (k-1)+71}]][[All, 2 ;; -1]];   (* right rhodanine *)
btd2Atoms[k_] := atomIDTBR[[88 (k - 1) + 56 ;; 88 (k - 1) + 61]][[All, 2 ;; -1]]; (* right BT      *)
idt2Atoms[k_] := atomIDTBR[[{88 (k-1)+48, 88 (k-1)+50, 88 (k-1)+53,
                            88 (k-1)+54, 88 (k-1)+55}]][[All, 2 ;; -1]];   (* right Th of IDT*)

(*----------------------- Periodic boundary conditions ---------------------*)
pbc[r_]          := Mod[r, pbc3D];                       (* wrap into box      *)
dpbc[dr_]        := Mod[dr, pbc3D, -pbc3D/2];            (* minimum image      *)
nearestImage[r0_, r_] := r0 + dpbc[r - r0];              (* unwrap r about r0  *)

(*----------------------- Ring centre (minimum image) ----------------------*)
(* Centre of a ring = mean of its atoms after unwrapping them about the first  *)
(* atom, then wrapped back into the box.                                      *)
ringCenter[atomList_] := Module[{posns, newPosns},
  posns   = atomList;
  newPosns = nearestImage[posns[[1]], #] & /@ posns;
  pbc[Mean[newPosns]]
];

(*---- Centres of every thiophene (P3HT) and BT (IDTBR) ring ----*)
ccListP3HT = Array[ringCenter[p3htAtoms[#]] &, 128 * 24];
ccListBTD1 = Array[ringCenter[btd1Atoms[#]] &, 300];
ccListBTD2 = Array[ringCenter[btd2Atoms[#]] &, 300];
ccListBTD  = Join[ccListBTD1, ccListBTD2];   (* 600 BT centres: 1..300 left, 301..600 right *)

(*----------------------- Chain bookkeeping ---------------------------------*)
nChain[k_] := Ceiling[k/24];                            (* chain # from ring # *)
isEnd[k_]  := If[Mod[k, 24, 1] == 1 || Mod[k, 24, 1] == 24, True, False]; (* chain end? *)

(*----------------------- Close-contact search ------------------------------*)
(* Pairs (Th index, BT index) whose centre-to-centre distance <= max (A).     *)
findCloseContacts[k_, max_] := Module[{thiopheneCenter, closeBTDs, indexPairs},
  thiopheneCenter = ccListP3HT[[k]];
  closeBTDs  = Position[ccListBTD, #] & /@
               Select[ccListBTD, Norm[dpbc[# - thiopheneCenter]] <= max &] // Flatten;
  indexPairs = {k, #} & /@ closeBTDs;
  indexPairs
];

(* Cutoff 5 A (paper, Sec. 2.3). *)
ThioBTDcloseContacts = Array[findCloseContacts[#, 5.] &, Length[ccListP3HT]];
ThioBTDcloseContacts = DeleteCases[Flatten[ThioBTDcloseContacts, 1], {}];

(* Exclude contacts near chain ends (thiophene within 3 monomers of an end),  *)
(* leaving 237 donor-acceptor contacts.                                       *)
ThioBTDpairList = Select[ThioBTDcloseContacts,
     ! (isEnd[#[[1]]] || isEnd[#[[1]] - 2] || isEnd[#[[1]] + 2]) &];

(*----------------------- Contact-pair helpers ------------------------------*)
(* GetCloseContactRingPair[i] -> {thiophene index, IDTBR molecule index}.     *)
GetCloseContactRingPair[i_] := Module[{ThioNo, BTDno, IDTBRno},
  {ThioNo, BTDno} = ThioBTDpairList[[i]];
  IDTBRno = If[BTDno > 300, BTDno - 300, BTDno];   (* 301..600 right BT -> molecule 1..300 *)
  {ThioNo, IDTBRno}
];

(* btdAtoms[i] / thioAtoms[i] : ring atoms of the contacting BT / thiophene.  *)
btdAtoms[i_] := Module[{ThioNo, BTDno},
  {ThioNo, BTDno} = ThioBTDpairList[[i]];
  If[BTDno > 300, btd2Atoms[BTDno - 300], btd1Atoms[BTDno]]
];
thioAtoms[i_] := p3htAtoms[ThioBTDpairList[[i, 1]]];

(*----------------------- (Optional) visualization --------------------------*)
(* Scatter of all close-contact pairs (red = thiophene, blue = BT, green =    *)
(* connecting line under minimum image). Uncomment to render.                 *)
(*
plot3D = Graphics3D[{
   Style[Point[ccListP3HT], PointSize[Small],  Red],
   Style[Point[ccListBTD],  PointSize[Medium], Blue]},
  Boxed -> True, Axes -> True, AxesLabel -> {"X","Y","Z"},
  PlotRange -> {{0, pbc3D[[1]]}, {0, pbc3D[[2]]}, {0, pbc3D[[3]]}}];

drawInterfaceContact[n_] := Module[{mid1, mid2},
  {mid1, mid2} = GetCloseContactRingPair[n];
  ringCenter[p3htAtoms[mid1]];
  ListPointPlot3D[{
    Table[nearestImage[ringCenter[idbAtoms[mid2]], #] & /@ p3htAtoms[i],
          {i, mid1 - 3, mid1 + 3}] // Flatten[#, 1] &,
    rh1Atoms[mid2], btd1Atoms[mid2], idt1Atoms[mid2], idbAtoms[mid2],
    idt2Atoms[mid2], btd2Atoms[mid2], rh2Atoms[mid2]},
   BoxRatios -> Automatic,
   PlotStyle -> {Orange, Cyan, Green, Cyan, Cyan, Cyan, Green, Cyan}]
];
*)
