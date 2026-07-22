(* ::Package:: *)

(*============================================================================*)
(*  04_interchain_hopping_ML.wl                                               *)
(*  Interchain hopping matrix element t_ab between the contacting BT (IDTBR)  *)
(*  and thiophene (P3HT) rings, predicted by a machine-learning model.        *)
(*                                                                            *)
(*  Paper ref:  Section 2.3 "Hopping matrix parameters", Fig. 6.             *)
(*  The ML model estimates the HOMO-HOMO overlap between two non-bonded rings *)
(*  from a small set of geometric descriptors (inter-ring separations and the *)
(*  relative orientation angle). For simplicity the hole and electron         *)
(*  interchain hoppings are taken to be equal (paper, Sec. 2.3).              *)
(*                                                                            *)
(*  This module:                                                              *)
(*    1. computes the geometric descriptors for all 237 contacts,            *)
(*    2. exports them as the ML input file,                                   *)
(*    3. reads the ML-predicted t_ab values (predictedt600.csv).             *)
(*                                                                            *)
(*  Inputs : data/predictedt600.csv  (ML output: one t_ab per contact)        *)
(*  Output: tInterchain  (length-237 list of t_ab in eV)                     *)
(*============================================================================ *)

Clear[thioCOM, btdCOM, thioCC1, thioCC2, btdCC1, btdCC2,
      xyzRing, xyzThio, xyzBTD, rotationmatrix];

(*----------------------- Centre-of-mass and sub-centres --------------------*)
(* thioCOM/btdCOM : centre of the contacting thiophene / BT ring (equal atom  *)
(* weights 0.2 on the 5 ring atoms).                                         *)
(* thioCC1/2, btdCC1/2 : midpoints of the two halves of each ring, used to    *)
(* define the inter-ring geometry descriptors.                                *)
thioCOM[k_] := Module[{posns, newPosns},
  posns = thioAtoms[k];
  newPosns = nearestImage[posns[[1]], #] & /@ posns;
  pbc[Sum[0.2 * newPosns[[i]], {i, 1, 5}]]
];
btdCOM[k_] := Module[{posns, newPosns},
  posns = btdAtoms[k];
  newPosns = nearestImage[posns[[1]], #] & /@ posns;
  pbc[Sum[0.2 * newPosns[[i]], {i, 1, 5}]]
];
thioCC1[k_] := Module[{posns, newPosns},
  posns = thioAtoms[k][[1 ;; 2]];
  newPosns = nearestImage[posns[[1]], #] & /@ posns;
  pbc[Sum[0.5 * newPosns[[i]], {i, 1, 2}]]
];
thioCC2[k_] := Module[{posns, newPosns},
  posns = thioAtoms[k][[3 ;; 4]];
  newPosns = nearestImage[posns[[1]], #] & /@ posns;
  pbc[Sum[0.5 * newPosns[[i]], {i, 1, 2}]]
];
btdCC1[k_] := Module[{posns, newPosns},
  posns = btdAtoms[k][[1 ;; 2]];
  newPosns = nearestImage[posns[[1]], #] & /@ posns;
  pbc[Sum[0.5 * newPosns[[i]], {i, 1, 2}]]
];
btdCC2[k_] := Module[{posns, newPosns},
  posns = btdAtoms[k][[3 ;; 4]];
  newPosns = nearestImage[posns[[1]], #] & /@ posns;
  pbc[Sum[0.5 * newPosns[[i]], {i, 1, 2}]]
];

(*---- Build centre lists for all contacts ----*)
centerposn1Th = Table[thioCC1[i], {i, Length[ThioBTDpairList]}];
centerposn2Th = Table[thioCC2[i], {i, Length[ThioBTDpairList]}];
centerposn1BT = Table[btdCC1[i], {i, Length[ThioBTDpairList]}];
centerposn2BT = Table[btdCC2[i], {i, Length[ThioBTDpairList]}];
centerposnTh  = Table[thioCOM[i], {i, Length[ThioBTDpairList]}];
centerposnBT  = Table[btdCOM[i], {i, Length[ThioBTDpairList]}];

(*----------------------- Distance descriptors -----------------------------*)
distcom12 = Table[Norm[dpbc[centerposnTh[[i]] - centerposnBT[[i]]]],
                  {i, Length[ThioBTDpairList]}];          (* ring COM separation *)
distcc11  = Table[Norm[dpbc[centerposn1Th[[i]] - centerposn1BT[[i]]]],
                  {i, Length[ThioBTDpairList]}];
distcc12  = Table[Norm[dpbc[centerposn1Th[[i]] - centerposn2BT[[i]]]],
                  {i, Length[ThioBTDpairList]}];
distcc21  = Table[Norm[dpbc[centerposn2Th[[i]] - centerposn1BT[[i]]]],
                  {i, Length[ThioBTDpairList]}];
distcc22  = Table[Norm[dpbc[centerposn2Th[[i]] - centerposn2BT[[i]]]],
                  {i, Length[ThioBTDpairList]}];

distctc   = distcom12;
distc1c2  = Table[{distcc11[[i]], distcc12[[i]], distcc21[[i]], distcc22[[i]]},
                  {i, Length[distcc11]}];

(* Reorder the four sub-centre distances so that the closest pair comes      *)
(* first; its complementary pairs follow (used as ML features).               *)
distc1c2mod = Table[
  Module[{p1, p2, p3, p4},
    p1 = Flatten[Position[distc1c2[[i]], Min[distc1c2[[i]]]]];
    p2 = 5 - # & /@ p1;
    p3 = If[OddQ[p1] == {True}, p1 + 1, p1 - 1];
    p4 = If[OddQ[p2] == {True}, p2 + 1, p2 - 1];
    distc1c2[[i, Join[p1, p2, p3, p4]]]],
  {i, 1, Length[distc1c2]}];

(*----------------------- Orientation (yaw-pitch-roll) ---------------------*)
(* xyzRing builds a local orthonormal frame (normal, in-plane axis, third)   *)
(* from three ring atoms; the relative frame rotation gives the contact      *)
(* angle used as an ML feature.                                               *)
xyzRing[a_, b_, c_] := {
  Cross[b - a, c - a]/Norm[Cross[b - a, c - a]],
  (b - a)/Norm[b - a],
  Cross[Cross[b - a, c - a]/Norm[Cross[b - a, c - a]], (b - a)/Norm[b - a]]};

xyzThio[k_] := Module[{a, b, c, posns, newPosns},
  posns = thioAtoms[k];
  newPosns = nearestImage[posns[[1]], #] & /@ posns;
  {a, b, c} = newPosns[[{1, 4, 5}]];
  {Cross[b - a, c - a]/Norm[Cross[b - a, c - a]],
   (b - a)/Norm[b - a],
   Cross[Cross[b - a, c - a]/Norm[Cross[b - a, c - a]], (b - a)/Norm[b - a]]}
];
xyzBTD[k_] := Module[{a, b, c, posns, newPosns},
  posns = btdAtoms[k];
  newPosns = nearestImage[posns[[1]], #] & /@ posns;
  {a, b, c} = {newPosns[[1]], newPosns[[4]], Mean[{newPosns[[5]], newPosns[[6]]}]};
  {Cross[b - a, c - a]/Norm[Cross[b - a, c - a]],
   (b - a)/Norm[b - a],
   Cross[Cross[b - a, c - a]/Norm[Cross[b - a, c - a]], (b - a)/Norm[b - a]]}
];
xyzvecTh = Table[xyzThio[k], {k, Length[ThioBTDpairList]}];
xyzvecBT = Table[xyzBTD[k],  {k, Length[ThioBTDpairList]}];

rotationmatrix[{a_, b_, c_}, {x_, y_, z_}] := {x, y, z} . Inverse[{a, b, c}];
eulerangle = Table[RollPitchYawAngles[rotationmatrix[xyzvecTh[[i]], xyzvecBT[[i]]]],
                   {i, 1, Length[ThioBTDpairList]}];

(*----------------------- ML input / output ---------------------------------*)
(* Six features per contact: the four reordered sub-centre distances, the     *)
(* ring COM distance (all scaled by 0.1 nm), and |cos(theta)| of the relative *)
(* orientation. Exported for the external ML model.                          *)
InputToML = Table[Flatten[{0.1 distc1c2mod[[i]], 0.1 distctc[[i]],
                          Abs[Cos[eulerangle[[i, 3]]]]}],
                  {i, 1, Length[ThioBTDpairList]}];

Export[FileNameJoin[{resultsDir, "SixParametersCTexciton1.txt"}], InputToML, "Table"];

(* Read the ML-predicted interchain hopping t_ab (eV) for every contact.      *)
tInterchain = Import[FileNameJoin[{dataDir, "predictedt600.csv"}]] // Flatten;
