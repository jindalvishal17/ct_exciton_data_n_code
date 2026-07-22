(* ::Package:: *)

(*============================================================================*)
(*  07_polaron_initial_condition.wl                                           *)
(*  Polaron (single-carrier) states on the isolated P3HT segment and IDTBR   *)
(*  molecule, used as initial guesses for the CT-exciton minimization.       *)
(*                                                                            *)
(*  Paper ref:  Section 2.2 and the polaron model of refs 39,43,44.          *)
(*  A polaron minimizes the one-carrier energy                                *)
(*                                                                            *)
(*      E_pol = -2 Sum_k t_k psi_k psi_{k+1}                                  *)
(*              - (1/2)(1 - 1/eps) Sum_{i,j} |psi_i|^2 U_ij |psi_j|^2         *)
(*                                                                            *)
(*  i.e. the tight-binding kinetic energy plus the dielectric polarization    *)
(*  self-stabilization of the carrier's own charge distribution. The CT       *)
(*  exciton (module 09) is then minimized starting from these polaron        *)
(*  wavefunctions.                                                           *)
(*============================================================================ *)

Clear[uCs, getIntP3HT, getIntIDTBR, ePol, ePolMin, holePolP3HT, elecPolIDTBR,
      psi0, wv7mer, zero7];

(*----------------------- Smeared Coulomb kernel (single carrier) ----------*)
(* A single fixed smearing length sig is used for the isolated polaron        *)
(* calculation (contrast with the per-site deltaCT of module 03).             *)
sig = 2.0;
uCs[r_]   := Erf[r/(2 sig)]/r;             (* offsite Gaussian Coulomb, Eq. 6 *)
uCs[0]    := 1/(Sqrt[Pi] sig);             (* onsite limit                   *)
uCs[0.]   := 1/(Sqrt[Pi] sig);

(* 7x7 Coulomb kernel for the isolated P3HT segment (sites 1..7 of a contact).*)
getIntP3HT[nThpair_] := Module[{seg},
  seg = getCloseContactPosn[nThpair][[1 ;; 7]];
  eFac * Outer[uCs[Norm[dpbc[#1 - #2]]] &, seg, seg, 1]
];

(* 7x7 Coulomb kernel for the isolated IDTBR molecule (sites 8..14), using    *)
(* the per-monomer deltaIDTBR smearing lengths.                                *)
getIntIDTBR[nThpair_] := Module[{mat, posns, i, j, R},
  posns = getCloseContactPosn[nThpair][[8 ;; 14]];
  mat = eFac * Outer[0 &, posns, posns, 1];
  Do[
    If[i != j,
       R = Norm[dpbc[posns[[i]] - posns[[j]]]];
       mat[[i, j]] = (eFac / R) * Erf[R / ((deltaIDTBR[[i]] + deltaIDTBR[[j]]) / 2)]],
    {i, 7}, {j, 7}];
  Do[mat[[i, i]] = eFac / (Sqrt[Pi] (deltaIDTBR[[i]] / 2)), {i, 7}];
  mat
];

(*----------------------- Polaron energy and minimizer --------------------*)
(* ePol[psi, t, int] : polaron energy for amplitudes psi, hopping list t and  *)
(* Coulomb kernel int (Eq. for E_pol above). psi is normalized internally.    *)
ePol[psiL_, tL_, int_] := Module[{psi, qVals},
  psi   = psiL/Norm[psiL];
  qVals = psi^2;
  -2 Total[tL * psi[[1 ;; -2]] * psi[[2 ;; -1]]]
    - (1/2) (1 - 1/eps) qVals . int . qVals
];

(* Starting guess for a localized polaron on a 7-site chain.                  *)
psi0 = {0.02, 0.164, 0.4468, 0.7359, 0.4468, 0.164, 0.02};

(* ePolMin minimizes the polaron energy subject to normalization, returning    *)
(* {energy, rules} from FindMinimum.                                          *)
ePolMin[psiL_, tL_, int_] := Module[{arg, guess},
  psi0  = psiL/Norm[psiL];
  arg   = Array[psi, 7];
  guess = Transpose[{arg, psi0}];
  FindMinimum[{ePol[arg, tL, int], arg . arg == 1}, guess]
];

(*----------------------- Per-contact polaron solvers ---------------------*)
(* Hole polaron on the P3HT 7-thiophene window of contact n.                  *)
holePolP3HT[nThpair_] := Module[{p3ht, idtbr},
  {p3ht, idtbr} = GetCloseContactRingPair[nThpair];
  ePolMin[psi0, tH[p3ht], getIntP3HT[nThpair]]
];

wv7mer = Table[psi[i], {i, 7}];   (* symbolic 7-site wavevector            *)
zero7  = Table[0, {7}];           (* zero padding for embedding on 14 sites *)

(* Electron polaron on the IDTBR molecule of contact n.                       *)
elecPolIDTBR[nThpair_] := Module[{p3ht, idtbr},
  {p3ht, idtbr} = GetCloseContactRingPair[nThpair];
  ePolMin[psi0, tLidtbr[idtbr], getIntIDTBR[nThpair]]
];
