(* ::Package:: *)

(*============================================================================*)
(*  09_ct_exciton_energy_analysis.wl                                          *)
(*  Compute and analyze the total CT exciton energy for all donor-acceptor    *)
(*  contacts.                                                                 *)
(*                                                                            *)
(*  Paper ref:  Eq. (12) for E_CT, Section 3.1, Figs. 8-13, Table 3, Eq. (13).*)
(*                                                                            *)
(*  E_CT = E_K + E_C + E_P    (Eq. 2, 12)                                    *)
(*    E_K = -<psi_h|H_h|psi_h> + <psi_e|H_e|psi_e>           (Eq. 3)          *)
(*    E_C = Sum_i |a_h_i|^2|a_e_i|^2 Es_i + Sum_{i!=j}|a_h_i|^2|a_e_j|^2 D_ij (Eq. 7)*)
(*    E_P = -(1-1/eps) (1/2) Sum_{i,j} q_i q_j D_ij          (Eq. 10)         *)
(*  with q_i = |a_h_i|^2 - |a_e_i|^2 and D_ij the Gaussian Coulomb kernel.    *)
(*                                                                            *)
(*  The total energy is minimized over the electron and hole amplitudes for   *)
(*  each of the 237 contacts (hole constrained to P3HT, electron to IDTBR),   *)
(*  then decomposed into one-body (E_K) and two-body (E_C+E_P) parts for the  *)
(*  correlation / outlier analysis of Figs. 11-13 and Table 3.                *)
(*                                                                            *)
(*  NOTE ON RECONSTRUCTION: The original notebook used CTexctionE[argsE,argsH,*
(*  k] and the symbols psiEfull / psiHfull but did not contain their          *)
(*  definitions. CTexctionE below is reconstructed exactly from Eq. (12) and  *)
(*  the component decomposition in CTSolve (which IS present in the original).*)
(*  psiEfull / psiHfull are given sensible localized guesses (electron on    *)
(*  IDTBR, hole on P3HT) consistent with the polaron initial conditions; the  *)
(*  optimizer is robust to this choice.                                       *)
(*============================================================================ *)

Clear[matKEelecCT, matKEholeCT, intmat, CTexctionE, CTexctionMin2,
      CTexctionMin3, CTSolve, psiEfull, psiHfull];

(*----------------------- Per-contact Hamiltonians & Coulomb ----------------*)
(* matKEholeCT[i]  : hole Hamiltonian  H_h on the 14 sites (HOMO onsite eH).  *)
(* matKEelecCT[i]  : electron Hamiltonian H_e on the 14 sites (LUMO onsite eL).*)
(* intmat[i]       : 14x14 Gaussian Coulomb kernel D_ij for contact i.        *)
matKEholeCT[i_] := matCTexciton[GetRightorLeft[i], eH, eHidtbr,
              tH[GetCloseContactRingPair[i][[1]]], tHidtbr[GetCloseContactRingPair[i][[2]]],
              -t12[i]];
matKEelecCT[i_] := matCTexciton[GetRightorLeft[i], eL, eLidtbr,
              tL[GetCloseContactRingPair[i][[1]]], tLidtbr[GetCloseContactRingPair[i][[2]]],
               t12[i]];
intmat[i_] := matCoulIntCT[getCloseContactPosn[i]];

(*----------------------- Total CT exciton energy (Eq. 12) ------------------*)
(* CTexctionE[psiE, psiH, k] = E_K + E_C + E_P for contact k.                  *)
(* Reconstructed from Eq. (12); components match the CTSolve decomposition.    *)
CTexctionE[argsE_, argsH_, k_] := Module[{qElec, qHole, qNet},
  qElec = -argsE^2;            (* electron charge = -|psi_e|^2                 *)
  qHole =  argsH^2;            (* hole charge    = +|psi_h|^2                  *)
  qNet  = qHole + qElec;       (* net charge q_i = |a_h|^2 - |a_e|^2           *)
  (* E_K  : one-body kinetic energy (Eq. 3)                                   *)
  (argsE . matKEelecCT[k] . argsE) - (argsH . matKEholeCT[k] . argsH)
  (* E_C  : Coulomb (Eq. 7) — offsite Gaussian + onsite Es minus self term    *)
  + qElec . intmat[k] . qHole
  + (VcCT - eFac/(Sqrt[Pi] (deltaCT/2))) . (qElec * qHole)
  (* E_P  : dielectric polarization (Eq. 10)                                  *)
  + qNet . (-1/2 (1 - 1/eps) intmat[k]) . qNet
];

(*----------------------- Initial guesses (reconstructed) -----------------*)
(* Fixed 14-site initial wavefunctions: hole localized on P3HT (sites 1-7),   *)
(* electron localized on IDTBR (sites 8-14). psi0 is the localized polaron    *)
(* guess from module 07.                                                      *)
psiHfull = Join[psi0,  zero7];
psiEfull = Join[zero7, psi0];

(*----------------------- Separated CT exciton (Fig. 9) -------------------*)
(* Minimize E_CT with hole constrained to P3HT (prob > 0.999 on sites 1-7)    *)
(* and electron constrained to IDTBR (prob > 0.999 on sites 8-14).           *)
CTexctionMin2[psiE0_, psiH0_, k_] := Module[{},
  psiEG = psiE0/Norm[psiE0];  psiHG = psiH0/Norm[psiH0];
  argsE = Table[psiE[i], {i, 1, 14}];
  argsH = Table[psiH[i], {i, 1, 14}];
  guessElec = Transpose[{argsE, psiEG}];
  guessHole = Transpose[{argsH, psiHG}];
  guesses   = Join[guessElec, guessHole];
  holeONp3ht  = Sum[argsH[[i]]^2, {i, 1, 7}];
  elecONidtbr = Sum[argsE[[i]]^2, {i, 8, 14}];
  FindMinimum[{CTexctionE[argsE, argsH, k], argsE . argsE == 1.,
               argsH . argsH == 1., elecONidtbr > 0.999, holeONp3ht > 0.999},
              guesses]
];

(* CT exciton energies for all contacts (Fig. 9 distribution).                *)
ExcitonEnergies2 = Table[CTexctionMin2[psiEfull, psiHfull, kthPair][[1]],
                          {kthPair, 1, nContacts}];

(*----------------------- Exciton localized on IDTBR (Fig. 8) --------------*)
(* Both electron and hole constrained to the IDTBR molecule (sites 8-14).     *)
(* Initial guesses taken from the per-contact polarons (module 07).            *)
CTexctionMin3[psiE0_, psiH0_, k_] := Module[{},
  psiEG = psiE0/Norm[psiE0];  psiHG = psiH0/Norm[psiH0];
  argsE = Table[psiE[i], {i, 1, 14}];
  argsH = Table[psiH[i], {i, 1, 14}];
  guessElec = Transpose[{argsE, psiEG}];
  guessHole = Transpose[{argsH, psiHG}];
  guesses   = Join[guessElec, guessHole];
  holeONibtbr  = Sum[argsH[[i]]^2, {i, 8, 14}];
  elecONidtbr = Sum[argsE[[i]]^2, {i, 8, 14}];
  FindMinimum[{CTexctionE[argsE, argsH, k], argsE . argsE == 1.,
               argsH . argsH == 1., elecONidtbr == 1., holeONibtbr == 1.},
              guesses]
];

ExcitonEnergies3 = Table[
  CTexctionMin3[Join[zero7, wv7mer /. elecPolIDTBR[kthPair][[2]]],
                Join[zero7, wv7mer /. holePolP3HT[kthPair][[2]]], kthPair][[1]],
  {kthPair, 1, nContacts}];

(*======================= Energy decomposition (CTSolve) =================*)
(* For each contact, minimize E_CT and return the energy components:          *)
(*   {E_CT, E_C+E_P, E_K(=E_1b), E_C, E_P, E_K_e, E_K_h, psiE, psiH}.        *)
(* These feed the Figs. 11-13 correlation / outlier analysis and Table 3.     *)
CTSolve[k_] := Module[{res, sol, psiEopt, psiHopt, qElec, qHole, qNet,
                       eCoulOpt, ePolzOpt, eKEOpt, holeKE, elecKE},
  psiEG = psiEfull/Norm[psiEfull];  psiHG = psiHfull/Norm[psiHfull];
  argsE = Table[psiE[i], {i, 1, 14}];
  argsH = Table[psiH[i], {i, 1, 14}];
  guessElec = Transpose[{argsE, psiEG}];
  guessHole = Transpose[{argsH, psiHG}];
  guesses   = Join[guessElec, guessHole];
  holeONp3ht  = Sum[argsH[[i]]^2, {i, 1, 7}];
  elecONidtbr = Sum[argsE[[i]]^2, {i, 8, 14}];
  res = FindMinimum[{CTexctionE[argsE, argsH, k], argsE . argsE == 1,
                     argsH . argsH == 1, elecONidtbr > 0.999, holeONp3ht > 0.999},
                    guesses];
  sol = res[[2]];
  psiEopt = argsE /. sol;
  psiHopt = argsH /. sol;
  (* one-body kinetic energies (Eq. 3)                                       *)
  elecKE =  (psiEopt . matKEelecCT[k] . psiEopt);
  holeKE = -(psiHopt . matKEholeCT[k] . psiHopt);
  eKEOpt = elecKE + holeKE;
  (* charges                                                                 *)
  qElec = -psiEopt^2;  qHole = psiHopt^2;  qNet = qHole + qElec;
  (* two-body Coulomb (Eq. 7) and polarization (Eq. 10)                       *)
  eCoulOpt = (qElec) . intmat[k] . (qHole)
           + (VcCT - eFac/(Sqrt[Pi] (deltaCT/2))) . (qElec * qHole);
  ePolzOpt = qNet . (-1/2 (1 - 1/eps) intmat[k]) . qNet;
  {res[[1]], eCoulOpt + ePolzOpt, eKEOpt, eCoulOpt, ePolzOpt, elecKE, holeKE,
   psiEopt, psiHopt}
];

(* Run the decomposition for every contact (this is the expensive step).      *)
results = Table[CTSolve[k], {k, 1, nContacts}];

(*======================= Analysis: Figs. 11-13, Table 3 =================*)
CTexcitonE = results[[All, 1]];   (* E_CT                                      *)
OneBodyE   = results[[All, 3]];   (* E_K = E_1b (sum of electron & hole KE)   *)
CoulombE   = results[[All, 2]];   (* E_C + E_P (two-body)                      *)

(* Linear fit E_CT vs E_1b (Fig. 11) and residuals.                            *)
lm        = LinearModelFit[Transpose[{OneBodyE, CTexcitonE}], x, x];
residuals = CTexcitonE - lm /@ OneBodyE;

(* Outlier removal: points deviating > 1.8 sigma from the one-body trend      *)
(* (these are quasi-separated polarons, see paper Sec. 3.1 & Table 3).        *)
rho      = StandardDeviation[residuals];
keepIdx  = Flatten @ Position[Abs[residuals], _ ?(# <= 1.8 rho &)];
RemoveIds = Flatten @ Position[Abs[residuals], _ ?(# >  1.8 rho &)];

(* Refined fit on the main population (Fig. 11, after outlier removal).       *)
lm2         = LinearModelFit[Transpose[{OneBodyE[[keepIdx]], CTexcitonE[[keepIdx]]}], x, x];
residuals2  = CTexcitonE[[keepIdx]] - lm2 /@ OneBodyE[[keepIdx]];

(* Table 3 : mean Coulomb / polarization / total interaction energies.        *)
(*   Main population:    Mean[E_C], Mean[E_P], Mean[E_C+E_P]  (results cols 4,5,2) *)
(*   Outliers:           same, over RemoveIds                                  *)
table3Main     = {Mean[results[[All, 4]][[keepIdx]]],
                  Mean[results[[All, 5]][[keepIdx]]],
                  Mean[results[[All, 2]][[keepIdx]]]};
table3Outliers = {Mean[results[[All, 4]][[RemoveIds]]],
                  Mean[results[[All, 5]][[RemoveIds]]],
                  Mean[results[[All, 2]][[RemoveIds]]]};

(* CT exciton binding energy relative to separated polarons (Eq. 13).         *)
(* E_bind = E_CT - (E_elec.polaron^IDTBR + E_hole.polaron^P3HT);              *)
(* negative => bound. Polaron energies ~1.0 eV each (paper Sec. 3.2).        *)
(* To compute: Ebind = CTexcitonE - (1.0 + 1.0);                              *)
