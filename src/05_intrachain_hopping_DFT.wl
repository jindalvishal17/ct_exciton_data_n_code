(* ::Package:: *)

(*============================================================================*)
(*  05_intrachain_hopping_DFT.wl                                              *)
(*  Intrachain hopping matrix elements from DFT (energy-splitting method).    *)
(*                                                                            *)
(*  Paper ref:  Section 2.3 "Hopping matrix parameters", Figs. 4 & 5,         *)
(*              Table 2.  DFT at B3LYP/6-311g(d).                             *)
(*                                                                            *)
(*  Hopping t between two bonded monomers is obtained from the DFT energies  *)
(*  of the anion, cation and neutral ground state of each monomer and of the *)
(*  co-dimer, using the energy-splitting relation                             *)
(*                                                                            *)
(*        t = (1/2) Sqrt[ dE^2 - (eA - eB)^2 ],                               *)
(*        dE = eA + eB - 2 eAB,                                               *)
(*                                                                            *)
(*  where eA, eB are the monomer cation (hole) or anion (electron) energies  *)
(*  relative to the neutral ground state, and eAB is the corresponding dimer *)
(*  energy.  (tDist[] below.)                                                 *)
(*                                                                            *)
(*  Dimer types needed (all per-contact, from MD geometries):                 *)
(*    ThDimer  : thiophene-thiophene along the P3HT backbone                  *)
(*    RhBT1/2  : rhodanine-BT        (left/right of the IDTBR)               *)
(*    ThBT1/2  : thiophene-BT         (left/right of the IDTBR)               *)
(*  plus the constant phenylene-thiophene (Th-Ben) coupling of the rigid IDT  *)
(*  core, computed once from planar optimized geometries.                    *)
(*                                                                            *)
(*  Input files: data/energies_extracted/*.txt  (one energy per row).        *)
(*============================================================================ *)

Clear[tDist, ThDimerGS, ThDimerCation, ThDimerAnion, RhBT1GS, RhBT1Cation,
      RhBT1Anion, RhBT2GS, RhBT2Cation, RhBT2Anion, ThBT1GS, ThBT1Cation,
      ThBT1Anion, ThBT2GS, ThBT2Cation, ThBT2Anion];

eDir = FileNameJoin[{dataDir, "energies_extracted"}];

(*---- Import DFT energies (Hartree) for every dimer configuration ----*)
ThDimerGS     = Import[FileNameJoin[{eDir, "ThDimer_gsE.txt"}], "Table"];
ThDimerCation = Import[FileNameJoin[{eDir, "Unrestricted", "Cation_ThDimer.txt"}], "Table"];
ThDimerAnion  = Import[FileNameJoin[{eDir, "Unrestricted", "Anion_ThDimer.txt"}], "Table"];

RhBT1GS     = Import[FileNameJoin[{eDir, "RhBT1_gsE.txt"}], "Table"];
RhBT1Cation = Import[FileNameJoin[{eDir, "Unrestricted", "cations", "RhBT1_cationU.txt"}], "Table"];
RhBT1Anion  = Import[FileNameJoin[{eDir, "Unrestricted", "anions",  "RhBT1_anionU.txt"}], "Table"];

RhBT2GS     = Import[FileNameJoin[{eDir, "RhBT2_gsE.txt"}], "Table"];
RhBT2Cation = Import[FileNameJoin[{eDir, "Unrestricted", "cations", "RhBT2_cationU.txt"}], "Table"];
RhBT2Anion  = Import[FileNameJoin[{eDir, "Unrestricted", "anions",  "RhBT2_anionU.txt"}], "Table"];

ThBT1GS     = Import[FileNameJoin[{eDir, "ThBT1_gsE.txt"}], "Table"];
ThBT1Cation = Import[FileNameJoin[{eDir, "Unrestricted", "cations", "ThBTD1_cationU.txt"}], "Table"];
ThBT1Anion  = Import[FileNameJoin[{eDir, "Unrestricted", "anions",  "ThBTD1_anionU.txt"}], "Table"];

ThBT2GS     = Import[FileNameJoin[{eDir, "ThBT2_gsE.txt"}], "Table"];
ThBT2Cation = Import[FileNameJoin[{eDir, "Unrestricted", "cations", "ThBTD2_cationU.txt"}], "Table"];
ThBT2Anion  = Import[FileNameJoin[{eDir, "Unrestricted", "anions",  "ThBTD2_anionU.txt"}], "Table"];

(*---- Formatting: flatten each file to a 1-D energy list ----*)
{ThDimerGS     = ThDimerGS[[All, 1]],
 ThDimerCation = ThDimerCation[[All, 1]],
 ThDimerAnion  = ThDimerAnion[[All, 1]]};

{ThBT1GS     = ThBT1GS[[All, 1]],
 ThBT1Cation = ThBT1Cation[[All, 1]],
 ThBT1Anion  = ThBT1Anion[[All, 1]],
 ThBT2GS     = ThBT2GS[[All, 1]],
 ThBT2Cation = ThBT2Cation[[All, 1]],
 ThBT2Anion  = ThBT2Anion[[All, 1]]};

{RhBT1GS     = RhBT1GS[[All, 1]],
 RhBT1Cation = RhBT1Cation[[All, 1]],
 RhBT1Anion  = RhBT1Anion[[All, 1]],
 RhBT2GS     = RhBT2GS[[All, 1]],
 RhBT2Cation = RhBT2Cation[[All, 1]],
 RhBT2Anion  = RhBT2Anion[[All, 1]]};

(* P3HT backbone indices that are interior (exclude the last monomer of each   *)
(* 24-mer chain, which has no forward neighbour).                            *)
ringNums = Select[Range[128 * 24], (Mod[#, 24] != 0) &];

(*======================= Monomer energies (planar, optimized) =============*)
(* B3LYP optimized-geometry cation / anion / ground-state energies (Hartree)  *)
(* for the four monomer types: benzene(phenylene), thiophene, BT, rhodanine. *)
{{BenAnion, BenCation, BenGS} = {-232.227630711, -231.959359844, -232.297902709};
 {ThAnion,  ThCation,  ThGS}  = {-552.999105696, -552.735570737, -553.061979360};
 {BTDAnion, BTDCation, BTDGS} = {-738.836969934, -738.491016167, -738.814871833};
 {RhAnion,  RhCation,  RhGS}  = {-1159.38182101,  -1159.02889428, -1159.34861317};}

(* Monomer cation (hole) / anion (electron) energies relative to neutral,    *)
(* converted to eV (Table 1 onsite energies epsilon_k).                      *)
{{eCationBen, eAnionBen} = fac * {BenCation - BenGS, BenAnion - BenGS};
 {eCationThio, eAnionThio} = fac * {ThCation - ThGS, ThAnion - ThGS};
 {eCationBTD, eAnionBTD} = fac * {BTDCation - BTDGS, BTDAnion - BTDGS};
 {eCationRh, eAnionRh} = fac * {RhCation - RhGS, RhAnion - RhGS};}

(*======================= Energy-splitting hopping function ===============*)
(* tDist[eA, eB, eAB] -> hopping t (eV) from the splitting relation.          *)
tDist[eA_, eB_, eABdist_] := Module[{dE},
  dE = eA + eB - 2 eABdist;
  Sqrt[dE^2 - (eA - eB)^2]/2
];

(*======================= Phenylene-thiophene (rigid IDT core) ============*)
(* The Th-Ph dihedral is fixed, so this coupling is constant (paper Sec.2.3). *)
pOptAnion  = fac * (-784.147882137 + 784.166171694);
pOptCation = fac * (-783.881088822 + 784.166171694);
tHoleThBen = tDist[eCationThio, eCationBen, pOptCation];   (* ~1.28 eV, Table 2 *)
tElecThBen = tDist[eAnionThio, eAnionBen, pOptAnion];     (* ~1.31 eV, Table 2 *)

(*======================= Dimer cation/anion energies (eV) ================*)
eCationP3HT = fac * (ThDimerCation - ThDimerGS);           (* P3HT Th-Th dimers *)
eAnionP3HT  = fac * (ThDimerAnion  - ThDimerGS);

eCationThBT = fac * Join[ThBT1Cation - ThBT1GS, ThBT2Cation - ThBT2GS];
eAnionThBT  = fac * Join[ThBT1Anion  - ThBT1GS, ThBT2Anion  - ThBT2GS];

eCationRhBT = fac * Join[RhBT1Cation - RhBT1GS, RhBT2Cation - RhBT2GS];
eAnionRhBT  = fac * Join[RhBT1Anion  - RhBT1GS, RhBT2Anion  - RhBT2GS];

(*======================= Intrachain hopping distributions ===============*)
(* tHole/tElec for each bond family (eV).                                       *)
tHoleRhBT = tDist[eCationRh, eCationBTD, eCationRhBT];
tElecRhBT = tDist[eAnionRh,  eAnionBTD, eAnionRhBT];
tElecThBT = tDist[eAnionThio, eAnionBTD, eAnionThBT];
tHoleThBT = tDist[eCationThio, eCationBTD, eCationThBT];

(* P3HT backbone Th-Th hopping (broad distribution from dihedral disorder,    *)
(* Figs. 4a/b).                                                              *)
tHoleP3HT = tDist[eCationThio, eCationThio, eCationP3HT];
tElecP3HT = tDist[eAnionThio,  eAnionThio,  eAnionP3HT];

(* Median values (compare Table 2: th_Th = 1.19, te_Th = 1.47 eV).            *)
(* Median[tHoleP3HT[[ringNums]]]   -> ~1.25 eV                                *)
(* Median[tElecP3HT[[ringNums]]]   -> ~1.51 eV                                *)
