(* ::Package:: *)

(* ::Section:: *)
(*Package setup*)


BeginPackage["EarthShadow`"]


VelDistFree::usage = 
"VelDistFree[v, \[Theta], \[Phi], \[Gamma]] calculates the free SHM velocity distribution (in [km/s]^-3).

The parameters describing the SHM can be changed in the EarthShadow.m package file.
If you do make changes, don't forget to delete and recalculate PhiIntegrals.dat.

Input:
	-v: DM speed measured in km/s
	-\[Theta], \[Phi]: angular coordinates (in Radians) for the DM direction
		    in a coordinate system with the detector along the z-axis
		    and the average DM velocity x-z plane
	-\[Gamma]: angle between z-axis and average DM velocity (see Fig. 1)
";


SpeedDistFree::usage = 
"SpeedDistFree[v] calculates the free SHM speed distribution (in [km/s]^-1).

This corresponds to the velocity distribution integrated over all directions.
The parameters describing the SHM can be changed in the EarthShadow.m package file.
If you do make changes, don't forget to delete and recalculate PhiIntegrals.dat.

Input:
	-v:   DM speed measured in km/s
";


VelDistAttenuation::usage = 
"VelDistAttenuation[v, \[Theta], \[Phi], \[Gamma], mx, op, c] calculates difference between the free and attenuated DM velocity distribution (in [km/s]^-3).

NB: VelDistAttenuation calculates the attenuated distribution MINUS the free distribution.
So you have to sum up VelDistFree + VelDistAttenuation + VelDistDeflection to get the full velocity distribution after scattering.

Input: 
	-v: DM speed measured in km/s
	-\[Theta], \[Phi]: angular coordinates (in Radians) for the DM direction
		    in a coordinate system with the detector along the z-axis
		    and the average DM velocity x-z plane
	-\[Gamma]: angle between z-axis and average DM velocity (see Fig. 1)
	-mx: DM mass in GeV (in the range [0.1, 300])
    -op: DM-nucleon operator (these are listed in Table 1 of the paper)
	-c: DM-nucleon isoscalar coupling (in units of 1/GeV^2)
";


VelDistDeflection::usage = 
"VelDistDeflection[v, \[Theta], \[Phi], \[Gamma], mx, op, c] calculates the contribution to the DM velocity distribution (in [km/s]^-3) from deflection.

Input:
	-v: DM speed measured in km/s
	-\[Theta], \[Phi]: angular coordinates (in Radians) for the DM direction
		    in a coordinate system with the detector along the z-axis
		    and the average DM velocity x-z plane
	-\[Gamma]: angle between z-axis and average DM velocity (see Fig. 1)
	-mx: DM mass in GeV (in the range [0.1, 300])
    -op: DM-nucleon operator (these are listed in Table 1 of the paper)
	-c: DM-nucleon isoscalar coupling (in units of 1/GeV^2)
";


SpeedDistAttenuation::usage = 
"SpeedDistAttenuation[v, \[Gamma], mx, op, c] calculates the difference between the free and attenuated DM speed distribution (in [km/s]^-1) 

NB: SpeedDistAttenuation calculates the attenuated distribution MINUS the free distribution.
So you have to sum up SpeedDistFree + SpeedDistAttenuation + SpeedDistDeflection to get the full speed distribution after scattering.

Input: 
	-v: DM speed measured in km/s
	-\[Gamma]: angle between z-axis and average DM velocity (see Fig. 1)
	-mx: DM mass in GeV (in the range [0.1, 300])
    -op: DM-nucleon operator (these are listed in Table 1 of the paper)
	-c: DM-nucleon isoscalar coupling (in units of 1/GeV^2)

";


SpeedDistDeflection::usage = 
"SpeedDistDeflection[v, \[Gamma], mx, op, c] calculates the contribution to the DM speed distribution (in [km/s]^-1) from deflection.

Input: 
	-v: DM speed measured in km/s
	-\[Gamma]: angle between z-axis and average DM velocity (see Fig. 1)
	-mx: DM mass in GeV (in the range [0.1, 300])
    -op: DM-nucleon operator (these are listed in Table 1 of the paper)
	-c: DM-nucleon isoscalar coupling (in units of 1/GeV^2)
";



Pc\[Alpha]plus::usage = 
"Pc\[Alpha]plus[c\[Alpha], mx, iso, op, v] calculates the probability of DM deflecting through an angle \[Alpha] (for the positive kinematic solution in Eq. 3.17 of the paper).

Input:
	-c\[Alpha]: cosine of the deflection angle \[Alpha]
	-mx: DM mass in GeV (in the range [0.1, 300])
	-iso: index of the Earth isotope (these are listed in Isotopes.txt in the data folder
	-op: DM-nucleon operator (these are listed in Table 1 of the paper)
	-v: DM speed measured in km/s
";

Pc\[Alpha]minus::usage = 
"Pc\[Alpha]minus[c\[Alpha], mx, iso, op, v] calculates the probability of DM deflecting through an angle \[Alpha] (for the negative kinematic solution in Eq. 3.17 of the paper).

Input:
	-c\[Alpha]: cosine of the deflection angle \[Alpha]
	-mx: DM mass in GeV (in the range [0.1, 300])
	-iso: index of the Earth isotope (these are listed in Isotopes.txt in the data folder
	-op: DM-nucleon operator (these are listed in Table 1 of the paper and in Operaters.txt)
	-v: DM speed measured in km/s";

pscat::usage = "pscat[mx, op, c] calculates the average probability of DM scattering in the Earth.

Input:
	-mx: DM mass in GeV (in the range [0.1, 300])
	-op: DM-nucleon operator (these are listed in Table 1 of the paper)
	-c: DM-nucleon isoscalar coupling (in units of 1/GeV^2)
";


Begin["`Private`"]


(* ::Section:: *)
(*Initialisation*)


(* ::Subsection:: *)
(*Some global parameters*)


(*Earth parameters*)
Subscript[R, E] = 6378.0;
Subscript[l, D] = 0.0;

(*Velocity distribution*)
ve = 220.0;
\[Sigma]v = 156.0;
vesc = 533.0;

(*Normalisation constant for the velocity distribution*)
Nesc = Erf[vesc/(Sqrt[2]*\[Sigma]v)] - Sqrt[2.0/\[Pi]]*(vesc/\[Sigma]v)*Exp[-vesc^2/(2\[Sigma]v^2)];
N1 = 1.0/(\[Sigma]v^3 Sqrt[2\[Pi]]);
N1 = N1/Nesc;

(*Location of data files*)
DataDirectory = NotebookDirectory[] <> "../data/";
PrintTemporary["Loading data files from: " <> DataDirectory];
PrintTemporary["This may take a couple of minutes..."];


(* ::Subsection:: *)
(*Load data tables*)


(* ::Subsubsection:: *)
(*Load data for isotopes*)


PrintTemporary["   Loading isotopes and density profiles..."];


(* ::Text:: *)
(*Load in list of isotopes*)


(* ::Text::Bold:: *)
(*NB: we fix the number of isotopes at 8 (i.e. we read data for the first 8 in the list). These are the ones relevant for Earth-scattering. The 9th in the list (iso = 12) is Xenon, which might be relevant for event rate calculations, but isn't needed here...*)


ClearAll[Avals]
isodata =  Import[DataDirectory <> "Isotopes.txt", "Table"];
isovals = isodata[[2;;9,1]]; (*Read only the first 8 isotopes...*) 
Table[Avals[isovals[[i]]] = isodata[[1+i,2]];, {i,Length[isovals]}];
niso = Length[isovals];


(* ::Text:: *)
(*Load in nbar values*)


nbardata =  Import[DataDirectory <> "/dens_profiles/nbar.dat", "Table"][[;;,2]];
Table[nbar[isovals[[i]]] = nbardata[[i]], {i, niso}];


(* ::Text:: *)
(*Load in tables for density correction factor D (cos\[Theta]) and define interpolation functions*)


Dtab = Table[ Import[DataDirectory <> "/dens_corrections/D_" <> ToString[isovals[[i]]] <>".dat", "Table"], {i, niso}] ;

Table[Dinterp[isovals[[i]]] = Interpolation[Dtab[[i]], InterpolationOrder->1], {i, niso}];


(* ::Subsubsection:: *)
(*Load in data for NR operators*)


PrintTemporary["   Loading operators and cross sections..."];


(* ::Text:: *)
(*Load in list of NR operators *)


opvals =  Import[DataDirectory <> "/Operators.txt", "Table"][[1]];
nop = Length[opvals];


(* ::Text:: *)
(*Load in tables of NR operators*)


(*Differential cross sections*)
ClearAll[xsecTable]
Table[xsecTable[i][j] = Import[DataDirectory <> "/dcross_sections/dsigmadE_" <>  ToString[j] <> "_" <> ToString[i] <>".dat", "Table"] ,{i,isovals},{j, opvals}] ;

(*Total cross sections*)
ClearAll[totxsectab]
Table[ totxsectab[i][j] = Import[DataDirectory <> "/totcross_sections/sigma_" <> ToString[j]  <> "_" <> ToString[i] <>".dat", "Table"], {i, isovals}, {j, opvals}] ;


(* ::Text:: *)
(*Define cross section interpolation tables*)


(*The two functions are the 1/v^2 and v_T^2/v^2 pieces of the xsection*)
Table[xsecinterp1[i][j] = Interpolation[xsecTable[i][j][[All,{1,2,3}]], InterpolationOrder->{1,1}];
xsecinterp2[i][j] = Interpolation[xsecTable[i][j][[All,{1,2,4}]], InterpolationOrder->{1,1}];, {i, isovals}, {j, opvals}];

(*Interpolation for total cross section as a function of m\[Chi], v*)
Table[totxsecinterp[i][j] = Interpolation[totxsectab[i][j][[All,{1,2,3}]],InterpolationOrder->{1,1}],{i,isovals}, {j,opvals}];


(* ::Subsubsection:: *)
(*Import/calculate tabulated integrals over \[Phi]'*)


PrintTemporary["   Loading PhiIntegrals.dat..."];


(* ::Text:: *)
(*Tabulate or import the integral of f(v', \[Theta]', \[Phi]') over \[Phi]', depending on whether it exists...*)


filen = DataDirectory <> "/PhiIntegrals.dat";
If[FileExistsQ[filen],
integralTab= Import[filen];,
integralTab = Flatten[Table[{co, \[Phi]max,NIntegrate[Exp[co Cos[\[Phi]]], {\[Phi], 0, \[Phi]max}]},{co, -7,7,0.05},{\[Phi]max, 0, \[Pi]+0.1, 0.05}],1];
Export[filen, integralTab, "Table"];];
IntFun\[Phi] = Interpolation[integralTab,InterpolationOrder->1];


(* ::Section:: *)
(*Auxiliary functions*)


(* ::Subsubsection::Closed:: *)
(*Distance to surface (d)*)


d1[\[Theta]_]:=
(
\[CapitalDelta] = 2*Subscript[R, E]*Subscript[l, D] - Subscript[l, D]^2 + (Subscript[R, E]-Subscript[l, D])^2*Cos[\[Theta]]^2;
d = (Subscript[R, E]- Subscript[l, D])Cos[\[Theta]] + Sqrt[\[CapitalDelta]]
)


(* ::Subsubsection::Closed:: *)
(*Effective distance to surface (deff)*)


(*Distance to surface, taking into account density correction factor*)
deff[\[Theta]_, iso_]:= d1[\[Theta]]*DensityCorrection[\[Theta], iso];


(* ::Subsubsection::Closed:: *)
(*Density correction function *)


(*Need to read in tables from /dens_corrections/ first...*)
ClearAll[DensityCorrection]
DensityCorrection[\[Theta]2_, iso_]:=
(
Piecewise[{{ 0, \[Theta]2 >= \[Pi]/2},
{ Dinterp[iso][\[Theta]2], \[Theta]2 < \[Pi]/2}}]
)


(* ::Subsubsection::Closed:: *)
(*Distance of interaction point from centre of Earth*)


rEarth[r_, \[Theta]2_]:= Sqrt[r^2 + Subscript[R, E]^2 - 2 r Subscript[R, E]Cos[\[Theta]2]];


(* ::Subsubsection::Closed:: *)
(*Calculating cos\[Alpha] (from angular coordinates)*)


cosine\[Alpha][\[Theta]2_, \[Phi]2_, \[Theta]1_, \[Phi]1_]:=
(
Sin[\[Theta]1]*Sin[\[Theta]2]*Cos[\[Phi]1 - \[Phi]2] + Cos[\[Theta]1]*Cos[\[Theta]2]
)


(* ::Subsubsection::Closed:: *)
(*Calculating ER from cos\[Alpha]*)


ERplus[c\[Alpha]_, mx_, mA_, v_]:=
(
vrel = v/3*^5;
1*^6*(mx*(1-c\[Alpha]^2) + mA - c\[Alpha]*Sqrt[mA^2 - mx^2*(1-c\[Alpha]^2)])*mx^2*vrel^2/(mx+mA)^2
)

ERminus[c\[Alpha]_, mx_, mA_, v_]:=
(
vrel = v/3*^5;
1*^6*(mx*(1-c\[Alpha]^2) + mA + c\[Alpha]*Sqrt[mA^2 - mx^2*(1-c\[Alpha]^2)])*mx^2*vrel^2/(mx+mA)^2
)


(* ::Subsubsection::Closed:: *)
(*Calculating dER/dcos\[Alpha] (for changing variables)*)


dEdc\[Alpha]plus[c\[Alpha]_, mx_, mA_, v_]:=
(
vrel = v/3*^5;
1*^6*(2*mx*c\[Alpha]*Sqrt[mA^2 - mx^2(1-c\[Alpha]^2)] + (mx^2*(2c\[Alpha]^2 -1) + mA^2))*mx^2*vrel^2/((mx+mA)^2*Sqrt[mA^2 - mx^2(1-c\[Alpha]^2)])
)

dEdc\[Alpha]minus[c\[Alpha]_, mx_, mA_, v_]:=
(
vrel = v/3*^5;
-1*^6*(2*mx*c\[Alpha]*Sqrt[mA^2 - mx^2(1-c\[Alpha]^2)] -(mx^2*(2c\[Alpha]^2 -1) + mA^2))*mx^2*vrel^2/((mx+mA)^2*Sqrt[mA^2 - mx^2(1-c\[Alpha]^2)])
)


(* ::Subsubsection::Closed:: *)
(*Minimum value of cos\[Alpha]*)


(*Set by kinematics*)
Minc\[Alpha][mx_, mA_]:=
(
If[mA > mx, -1, Sqrt[1-mA^2/mx^2]]
)

(*If you're reading this, then you decided to actually look at the code. 
You're not supposed to do that. We're supposed to talk about making code
publicly available, but you're not supposed to actually LOOK AT IT or USE
IT! Shame on you...*)


(* ::Subsubsection::Closed:: *)
(*Velocity ratio \[Kappa]*)


(* \[Kappa] = v/v' (fixed by kinematics) *)
\[Kappa]plus[c\[Alpha]_,mx_, mA_] :=
(
(mx + mA)/(mx*c\[Alpha] + Sqrt[mx^2(c\[Alpha]^2 - 1.0) + mA^2])
)

\[Kappa]minus[c\[Alpha]_,mx_, mA_] :=
(
(mx + mA)/(mx*c\[Alpha] -Sqrt[mx^2(c\[Alpha]^2 - 1.0) + mA^2])
)


(* ::Subsubsection::Closed:: *)
(*Transverse velocity*)


(*vT in natural units*)
vTsquared[mx_, mA_,v_, Er_]:=
(
\[Mu] = (mx*mA)/(mx + mA);
v^2  - (1*^-6)*2.0*mA*Er/(4 \[Mu]^2)
)


(* ::Subsubsection::Closed:: *)
(*Maximum recoil energy*)


(*Max recoil energy set by kinematics...*)
ERmax[mx_,iso_,v_]:=
(
mA = 0.9312*Avals[iso];
1*^6*2.0*(v/3*^5)^2*(mA*mx^2)/(mx+mA)^2
)


(* ::Subsubsection::Closed:: *)
(*Differential cross section*)


d\[Sigma]dE[mx_,iso_,op_,v_,Er_]:=
(
(*If[Er > ERmax[mx, iso, v], Return[0];];*)
vrel = v/(3*^5);
mA = 0.9312*Avals[iso];
(xsecinterp1[iso][op][mx, Er] + vTsquared[mx, mA, vrel, Er]*xsecinterp2[iso][op][mx, Er])/(vrel^2)
)


(* ::Subsubsection::Closed:: *)
(*Total cross section*)


Clear[\[Sigma]tot]
\[Sigma]tot[mx_,iso_,op_,v_] :=
(
vrel = v/(3*^5);
totxsecinterp[iso][op][mx, vrel]
)


(* ::Subsubsection::Closed:: *)
(*Mean free path (for each element)*)


(*Coupling c in units of GeV^-2*)
(*Returns 1/MFP in km*)
MFPinv[op_, c_, v_, mx_, iso_]:=
(
cmtokm = 10^-5;
cscaled = 0.5*c*246.2^2;
(cscaled^2)*(nbar[iso]*\[Sigma]tot[mx,iso,op,v])/cmtokm
)


(* ::Subsubsection::Closed:: *)
(*Probability distribution of cos\[Alpha]*)


(* Probability of scattering through angle cos\[Alpha] = c\[Alpha] *)
(* Normalised to 1 over cos\[Alpha] \[Element] [Minc\[Alpha], 1] *)
Pc\[Alpha]plus[c\[Alpha]_, mx_, iso_, op_, v_]:=
(
mA =  0.9312*Avals[iso];
Abs[d\[Sigma]dE[mx,iso,op,v,ERplus[c\[Alpha], mx, mA, v]]*dEdc\[Alpha]plus[c\[Alpha], mx, mA, v]]/\[Sigma]tot[mx, iso, op, v]
)

Pc\[Alpha]minus[c\[Alpha]_, mx_, iso_,op_, v_]:=
(
mA =  0.9312*Avals[iso];
Abs[d\[Sigma]dE[mx,iso,op,v,ERminus[c\[Alpha], mx, mA, v]]*dEdc\[Alpha]minus[c\[Alpha], mx, mA, v]]/\[Sigma]tot[mx, iso, op, v]
)


(* ::Subsubsection::Closed:: *)
(*Mean scattering probability pscat*)


pscat[mx_, op_, c_?NumericQ]:=
(
(*vavg is the average DM speed, calculated elsewhere...*)
(2.0/vavg)*NIntegrate[ v*SpeedDistFree[v]*(1-Exp[-Sum[MFPinv[op, c , v, mx, iso]*deff[\[Theta], iso],{iso, isovals}]])Cos[\[Theta]]Sin[\[Theta]], {v,0,753},{\[Theta], 0, \[Pi]/2},Method->{"QuasiMonteCarlo", MaxPoints->10000}]  
)


(* ::Subsubsection::Closed:: *)
(*Integrands for calculating deflection distribution*)


(* ::Text:: *)
(*fDintegrand1 is used for calculating the deflected velocity distribution*)


ClearAll[fDintegrand1]
fDintegrand1[v2_, \[Theta]2_, \[Phi]2_, \[Theta]1_, \[Phi]1_, \[Gamma]_, mx_, iso_, op_, cint_]:=
(
c\[Alpha] = cosine\[Alpha][\[Theta]2, \[Phi]2, \[Theta]1, \[Phi]1];
mA =  0.9312*Avals[iso];

Piecewise[{{
kp = \[Kappa]plus[c\[Alpha],mx, mA];
aplus = (kp^4)*(Pc\[Alpha]plus[c\[Alpha],mx, iso, op, v2*kp]/(6.283185))*VelDistFree[v2*kp, \[Theta]1, \[Phi]1, \[Gamma]]*MFPinv[op, cint, v2*kp, mx, iso];
aminus = 0;
aminus = Piecewise[{{km = \[Kappa]minus[c\[Alpha],mx, mA];
(km^4)*(Pc\[Alpha]minus[c\[Alpha],mx, iso, op, v2*km]/(6.283185))*VelDistFree[v2*km, \[Theta]1, \[Phi]1, \[Gamma]]*MFPinv[op, cint, v2*km, mx, iso],  mx > mA}, {0,  mx <= mA}}];
deff[\[Theta]2, iso]*(aplus + aminus), c\[Alpha] > Minc\[Alpha][mx, mA]},{0, c\[Alpha] < Minc\[Alpha][mx, mA]}}]
)


(* ::Text:: *)
(*fDintegrand2 is used for calculating the deflected speed distribution*)


ClearAll[fDintegrand2]
fDintegrand2[v2_, \[Theta]2_, \[Phi]2_, \[Theta]1_, \[Gamma]_, mx_, iso_, op_, cint_]:=
(
c\[Alpha] = cosine\[Alpha][\[Theta]2, \[Phi]2, \[Theta]1, 0];
mA =  0.9312*Avals[iso];

Piecewise[{{
kp = \[Kappa]plus[c\[Alpha],mx, mA];
aplus = (kp^4)*(Pc\[Alpha]plus[c\[Alpha],mx, iso, op, v2*kp]/(6.283185))*VelDistInt[v2*kp, \[Theta]1, \[Gamma]]*MFPinv[op, cint, v2*kp, mx, iso];
aminus = 0;
aminus = Piecewise[{{km = \[Kappa]minus[c\[Alpha],mx, mA];
(km^4)*(Pc\[Alpha]minus[c\[Alpha],mx, iso, op, v2*km]/(6.283185))*VelDistInt[v2*km, \[Theta]1, \[Gamma]]*MFPinv[op, cint, v2*km, mx, iso],  mx > mA}, {0,  mx <= mA}}];
deff[\[Theta]2, iso]*(aplus + aminus), c\[Alpha] > Minc\[Alpha][mx, mA]},{0, c\[Alpha] < Minc\[Alpha][mx, mA]}}]
)


(* ::Section:: *)
(*Calculating free distributions*)


(* ::Subsubsection::Closed:: *)
(*Integral over \[Phi]*)


ClearAll[Integralover\[Phi]]
Integralover\[Phi][x_, \[Phi]max_]:=
(
Piecewise[{{0, \[Phi]max  <= 0},{2\[Pi] BesselI[0, x], \[Phi]max >= \[Pi]}}, 2IntFun\[Phi][x, \[Phi]max]]
)


(* ::Subsubsection::Closed:: *)
(*SHM velocity distribution*)


VelDistFree[v_,\[Theta]_,\[Phi]_, \[Gamma]_]:=
(
c\[Delta] = Sin[\[Gamma]]*Sin[\[Theta]]*Cos[\[Phi]] + Cos[\[Gamma]]Cos[\[Theta]];
\[CapitalDelta]sq = v^2 - 2v ve c\[Delta] + ve^2; 
Piecewise[{{0, \[CapitalDelta]sq > vesc^2}, {(N1/(2\[Pi] ))*Exp[-(\[CapitalDelta]sq)/(2\[Sigma]v ^2)], \[CapitalDelta]sq <= vesc^2}}]
)


(* ::Subsubsection::Closed:: *)
(*SHM velocity distribution (integrated over \[Phi])*)


VelDistInt[v_?NumericQ,\[Theta]_?NumericQ, \[Gamma]_?NumericQ]:=
(
(*In principle, don't need to worry about dividing by zero. I should still work in the case that cosmin \[Rule] \[Infinity]*)
\[CapitalDelta]sq = v^2 - 2v ve Cos[\[Gamma]]Cos[\[Theta]] + ve^2;

If[Sin[\[Gamma]]Sin[\[Theta]] == 0, 2 \[Pi] VelDistFree[v, \[Theta], 0, \[Gamma]],

cosmin = (v^2+ve^2-vesc^2)/(2 v ve Sin[\[Gamma]]Sin[\[Theta]])-(Cos[\[Gamma]]Cos[\[Theta]])/(Sin[\[Gamma]]Sin[\[Theta]]);
(*fudge = ArcCos[Max[tester,-1]]/(\[Pi]);*)
Integralover\[Phi][Sin[\[Theta]]Sin[\[Gamma]] v ve /(\[Sigma]v^2),Re[ArcCos[Min[Max[cosmin,-1.0],1.0]]]](N1/(2\[Pi] ))*Exp[-(\[CapitalDelta]sq)/(2\[Sigma]v ^2)]]
)


(* ::Subsubsection::Closed:: *)
(*SHM speed distribution*)


SpeedDistFree[v_?NumericQ]:=
(
Piecewise[{{0,v >= vesc + ve }, 
{ \[Beta] = v*ve/\[Sigma]v^2;
(N1/\[Beta])*v^2*Exp[-(v^2 + ve^2)/(2.0\[Sigma]v^2)]*(Exp[\[Beta]] - Exp[(v^2 + ve^2 - vesc^2)/(2.0 \[Sigma]v^2)]), (v < vesc + ve)&&(v > vesc - ve)},
 { \[Beta] = v*ve/\[Sigma]v^2;(N1/\[Beta])*v^2*Exp[-(v^2 + ve^2)/(2.0\[Sigma]v^2)]*(Exp[\[Beta]] - Exp[-\[Beta]]) ,v <= vesc - ve}}]
)


(* ::Subsubsection::Closed:: *)
(*Average DM speed*)


vavg = NIntegrate[v SpeedDistFree[v], {v, 0, 753}];


(* ::Section:: *)
(*Calculating perturbed distributions*)


(* ::Subsection:: *)
(*Velocity distributions*)


(* ::Subsubsection:: *)
(*Attenuated distribution*)


(* ::Text:: *)
(*NB: We actually calculate and return fA - f0*)


VelDistAttenuation[v2_, \[Theta]2_, \[Phi]2_, \[Gamma]_, mx_, op_, cint_] :=(Exp[-Sum[deff[\[Theta]2,iso]*MFPinv[op, cint, v2*1.0, mx, iso],{iso, isovals}]] - 1)*VelDistFree[v2,\[Theta]2, \[Phi]2, \[Gamma]];


(* ::Subsubsection:: *)
(*Deflected distribution*)


(*Velocity distribution of particles which scatter once*)
VelDistDeflection[v2_, \[Theta]2_, \[Phi]2_,\[Gamma]_, mx_, op_, cint_] := NIntegrate[Sin[\[Theta]1]*Sum[fDintegrand1[v2, \[Theta]2, \[Phi]2, \[Theta]1, \[Phi]1,\[Gamma], mx, iso, op, cint], {iso, isovals}], {\[Theta]1, 0, \[Pi]}, {\[Phi]1, 0, 2\[Pi]}, Method->{"QuasiMonteCarlo", MaxPoints->20000}];


(* ::Subsection:: *)
(*Speed distributions*)


(* ::Subsubsection:: *)
(*Attenuated distribution*)


SpeedDistAttenuation[v_, \[Gamma]_, mx_, op_, cint_]:=NIntegrate[v^2*Sin[\[Theta]]*(Exp[-Sum[deff[\[Theta],iso]*MFPinv[op, cint, v*1.0, mx, iso],{iso, isovals}]] -1)*VelDistInt[v,\[Theta],\[Gamma]], {\[Theta], 0, 3.141593},   Method->{"QuasiMonteCarlo", MaxPoints->10000}];


(* ::Subsubsection:: *)
(*Deflected distribution*)


SpeedDistDeflection[v_, \[Gamma]_,mx_,op_,cint_] :=
(
Np = Round[If[mx <= 15.0, 10000, 10000+(mx-15.0)*500],10];
Sum[NIntegrate[v^2*Sin[\[Theta]2]*Sin[\[Theta]1]*fDintegrand2[v, \[Theta]2, \[Phi]2, \[Theta]1,\[Gamma], mx, iso, op, cint],{\[Theta]2, 0, 3.14159265}, {\[Theta]1, 0, 3.14159265}, {\[Phi]2, 0, 6.283185}, Method->{"QuasiMonteCarlo", MaxPoints->Np}], {iso, isovals}])


(* ::Section:: *)
(*Package post-script*)


End[]


Print[
"The EarthShadow package allows you to calculate the effects of DM-scattering in the Earth.
Send corrections/suggestions to bradkav@gmail.com.
To get a list of available functions, type: ?EarthShadow`*"]
EndPackage[]
