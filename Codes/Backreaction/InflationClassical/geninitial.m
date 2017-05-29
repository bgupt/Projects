(* ::Package:: *)

(*SetDirectory[NotebookDirectory[]]*)


(* ::Subtitle:: *)
(*Pre-processing:*)


(* ::Text::GrayLevel[0]:: *)
(*Call the input parameter file named mathinitial.m. This file contains the initial parameter for both the background *)
(*and the perturbations.*)


<<"mathinitial.m"


(* ::Text:: *)
(*Call Numerical package. This may not be necessary on most of the machines. Does not hurt to call anyway.*)


Needs["NumericalCalculus`"]


(* ::Text:: *)
(*Set HistoryLength to zero. This is extremely important when running on HPC for large number of k's. If not set to zero, the default value is Infinity, which means that Mathematica stores the output of all the calculations even the temporaryones. This will result in very high memory requirement and the code will slow down a lot. *)


$HistoryLength=0;


(* ::Subtitle:: *)
(*Evolution equations:*)


V[x_]:=(m^2 x^2)/2;


evolEQN[t_,k_]:=Derivative[2][s][t]+(3 H[t]) Derivative[1][s][t]+((k^2) s[t])/a[t]^2;


SQ[t_,k_]:=Derivative[2][s][t]+3 Hub[t] Derivative[1][s][t]+((k^2) s[t])/a[t]^2


\[Rho]pert[t_,k_]:=(1/(2 \[Pi]))^3 (1/2 ((Derivative[1][si][t]^2+Derivative[1][sr][t]^2)+((Ut[t]+k^2) (si[t]^2+sr[t]^2))/a[t]^2)-C\[Rho]q[t,k]);
Ppert[t_,k_]:=(1/(2 \[Pi]))^3 (1/2 (Derivative[1][si][t]^2+Derivative[1][sr][t]^2-((Ut[t]+k^2) (si[t]^2+sr[t]^2))/(3 a[t]^2))-Cpq[t,k]);


\[Rho]pertS[t_,k_]:=(1/(2 \[Pi]))^3 (1/2 (Abs[Derivative[1][s][t]]^2+((Ut[t]+k^2) Abs[s[t]]^2)/a[t]^2)-C\[Rho]q[t,k]);
PpertS[t_,k_]:=(1/(2 \[Pi]))^3 (1/2 (Abs[Derivative[1][s][t]]^2-((Ut[t]+k^2) Abs[s[t]]^2)/(3 a[t]^2))-Cpq[t,k]);


(*eq1:=Derivative[2][\[Phi]][t]+(Derivative[1][v][t] Derivative[1][\[Phi]][t])/v[t]+Derivative[1][V][\[Phi][t]]==0;
eq20:=Derivative[2][v][t]==3 v[t] (-4 \[Pi] Derivative[1][\[Phi]][t]^2+8 \[Pi] (1/2 Derivative[1][\[Phi]][t]^2+V[\[Phi][t]]));
eq2[t_,k_]:=Derivative[2][v][t]==3 v[t] (-4 \[Pi] (Derivative[1][\[Phi]][t]^2+\[Rho]pert[t,k]+Ppert[t,k])+8 \[Pi] (1/2 Derivative[1][\[Phi]][t]^2+V[\[Phi][t]]+\[Rho]pert[t,k]));
eq3:=Derivative[2][v][t]==3 v[t] (-4 \[Pi] (Derivative[1][\[Phi]][t]^2+\[Rho]Q[t]+pQ[t])+8 \[Pi] (1/2 Derivative[1][\[Phi]][t]^2+V[\[Phi][t]]+\[Rho]Q[t]));
eq4:=Derivative[1][v][t]/(3 v[t])==Sqrt[8/3 \[Pi] (\[Rho]Q[t]+(1/2 Derivative[1][\[Phi]][t]^2+V[\[Phi][t]]))];*)


U[t_]:=0


(* ::Subtitle:: *)
(*Definitions of initial perturbations states for ANA vacuum:*)


(* ::Text:: *)
(*Define the adiabatic subtraction terms defined in eq. (A3) and (A4) of arXiv:1412.3524. Note that these definitions are in conformal time. We will convert them to the cosmic time later. After obtaining the full expression for the mode functions that is.*)


C\[Rho]q[t_, k_] := k/(2*a[t]^4) + ((4*a[t]^2*U[t]*k^2 + 4*k^2*Derivative[1][a][t]^2))/
        (16*a[t]^6*k^3) + (1/(256*a[t]^7*f[t]^7))*
      ((-16*a[t]^3*f[t]^4*U[t]^2 - 32*a[t]*f[t]^4*U[t]*Derivative[1][a][t]^2 - 
        96*a[t]^2*f[t]^3*U[t]*Derivative[1][a][t]*Derivative[1][f][t] - 
        40*a[t]^3*f[t]^2*U[t]*Derivative[1][f][t]^2 - 24*a[t]*f[t]^2*
         Derivative[1][a][t]^2*Derivative[1][f][t]^2 - 120*a[t]^2*f[t]*
         Derivative[1][a][t]*Derivative[1][f][t]^3 - 45*a[t]^3*Derivative[1][f][t]^4 + 
        32*a[t]^2*f[t]^4*Derivative[1][a][t]*Derivative[1][U][t] + 
        16*a[t]^3*f[t]^3*Derivative[1][f][t]*Derivative[1][U][t] + 
        64*f[t]^4*Derivative[1][a][t]^2*Derivative[2][a][t] + 
        112*a[t]*f[t]^3*Derivative[1][a][t]*Derivative[1][f][t]*Derivative[2][a][t] + 
        16*a[t]^2*f[t]^2*Derivative[1][f][t]^2*Derivative[2][a][t] + 
        16*a[t]*f[t]^4*Derivative[2][a][t]^2 + 16*a[t]*f[t]^3*Derivative[1][a][t]^2*
         Derivative[2][f][t] + 112*a[t]^2*f[t]^2*Derivative[1][a][t]*
         Derivative[1][f][t]*Derivative[2][f][t] + 40*a[t]^3*f[t]*Derivative[1][f][t]^2*
         Derivative[2][f][t] + 16*a[t]^2*f[t]^3*Derivative[2][a][t]*
         Derivative[2][f][t] + 4*a[t]^3*f[t]^2*Derivative[2][f][t]^2 - 
        32*a[t]*f[t]^4*Derivative[1][a][t]*Derivative[3][a][t] - 
        16*a[t]^2*f[t]^3*Derivative[1][f][t]*Derivative[3][a][t] - 
        16*a[t]^2*f[t]^3*Derivative[1][a][t]*Derivative[3][f][t] - 
        8*a[t]^3*f[t]^2*Derivative[1][f][t]*Derivative[3][f][t])) /. 
    {f[t] -> Sqrt[k^2 + a[t]^2*\[Mu]^2], 
       Derivative[1][f][t] -> D[Sqrt[k^2 + a[x]^2*\[Mu]^2], {x, 1}]/.x->t, 
     Derivative[2][f][t] -> D[Sqrt[k^2 + a[x]^2*\[Mu]^2], {x, 2}]/.x->t, 
     Derivative[3][f][t] -> D[Sqrt[k^2 + a[x]^2*\[Mu]^2], {x, 3}]/.x->t, 
     Derivative[4][f][t] -> D[Sqrt[k^2 + a[x]^2*\[Mu]^2], {x, 4}]/.x->t}; 


Cpq[t_, k_] := k/(6*a[t]^4) + ((a[t]^2*U[t] + 3*Derivative[1][a][t]^2 - 
        2*a[t]*Derivative[2][a][t])/(12*k*a[t]^6)) + 
     (1/(768*a[t]^7*f[t]^9))*(-48*\[Mu]^2*a[t]^5*f[t]^4*U[t]^2 - 
       16*a[t]^3*f[t]^6*U[t]^2 - 96*a[t]*f[t]^6*U[t]*Derivative[1][a][t]^2 - 
       288*a[t]^2*f[t]^5*U[t]*Derivative[1][a][t]*Derivative[1][f][t] - 
       200*\[Mu]^2*a[t]^5*f[t]^2*U[t]*Derivative[1][f][t]^2 - 424*a[t]^3*f[t]^4*U[t]*
        Derivative[1][f][t]^2 - 72*a[t]*f[t]^4*Derivative[1][a][t]^2*
        Derivative[1][f][t]^2 - 360*a[t]^2*f[t]^3*Derivative[1][a][t]*
        Derivative[1][f][t]^3 - 315*\[Mu]^2*a[t]^5*Derivative[1][f][t]^4 - 
       765*a[t]^3*f[t]^2*Derivative[1][f][t]^4 + 96*a[t]^2*f[t]^6*Derivative[1][a][t]*
        Derivative[1][U][t] + 80*\[Mu]^2*a[t]^5*f[t]^3*Derivative[1][f][t]*
        Derivative[1][U][t] + 208*a[t]^3*f[t]^5*Derivative[1][f][t]*
        Derivative[1][U][t] + 96*\[Mu]^2*a[t]^4*f[t]^4*U[t]*Derivative[2][a][t] + 
       64*a[t]^2*f[t]^6*U[t]*Derivative[2][a][t] + 32*\[Mu]^2*a[t]^2*f[t]^4*
        Derivative[1][a][t]^2*Derivative[2][a][t] + 256*f[t]^6*Derivative[1][a][t]^2*
        Derivative[2][a][t] + 80*\[Mu]^2*a[t]^3*f[t]^3*Derivative[1][a][t]*
        Derivative[1][f][t]*Derivative[2][a][t] + 496*a[t]*f[t]^5*Derivative[1][a][t]*
        Derivative[1][f][t]*Derivative[2][a][t] + 200*\[Mu]^2*a[t]^4*f[t]^2*
        Derivative[1][f][t]^2*Derivative[2][a][t] + 448*a[t]^2*f[t]^4*
        Derivative[1][f][t]^2*Derivative[2][a][t] - 64*\[Mu]^2*a[t]^3*f[t]^4*
        Derivative[2][a][t]^2 - 80*a[t]*f[t]^6*Derivative[2][a][t]^2 + 
       80*\[Mu]^2*a[t]^5*f[t]^3*U[t]*Derivative[2][f][t] + 96*a[t]^3*f[t]^5*U[t]*
        Derivative[2][f][t] + 48*a[t]*f[t]^5*Derivative[1][a][t]^2*
        Derivative[2][f][t] + 336*a[t]^2*f[t]^4*Derivative[1][a][t]*Derivative[1][f][t]*
        Derivative[2][f][t] + 420*\[Mu]^2*a[t]^5*f[t]*Derivative[1][f][t]^2*
        Derivative[2][f][t] + 960*a[t]^3*f[t]^3*Derivative[1][f][t]^2*
        Derivative[2][f][t] - 80*\[Mu]^2*a[t]^4*f[t]^3*Derivative[2][a][t]*
        Derivative[2][f][t] - 112*a[t]^2*f[t]^5*Derivative[2][a][t]*
        Derivative[2][f][t] - 60*\[Mu]^2*a[t]^5*f[t]^2*Derivative[2][f][t]^2 - 
       108*a[t]^3*f[t]^4*Derivative[2][f][t]^2 - 16*\[Mu]^2*a[t]^5*f[t]^4*
        Derivative[2][U][t] - 32*a[t]^3*f[t]^6*Derivative[2][U][t] - 
       32*\[Mu]^2*a[t]^3*f[t]^4*Derivative[1][a][t]*Derivative[3][a][t] - 
       160*a[t]*f[t]^6*Derivative[1][a][t]*Derivative[3][a][t] - 
       80*\[Mu]^2*a[t]^4*f[t]^3*Derivative[1][f][t]*Derivative[3][a][t] - 
       208*a[t]^2*f[t]^5*Derivative[1][f][t]*Derivative[3][a][t] - 
       48*a[t]^2*f[t]^5*Derivative[1][a][t]*Derivative[3][f][t] - 
       80*\[Mu]^2*a[t]^5*f[t]^2*Derivative[1][f][t]*Derivative[3][f][t] - 
       184*a[t]^3*f[t]^4*Derivative[1][f][t]*Derivative[3][f][t] + 
       16*\[Mu]^2*a[t]^4*f[t]^4*Derivative[4][a][t] + 32*a[t]^2*f[t]^6*
        Derivative[4][a][t] + 8*\[Mu]^2*a[t]^5*f[t]^3*Derivative[4][f][t] + 
       16*a[t]^3*f[t]^5*Derivative[4][f][t]) /. {f[t] -> Sqrt[k^2 + a[t]^2*\[Mu]^2], 
     Derivative[1][f][t] -> D[Sqrt[k^2 + a[x]^2*\[Mu]^2], {x, 1}]/.x->t, 
     Derivative[2][f][t] -> D[Sqrt[k^2 + a[x]^2*\[Mu]^2], {x, 2}]/.x->t, 
     Derivative[3][f][t] -> D[Sqrt[k^2 + a[x]^2*\[Mu]^2], {x, 3}]/.x->t, 
     Derivative[4][f][t] -> D[Sqrt[k^2 + a[x]^2*\[Mu]^2], {x, 4}]/.x->t}; 


(* ::Text:: *)
(*The adiabatic frequence \[CapitalOmega]k defined in (3.2) of arXiv : 1412.3524 is called Wq here. The following is obtained by substituting m=0 in (3.2), as hte cosmological perturbations are treated as massless scalar fields.*)


Wq[t_, k_] :=-((k^2 + U[t])/(3*a[t]^4*(Cpq[t, k] - C\[Rho]q[t, k]))); 


(* ::Text:: *)
(*Define Vk of eq. (3.3) of arXiv : 1412.3524. Vq1 (Vq2) is used when the initial conditions are provided in the expanding (contracting) branch.*)


Vq1[t_, k_] := -2*(Sqrt[-(k^2 + U[t] - 4*a[t]^4*C\[Rho]q[t, k]*Wq[t,k] + Wq[t,k]^2)] - 
     (Derivative[1][a][t]/a[t])); 
Vq2[t_, k_] := 2*(Sqrt[-(k^2 + U[t] - 4*a[t]^4*C\[Rho]q[t, k]*Wq[t,k] + Wq[t,k]^2)] + 
     (Derivative[1][a][t]/a[t])); 


(* ::Text:: *)
(*In all the simulations t=0 is chosen to be the bounce time. So tmin>0 (<0) means we are giving the initial conditions for the perturbations in the expanding (contracting) branch .*)


Vq[t_,k_] := If[tmin>0, Vq1[t,k], Vq2[t,k]];


(* ::Text:: *)
(*Having the adiabatic frequency and Vk, we define the mode function and its time derivative as given in eq. (3.1) of arXiv :1412.3524.*)


qvac[t_, k_] := 1/(a[t] Sqrt[2*Wq[t,k]]); 
qdvac[t_, k_] := ((-I)*Wq[t,k] +  Vq[t,k]/2 - Derivative[1][a][t]/a[t])*qvac[t, k]; 
qdvacconj[t_, k_] := (I*Wq[t,k] + Vq[t,k]/2 - Derivative[1][a][t]/a[t])*qvac[t, k]


Print["\[Psi](eta), \[CapitalOmega]4(eta) and Vq of eq(3.1) of arXiv:1412.3524 defined. So far all the
 definition are in conformal time."];


(* ::Text:: *)
(*Note that all the definitions are still in conformal time. We now convert them to cosmic time as follow. \[Psi]0 is the mode functions and \[Psi]0p is its tiem derivative.*)


\[Psi]0[t_, k_]= qvac[\[Tau], k] /. {a[\[Tau]] -> a[t], Derivative[1][a][\[Tau]] -> 
      Derivative[1][a][t]*a[t], Derivative[2][a][\[Tau]] -> 
      a[t]*(Derivative[1][a][t]^2 + a[t]*Derivative[2][a][t]), 
     Derivative[3][a][\[Tau]] -> a[t]*(Derivative[1][a][t]^3 + 4*a[t]*Derivative[1][a][t]*
         Derivative[2][a][t] + a[t]^2*Derivative[3][a][t]), 
     Derivative[4][a][\[Tau]] -> a[t]*(Derivative[1][a][t]^4 + 
       12*a[t]*Derivative[1][a][t]^2*Derivative[2][a][t] + 
       4*a[t]^2*Derivative[2][a][t]^2 + 7*a[t]^2*Derivative[1][a][t]*
       Derivative[3][a][t] + a[t]^3*Derivative[4][a][t]),f[\[Tau]]->f[t], 
       U[\[Tau]]->Ut[t], 
       Derivative[1][U][\[Tau]]->a[t] Derivative[1][Ut][t],
       Derivative[2][U][\[Tau]]-> a[t] (Derivative[1][a][t] Derivative[1][Ut][t]+a[t] Derivative[2][Ut][t])(*,
       g[\[Tau]] -> G[t],  Derivative[1][g][\[Tau]] -> Derivative[1][G][t]*a[t], Derivative[2][g][\[Tau]] -> 
      a[t]*(Derivative[1][a][t]*Derivative[1][G][t] + a[t]*Derivative[2][G][t])*), 
     \[Delta] -> 1}; 


\[Psi]0p[t_, k_] = qdvac[\[Tau], k]/a[\[Tau]] /. {(*\[Phi][\[Tau]]->\[Phi]t[t],*)a[\[Tau]] -> a[t], Derivative[1][a][\[Tau]] -> 
      Derivative[1][a][t]*a[t], Derivative[2][a][\[Tau]] -> 
      a[t]*(Derivative[1][a][t]^2 + a[t]*Derivative[2][a][t]), 
     Derivative[3][a][\[Tau]] -> a[t]*(Derivative[1][a][t]^3 + 4*a[t]*Derivative[1][a][t]*
         Derivative[2][a][t] + a[t]^2*Derivative[3][a][t]), 
     Derivative[4][a][\[Tau]] -> a[t]*(Derivative[1][a][t]^4 + 
        12*a[t]*Derivative[1][a][t]^2*Derivative[2][a][t] + 
        4*a[t]^2*Derivative[2][a][t]^2 + 7*a[t]^2*Derivative[1][a][t]*
         Derivative[3][a][t] + a[t]^3*Derivative[4][a][t]), 
     Derivative[5][a][\[Tau]] -> a[t]*(Derivative[1][a][t]^5 + 
        28*a[t]*Derivative[1][a][t]^3*Derivative[2][a][t] + 
        36*a[t]^2*Derivative[1][a][t]*Derivative[2][a][t]^2 + 
        33*a[t]^2*Derivative[1][a][t]^2*Derivative[3][a][t] + 
        15*a[t]^3*Derivative[2][a][t]*Derivative[3][a][t] + 
        11*a[t]^3*Derivative[1][a][t]*Derivative[4][a][t] + 
        a[t]^4*Derivative[5][a][t]),f[\[Tau]]->f[t], U[\[Tau]]->Ut[t], 
       Derivative[1][U][\[Tau]]->a[t] Derivative[1][Ut][t],
       Derivative[2][U][\[Tau]]-> a[t] (Derivative[1][a][t] Derivative[1][Ut][t]+a[t] Derivative[2][Ut][t])(*, g[\[Tau]] -> G[t], 
       Derivative[1][g][\[Tau]] -> Derivative[1][G][t]*a[t], 
       Derivative[2][g][\[Tau]] -> a[t]*(Derivative[1][a][t]*Derivative[1][G][t] + a[t]*Derivative[2][G][t]), 
     Derivative[3][g][\[Tau]] -> a[t]*(Derivative[1][a][t]^2*Derivative[1][G][t] + 
        a[t]*Derivative[2][a][t]*Derivative[1][G][t] + 3*a[t]*Derivative[1][a][t]*
         Derivative[2][G][t] + a[t]^2*Derivative[3][G][t])*), \[Delta] -> 1}; 


Print["Initial conditions \[Psi]0 and \[Psi]0dot now defined in cosmic time."];


(* ::Text:: *)
(*rt[t] is the fractional KE, and Ut[t] is the effective potential felt by the scalar perturbations. Note that rt and Ut only depend on the background quantities.*)


(*rt[t_] := (3*8*Pi*(Derivative[1][\[Phi]][t] /. q)^2)/\[Rho][t]; 
Ut[t_] := If[perttype == scalar,((1/2)*M2*(\[Phi][t]/. q)^2*rt[t] - 2*M2*(\[Phi][t] /. q)*Sqrt[rt[t]] + M2)*a[t]^2, 0]; *)


(*R[t_] := (3*8*Pi*(Derivative[1][\[Phi]][t] /. q)^2)/\[Rho][t]; 
G[t_] := If[perttype == scalar,((1/2)*M2*(\[Phi][t]/. q)^2*rt[t] - 2*M2*(\[Phi][t] /. q)*Sqrt[rt[t]] + M2)*a[t]^2, 0];  
(*Print["RhoQ and pQ defined!"];*)*)


G[t_]:=0;
Ut[t_]:=0;
U[t_]:=0;


(* ::Subtitle:: *)
(*Background Evolution*)


(* ::Text:: *)
(*Now evaluate the background for the provided initial conditions.*)


(* ::Text:: *)
(*Initialize background quantities:*)


\[Rho]ratio = Rationalize[rhoratio,0]; 
M2 = Rationalize[mass^2,0]; 
\[Gamma] = 2375/10000; 
\[Lambda] = Sqrt[4*Sqrt[3]*Pi*\[Gamma]];
\[Rho]crit = 3/(8*Pi*\[Gamma]^2*\[Lambda]^2);
\[Rho]c = \[Rho]crit/\[Rho]ratio; 
Print["rho_crit = ", N[\[Rho]crit]]; 
rho[t_] := (1/2)*Derivative[1][\[Phi]][t]^2 + (1/2)*M2*\[Phi][t]^2; 
press[t_] := (1/2)*Derivative[1][\[Phi]][t]^2 + (1/2)*M2*\[Phi][t]^2; 
\[Phi]0 = Rationalize[phi0,0]; 


vppeq = Derivative[2][v][t] == -4*Pi*rho[t]*(1 - 4*(rho[t]/\[Rho]c)) + 
     16*Pi*rho[t]*(1 - rho[t]/\[Rho]c) - 12*Pi*press[t]*(1 - 2*(rho[t]/\[Rho]c)); 
d = NDSolve[{Derivative[2][\[Phi]][t] + (Derivative[1][v][t]*Derivative[1][\[Phi]][t])/v[t] + 
        M2*\[Phi][t] == 0, vppeq, \[Phi][0] == \[Phi]0, Derivative[1][\[Phi]][0] == 
       Sqrt[2*(\[Rho]c - (1/2)*M2*\[Phi]0^2)], v[0] == 1, Derivative[1][v][0] == 0}, 
     {\[Phi], v}, {t, 0, 1}, Method -> "StiffnessSwitching", AccuracyGoal -> 35, 
     PrecisionGoal -> 35, WorkingPrecision -> 90, MaxSteps -> Infinity][[1]]; 
vpp = Rationalize[Derivative[2][v][0] /. d, 0]; 
Print["Double derivative of v at t=0 is", N[vpp]]; 


Print["Calculating background ..."]; 


Print["Mass is ", N[M2]]; 
q = Flatten[NDSolve[{Derivative[2][\[Phi]][t] + (Derivative[1][v][t]*
          Derivative[1][\[Phi]][t])/v[t] + M2*\[Phi][t] == 0, Derivative[3][v][t] == 
       (Derivative[1][v][t]*(24*Pi*(((1/2)*Derivative[1][\[Phi]][t]^2)^2 + 
            (1/2)*M2*\[Phi][t]^2*(\[Rho]c - (1/2)*M2*\[Phi][t]^2))))/\[Rho]c + 
        (24*Pi*v[t]*(Derivative[1][\[Phi]][t]^3*(-((Derivative[1][v][t]*Derivative[1][\[Phi]][
                 t])/v[t]) - M2*\[Phi][t]^2) + M2*\[Phi][t]*Derivative[1][\[Phi]][t]*
            (\[Rho]c - (1/2)*M2*\[Phi][t]^2) - (1/2)*M2^2*\[Phi][t]^3*Derivative[1][\[Phi]][t]))/
         \[Rho]c, \[Phi][0] == \[Phi]0, Derivative[1][\[Phi]][0] == 
       Sqrt[2*(\[Rho]c - (1/2)*M2*\[Phi][0]^2)], v[0] == 1, Derivative[1][v][0] == 0, 
      Derivative[2][v][0] == vpp}, {\[Phi], v}, {t, -2000, 6 10^6}, 
     Method -> "StiffnessSwitching", AccuracyGoal -> 35, PrecisionGoal -> 35, 
     WorkingPrecision -> 90, MaxSteps -> Infinity]]; 

Print["Background calculated!"];

Print["a, H, \[Epsilon], z, \[Rho] and Pressure defined in terms of the background solution!"];


(* ::Text:: *)
(*Define background quantities in terms of the numerical solution.*)


a[t_] := v[t]^(1/3) /. q;  (*Scale factor*)
\[Phi]t[t_]:=\[Phi][t]/.q;     (*Scalar field*)
p\[Phi][t_]:= Derivative[1][\[Phi]][t] v[t]/.q; (*Conjugate momentum to the scalar field*)
H[t_] := Derivative[1][v][t]/(3*v[t]) /. q; (*Hubble rate*)
\[Epsilon][t_] := -(Derivative[1][H][t]/H[t]^2) /. q; (*First slow-roll parameter*)
z[t_] := (a[t]*(Derivative[1][\[Phi]][t] /. q))/H[t]; (*Red shift*)
\[Eta][t_] := Derivative[2][H][t]/(2*Derivative[1][H][t]*H[t]); (*Second slow-roll parameter*)
\[Rho][t_] := (1/2)*M2*(\[Phi][t] /. q)^2 + (1/2)*(Derivative[1][\[Phi]][t] /. q)^2; (*Background energy density*)
P[t_] := (1/2)*(Derivative[1][\[Phi]][t] /. q)^2 - (1/2)*M2*(\[Phi][t] /. q)^2; (**)


(* ::Text:: *)
(*Find the time for the onset of slow - roll.*)


(*tonsrNear=150000;


tonsr=Rationalize[t/.FindRoot[H[t]-Hs,{t,tonsrNear}],0];
*)

(* ::Section:: *)
(*Solving the classical EOMs*)


(* ::Subtitle:: *)
(*Initial conditions for perturbations:*)


(* ::Text:: *)
(*Compute and output the initial conditions for the mode functions.*)


(*numk=IntegerPart[(kmax-kmin)/(1)+1];
*)
stepk=Rationalize[(kmax-kmin)/(numk-1),0];
karray=Rationalize[Table[kmin+(kmax-kmin)/(numk-1) (i-1) ,{i,1,numk}],0];
\[Mu]=\[Mu]var;
(*tini = tmin;*)
tini = Rationalize[tmin,0];
tfinal = Rationalize[tmax,0]; 
len = Length[karray]; 
sinitmp = {}; 
dsinitmp = {}; 
For[i = 1, i < len + 1, i++, 
   If[ inidatatype == 1,
   stmp = \[Psi]0[tini, karray[[i]]];
   dstmp = \[Psi]0p[tini, karray[[i]]],
   stmp = 1/Sqrt[2*karray[[i]]]/a[tini];
   dstmp = -I*Sqrt[karray[[i]]/2]/a[tini]^2-(a'[tini]/a[tini]^2) stmp
   ]
   AppendTo[sinitmp, stmp]; 
   AppendTo[dsinitmp, dstmp]; 
];


(*\[Beta]0 = Rationalize[1.0,0];*)
\[Beta]array = \[Beta]0 Exp[- (karray/kbeta)^16];
\[Alpha]array = Sqrt[1+Abs[\[Beta]array]^2];


sini=Rationalize[\[Alpha]array sinitmp + Conjugate[\[Beta]array sinitmp],0];
dsini=Rationalize[\[Alpha]array dsinitmp + Conjugate[\[Beta]array dsinitmp],0];

(*sini = Rationalize[sinitmp, 0];
dsini = Rationalize[dsinitmp, 0];*)


(* ::Subsection:: *)
(*Initial conditions for the background*)


\[Phi]i = \[Phi]t[tini];
\[Phi]di = \[Phi]t'[tini];
vi = a[tini]^3;
p\[Phi]i= vi \[Phi]t'[tini];


out1=Transpose[{N[karray,50], N[Re[sini],50], N[Im[sini],50], N[Re[dsini],50], N[Im[dsini],50]}]//TableForm;


Export["./initialdata.dat",out1];


initial={{"&input"},{"mass=",N[mass,50]},
{"vol0=",N[a[tini]^3,50]},
{"H0=",N[H[tini],50]},
{"phi0=",N[\[Phi]i,50]},
{"rhoratio=",1},
{"rhofactor=",N[1/(5 10^4), 50]},
{"maxiter=",maxiter},
{"pphi0=",p\[Phi]i},
{"pert_type=",perttype},
{"potential_type=",potentialtype},
{"tmax=",N[tmax,50]},
{"tmin=",N[tini,50]},
{"numk=",numk},
{"dk=",N[karray[[2]]-karray[[1]],50]},
{"outevery=",outevery},
{"outbgrndevery=",outbgrndevery},
{"nsteps=",nsteps},
{"direc=",N[-1,50]},
{"inidatatype=",inidatatype},
{"mu=","1.0d0"},
{"beta0=",N[\[Beta]0,50]},
{"kmin=",kmin},
{"kmax=",kmax},
{"/"},
{"NOTE: Please do not modify this file. If you wish to change the initial condition, 
   please do so by modifying 'mathinitial.m' which is a mathematica script to provide 
   the initial data."}};


Export["input.dat",initial//TableForm];

