(* ::Package:: *)

lqc=Import["./powerspec.dat"];


bd=Import["./bdtonsr.dat"];


karray=Transpose[lqc][[2]];


atonsr=Transpose[bd][[6]];


qlqc= Transpose[lqc][[7]]+ I  Transpose[lqc][[8]];
qdlqc= Transpose[lqc][[9]]+ I  Transpose[lqc][[10]];


qbd= Transpose[bd][[2]]+ I  Transpose[bd][[3]];
qdbd= Transpose[bd][[4]]+ I  Transpose[bd][[5]];


\[Alpha] = (qdlqc Conjugate[qbd] - qlqc Conjugate[qdbd]) I atonsr^3;


\[Beta] = -(qdlqc qbd - qlqc qdbd) I atonsr^3;

(*bdratio = Abs[\[Alpha]]^2 + \[Beta]^2 + 2 Abs[\[Alpha]] \[Beta] Cos[Arg[\[Beta]]+2 Arg[qbd]]*)
bdratio = Abs[qlqc]^2/Abs[qbd]^2;

out=Transpose[N[{karray,Abs[\[Alpha]]^2+Abs[\[Beta]]^2,Abs[\[Alpha]]^2-Abs[\[Beta]]^2, Abs[\[Alpha]+\[Beta]]^2, Re[\[Alpha]], Im[\[Alpha]], Re[\[Beta]], Im[\[Beta]], bdratio, atonsr, Re[qlqc], Im[qlqc], Re[qdlqc], Im[qdlqc], Re[qbd], Im[qbd], Re[qdbd], Im[qdbd]},20]]//TableForm;


Export["./tonsrdata.dat",out]

Print["Exiting Mathematica"];
Exit[];



