(* ::Package:: *)

(* ::Section:: *)
(*Package Header*)


BeginPackage["Taggar`GERF`"];


(* ::Text:: *)
(*Declare your public symbols here:*)


GERFSolve::usage="GERFSolve[lhs == rhs, u[vars]] solves the NLPDE given by lhs==rhs using GERF method.";
GERFSolve::arglen="Exponents and coefficients lists must have the same length.";
GERFSolve::wavelen="WaveConstants list length must match the number of independent variables.";
\[ScriptCapitalC];
\[ScriptCapitalA];


Begin["`Private`"];


(* ::Section:: *)
(*Definitions*)


(* ::Text:: *)
(*Define your public and private symbols here:*)


ClearAll[GERFSolve];
Options[GERFSolve] =
	{"BalanceConstant" -> 1,
	"AnsatzExponents" -> {I, -I, I, -I},
	"AnsatzCoefficients" -> {1, 1, 1, -1},
	"WaveConstants" -> Automatic,
	"SolveFunction" -> Solve,
	"OutputMode" -> "Solutions"};
GERFSolve[eqn_, u_[vars__], opts : OptionsPattern[]]:=
	Module[
		
		{waveConsts, c, a, Y, eta, U, etaExpr, rTrial,
		uTrial, numExpanded, yTerms, system, sol, formattingRule,
		finalFuncs},
		
		(* coeffs and exps in ansatz to be the same in number,
		validating that: *)
		If[Length[OptionValue["AnsatzExponents"]] != Length[OptionValue["AnsatzCoefficients"]],
			Message[GERFSolve::arglen]; Return[$Failed]];
			
		(* provision for custom wave constants: *)
		waveConsts = OptionValue["WaveConstants"];
		If[waveConsts === Automatic,
			waveConsts = Array[c, Length[{vars}]],
			If[Length[waveConsts] != Length[{vars}],
				Message[GERFSolve::wavelen]; Return[$Failed]]];
		
		(* wave transformation: *)
		etaExpr = Sum[waveConsts[[i]] * {vars}[[i]], {i, Length[{vars}]}];
		
		(* trial solution: *)
		rTrial[eta_] := Sum[OptionValue["AnsatzCoefficients"][[i]] * Exp[OptionValue["AnsatzExponents"][[i]] * eta],
		{i, Length[OptionValue["AnsatzExponents"]]/2}] /
			Sum[OptionValue["AnsatzCoefficients"][[i]] * Exp[OptionValue["AnsatzExponents"][[i]] * eta],
				{i, Length[OptionValue["AnsatzExponents"]]/2+1, Length[OptionValue["AnsatzExponents"]]}];
		uTrial[eta_] := Sum[a[k] * rTrial[eta]^k, {k, -OptionValue["BalanceConstant"], OptionValue["BalanceConstant"]}];
		
		(* evaluate ODE, log-transform, and extract numerator in one sweep: *)
		numExpanded = Expand[Numerator[Together[
			If[Head[#] === Equal, #[[1]] - #[[2]], #]& [eqn /. u -> Function[Evaluate[{vars}], Evaluate[U[etaExpr]]]] /. 
				etaExpr -> eta /. U -> Function[{eta}, Evaluate[uTrial[eta]]] /. eta -> Log[Y]]]];
		
		(* system formulation: *)
		yTerms = DeleteDuplicates[Cases[numExpanded, Y^p_ :> Y^p, {0, Infinity}]];
		If[!FreeQ[numExpanded, Y, {0,0}] &&! MemberQ[yTerms, Y],
			AppendTo[yTerms, Y]];
		system = Map[Coefficient[numExpanded, #] == 0 &, yTerms];
		If[(numExpanded /. (Y | Y^_ -> 0)) =!= 0, 
			AppendTo[system, (numExpanded /. (Y | Y^_ -> 0)) == 0]];
			
		(* solve and format: *)
		sol = OptionValue["SolveFunction"][system];
		sol = Quiet[Select[sol, FreeQ[uTrial[etaExpr] /. #,
			ComplexInfinity | Indeterminate | DirectedInfinity] &]];
		sol = DeleteDuplicatesBy[sol, Expand[uTrial[etaExpr] /. #] &];
		formattingRule = {c[i_] :> Subscript[\[ScriptCapitalC], i], a[i_] :> Subscript[\[ScriptCapitalA], i]};
		finalFuncs = Map[(u @@ {vars}) -> (uTrial[etaExpr] /. #) &, sol] /. formattingRule;
		
		Switch[OptionValue["OutputMode"],
			"SolutionSets", sol /. formattingRule,
			"Solutions", finalFuncs,
			"Both", MapThread[{#1, #2} &, {sol /. formattingRule, finalFuncs}],
			_, sol /. formattingRule]]


(* ::Section::Closed:: *)
(*Package Footer*)


End[];
EndPackage[];
