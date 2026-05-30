(* ::Package:: *)

(* ::Section:: *)
(*Package Header*)


BeginPackage["Taggar`GERF`"];


(* ::Text:: *)
(*Declare your public symbols here:*)


GERFSolve::usage="GERFSolve[lhs == rhs, u[vars]] solves the NLPDE given by lhs==rhs using GERF method.";


Begin["`Private`"];


(* ::Section:: *)
(*Definitions*)


(* ::Text:: *)
(*Define your public and private symbols here:*)


(* ::Text:: *)
(*Helper function to calculate the balance constant of given equation.*)


BalanceConstant[eqn_, u_[vars__]] :=
	Module[
		{expr, terms, ld, nld, m, bal},
		
		expr = Expand @ If[Head[eqn] === Equal, Subtract @@ eqn, eqn];
		terms = If[Head[expr] === Plus, List @@ expr, {expr}];
		
		(* extract degrees of linear and nonlinear terms *)
		ld = Simplify[GetDegreeofTerm[#, u[vars], m] & /@
			Select[terms, LinearQ[#, u[vars]] &]];
		nld = Simplify[GetDegreeofTerm[#, u[vars], m] & /@
			Select[terms, ! LinearQ[#, u[vars]] &]];
		
		(* if either is zero, error out- *)
		If[Length[ld] == 0 || Length[nld] == 0,
			Message[GERFSolve::GERFPackageError, 
				"Balance constant could not be calculated. Please provide one using the 'BalanceConstant' option. Probable reason: equation does not contain both linear and nonlinear terms."];
			Return[$Failed]];
		
		(* solve and appropriately return the balance constant *)
		bal = Solve[
			Last[SortBy[ld, # /. m -> 100 &]] ==
				Last[SortBy[nld, # /. m -> 100 &]],
			m];
		
		If[Length[bal] > 0,
			bal = m /. bal[[1]];
			If[bal > 0, bal,
				Message[GERFSolve::GERFPackageError,
					"Extracted balance constant is non-positive or invalid."];
				$Failed],
			(* if everything fails, user may provide a balance constant as input *)
			Message[GERFSolve::GERFPackageError,
				"Balance constant could not be calculated. Please consider providing one using the 'BalanceConstant' option"];
				$Failed]]


(* ::Text:: *)
(*Helper function to check if the given equation is an NLPDE or not, and some other related functions.*)


NonlinearQ[eqn_, u_[vars__]] :=
	Module[
		{expr, pat},
		
		(* has derivative terms or not: *)
		If[Not[Length[{vars}] > 1 &&
			Not @ FreeQ[eqn, Derivative | FractionalD | CaputoD]],
				Throw @ False];
		(* lhs == rhs -> lhs-rhs == 0 *)
		expr = Expand @ If[Head[eqn] === Equal, Subtract @@ eqn, eqn];
		(* testing the pattern *) pat = u[vars] | Derivative[__][u][vars] |
			FractionalD[__] | CaputoD[__];
		(* against expr *)
		Throw @ Not @ FreeQ[expr, Times[a_, b_] /;
			MatchQ[a, pat] && MatchQ[b, pat]]]


FractionalPDEQ[eqn_] := Not @ FreeQ[eqn, FractionalD | CaputoD]


FractionalOrders[eqn_] :=
	Return @ DeleteDuplicates @
		Cases[eqn, (FractionalD | CaputoD)[_, {_, alpha_}] :> alpha, Infinity]


(* ::Text:: *)
(*For extracting balance constants:*)


GetDegreeofTerm[term_, u_[vars__], m_] :=
	Module[
		{k, subbed},
		
		subbed = term /. {
			u[vars] -> Power[k, m],
			Derivative[orders__][u][vars] :> Power[k, m + Total[{orders}]]};
		
		Simplify[Exponent[subbed, k]]]


LinearQ[term_, u_[vars__]] :=
	Module[
		{k, subbed},subbed = term /. {
			u[vars] -> k,
			Derivative[__][u][vars] :> k,
			(FractionalD | CaputoD)[__] :> k};
		Simplify[Exponent[subbed, k]] == 1]


(* ::Text:: *)
(*Reducing the PDE to ODE:*)


ReducetoODE[eqn_, u_[vars__], U_, eta_, c_] :=
	Module[
		{dr, waveConstants},
		
		(* default wave constants *)
		waveConstants = Table[c[k], {k, {vars}}];
		
		(* applying the derivative rule given as: *)
			dr = {(FractionalD | CaputoD)[u[vars], {var_, alpha_}] :>
				c[var] * Derivative[1][U][eta],
				(* only one FractionalD at a time please *)
				(* D^alpha u_x = x^(1-alpha) U'(eta) deta/dx
					= x^(1-alpha) U'(eta) a x^(alpha-1)
					= a U'(eta) *)
				Derivative[orders__][u][vars] :>
					Times @@ Power[waveConstants, {orders}] *
						Derivative[Total[{orders}]][U][eta]};
				(* simply using this rule, for example u_x converts into 
					a*U', or u_xxy -> a^2bU''', and so on *)
		(* to the equation and return *)
			Return[eqn /. dr /. u[vars] :> U[eta]]]


(* ::Text:: *)
(*Constructing the ansatz:*)


TrialSolution[m_, eta_,
	{
		A_ (* coefficient in trial solution *),
		R_ (* method dependent ansatz *)
	}] :=
		Sum[A[k] Power[R[eta], k], {k, -m, m}]


(* ::Text:: *)
(*Algebraic system solver*)


SolveAlgebraicSystem[
	{ts_, poly_, R_[eta_], eqn_, vars_},
	{R0_, {A_, c_, \[Lambda]_, \[Mu]_}},
	{Solver_, out_}] :=
	Module[
		{expr, u, v, sys, sol, pairs, valid, dummy, rad},
		
		(* symbolic R must be replaced with the actual corresponding
			form before any further evaluation to ensure proper
			solutions appear *)
		sys = Thread[CoefficientList[
					Expand @ Numerator @ Together @ TrigToExp @ First @
						(poly /. R -> R0) /.
						Exp[d_. * eta] :> u^Re[d] v^Im[d], {u, v}] == 0];
				
		sol = Select[
			Solver[sys, Table[A[i], {i, -3, 3}]],
			Length @ DeleteDuplicates @ #[[All,2]] > 1 &];
		
		pairs =
			DeleteDuplicates @ Table[
				With[{f = ts /. R -> R0 /. eta -> 
					FormWaveTransformation[eqn, vars, c] /. v},
					{v, u@@vars -> f}], {v, sol}];
		
		valid = Select[pairs,
			FreeQ[Last[#], Indeterminate | ComplexInfinity | Infinity | DirectedInfinity] &];

		CleanSymbols[
			Switch[out,
				"SolutionSets", First /@ valid,
				All, Transpose[{First /@ valid, Last /@ valid}],
				_, Last /@ valid],
		  A, c]]


(* ::Text:: *)
(*Make solution forms*)


FormWaveTransformation[eqn_, vars_, c_] :=
	Module[
		{raw, fds, defaulters},
		
		(* all cases of fractional orders: *)
		raw = DeleteDuplicates @ Cases[eqn,
			(CaputoD | FractionalD)[_, {x_, alpha_}] :> (x -> alpha), 
			Infinity];
		
		(* if at all, fractional orders must appear uniquely for
			each dimension. hence, the following are troublesome: *)
		defaulters = Select[
			GroupBy[raw, First -> Last, DeleteDuplicates],
				Length @ # > 1 &];
		If[Length @ defaulters > 1,
			Message[GERFSolve::InvalidFractionalDerivatives, defaulters];
			Throw @ $Failed]; (* throw to top level *)
		
		(* finally, when all is well, we now formulate the transformation *)
		fds = Association[raw];
		Sum[c[k] * k^fds[k] / fds[k] , {k, Keys @ fds}] +
			Sum[c[k] * k, {k, Complement[vars, Keys @ fds]}]]


CleanSymbols[expr_, A_, c_] :=
	Module[
		{\[ScriptCapitalA], \[ScriptC]},
		expr /. {A[n_] :> Subscript[\[ScriptCapitalA], n], c[x_] :> Subscript[\[ScriptW], x]} /.
			sym_Symbol /; StringMatchQ[Context[sym], "*Private*"] :> 
				Symbol[StringSplit[SymbolName[sym], "$"][[1]]]]


(* ::Text:: *)
(*Main solver*)


ClearAll[GERFSolve];
GERFSolve::GERFPackageError = "`1`";
GERFSolve::InvalidFractionalDerivatives = "Multiple fractional orders received for these dimensions: `1`. This is not allowed.";

Options[GERFSolve] = {
	"RationalFunction" -> Function[n, 1/(1+E^n)],
	"WaveConstants" -> Automatic,
	"OutputMode" -> "Solutions",
	"BalanceConstant" -> Automatic
};
GERFSolve[ieqn_, u_[vars__], opts : OptionsPattern[]] :=
	Module[
		{eqn, ode, ts, U, eta, bal, A, R, poly, sol,
			c, \[Lambda], \[Mu], frac, options, R0, out},
		
		eqn = ((ieqn /. Equal -> Subtract) == 0);
		
		out = OptionValue["OutputMode"];
		
		(* Provision for custom wave constants *)
		If[
			OptionValue["WaveConstants"] =!= Automatic,
				If[Length @ OptionValue["WaveConstants"] != Length @ {vars},
					Message[GERFSolve::ConstantsLengthMismatch];
					Throw[$Failed]];
				Do[
					c[{vars}[[i]]] = OptionValue["WaveConstants"][[i]],
						{i, Length @ {vars}}]];
		
		R0 = OptionValue["RationalFunction"];

		Catch[
			(* convert eqn to ODE using wave transformation *)
			ode = ReducetoODE[eqn, u[vars], U, eta, c];
			bal = If[
				OptionValue["BalanceConstant"] === Automatic,
				BalanceConstant[eqn, u[vars]],
				OptionValue["BalanceConstant"]];
			If[!IntegerQ[bal] || bal <= 0, Throw[$Failed]];
			(* set up the trial solution *)
			ts = TrialSolution[bal, eta, {A, R}];
			poly = ExpandAll[ode /. U -> Function[{ieta},
				Evaluate[ts /. eta -> ieta]]];
			(* to solve the polynomial in R[eta] *)
			Throw[SolveAlgebraicSystem[
				(* feed everything to SolveAlgebraicSystem *)
				{ts, poly, R[eta], eqn, {vars}},
				{R0, {A, c, \[Lambda], \[Mu]}},
				{Solve @* Reduce, out}]]]]


(* ::Section::Closed:: *)
(*Package Footer*)


End[];
EndPackage[];
