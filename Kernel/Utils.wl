(* ::Package:: *)

(* ::Section:: *)
(*Package Header*)


Begin["Taggar`GERF`Private`"];


(* ::Section:: *)
(*Definitions*)


(* ::Text:: *)
(*For initial checks of input equation:*)


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


GetDegreeofTerm[term_, m_] :=
	Module[{K, subbed, drules},
		
		(* Substitute all functions u[k] with K^m[k] *)
		drules = Flatten @ Table[With[{k = i}, {
			$GERFState["Function"][k] @@ $GERFState["Variables"] -> Power[K, m[k]],
			Derivative[orders__][$GERFState["Function"][k]] @@ $GERFState["Variables"] :> Power[K, m[k] + Total[{orders}]],
			(FractionalD | CaputoD)[$GERFState["Function"][k] @@ $GERFState["Variables"], {_, alpha_}] :> Power[K, m[k] + alpha]
		}], {i, 1, $GERFState["Length"]}];
		
		subbed = term /. drules;
		Simplify[Exponent[subbed, K]]]


LinearQ[term_] :=
	Module[{K, subbed, lrules},
		
		(* Substitute all functions u[k] with K to check for linearity (degree 1) *)
		lrules = Flatten @ Table[{
			$GERFState["Function"][k] @@ $GERFState["Variables"] -> K,
			Derivative[__][$GERFState["Function"][k]] @@ $GERFState["Variables"] :> K,
			(FractionalD | CaputoD)[$GERFState["Function"][k] @@ $GERFState["Variables"], __] :> K
		}, {k, 1, $GERFState["Length"]}];
		
		subbed = term /. lrules;
		Simplify[Exponent[subbed, K]] == 1]


(* ::Text:: *)
(*Function to extract the balance constant of equation when "BalanceConstant" option is set to Automatic:*)


BalanceConstant[] :=
	Module[
		{sys = {}, k, m, sols, ld, nld, expr, terms},
		
		For[k = 1, k <= $GERFState["Length"], k++,
			expr = Expand @ If[Head[$GERFState["Equation"][[k]]] === Equal,
				Subtract @@ $GERFState["Equation"][[k]], 
				$GERFState["Equation"][[k]]];
			terms = If[Head[expr] === Plus, List @@ expr, {expr}];
			
			(* Extract degrees using the new system-aware helpers *)
			ld = Simplify[GetDegreeofTerm[#, m] & /@ Select[terms, LinearQ]];
			nld = Simplify[GetDegreeofTerm[#, m] & /@ Select[terms, Not @* LinearQ]];
			
			If[Length[ld] == 0 || Length[nld] == 0,
				Message[GERFSolve::GERFPackageError, "Balance constant could not be calculated for equation " <> ToString[k] <> "."];
				Return[$Failed]];
			
			(* Build the system by equating the highest degree terms *)
			AppendTo[
				sys, 
				Last[SortBy[ld, # /. {m[_] :> 100} &]] == Last[SortBy[nld, # /. {m[_] :> 100} &]]
			];
		];
		
		Print[sys];
		(* Solve the simultaneous system for all m[k] *)
		sols = Solve[sys, Table[m[i], {i, 1, $GERFState["Length"]}]];
		
		If[Length[sols] > 0,
			(* Grab the first solution and format as an Association for $GERFState *)
			Last /@ sols[[1]],
			Message[GERFSolve::GERFPackageError, "No valid balance constants found for the system."];
			$Failed
		]]


(* ::Text:: *)
(*Reducing the PDE to ODE:*)


ReducetoODE[] :=
	Module[
		{dr, U, eta},
		
		$GERFState @ "AuxiliaryFunction" = AssociationMap[U, Range @ $GERFState @ "Length"];
		$GERFState @ "Eta" = eta;
		
		(* applying the derivative rule given as: *)
			dr = Flatten @ Table[
				With[{k = i}, {(FractionalD | CaputoD)[$GERFState["Function"][k]@@$GERFState["Variables"], {var_, alpha_}] :>
					$GERFState["WaveConstant"] @ var * Derivative[1][U[k]][eta],
					(* only one FractionalD at a time please *)
					(* D^alpha u_x = x^(1-alpha) U'(eta) deta/dx
						= x^(1-alpha) U'(eta) a x^(alpha-1)
						= a U'(eta) *)
					Derivative[orders__][$GERFState["Function"][k]]@@$GERFState["Variables"] :>
						Times @@ Power[$GERFState["WaveConstant"] /@ $GERFState["Variables"], {orders}] *
							Derivative[Total[{orders}]][U[k]][eta]}],
					(* simply using this rule, for example u_x converts into 
						a*U', or u_xxy -> a^2bU''', and so on *)
					{i, 1, $GERFState @ "Length"}];
		(* to the equation and return *)
			$GERFState["Equation"] /. dr /. Table[With[{k = i},
				$GERFState["Function"][k]@@$GERFState["Variables"] :> U[k][eta]],
				{i, 1, $GERFState @ "Length"}]]


(* ::Text:: *)
(*Constructing the ansatz for Subscript[U, k]:*)


TrialSolution[k_] :=
	Module[
		{A, R, eta},
		
		$GERFState @ "TrialSolutionCoefficient" = A;
		$GERFState @ "SymbolicRationalHead" = R;
		
		(* a function of eta *)
		eta = $GERFState @ "Eta";
		Return @ Function[
			eta,
			Sum[
				A[i, k] R[eta]^i, 
				{i, -$GERFState["BalanceConstant"] @ k, $GERFState["BalanceConstant"] @ k}]]]


(* ::Text:: *)
(*Construct the auxiliary polynomials:*)


AuxiliaryPolynomial[] :=
	ExpandAll[
		ReplaceAll[$GERFState @ "ODE",
			Table[With[{k = i},
			$GERFState["AuxiliaryFunction"][k] :> $GERFState["TrialSolution"][k]],
			{i, 1, $GERFState @ "Length"}]]]


(* ::Text:: *)
(*Algebraic system solver*)


SolveAuxiliaryPolynomial[] :=
	Module[
		{tosolve, sys, u, v, sols, for, pairs},
		
		(* symbolic R must be replaced with the actual corresponding
			form before any further evaluation to ensure proper
			solutions appear *)
		tosolve = ReplaceAll[
			$GERFState @ "AuxiliaryPolynomial",
			$GERFState @ "SymbolicRationalHead" -> $GERFState["Options"] @ "RationalFunction"];

		sys = Thread[
			CoefficientList[
				Expand @ Numerator @ Together @ TrigToExp @ (First /@
					tosolve) /. Exp[d_. * $GERFState @ "Eta"] :> u^Re[d] v^Im[d],
					{u, v}] == 0];
		
		(* solve the system usiing Solve@*Reduce *)
		for = Flatten[Table[
					Table[
						$GERFState["TrialSolutionCoefficient"][i, k],
						{i, -$GERFState["BalanceConstant"][k], $GERFState["BalanceConstant"][k]}],
					{k, 1, $GERFState @ "Length"}]];
		sols = Select[
			Solve @
				Reduce[sys, for],
			Length @ DeleteDuplicates @ #[[All, 2]] > 1 &];
		
		pairs = Select[
			DeleteDuplicates @
				Table[
					With[
						{form = Table[
							$GERFState["TrialSolution"][k][FormWaveTransformation[]],
							{k, 1, $GERFState["Length"]}] /.
								$GERFState @ "SymbolicRationalHead" ->
									$GERFState["Options"] @ "RationalFunction" /.
										sol},
						{sol, Thread[Through[($GERFState["Function"]/@Range[$GERFState["Length"]])@@$GERFState["Variables"]]
							-> form]}],
					{sol, sols}],
			FreeQ[Last[#], Indeterminate | ComplexInfinity | Infinity | DirectedInfinity] &];
		
		Throw @ CleanSymbols[
			Switch[$GERFState["Options"] @ "OutputMode",
				"SolutionSets", First /@ pairs,
				All, Transpose[{First /@ pairs, Last /@ pairs}],
				_, Last /@ pairs],
		  $GERFState @ "TrialSolutionCoefficient", $GERFState @ "WCH"]]


(* ::Text:: *)
(*Make solution forms*)


FormWaveTransformation[] :=
	Module[
		{raw, fds, defaulters},
		
		(* all cases of fractional orders: *)
		raw = DeleteDuplicates @ Cases[$GERFState @ "Equation",
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
		Sum[$GERFState["WaveConstant"][k] * k^fds[k] / fds[k] , {k, Keys @ fds}] +
			Sum[$GERFState["WaveConstant"][k] * k, {k, Complement[$GERFState @ "Variables", Keys @ fds]}]]


CleanSymbols[expr_, A_, w_] :=
	Module[
		{\[ScriptCapitalA], \[ScriptW]},
		expr /. {A[n_, k_] :> 
			If[$GERFState["Length"] == 1, Subscript[\[ScriptCapitalA], n], Subsuperscript[\[ScriptCapitalA], n, k]],
			w[x_] :> Subscript[\[ScriptW], x]} /.
			sym_Symbol /; StringMatchQ[Context[sym], "*Private*"] :> 
				Symbol[StringSplit[SymbolName[sym], "$"][[1]]]]


(* ::Section:: *)
(*Package Footer*)


End[];
