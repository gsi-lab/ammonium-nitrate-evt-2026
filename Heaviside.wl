(* ::Package:: *)

(* Define symbolic variables *)
Clear[x, alpha, beta, L, H]

(* Define the PDF f(x) for the Pareto-Lomax distribution *)
f[x_, alpha_, beta_] := (alpha / beta) (1 + x / beta)^(-(alpha + 1))

(* Define the lower and upper bounds as variables *)
L = 10; (* Replace with your lower bound *)
H = 581; (* Replace with your upper bound *)

(* Calculate the integral of x*f(x) over the specified bounds *)
int_f1 = Integrate[x * f[x, alpha, beta], {x, L, H}]
int_f2 = Integrate[f[x, alpha, beta], {x, L, H}]
hv_mean = int_f1 / int_f2

(* Display the result *)
Print["The integral of x*f(x) over [L, H] for the Pareto-Lomax distribution is:"]
simplify(int_f1)
simplify(int_f2)
simnplify(hv_mean)





