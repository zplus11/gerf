<span style="float:left;padding-right:10"><img src="img/gerf.png" width="200"></span> is a Mathematica paclet that solves nonlinear integral order partial differential equations using the GERF (generalized exponential rational function) expansion technique. For more info on the method, see [10.1140/epjp/i2018-11984-1](https://doi.org/10.1140/epjp/i2018-11984-1). To learn more about the paclet, visit the homepage at [Wolfram repository](https://resources.wolframcloud.com/PacletRepository/resources/Taggar/GERF/).

# Installation

is as easy as running

```mathematica
PacletInstall["Taggar/GERF"]
```

in your local or cloud notebook. Then, load it with

```mathematica
<<Taggar`GERF`
```

# Usage

Let

```mathematica
burgers = D[u[x, t], t] + u[x, t] * D[u[x, t], x] - \[Nu] * D[u[x, t], x, x] == 0
```

be the given equation (this is the Burgers' equation in (1+1)-dimensions). Then, use GERFSolve as follows:

```mathematica
sol = GERFSolve[burgers, u[x, t]]
```

which returns the following output:

![Solutions](img/sols.png)

Pick any of these and plot it for appropriate values:

```mathematica
Plot3D[
    u[x, t] /. sol[[3]] /. {(* parameters *)},
    {x, -4, 4}, {t, 0, 4},
    PlotRange -> All
]
```
![Graph](img/plot.png)

A full demonstration is available on the Wolfram repository page of this paclet, see [GERFSolve.html](https://resources.wolframcloud.com/PacletRepository/resources/Taggar/GERF/ref/GERFSolve.html).

# Version log

**Version 1.0.0,** *on 31 March, 2026* &mdash; initial upload.

**Version 1.1.0,** *on 30 May, 2026* &mdash; support for fractional ordered equations.

**Version 1.2.0,** *15 June, 2026* &mdash; support for systems of equations.

# Contributions

are welcome.