# Computation of Bounds for Polynomial Dynamic Systems 

This project contains the source files for the computation of bounds on Lyapunov functions to approximate positive invariant sets such as attractors. The definiteness conditions are verified using quantifier elimination.

Röbenack, K., Gerbet, D.: 
*Computation of Bounds for Polynomial Dynamic Systems*    
Submitted 2025 to [Algorithms](https://www.mdpi.com/journal/algorithms).

## Prerequisites

### Quantifier elimination

You need to install the computer algebra system REDUCE. The REDUCE distribution is available for several operating systems:

http://reduce-algebra.sourceforge.net/

There are two main versions of REDUCE available, which are based on different Lips libraries:

Command | Library version 
:--- | :--- 
`redcsl`   | Codemist Standard Lisp (CSL) 
`redpsl`   | Portable Standard Lisp (PSL) 

The quantifier elimination is carried out using the package REDLOG. The program code is already part of REDUCE. The documentation can be found on the REDLOG website:

http://www.redlog.eu/

### Simplification

SLFQ (Simplifying Large Formulas with QEPCAD B) is a programm to simplify very large quantifier-free formulas.

https://www.usna.edu/Users/cs/wcbrown/qepcad/SLFQ/Home.html

The program requires QEPCAD B (Quantifier Elimination by Partial Cylindrical Algebraic Decomposition):

https://www.usna.edu/Users/cs/wcbrown/qepcad/B/QEPCAD.html

## Computation

### Structure of the Lyapunov candidate function

The necessary conditions for the time derivative of the Lyapunov candidate funcion along the system dynamics to be negative semi-definite are computed with REDUCE:

```
redcsl <chua-struture.red
```

The computation results in polynomial in Boolean combination of polynomial inequalities stored in the file `c.ineq` (1.8K). These formulas can be simplified with SLFQ:

```
slfq <c.ineq -a "b>0 /\ g>0 /\ c1>0 /\ c2>0 /\ a<0 /\ l>0"
```

The simplification yields the following formula:

```
[ p13 = 0 /\ p12 - p13 = 0 /\ p11 >= 0 ]
```

### Parametrization of the Lyapunov candidate function

Further conditions on the parameters of the Lyapunov canidate function using specific parameter values for the Chua system are calculated as follows:

```
redcsl <chua-parameters.red
```

The result is stored in the file `c0.ineq`. A simplification with SLFQ

```
slfq <c0.ineq
```

yields the following formula:

```
[ p23 < 0 /\ p22 + 10 p23 > 0 /\ p22 - p23^2 > 0 /\ 70 p33 - 10 p22 + 7 p23 = 0 ]
```


## Licence

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.