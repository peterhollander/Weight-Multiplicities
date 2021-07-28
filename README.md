# Weight-Multiplicities

This is a group of programs in addendum to the article [*On Kostant's Weight $q$-Multiplicity Formula for $\mathfrak{sp}_6(\mathbb{C})$*]("Will have a link when we have a link") written by the authors of this repository. These programs perform various computations related to representations of the Lie algebra of type $C_3$ ($\mathfrak{sp}_6(\mathbb{C})$).

## Programs Pertaining Solely to the Article

### sigma calculations.nb

In the article we need to calculate $\sigma(\lambda + \rho) - \mu - \rho$ for every $\sigma$ in the Weyl group, $W$. To use this program, first run the first cell to send the definitions to the kernel. The parameter `Sigma` is a $\sigma \in W$ written with the $s_i$'s separated by commas (ie. $s_1s_2s_3$ becomes s1, s2, s3).

```Mathematica
Simplify[Composition[Simplify, Sigma][
   Distribute[(m + n + k + 3)*alpha1] + 
    Distribute[(m + 2 n + 2 k + 5)*alpha2] + 
    Distribute[(m/2 + n + (3/2)*k + 3)*alpha3]] + 
  Distribute[-(x + y + z + 3)*alpha1] + 
  Distribute[-(x + 2 y + 2 z + 5)*alpha2] + 
  Distribute[-(x/2 + y + (3/2)*z + 3)*alpha3]]
  ```

### remove-contradictions.py (I, II, and I and II)

In the article, we define 2 types of contradictions which arise forbidding certain subsets of the Weyl group from appearing as Weyl alternation sets. updated-alt-sets.py contains the lists returned by remove-contradictions-type-I.py (`currentAltSets1`), remove-contradictions-type-II.py when `currentAltSets1` is passed (`currentAltSets2`), and remove-contradictions-type-I-and-II.py when `currentAltSets2` is passed (`currentAltSets3`).

- remove-contradictions-type-I.py returns the subsets of the Weyl group which do not cause a contradiction of type I. 
- remove-contradictions-type-II.py, when given a list of subsets of the Weyl group, returns the subsets which do not cause a contradiction of type II. The parameter `sets` is a list of sets (the sets being subsets of the Weyl group). `currentAltSets1` is the list returned by remove-contradictions-type-I.py
  ```Python
  print(removeSubsets(sets))
  ```
- remove-contradictions-type-I-and-II.py, when given a list of subsets of the Weyl group, returns the subsets which do not cause a contradiction of type I or II or both together. The parameter `sets` is a list of sets (the sets being subsets of the Weyl group). `currentAltSets2` is the list returned by remove-contradictions-type-II.py
  ```Python
  print(removeSubsets(sets))
  ```

## Other Programs

### Evaluating the q-analog of Kostant's partition function.nb

Computes the $q$-analog of Kostant's partition function for a weight. Parameters ``m``, ``n``, and ``k`` are the triple representing a linear combination of the simple roots $m\alpha_1+n\alpha_2+k\alpha_3$. 

```Mathematica
P[m, n, k]
```

### q-mult-formula.py

Calculates and Kostant's Weight $q$-multiplicity formula in terms of the capital letters $A$ through $Q$ (defined in the article). To run, type python3 q-mult-formula.py, you will be prompted to enter lambda, then mu as triples of the form `(x,y,z)` representing linear combinations of the fundamental weights $\omega_1,\omega_2,\omega_3$.

### calculate-alt-sets.sage

- Compute Weyl alternation sets. Parameters ``lam`` and ``mu`` are triples representing linear combinations of the fundamental weights $\omega_1,\omega_2,\omega_3$. 
```Python
weyl_alternation_set(lam,mu=(0,0,0))
weyl_alternation_set_17(lam, mu=(0,0,0))
```

- Gives the Weyl alternation sets, together with a $\lambda$ and $\mu$ pair (printed as triples) which induce the alternation set, which appear for $\lambda$ with coefficients between 0 and `lamnum` and $\mu$ with coefficients between 0 and `munum`.
```Python
alt_sets_with_mus(lamnum,munum)
```

- Gives the Weyl alternation sets, together with a $\lambda$ and $\mu$ pair (printed as triples) which induce the alternation set, which appear for $\lambda$ with coefficients between 0 and `lamnum` and $\mu$ with coefficients between 0 and `munum` for $m+k+x+z$ divisible by 2 (when $\lambda$ equals the triple `(m,n,k)` and $\mu$ equals the triple `(x,y,z)`).
```Python
alt_sets_new_lattice(lamnum, munum)
```

### C3_weight_multiplicities.sage

- Compute Kostant's partition function and its $q$-analog for a weight. Parameter ``xi`` is a triple ``(m,n,k)`` representing a linear combination of the simple roots $m\alpha_1+n\alpha_2+k\alpha_3$. 
```Python
kostant_partition_function(xi, q_analog=True)
```
- Compute Kostant's weight multiplicity formula and its $q$-analog. Parameters ``lam`` and ``mu`` are triples representing linear combinations of the fundamental weights $\omega_1,\omega_2,\omega_3$.
```Python
kostant_weight_multiplicity(lam, mu=(0,0,0), q_analog=True)
```
 - Compute Weyl alternation sets. Parameters ``lam`` and ``mu`` are triples representing linear combinations of the fundamental weights $\omega_1,\omega_2,\omega_3$.
 ```Python
 weyl_alternation_set(lam,mu=(0,0,0))
 ```
 - Compute and display empty regions (weights lambda for which no Weyl group elements contribute to the weight multiplicity m(lambda, mu)). Parameter ``mu`` is a triple representing a linear combination of the fundamental weights $\omega_1,\omega_2,\omega_3$. Parameter ``dist`` specifies the maximum magnitude of the components of lambda the program should check and plot. Parameters ``color`` and ``size`` specify the color and size of the points in the plot
 ```Python
 plot_empty_region(mu, dist=15, color='red', size=25)
 ```
 - Plot Weyl alternation diagrams for a given weight mu and a list of Weyl group elements. If not restricted, plots weights lambda such that the alternation set for lambda, mu is contains ``sigmas``. If restricted, plots weights lambda such that the alternation set for lambda, mu is exactly ``sigmas``. 
 ```Python
 plot_alternation_diagram(sigmas, mu=(0,0,0), restricted=False, color='red', dist=20, size=15)
 ```

## Prerequisites

##### For Sage programs

- A working installation of SageMath (Tested on Version 8.8)
- Java Runtime Environment that can run jmol (for empty region plots)

##### For Python programs

- A working installation of Python (Tested on Version 3)

##### For Mathematica Programs

- A working installation of Mathematica (Tested on Version 12.0.0 and higher)

## Usage

##### For Sage programs

- Open Sage in the project directory
- Run the following code to begin (for filename.sage the name of the file you desire to use):
  ```Python 
  load('filename.sage')
  ```
- Use functions as specified above

##### For Python programs

- Run the following code to the command line (where filename.py is the file you'd like to use, functionname is the functions you'd like to use, and parameters are the parameters of the function):
  ```shell
  python3 filename.py functionname(parameters)
  ```
- Alternatively, open the program in Python and use the functions as specified above

##### For Mathematica programs

- Open filename.nb in Mathematica (for filename.nb the file you'd like to use)
- Use functions as specified above

## Examples

First we compute the Weyl alternation set for $\lambda=2\omega_1+3\omega_2+4\omega_3$ and $\mu=0$. We then compute the weight $q$-multiplicity for $\lambda$ and $\mu$.
```Python
> load('C3_weight_multiplicities.sage')
> weyl_alternation_set(lam=(2,3,4), mu=(0,0,0))
{s1*s2*s1, s3*s1*s2, s1*s2, s2*s1, s3*s1, s1, s2*s3, s3*s2, s2, s3, 1}
> kostant_weight_multiplicity(lam=(2,3,4), mu=(0,0,0), q_analog=True)
q^35 + 2*q^34 + 5*q^33 + 8*q^32 + 14*q^31 + 19*q^30 + 28*q^29 + 34*q^28 + 45*q^27 + 50*q^26 + 61*q^25 + 63*q^24 + 72*q^23 + 69*q^22 + 74*q^21 + 64*q^20 + 64*q^19 + 48*q^18 + 44*q^17 + 27*q^16 + 22*q^15 + 10*q^14 + 7*q^13 + 2*q^12 + q^11
```

Here, we plot the restricted Weyl alternation diagram for {1, s1} and ``mu=0``:
```Python
> plot_alternation_diagram([e,s1], mu=(0,0,0), restricted=False, size=40)
Launched jmol viewer for Graphics3d Object
```
Output:

![alternation diagram](1s1.png?raw=true)
