# Updated Read Me In Progress

# Weight-Multiplicities

This is a group of programs in addendum to the article [*On Kostant's Weight $q$-Multiplicity Formula for $\mathfrak{sp}_6(\mathbb{C})$*](http://https://github.com/21mdr1/Weight-Multiplicities "Will have a link when we have a link") written by the authors of this repository. These programs perform various computations related to representations of the Lie algebra of type $C_3$ ($\mathfrak{sp}_6(\mathbb{C})$).

## Programs Pretaining Solely to the Article

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

### q-mult-formula.py

### calculat-alt-sets.ipynb

### C3_weight_multiplicities.sage

- Compute Kostant's partition function and its $q$-analog for a weight. Parameter ``xi`` is a triple ``(m,n,k)`` representing a linear combination of the simple roots $m\alpha_1+n\alpha_2+k\alpha_3$. 
```Python
kostant_partition_function(xi, q_analog=True, triple_sum=False)
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

## Prerequisites - will need to add subsections with Prerequisites for each program language, etc.

- A working installation of SageMath (Tested on Version 8.8)
- Java Runtime Environment that can run jmol (for empty region plots)

## Usage - do we even need this? if yes need to add subsections for each program - or just straight up add it to the program sections above

- Open Sage in the project directory
- Run the following code to begin:
```Python 
load('C3_weight_multiplicities.sage')
```
- Use functions as specified above

## Examples - may need more / new ones, these are for A3

First we compute the Weyl alternation set for $\lambda=2\omega_1+3\omega_2+4\omega_3$ and $\mu=0$. We then compute the weight $q$-multiplicity for $\lambda$ and $\mu$.
```Python
> load('C3_weight_multiplicities.sage')
> weyl_alternation_set(lam=(2,3,4), mu=(0,0,0))
{s3*s1, s1, s2, s3, 1}
> kostant_weight_multiplicity(lam=(2,3,4), mu=(0,0,0), q_analog=True)
q^15 + 2*q^14 + 4*q^13 + 5*q^12 + 6*q^11 + 5*q^10 + 4*q^9 + 2*q^8 + q^7
```
Here we demonstrate and compare runtimes for Kostant's partition function using the triple summation formula and our closed formulas:
```Python
> kostant_partition_function(xi=(2,5,4), triple_sum=False)
q^11 + 2*q^10 + 4*q^9 + 5*q^8 + 6*q^7 + 5*q^6 + 2*q^5
> kostant_partition_function(xi=(2,5,4), triple_sum=True)
q^11 + 2*q^10 + 4*q^9 + 5*q^8 + 6*q^7 + 5*q^6 + 2*q^5
> %timeit kostant_partition_function(xi=(50,70,60), q_analog=True, triple_sum=False)
1000 loops, best of 3: 1.54 ms per loop
> %timeit kostant_partition_function(xi=(50,70,60), q_analog=True, triple_sum=True)
1 loop, best of 3: 398 ms per loop
```
Here, we plot the restricted Weyl alternation diagram for {1, s1, s2} and ``mu=0``:
```Python
> plot_alternation_diagram([e,s1,s2], mu=(0,0,0), restricted=True, size=40)
Launched jmol viewer for Graphics3d Object
```
Output:

![alternation diagram](1s1s2.png?raw=true)

