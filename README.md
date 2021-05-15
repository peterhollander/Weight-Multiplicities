# Weight-Multiplicities
A Sage program to preform various computations related to representations of the Lie algebra of type $A_3$

## Features
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

## Prerequisites
- A working installation of SageMath (Tested on Version 8.8)
- Java Runtime Environment that can run jmol (for empty region plots)
## Usage
- Open Sage in the project directory
- Run the following code to begin:
```Python 
load('A3_weight_multiplicities.sage')
```
- Use functions as specified above
## Examples
First we compute the Weyl alternation set for $\lambda=2\omega_1+3\omega_2+4\omega_3$ and $\mu=0$. We then compute the weight $q$-multiplicity for $\lambda$ and $\mu$.
```Python
> load('A3_weight_multiplicities.sage')
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
