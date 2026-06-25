# ON-DIAGONALIZABLE-SYSTEMS-WITH-RANDOM-STRUCTURE_Simulations
MATLAB simulation code for studying the probability of structural diagonalizability in the random directed graph models G(n,p) and G(n,p,q), including Monte Carlo experiments, theoretical bounds, and reproduction of the numerical figures.

# Representative MATLAB Codes for Structural Diagonalizability in Random Directed Graphs

This repository provides representative MATLAB implementations for the numerical experiments on structural diagonalizability in random directed graphs. The codes evaluate the probability of structural diagonalizability for the random graph models \mathcal{G}(n,p) and \mathcal{G}(n,p,q), compare Monte Carlo estimates with the corresponding theoretical upper and lower bounds, and investigate finite-size behavior as the number of vertices increases.

The minimum-cost perfect matching problems appearing in the structural diagonalizability test are solved using the Munkres algorithm implemented in `munkres.m`.

## Repository Structure

```text
.
├── gnp_probability_vs_c.m
├── gnpq_probability_vs_q.m
├── gnp_finite_size_logscale.m
├── gnpq_finite_size_logscale.m
├── munkres.m
└── README.md
```

The suggested file names correspond to the following scripts.

| File | Description |
|---|---|
| `gnp_probability_vs_c.m` | Fixes n, varies the transition parameter c, and compares the empirical structural diagonalizability probability of \mathcal{G}(n,p) with the theoretical bounds. |
| `gnpq_probability_vs_q.m` | Fixes n and c, varies the self-loop probability q, and compares the empirical probability for \mathcal{G}(n,p,q) with the theoretical bounds. |
| `gnp_finite_size_logscale.m` | Studies finite-size behavior of \mathcal{G}(n,p) for several representative values of c and several graph sizes n. The horizontal axis uses a true logarithmic scale. |
| `gnpq_finite_size_logscale.m` | Studies finite-size behavior of \mathcal{G}(n,p,q) for several representative values of q, with c fixed, and uses a true logarithmic horizontal axis. |
| `munkres.m` | Implements the Munkres/Hungarian algorithm used to solve the minimum-cost perfect matching problem required by the structural diagonalizability test. |

## Random Graph Models

### The model \mathcal{G}(n,p)

For a directed graph with n vertices, every possible directed edge, including every self-loop, is independently present with probability

p=\frac{\log n+c}{n},


where c controls the location within the transition region.

The theoretical lower and upper bounds used in the codes are

exp(-2.*exp(-c_values)).*(1+2.*exp(-c_values)+exp(-2.*c_values)),

and

1-exp(-2*exp(-c_values)).*exp(-2*c_values).


### The model \(\mathcal{G}(n,p,q)\)

For \mathcal{G}(n,p,q), every non-self-loop directed edge is independently present with probability


p=\frac{\log n+c}{n},


whereas every self-loop is independently present with probability q.

The theoretical lower and upper bounds used in the codes are

 exp(-2 .* (1 - q) .* exp(-c)) .* (1 + 2 .* (1 - q) .* exp(-c) + (1 - q).^2 .* exp(-2 * c)),

and

1 - (1 - q).^2 .* exp(-2 * c) .* exp(-2 .* (1 - q) .* exp(-c)).


## Structural Diagonalizability Test

For each randomly generated graph, the adjacency pattern is first converted to the system-matrix convention by transposition. Its generic rank is computed using `sprank` in the large-scale scripts.

A cost matrix is then constructed according to the following rule:

- an existing edge has cost \(0\);
- a missing diagonal edge has cost \(1\);
- a missing non-diagonal edge has a sufficiently large cost.

The minimum-cost perfect matching is computed by

```matlab
[assignment, min_totalcost] = munkres(A_cost);
```

The realization is counted as structurally diagonalizable when

```matlab
generic_rank == n - min_totalcost
```

holds. The empirical probability is the fraction of Monte Carlo samples satisfying this criterion.

## Description of the Representative Experiments

### 1. Transition with respect to \(c\) in \(\mathcal{G}(n,p)\)

Run

```matlab
gnp_probability_vs_c
```

The script fixes \(n=4000\), samples \(c\) over the interval \([-1,2]\), and performs Monte Carlo simulations for every value of \(c\). It plots:

- the theoretical lower bound;
- the theoretical upper bound;
- the empirical structural diagonalizability probability.

The numerical values of \(c\), \(p\), and the estimated probabilities are also printed in the MATLAB command window.

### 2. Dependence on \(q\) in \(\mathcal{G}(n,p,q)\)

Run

```matlab
gnpq_probability_vs_q
```

The script fixes \(n=4000\) and \(c=0\), varies \(q\) from \(0\) to \(1\), and compares the Monte Carlo estimates with the theoretical bounds. A 600-dpi figure is exported as

```text
gnpq_probability_vs_q_c0.png
```

to the Windows Desktop when it is available; otherwise, the figure is saved in the current working directory.

### 3. Finite-size study for \(\mathcal{G}(n,p)\)

Run

```matlab
gnp_finite_size_logscale
```

The script considers representative values

```matlab
selected_c = [-0.5, 0, 0.5, 1.0, 1.5, 2.0];
```

and graph sizes

```matlab
n_values = [100, 500, 1000, 3000, 6000, 10000];
```

A six-panel figure is generated. Each panel contains the empirical probability curve and the theoretical interval associated with one value of \(c\). The graph sizes are displayed at their actual positions on a logarithmic horizontal axis.

Intermediate results are saved in

```text
gnp_logscale_checkpoint.mat
```

and the final figure is exported as

```text
subplot_gnp_true_logscale.png
```

### 4. Finite-size study for \(\mathcal{G}(n,p,q)\)

Run

```matlab
gnpq_finite_size_logscale
```

The script fixes \(c=0\), considers

```matlab
selected_q = [0, 0.2, 0.4, 0.6, 0.8, 1.0];
```

and uses the same graph sizes as the \(\mathcal{G}(n,p)\) finite-size experiment. Each panel shows the empirical probability and the theoretical interval corresponding to one value of \(q\).

Intermediate results are saved in

```text
gnpq_logscale_checkpoint.mat
```

and the final figure is exported as

```text
subplot_gnpq_true_logscale.png
```

## The `munkres.m` Function

The function

```matlab
[assignment, cost] = munkres(costMat)
```

solves a linear assignment problem, equivalently a minimum-cost perfect matching problem on a weighted bipartite graph.

### Input

- `costMat`: a real-valued cost matrix. Entry `costMat(i,j)` is the cost of assigning row \(i\) to column \(j\).
- `NaN` and `Inf` entries are treated as invalid assignments and internally replaced by a sufficiently large penalty.
- Rectangular matrices are supported by internally padding the problem to a square matrix.

### Output

- `assignment`: an assignment vector. When `assignment(i)=j`, row \(i\) is assigned to column \(j\). An unassigned or invalid row is represented by zero.
- `cost`: the total cost of the returned valid assignment.

### Main computational steps

The implementation follows the standard Munkres/Hungarian procedure:

1. invalid entries are identified and replaced by a large penalty;
2. rows and columns containing no valid entries are removed;
3. a rectangular problem is padded to a square problem;
4. row and column reductions are performed;
5. independent zeros are starred to construct an initial matching;
6. uncovered zeros are primed and alternating augmenting paths are generated;
7. the reduced-cost matrix is updated until a complete optimal assignment is obtained.

The local helper function

```matlab
[minval, rIdx, cIdx] = outerplus(M, x, y)
```

computes reduced costs relative to the current row and column potentials, finds the minimum uncovered reduced cost, and returns the locations at which this minimum is attained. This update is used when no uncovered zero is available.

In the present repository, `munkres.m` is the key optimization routine that converts the matching-based structural diagonalizability criterion into a computable minimum-cost assignment problem.

## Requirements

- MATLAB R2020a or later is recommended.
- `munkres.m` must be located in the current folder or on the MATLAB search path.
- The scripts use standard MATLAB functions including `sprank`, `sparse`, `tiledlayout`, and `exportgraphics`.

To verify that MATLAB can locate the matching routine, run

```matlab
which munkres
```

## Quick Test

The default settings are intended for the numerical experiments reported in the paper and can be computationally demanding. For a quick test, reduce the graph sizes and the number of Monte Carlo samples, for example:

```matlab
n = 200;
num_samples = 20;
```

or, for the finite-size scripts,

```matlab
n_values = [100, 200, 500];
num_samples = 20;
```

After confirming that the scripts run correctly, restore the original parameters to reproduce the reported experiments.

## Computational Notes

The assignment problem is solved on an \(n\times n\) cost matrix. The classical Munkres algorithm has cubic worst-case time complexity and the current implementation stores dense matrices. Consequently, experiments with large \(n\), especially \(n=6000\) or \(n=10000\), require substantial memory and computation time.

For the full experiments, a high-memory workstation or computing server is recommended. The checkpoint files produced by the finite-size scripts preserve completed numerical results after each parameter pair.

The Monte Carlo estimates may vary slightly between runs. Scripts containing

```matlab
rng(1);
```

use a fixed random seed for reproducibility. The same line may be added near the beginning of any script for which an exactly repeatable random sequence is desired.

## Reproducing the Figures

1. Place all `.m` files in the same directory.
2. Start MATLAB and set that directory as the current folder.
3. Confirm that `munkres.m` is visible by running `which munkres`.
4. Run one of the four main scripts.
5. Inspect the printed numerical summary and the generated figure.
6. For the finite-size experiments, retain the generated `.mat` checkpoint files.

## Code Availability

Representative MATLAB implementations used in the numerical experiments are provided in this repository. They include simulations for \(\mathcal{G}(n,p)\) and \(\mathcal{G}(n,p,q)\), comparisons with the theoretical probability bounds, finite-size studies, and the minimum-cost matching routine used in the structural diagonalizability test.
