# R for Frugal Hawkes Process

1. **Simulation_Preparation.R** prepares necessary functions for simulation of the frugal Hawkes processes.
2. **Inference_Preparation.R** prepares necessary functions for inference of the frugal Hawkes processes.
3. **Simulation_Methods.R** simulates frugal Hawkes processes with copula/thinning/marginal modelling based methods.
4. **Simulation_Test.R** compares the performance of different simulation algorithms with ks/cvm/ad tests via random time change theorem.
5. **Inference_MLE.R** implements the inference algorithm for frugal Hawkes processes.
6. **Seismic_Data_Example.R** models seismic activities of two regions with the frugal Hawkes process.

# Errata for Frugal Hawkes Process

```math
begin{bmatrix}X\\Y\end{bmatrix}
```
1. P9: In Assumption 13, replace “$N^{P_i}$” with “$ \widetilde{N}^{P_i} $” (also in P11 and the third-to-last line in P20). 
2. P10: In the remark, correct “$ C’(m, n) \neq C’(m, n) $” with “$ C’(m, n) \neq C’(n, m) $”. 
3. P18: In definition 18, correct “$ \nu(D) $” with “$ \nu(Q) $” and “$ M:[a, b] \rightarrow \mathbb{R}^{n \times n} $” with “$ A:[a, b] \rightarrow \mathbb{R}^{n \times n} $”. 
4. P23: In the first paragraph, clarify  "it is convenient to write the original distribution
$ dL(O) $ for the general process $ N $". 
5. P30: In the first paragraph, change “$ N^k $” with “$ N^K $”. 
5. 	P32: Replace the first equation with 
$$\text{"} \lambda_{t}^{i,i} = \nu_{i} + \int_{0}^{t} \phi_{i}(u) d \widetilde{N}_{t-u}^{i} = \nu_{i} + \int_{0}^{t} \alpha_{i} e^{-\beta_{i} u} d \widetilde{N}_{t-u}^{i} \text{"}$$ 
and "$ \lambda_s^{i}$" in 5.3 with "$ \lambda_s^{i,i} $".
