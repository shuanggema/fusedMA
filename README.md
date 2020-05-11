Algorithm
-------
fusedMA: Heterogeneity analysis via penalized fusion with model averaging


Usage
-------
f_group_T (groupTMA.R) is the main function to obtain the subgrouping structure and estimation.

f_group_T(M,n,Y,Z,Z_matrix,Z_inv,A,
          d,h,N,B,sample_z,Lambda1,
          Lambda2,kappa,gamma,iter_admm,eps,penalty_index,
         iter_fw,Intial_theta)
         

Arguments
-------  
M: the size of machines if you consider dividing the samples, here we consider "M=1", which indicates considering all the samples.

n: the sample size in each machine if you consider dividing samples.

Y: is the response vector.

Z: is the covariates matrix.

Z_matrix: is the diagonal matrix of the covariates matrix.

Z_inv: is the inverse item of equation (9).

A: is the matrix to formulate the difference between coefficients into matrix.

d: is the dimension of covariates.

h: is the dimension of covariates in each candidate model.

N: is the total sample size.

B: is the size of candidate model.

sample_z: is the index set for constructing the candidate model.

Lambda1: is the tuning parameter for sparsity, it can either be a value or a vector.

Lambda2: is the tuning parameter for subgrouping, it can either be a value or a vector.

kappa,gamma:  is the parameter for penalty "SCAD" and "MCP".

iter_admm:is the maximum iteration for ADMM algorithm.

eps: is the tolerance values for ADMM and greedy algorithm.

penalty_index: indicates which penalty would be adopted, "1" denotes "LASSO", "2" denotes "MCP"
and "3" denotes "SCAD".

iter_fw: is the maximum iterations for greedy algorithm.

Intial_theta: is the initial values for individual coefficients.



Output
-------  
a list containing the coefficients for each samples and subgroups, 
the size of estimated subgroup and the subgrouping structure and 
the EBIC values for different tuning parameter.


Example
-------  
Please refer to codes in main.R