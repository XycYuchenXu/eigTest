# eigTest
Eigenvector testing.

To install, use the code:
`devtools::install_github('XycYuchenXu/eigTest')`

The package has the following functionalities:

1. Check the means (<img src="https://render.githubusercontent.com/render/math?math=\mu_1, \mu_2">) of two random matrices (<img src="https://render.githubusercontent.com/render/math?math=A_1, A_2">) are simultaneously diagonalizable, either considering the commutator of the two matrices (use function `commutatorTest`), or using the log-likelihood ratio test framework (`projTest`). The random matrices <img src="https://render.githubusercontent.com/render/math?math=A_1, A_2"> are supposed to be consistent estimators of the means with asymptotic normal distribution with limiting covariance matrices <img src="https://render.githubusercontent.com/render/math?math=\Sigma_1, \Sigma_2">. For both functions, one needs to input a list of two matrices, a list of two limiting covariance matrices, and the convergence rate.

2. Given square matrices <img src="https://render.githubusercontent.com/render/math?math=A_1, \dots, A_p">, estimate the matrix <img src="https://render.githubusercontent.com/render/math?math=V"> such that for every <img src="https://render.githubusercontent.com/render/math?math=V">, <img src="https://render.githubusercontent.com/render/math?math=V^{-1} A_i V = D_i"> with diagonal matrices <img src="https://render.githubusercontent.com/render/math?math=D_i"> (by supplying the list of matrices <img src="https://render.githubusercontent.com/render/math?math=A_1, \dots, A_p"> for function `JDTE`).

2. Test the means (<img src="https://render.githubusercontent.com/render/math?math=\mu_1, \dots, \mu_p">) corresponding to random matrices (<img src="https://render.githubusercontent.com/render/math?math=A_1, \dots, A_p">) are simultaneously diagonalizable, by checking whether the estimated <img src="https://render.githubusercontent.com/render/math?math=V"> from above is the common eigenvector matrix for all means. Use function `eigTest` with input of two lists and the convergence rate.
