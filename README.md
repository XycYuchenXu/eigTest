# eigTest
Eigenvector testing.

To install, use the code:
`devtools::install_github('XycYuchenXu/eigTest')`

The package has the following functionalities:

1. Check the means (<img src="https://render.githubusercontent.com/render/math?math=\mu_1, \mu_2">) of two random matrices ($$A_1, A_2$$) are simultaneously diagonalizable, either considering the commutators of the two matrices, or using the log-likelihood ratio test framework. The random matrices $$A_1$$ and $$A_2$$ are supposed to be consistent estimators of the means with asymptotic normal distribution.

2. Given square matrices $$A_1, \dots, A_p$$, estimate the matrix $$V$$ such that for every $$A_i$$, $$V^{-1} A_i V = D_i$$ with diagonal matrices $$D_i$$.

2. Test the means ($$\mu_i, ~ i = 1, \dots, p$$) corresponding to random matrices ($$A_i, ~ i = 1, \dots, p$$) are simultaneously diagonalizable, by checking whether the estimated $$V$$ from above is the common eigenvector matrix for all means.
