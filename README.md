# eigTest
Eigenvector testing.

To install, use the code:
`devtools::install_github('XycYuchenXu/eigTest')`

The package has the following functionalities:

1. Check whether the means of two random square matrices (`matrixList = list(matrixA, matrixB)`) with dimension `d`-by-`d` are simultaneously diagonalizable, either considering the commutator of the two matrices (use function `commutatorTest`), or using the log-likelihood ratio test framework (`projTest`). The random matrices are supposed to be consistent estimators of the means with asymptotic normal distribution with limiting covariance matrices (`covMatList = list(covMatA, covMatB)`) and a convergence rate (`cn`). For both functions, one needs to input a list of two matrices (`matrixList`), a list of two limiting covariance matrices (`covMatList`), and the convergence rate (`cn`).

2. Given a list of square random matrices (`matrixList`), estimate a common eigenvector matrix `V` (by supplying the list of matrices for function `JDTE`).

3. Test whether the list of means matrices (`matrixList`, list length `p` ≥ 2) for random matrices are simultaneously diagonalizable, by checking whether the eigenvector matrix estimated above (`V`) is indeed the common eigenvector matrix. Use function `eigTest` with input of two lists (`matrixList` and `covMatList`) and the convergence rate `cn`.

4. Test whether the list of mean matrices (`matrixList`, list length `p` ≥ 2) share `k` of the all `d` eigenvectors (`k` ≤ `d`), by checking whether the `k` estimated common eigenvectors from optimization function `expmPartSchur` (with input `matrixList` and parameter `k`) are indeed the common eigenvectors. Use function `partialTest` with input of two lists (`matrixList` and `covMatList`), the convergence rate `cn`, and number of common eigenvectors `k`.

5. Generate gaussian samples for simulations. The function `generateMeans` is used to generate `p` mean matrices with dimension `d` and noise level `snr`. Function `simuSamples` takes the output from `generateMeans` as input for mean matrices. For each mean matrix, it generates gaussian samples (with identity covariance matrix) of size `cn^2` to get a consistent mean estimators (`mu.bar`) and the covariance matrix estimator. It then combines those estimates into two lists, `mu.bar` for means, and `cov.bar` for covariances.
