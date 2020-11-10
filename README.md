# eigTest

## Intro & Setup
Package is developed for eigenvector testing.

To install, use the code:
`devtools::install_github('XycYuchenXu/eigTest')`

## Usage
The package has the following functionalities:

1. Check whether the means of two random square matrices (`matrixList = list(matrixA, matrixB)`) with dimension `d`-by-`d` are simultaneously diagonalizable, either considering the commutator of the two matrices (use function `commutatorTest`), or using the log-likelihood ratio test framework (`projTest`). The random matrices are supposed to be consistent estimators of the means with asymptotic normal distribution with limiting covariance matrices (`covMatList = list(covMatA, covMatB)`) and a convergence rate (`cn`). For both functions, one needs to input a list of two matrices (`matrixList`), a list of two limiting covariance matrices (`covMatList`), and the convergence rate (`cn`).

2. Given a list of square random matrices (`matrixList`), estimate a common eigenvector matrix `V` (by supplying the list of matrices `matrixList` for function `JDTE`).

3. Test whether the list of means matrices (`matrixList`, list length `p` ≥ 2) for random matrices are simultaneously diagonalizable, by checking whether the eigenvector matrix estimated above (`V`) is indeed the common eigenvector matrix. Use function `eigTest` with input of two lists (`matrixList` and `covMatList`) and the convergence rate `cn`.

4. Given a list of square random matrices (`matrixList`), estimate `k` (`k` ≤ `d`) common eigenvectors `Vk` (by supplying the list of matrices `matrixList` and parameter `k` for function `expmPartSchur`).

5. Test whether the list of the random estimates (`matrixList`, list length `p` ≥ 2) for mean matrices share `k` of the all `d` eigenvectors (`k` ≤ `d`), by checking whether the `k` estimated common eigenvectors `Vk` from optimization above are indeed the common eigenvectors. Use function `partialTest` with input of two lists (`matrixList` and `covMatList`), the convergence rate `cn`, and number of common eigenvectors `k`.

6. Generate gaussian samples for simulations. The function `generateMeans` is used to generate `p` mean matrices with dimension `d` and noise level `snr` of common eigenvector `V`. Function `simuSamples` takes the output from `generateMeans` as input for mean matrices. For each mean matrix, it generates `cn^2` gaussian random matrices (with identity covariance matrix), and compute the sample mean and the sample covariance matrix as estimators. It then combines those estimates into two lists, `mu.bar` for means, and `cov.bar` for covariances.

7. For more details, read the vignette: `browseVignettes('XycYuchenXu/eigTest')`
