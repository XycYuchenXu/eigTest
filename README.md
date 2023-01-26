# eigTest: Jointly Estimate and Test for Common Eigenvectors

## Intro & Setup
This `R` package is developed for testing simultaneous diagonalizability.

To install, use the code in `R`:
`devtools::install_github('XycYuchenXu/eigTest', force = T, build_vignettes = T)`

## Usage
The package has the following functionalities:

1. Check whether the means of two random square matrices (`matrixA` and `matrixB`) with dimension `d`-by-`d` are simultaneously diagonalizable, either considering the commutator of the two matrices (use function [`commutatorTest`](R/commutatorTest.R)), or using the log-likelihood ratio test framework ([`projTest`](R/projTest.R)), given asymptotic limiting covariance matrices (`covMatA` and `covMatB`) and a convergence rate (`cn`). For both functions, one needs to input an array of two matrices (`A` such that `A[1,,] = matrixA` and `A[2,,] = matrixB`), an array of two limiting covariance matrices (`cov.arr` such that `cov.arr[1,,] = covMatA` and `cov.arr[2,,] = covMatB`), and the convergence rate (`cn`).

2. Given an array of square random matrices (`A`), estimate a common eigenvector matrix `V` (by supplying the array of matrices `A` for function [`JDTE`](R/JDTE.R)).

3. Test whether the array of means matrices (`A` with dimension `p`-by-`d`-by-`d` and `p` ≥ 2) for random matrices are simultaneously diagonalizable, by checking whether the eigenvector matrix estimated above (`V`) is indeed the common eigenvector matrix. Use function [`eigTest`](R/eigTest.R) with input of two arrays (`A` and `cov.arr`) and the convergence rate `cn`.

4. Given an array of square random matrices (`A` with dimension `p`-by-`d`-by-`d`), estimate `k` (`k` ≤ `d`) common eigenvectors `Vk` (by supplying the array of matrices `A` and parameter `k` for function [`expmPartSchur`](R/expmPartSchur.R)).

5. Test whether the array of the random estimates (`A` with dimension `p`-by-`d`-by-`d` and `p` ≥ 2) for mean matrices share `k` of the all `d` eigenvectors (`k` ≤ `d`), by checking whether the `k` estimated common eigenvectors `Vk` from optimization above are indeed the common eigenvectors. Use function [`partialTest`](R/partialTest.R) with input of two arrays (`A` and `cov.arr`), the convergence rate `cn`, and number of common eigenvectors `k`. `Vk` is optional here as `partialTest` can call the optimization function.

6. Generate Gaussian samples for simulations. The function [`generateMeans`](R/generateMeans.R) is used to generate `p` mean matrices with dimension `d`-by-`d` and noise levels `snr` of common eigenvectors `V`. Function [`simuSamples`](R/simuSamples.R) takes the output from [`generateMeans`](R/generateMeans.R) as input for mean matrices. For each mean matrix, it generates `cn^2` Gaussian random matrices (with identity covariance matrix), and compute the sample mean and the sample covariance matrix as consistent estimators. It then combines those estimates into a list of sub-lists that include estimated mean `mu.bar` and covariance `cov.bar`.

7. For more details, read the vignette: `browseVignettes('eigTest')` and the paper ['Testing Simultaneous Diagonalizability'](https://arxiv.org/abs/2101.07776).
