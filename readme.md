# Efficient classical algorithms for simulating symmetric quantum systems

Classical simulation of Hamiltonians that are invariant under the permutation group. Here, we confirm the algorithm works by comparing to random examples which are small enough for exact diagonalization. The algorithm here is only used for confirming results and has not been optimized for runtime.

## Dependencies
The code is written completely in Python. We use the following packages:
* [numpy](https://numpy.org/)
* [scipy](https://scipy.org/)

All simulations were run on a laptop. Scaling to larger number of qubits may need a computer with more memory or faster cpu.

## Authors

* [Bobak Kiani](https://github.com/bkiani) (MIT) 
* Andreas Bauer
* [Eric Anschuetz](https://github.com/eanschuetz) (MIT)
