
# Efficient classical algorithms for simulating symmetric quantum systems

This code classically simulates Hamiltonians on qubit systems that are invariant under the permutation group. The code works by comparing to random examples which are small enough for exact diagonalization. 

*As a disclaimer, though this code runs in polynomial time in the number of qubits, it has not been optimized to have the actual scaling reported in the paper.* 

For further details, please refer to the corresponding paper: [https://arxiv.org/abs/2211.16998](https://arxiv.org/abs/2211.16998)

## Dependencies
The code is written completely in Python. We use the following packages:
* [numpy](https://numpy.org/)
* [scipy](https://scipy.org/)

All simulations were run on a laptop. Scaling to larger number of qubits may need a computer with more memory or faster cpu.

## Authors

* [Bobak Kiani](https://github.com/bkiani) (MIT) 
* Andreas Bauer
* [Eric Anschuetz](https://github.com/eanschuetz) (MIT)
