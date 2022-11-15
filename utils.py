import numpy as np
import itertools
from sympy.utilities.iterables import multiset_permutations
import math
from scipy.special import comb

data_type = np.complex128
I = np.array([[1,0],[0,1]]).astype(data_type)
X = np.array([[0,1],[1,0]]).astype(data_type)
Y = np.array([[0,-1j],[1j,0]]).astype(data_type)
Z = np.array([[1,0],[0,-1]]).astype(data_type)

pauli_dict = {'I':I, 'X':X, 'Y':Y, 'Z':Z}

fac_fun = np.math.factorial

def get_pauli_matrix(string):
	out = np.array([[1.]]).astype(data_type)
	for p in string:
		out = np.kron(out, pauli_dict[p])
	return out

def get_symmetric_pauli_term(n,nx,ny,nz, normalize = False):
	if nx+ny+nz > n:
		raise valueError('too many Paulis')
	ni = n - nx - ny - nz

	strings = ['I']*ni + ['X']*nx + ['Y']*ny + ['Z']*nz
	strings = list(multiset_permutations(strings))
	out = 0
	for st in strings:
		out += get_pauli_matrix(st)
	if normalize:
		out /= np.sqrt(len(strings))
	return out
	
def construct_matrix_large(n,paulis,coeffs):
	out = np.zeros((2**n,2**n)).astype(data_type)
	for p,c in zip(paulis,coeffs):
		out += c*get_symmetric_pauli_term(n,*p)
	return out



def sn_spin_transfer_term(n, mu, d1, d2, ix, iy, iz):
	i1 = n - ix - iy - iz

	def sn_spin_iterator():
		for f11 in range(i1+1):
			for g010 in range(i1+1):
				for g111 in range(i1+1):
					for fxx in range(ix+1):
						for g0x1 in range(ix+1):
							for g1x0 in range(ix+1):
								for fyy in range(iy+1):
									for g0y1 in range(iy+1):
										for g1y0 in range(iy+1):
											for fzz in range(iz+1):
												for g0z0 in range(iz+1):
													for g1z1 in range(iz+1):
														checks = []
														checks.append( (g010 + g0z0 + g0x1 + g0y1) == (n-2*mu-d1) )
														checks.append( (g010 + g0z0 + g1x0 + g1y0) == (n-2*mu-d2) )
														checks.append( (g111 + g1z1 + g1x0 + g1y0) == d1 )
														checks.append( (g111 + g1z1 + g0x1 + g0y1) == d2 )
														checks.append( (2*f11 + g010 + g111) == i1 )
														checks.append( (2*fxx + g0x1 + g1x0) == ix )
														checks.append( (2*fyy + g0y1 + g1y0) == iy )
														checks.append( (2*fzz + g0z0 + g1z1) == iz )
														passed = True
														for check in checks:
															if not check:
																passed = False
														if passed:
															yield f11,g010,g111,fxx,g0x1,g1x0,fyy,g0y1,g1y0,fzz,g0z0,g1z1

	iterator = sn_spin_iterator()
	term = 0
	for f11,g010,g111,fxx,g0x1,g1x0,fyy,g0y1,g1y0,fzz,g0z0,g1z1 in iterator:
		factor = fac_fun(mu)*fac_fun(n-2*mu)
		for value in [f11,g010,g111,fxx,g0x1,g1x0,fyy,g0y1,g1y0,fzz,g0z0,g1z1]:
			factor /= fac_fun(value)
		prefac = (2*fxx+2*fyy+2*fzz+2*g1z1-g0y1+g1y0) % 4# (2*f11+2*fxx+2*fyy+2*fzz+2*g1z1+g0y1-g1y0) % 4
		term += factor*(1j)**prefac
	return term/np.sqrt(comb(n-2*mu, d1)*comb(n-2*mu, d2))



def get_matrices(n,nx,ny,nz, c=1):
	n_blocks = n//2+1
	out = []
	for block_i in range(n_blocks):
		block_size = n-block_i*2+1
		mat = np.zeros((block_size, block_size)).astype(data_type)
		for i in range(block_size):
			for j in range(block_size):
				mat[i,j] = sn_spin_transfer_term(n, block_i, i, j, nx, ny, nz)
		out.append(mat*c)
	return out


def add_matrix_blocks(ms):
	out = ms[0]
	ms = ms[1:]
	for m in ms:
		for i in range(len(m)):
			out[i]+=m[i]
	return out

def construct_matrix_blocks(n,paulis,coeffs):
	ms = [get_matrices(n, *p,c=c) for p,c in zip(paulis,coeffs)]
	return add_matrix_blocks(ms)


def dicke(nn,dd):
	def binary_to_int(st):
		out = 0
		for i,s in enumerate(st):
			out += (2**i)*s
		return out

	strings = [0]*(nn-dd) + [1]*dd
	combos = list(multiset_permutations(strings))
	out = np.zeros(2**nn).astype(data_type)
	for s in combos:
		out[binary_to_int(s)] = 1.
	return out / np.sqrt(len(combos))

def state_constructor(n_sing,local_state):
	def singlet():
		return np.asarray([0,1,-1,0]).astype(data_type)/np.sqrt(2)
	def all_dicke(nn):
		out = np.zeros((2**nn,nn+1)).astype(data_type)
		for i in range(nn+1):
			out[:,i] = dicke(nn,i)
		return out

	state = np.ones((1)).astype(data_type)
	for _ in range(n_sing):
		state = np.kron(state, singlet())
	if len(local_state) == 1:
		return state
	else:
		# print('dickes')
		dickes = all_dicke(len(local_state)-1)
		# print(dickes)
		# print(local_state)
		local_state = dickes@local_state
		# print(local_state)
		return np.kron(state,local_state.reshape(-1))
