from pyquest import Register

import numpy as np
from qarsen_lite.singleparticlefunctions import gwp_1D, harmonic_osc_1D
from qarsen_lite.unitaries import *

from tqdm import tqdm

# define spatial grid
L = 20 # box size is 15
n_qubits = 6 # per spatial dimension!

dx = L/(2**n_qubits)
x_grid = np.array([dx*n for n in range(int(2**n_qubits))])

# define wavefunction in 1D
psi_x1 = gwp_1D(x_grid, xc=L*0.3, pc=0.1, alpha=0.8, gamma=0)* np.sqrt(dx)
psi_x2 = gwp_1D(x_grid, xc=L*0.7, pc=-0.1, alpha=0.8, gamma=0)* np.sqrt(dx)


reg_init = Register(n_qubits + n_qubits) # 2x5 qubits!
reg_init[:] = np.tensordot(psi_x2, psi_x1, axes=0).flatten()


def gen_propagators(n_qubits, dx, xc, dt, m, k):
    """
    Specific for 2 particles in 1D harmonic potential, second order trotter - edit as necessary to generalise
    """
    # qubit targets
    init_qbit_idx = 0
    targets = []
    for _ in range(2):
        targets.append(list(range(init_qbit_idx, init_qbit_idx+n_qubits)))
        init_qbit_idx += n_qubits
    
    # QFT circuits
    cQFT_circs = [gen_QFT_circ(targets[d][0], targets[d][-1]) for d in range(2)]
    icQFT_circs = [circ.inverse for circ in cQFT_circs]
    # Kinetic propagator
    dp = 2*np.pi/(2**n_qubits * dx)
    pc = 2**n_qubits * dp / 2 # calculate the centre of the momentum space
    coeff = -dt/2*m
    K_props = [gen_Quadratic_PhaseFunc(targets[d], dp, 0.5*coeff, pc) for d in range(2)]

    # Potential 
    # - external harmonic potential
    coeff = -k*dt/2
    V_props = [gen_Quadratic_PhaseFunc(targets[d], dx, coeff, xc) for d in range(2)]
    # - pairwise interaction potential
    Z = 1
    coeff = Z*dt
    V_props.append(gen_pairwise_Coulomb_PhaseFunc([targets[1], targets[0]], coeff, dx))
    

    return targets, cQFT_circs, icQFT_circs, K_props, V_props


def apply_single_step(reg, qft_circs, iqft_circs, K_props, V_props):
    """
    """
    for qft in qft_circs:
        reg.apply_circuit(qft)
    for k in K_props:
        reg.apply_operator(k)
    for iqft in reversed(iqft_circs):
        reg.apply_circuit(iqft)
    
    for v in V_props:
        reg.apply_operator(v)

    for qft in qft_circs:
        reg.apply_circuit(qft)
    for k in K_props:
        reg.apply_operator(k)
    for iqft in reversed(iqft_circs):
        reg.apply_circuit(iqft)


def U_RTE(reg_init, n_qubits, dx, t_tot, n_steps, xc, m, k, sample_every=10):
    """
    """
    dt = t_tot/n_steps # time resolution
    reg = Register(copy_reg=reg_init) # copy initial register
    targets, Qft, iQft, K, V = gen_propagators(n_qubits, dx, xc, dt, m, k) # generate propagators
    
    print(targets)

    autocorrelation = []
    particle_densities = [[], []]
    print("Running SO-QFT propagation...")
    for n in tqdm(range(n_steps)):
        apply_single_step(reg, Qft, iQft, K, V)

        # saving data
        if n%sample_every==0:
            autocorrelation.append(reg*reg_init)
            for i, t_reg in enumerate(targets):
                particle_densities[i].append(reg.prob_of_all_outcomes(t_reg))
    print("Done!")

    return autocorrelation, particle_densities



t_tot = 50
n_steps = 400
sample_every= 2
t_grid = np.linspace(0, t_tot, int(n_steps/sample_every))

k = 0.2

auto, particles = U_RTE(
    reg_init, 
    n_qubits, 
    dx, 
    t_tot=t_tot, 
    n_steps=n_steps, 
    xc=L/2, 
    m=1, 
    k=k, 
    sample_every=sample_every)