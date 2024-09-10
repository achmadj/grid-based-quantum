import numpy as np
from pyquest.unitaries import H, X, Phase, Swap
from pyquest import Circuit
from pyquest.operators import PhaseFunc
from typing import List, Tuple


def gen_QFT_circ(
    start_idx: int, 
    end_idx: int, 
    centred=True) -> Circuit:
    """ 
    Performing QFT circuit across qubit[start_idx] to qubit[end_idx]
    """
    gate_sequence = []
    
    if centred:
        gate_sequence.append(X(end_idx))
    # QFT
    for j in range(end_idx, start_idx - 1, -1):
        gate_sequence.append(H(j))
        for k in range(j - 1, start_idx - 1, -1):
            phase = np.pi / (2**(j - k))
            gate_sequence.append(Phase(k, phase, controls=[j]))
    # Swaps
    for i in range(start_idx,
                    start_idx + ((end_idx-start_idx)+1)//2):
        gate_sequence.append(Swap([i,
                                    end_idx + start_idx -i]))    
    if centred:
        gate_sequence.append(X(end_idx))
    
    return Circuit(gate_sequence)


def gen_Quadratic_PhaseFunc(
    qubit_list: List[int],
    dx: float,
    coeff: float,
    xc: float) -> PhaseFunc:
    """
    Quadratic phase function
    exp(i*coeff*(n*dx-xc)^2) = exp(i*coeff*(dx^2*n^2 - 2n*dx*xc + xc^2))
    coeff would absorb any constants e.g. dt, k, 1/2m etc
    """
    terms = [(coeff*dx**2, 0, 2),
             (-2*coeff*dx*xc, 0, 1),
             (coeff*xc**2, 0, 0)]
    return PhaseFunc(targets=qubit_list, terms=terms)


def gen_Coulomb_PhaseFunc(
    target_regs,
    dx: float,
    coeff: float,
    xc: List,
    override: float) -> PhaseFunc:
    """
    Scaled shifted norm phase function
    exp(i*coeff*1/sqrt((x-xc)^2 + (y-yc)^2 + ...))
    
    override: override value for singularity
    """
    assert len(target_regs) < 4
    assert len(target_regs)==len(xc)

    scale = coeff/dx
    shifts = [coord/dx for coord in xc]
    
    return PhaseFunc(
            targets=target_regs,
            func_type="scaled_inverse_shifted_norm",
            scaling=scale,
            shifts=shifts,
            divergence_override=override)


def gen_pairwise_Coulomb_PhaseFunc(
    target_regs,
    coeff: float,
    dx: float
    ) -> PhaseFunc:
    """ Pairwise Coulombic interaction between particle registers
    """
    scale = coeff/dx
    PhaseFunc(targets=target_regs,
              func_type="scaled_inverse_distance",
              scaling=scale
              )
    




