# symplectic_double_code.jl

"""
Symplectic Double Code implementation, based on https://arxiv.org/pdf/2509.15457.
"""
module SymplecticDoubleCode

export SymplecticDouble

import QECCore: AbstractCSSCode, code_n, parity_matrix_x, parity_matrix_z


struct SymplecticDouble <: AbstractCSSCode end

_Hx = Bool[
    1 0 0 0 0 0 0 0 0 1 0 1 0 1 0   
    0 0 0 0 0 0 0 1 1 1 1 1 0 0 1   
    0 1 0 0 0 0 0 1 1 1 1 0 0 1 0   
    0 0 0 0 0 0 1 1 1 1 1 1 1 0 1   
    0 0 1 0 0 0 0 1 0 0 1 1 1 1 0   
    0 0 0 0 0 0 1 0 1 1 1 1 1 1 1  
    0 0 0 1 0 0 0 1 0 1 0 1 0 0 0   
    0 0 0 0 0 0 1 0 0 1 1 1 1 1 0   
    0 0 0 0 1 0 0 0 1 0 1 0 1 0 0   
    0 0 0 0 0 0 0 1 0 0 1 1 1 1 1  
    0 0 0 0 0 1 0 1 1 0 0 1 1 0 1   
    0 0 0 0 0 0 1 1 1 0 0 1 1 1 0 
]

_Hz = Bool[
    0 0 0 0 0 0 0 1 1 1 1 1 0 0 1
    1 0 0 0 0 0 0 1 1 0 1 0 0 1 1
    0 0 0 0 0 0 1 1 1 1 1 1 1 0 1
    0 1 0 0 0 0 1 0 0 0 0 1 1 1 1
    0 0 0 0 0 0 1 0 1 1 1 1 1 1 1
    0 0 1 0 0 0 1 1 1 1 0 0 0 0 1
    0 0 0 0 0 0 1 0 0 1 1 1 1 1 0
    0 0 0 1 0 0 1 1 0 0 1 0 1 1 0
    0 0 0 0 0 0 0 1 0 0 1 1 1 1 1
    0 0 0 0 1 0 0 1 1 0 0 1 0 1 1
    0 0 0 0 0 0 1 1 1 0 0 1 1 1 0
    0 0 0 0 0 1 1 0 0 0 0 0 0 1 1
]

parity_matrix_x(c::SymplecticDouble) = hcat(_Hx, _Hz)
parity_matrix_z(c::SymplecticDouble) = hcat(_Hz, _Hx)
code_n(c::SymplecticDouble) = 30

end