
module SymplecticDoubleCode
"""
Markus Grassl.
"Bounds on the minimum distance of linear codes and quantum codes."
Online available at http://www.codetables.de.
Accessed on 2026-05-20.
"""

"""
Discussed in https://arxiv.org/pdf/2509.15457:
CSS code with logical ops known
"""

import QECCore: AbstractCSSCode, code_n, code_k, parity_matrix_x, parity_matrix_z, distance

export SymplecticDouble

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
#code_k(c::sd30_6_5) = 6

#distance(c::sd30_6_5) = 5


end