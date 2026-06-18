# trivariate_bicycle_code.jl

"""
Implementation of the trivariate bicycle code as discussed in https://arxiv.org/pdf/2406.19151 and https://arxiv.org/abs/2601.18879.
This implementation is based on the bivariate case, implemented in 
https://github.com/QuantumSavory/QuantumClifford.jl/blob/2b6cee1f2db6b2c781c4b0d168b9c07b10d47888/lib/QECCore/src/codes/quantum/generalized_circulant_bivariate_bicycle.jl.

We define three variates `x`, `y` and `z` based on parameters `l` and `m` as 
`x = S_l ⊗ I_m`, `y = I_l ⊗ S_m` and `z = S_l ⊗ S_m`, where `S_i` denotes the cyclic shift matrix of dimension `i`, 
and `I_i` is the `i×i` identity matrix. The matrices `A` and `B` are polynomials in those variables, from which we form the CSS code
stabiliser matrix with `H_X = [A|B]` and `H_Z = [B^t|A^t]`.
"""
module TrivariateBicycleCode

export TrivariateBicycleViaCirculantMat

using LinearAlgebra

import QECCore: AbstractCSSCode, code_n, parity_matrix_x, parity_matrix_z


struct TrivariateBicycleViaCirculantMat <: AbstractCSSCode
    """Dimension of cyclic shift matrix `Sₗ` where `x = Sₗ ⊗ Iₘ`"""
    l::Int
    """ Dimension of cyclic shift matrix `Sₘ` where `y = Iₗ ⊗ Sₘ`"""
    m::Int
    """Terms in matrix A, where each tuple is (:x or :y, power)"""
    A::Vector{Tuple{Symbol,Int}}
    """Terms in matrix B, where each tuple is (:x or :y, power)"""
    B::Vector{Tuple{Symbol,Int}}

    function TrivariateBicycleViaCirculantMat(l, m, A, B)
        (l >= 0 && m >= 0) || throw(ArgumentError("l and m must be non-negative"))
        (length(A) >= 1 && length(B) >= 1) || throw(ArgumentError("A and B must each have at least one entry"))
        z_period = lcm(l,m)
        for (mat, terms) in [(:A, A), (:B, B)]
            for (var, pow) in terms
                var ∈ [:x, :y, :z] || throw(ArgumentError("Matrix $mat contains invalid variable $var (must be :x or :y or :z)"))
                pow >= 0 || throw(ArgumentError("Matrix $mat contains negative power $pow"))
                max_pow = var == :x ? (l-1) : (var == :y ? (m-1) : (z_period-1))
                # x^l = 1, so x is l-dimensional and y is m-dimensional
                pow <= max_pow || throw(ArgumentError("Power $pow in matrix $mat exceeds maximum $max_pow for $var"))
            end
        end
        new(l, m, A, B)
    end
end

function parity_matrix_xz(c::TrivariateBicycleViaCirculantMat)
    Iₗ = Matrix{Bool}(I, c.l, c.l)
    Iₘ = Matrix{Bool}(I, c.m, c.m)
    xₚ = Dict(i => kron(circshift(Iₗ, (0,i)), Iₘ) for i in 0:(c.l-1))
    yₚ = Dict(i => kron(Iₗ, circshift(Iₘ, (0,i))) for i in 0:(c.m-1))
    z_period = lcm(c.l,c.m)
    zₚ = Dict(i => kron(circshift(Iₗ, (0,i)), circshift(Iₘ, (0,i))) for i in 0:(z_period-1))
    A = zeros(Bool, c.l*c.m, c.l*c.m)
    for (var, pow) in c.A
        mat = var == :x ? xₚ[pow] : (var ==:y ? yₚ[pow] : zₚ[pow] ) 
        A .+= mat
    end
    A = mod.(A, 2)
    B = zeros(Bool, c.l*c.m, c.l*c.m)
    for (var, pow) in c.B
        mat = var == :x ? xₚ[pow] : (var ==:y ? yₚ[pow] : zₚ[pow] ) 
        B .+= mat
    end
    B = mod.(B, 2)
    Hx = hcat(A, B)
    Hz = hcat(B', A')
    return Hx, Hz
end

code_n(c::TrivariateBicycleViaCirculantMat) = 2*c.l*c.m

parity_matrix_x(c::TrivariateBicycleViaCirculantMat) = parity_matrix_xz(c)[1]

parity_matrix_z(c::TrivariateBicycleViaCirculantMat) = parity_matrix_xz(c)[2]

end