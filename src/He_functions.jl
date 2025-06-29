module HeNucleus

using Enzyme

# Define first the Jastrow factors to then construct the trial state function.  par = [a, β, γ], r is a number
function Jastrow_factor(par::Vector{Float64}, r::Float64)
    exp( - par[3] * r^2) + par[1] * exp(- (par[2] + par[1]) * r^2)
end

# Define the distance between two vectors
function vecdistance(r1, r2)
    @assert(length(r1) == length(r2)) 
    sqrt(sum((r1[i] - r2[i])^2 for i in eachindex(r1)))
end


# The He⁴ nucleus model has 12 dof represented by the positon vector r = [r₁,r₂,r₃,r₄]; where every rᵢ is a 3D vector which represents the position of the i-nth nucleons
# Since internucleus distances are used multiple times, define a function which returns them
function internucleusd(r)
    r1 = r[1:3]
    r2 = r[4:6]
    r3 = r[7:9]
    r4 = r[10:12]

    return vecdistance(r1, r2), vecdistance(r1, r3), vecdistance(r1, r4), vecdistance(r2, r3), vecdistance(r2, r4), vecdistance(r3, r4)
end

# Define the trial state function. par = [a, β, γ], r is a 12 dimensional vector which represents the positions coordinates of the 4 nucleons
function Hetrialfunction(par::Vector{Float64}, r::Vector{Float64})
    # Verify the input are of the correct size
    @assert(length(par) == 3)
    @assert(length(r) == 12)
    # Evaluate the internucleons distances
    r_12, r_13, r_14, r_23, r_24, r_34 = internucleusd(r)

    (Jastrow_factor(par, r_12) * Jastrow_factor(par, r_13) 
    * Jastrow_factor(par, r_14) * Jastrow_factor(par, r_23) 
    * Jastrow_factor(par, r_24) * Jastrow_factor(par, r_34))

end

# Define a function which return the i-nth standard basis vector of dimension N
function stdbasis(n, i)
    e_i = zeros(n)
    e_i[i] = 1.
    e_i
end

# Define the potential.
# Only 2 bodies interaction are considered
# Consider a S3 potential dependant only on the internucleus distance
function S3Pot(r)
    (1000 * exp(-3  *  r^2) - 163.5 * exp(-1.05  *  r^2) - 21.5 * exp(-0.6  *  r^2)
    - 83 * exp(-0.8  *  r^2) - 11.5 * exp(-0.4  *  r^2))
end

# Define the potential for the He⁴ HeNucleus
function HeS3Pot(r)
    @assert(length(r) == 12)
    # Evaluate the internucleons distances
    distances = internucleusd(r)
    # The potential is the sum of the potential for each internucleons distance (counted only one time)
    res = 0.
    for d in distances
        res += S3Pot(d)
    end
    res
end


function Laplacian_Hetrialfunction(par::Vector{Float64}, r::Vector{Float64})
    # Verify the input are of the correct size
    @assert(length(par) == 3)
    @assert(length(r) == 12)
    # Evaluate the Laplacian using automatic differentiation from Enzyme.jl.
    # Evaluate each row of the Hessian and consider only the diagonal elements δ²f/δx²
    lapl = 0.
    row = Vector{Float64}(undef, length(r))
    for i in eachindex(r)
        # Evaluate the Hessian vector product with the i_nth standard basis vector.
        # So 
        hvp!(row, x -> Hetrialfunction(par, x), r, stdbasis(length(r), i))
        # Accumulate the Laplacian
        lapl += row[i]
    end
    lapl
end

function Hehamiltonian(par::Vector{Float64}, r::Vector{Float64})
    # Verify the input are of the correct size
    @assert(length(par) == 3)
    @assert(length(r) == 12)
    # Evaluate the kinetic kinetic part
    kinPart = - Laplacian_Hetrialfunction(par, r)
    potPart = HeS3Pot(r) * Hetrialfunction(par, r)
    kinPart + potPart
end

# Construct the local energy
function Helocalen(par, r)
    @assert(length(par) == 3)
    @assert(length(r) == 12)

    Hehamiltonian(par, r) / Hetrialfunction(par, r)
end

end
