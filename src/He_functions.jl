using DrWatson, Test
@quickactivate "VariationalMonteCarlo"
module HeNucleus

using Enzyme

# Define first the Jastrow factors to then construct the trial state function.  par = [a, β, γ], r is a number
function Jastrow_factor(par::Vector{Float64}, r2::Float64)
    exp(-par[3] * r2) + par[1] * exp(-(par[2] + par[3]) * r2)
end

# Define the distance between two vectors
function vecdistance2(r1, r2)
    @assert(length(r1) == length(r2))
    sum((r1[i] - r2[i])^2 for i in eachindex(r1))
end


# The He⁴ nucleus model has 12 dof represented by the positon vector r = [r₁,r₂,r₃,r₄]; where every rᵢ is a 3D vector which represents the position of the i-nth nucleons
# Since internucleus distances are used multiple times, define a function which returns them
function internucleusd(r)
    # Separate from r the 4 position vectors, using @view to not allocate
    r1 = @view r[1:3]
    r2 = @view r[4:6]
    r3 = @view r[7:9]
    r4 = @view r[10:12]

    return vecdistance2(r1, r2), vecdistance2(r1, r3), vecdistance2(r1, r4), vecdistance2(r2, r3), vecdistance2(r2, r4), vecdistance2(r3, r4)
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

function Heimpsampling(par::Vector{Float64})
    function inner(r::Vector{Float64})
        Hetrialfunction(par, r)^2
    end
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
function S3Pot(r2)
    (1000.0 * exp(-3.0 * r2) - 163.5 * exp(-1.05 * r2) - 21.5 * exp(-0.6 * r2)
     -
     83.0 * exp(-0.8 * r2) - 11.5 * exp(-0.4 * r2))
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
    # Evaluate each row of the Hessian and consider only the diagonal elements ∂²f/∂x² 
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
    h2d2m = 20.74
    kinPart = -h2d2m * Laplacian_Hetrialfunction(par, r)
    potPart = HeS3Pot(r) * Hetrialfunction(par, r)
    kinPart + potPart
end

# Construct the local energy
function Helocalen(par)
    @assert(length(par) == 3)
    function inner(r)
        @assert(length(r) == 12)
        Hehamiltonian(par, r) / Hetrialfunction(par, r)
    end
end

# Define a function which acts on r = (r1,r2,r3,r4) and return the centerofmass
function centerofmass(r)
    # Separate the 4 positions vectors
    r1 = @view r[1:3]
    r2 = @view r[4:6]
    r3 = @view r[7:9]
    r4 = @view r[10:12]
    cm = Vector{Float64}(undef, 3)
    for i in eachindex(cm)
        cm[i] = 0.25 * (r1[i] + r2[i] + r3[i] + r4[i])
    end
    cm
end

# Define a function which recenter r =(r1,r2,r3,r4) around its center of mass
function recenter!(r)
    # Separate the 4 positions vectors
    r1 = @view r[1:3]
    r2 = @view r[4:6]
    r3 = @view r[7:9]
    r4 = @view r[10:12]
    cm = centerofmass(r)
    for i in eachindex(cm)
        r1[i] -= cm[i]
        r2[i] -= cm[i]
        r3[i] -= cm[i]
        r4[i] -= cm[i]
    end
end


# Try with finite differences method to see if it gives the same results
function Laplacian_Hetrialfunctiondiff(par::Vector{Float64}, r::Vector{Float64}, h::Float64)
    # Verify the input are of the correct size
    @assert(length(par) == 3)
    @assert(length(r) == 12)
    #Evaluate the laplacian using finite difference differentiation
    #First define and allocate the two vectors rp and Random
    rp = deepcopy(r)
    rm = deepcopy(r)
    #Define the variable in which we accumulate the laplacian
    lapl = 0.
    for i in eachindex(r)
        #Calculate rp and rm for the i-th coordinate
        rp[i] = r[i] + h
        rm[i] = r[i] - h
        #Evaluate the second derivative and sum it to the laplacian
        second_deriv = (Hetrialfunction(par, rp) +  Hetrialfunction(par, rm) - 2*Hetrialfunction(par, r)) / (h^2)
        lapl += second_deriv
        #reset rp and rm to r
        rp[i] = r[i]
        rm[i] = r[i]
    end
    lapl
end

function Hehamiltoniandiff()
    
end

end
