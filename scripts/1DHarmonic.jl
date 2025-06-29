using DrWatson
@quickactivate "VariationalMonteCarlo"

# Include modules from src
includet(srcdir("MonteCarlo.jl"))
includet(srcdir("1DHarm_functions.jl"))

# Start the plotting module
using Plots
gr()

# Choose an interval for the parameter alpha to analyze
alphas = collect(range(0.9, 1.1, 21))

# Decide all parameters for MonteCarlo method
nMoves = 100000
nThermMoves = 100
metroStep = 1.0
startingPoint = 0.0

# Initialize the vector in which to store points
#points = Vector{Float64}(undef, nMoves)
# Initialize the vectors in  which to store the results
means = Vector{Float64}()
sigmas = Vector{Float64}()

for alpha in alphas
    # Generate points for the given alpha
    points, roa = MonteCarlo.metropolis(OneDHarm.harmimpsampling(alpha), nMoves, nThermMoves, metroStep,  startingPoint)
    # Evaluate the results
    mean, sigma = MonteCarlo.evaluate(points, OneDHarm.harmlocalen(alpha))
    println(mean, "  +/-  ", sigma)
    push!(means, mean)
    push!(sigmas, sigma)
end

scatter(alphas, means; yerror = sigmas)