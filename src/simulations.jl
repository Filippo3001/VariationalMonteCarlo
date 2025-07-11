using DrWatson
@quickactivate "VariationalMonteCarlo"

includet(srcdir("MonteCarlo.jl"))
includet(srcdir("1DHarm_functions.jl"))
includet(srcdir("He_functions.jl"))

using DataFrames

#Define a main-like function which takes a dictionary of parameters (data) and does a MonteCarlo simulation for the 1D harmonic oscillator
#It will also output all results in a file directory
function simulation_onedharm(data)
    #Unpack data in the various simulation parameters
    nMoves = data[:nMoves]
    nThermMoves = data[:nThermMoves]
    metroStep = data[:metroStep]
    startingPoint = data[:startingPoint]
    alpha = data[:alpha]

    #Generate the points
    points, roa = MonteCarlo.metropolis(OneDHarm.harmimpsampling(alpha), nMoves, nThermMoves, metroStep, startingPoint)
    # Evaluate the results
    mean, sigma = MonteCarlo.evaluate(points, OneDHarm.harmlocalen(alpha))
    println(mean, "  +/-  ", sigma)

    #points_savename = savename("onedharm_points", data, "jld2")
    results_savename = savename("onedharm_results", data, "jld2")

    #wsave(datadir("onedharm", "points", points_savename), @strdict alpha nMoves nThermMoves metroStep startingPoint points)
    wsave(datadir("onedharm", "results", results_savename), @strdict alpha nMoves nThermMoves metroStep startingPoint mean sigma roa)

end

function simulation_onedharm_corr(data)
    #Unpack data in the various simulation parameters
    nMoves = data[:nMoves]
    nThermMoves = data[:nThermMoves]
    metroStep = data[:metroStep]
    startingPoint = data[:startingPoint]
    starting_alpha = data[:starting_alpha]
    nPoints = data[:nPoints]
    alphastep = data[:alphastep]
    alphas = [starting_alpha + (i - nPoints / 2) *  alphastep for i in nMoves]
    
    #Generate the points
    points, roa = MonteCarlo.metropolis(OneDHarm.harmimpsampling(starting_alpha), nMoves, nThermMoves, metroStep, startingPoint)
    means = Vector{Float64}(undef, length(alphas))
    sigmas = Vector{Float64}(undef, length(alphas))
    for a in alphas
        mean, sigma = MonteCarlo.corr_evaluate(points, OneDHarm.harmimpsampling(alpha), OneDHarm.harmimpsampling(a), OneDHarm.harmlocalen(a))
        push!(means, mean)
        push!(sigmas, sigma)
    end

    results_savename = savename("onedharm_corr_results", data, "jld2")
    wsave(datadir("onedharm","corr_results", results_savename), @strdict starting_alpha nPoints alphastep nMoves nThermMoves metroStep startingPoint roa means sigmas)
end

#
function simulation_Henucleus(data)
    nMoves = data[:nMoves]
    nThermMoves = data[:nThermMoves]
    metroStep = data[:metroStep]
    startingPoint = data[:startingPoint]
    a = data[:a]
    β = data[:β]
    γ = data[:γ]
    par = [a, β, γ]

    points, roa = MonteCarlo.metropolis(HeNucleus.Heimpsampling(par), nMoves, nThermMoves, metroStep, startingPoint)

    mean, sigma = MonteCarlo.evaluate(points, HeNucleus.Helocalen(par))
    println(mean, "  +/-  ", sigma)

    #points_savename = savename("Henucleus_points", data, "jld2")
    results_savename = savename("Henucleus_results", data, "jld2")

    #wsave(datadir("Henucleus", "points", points_savename), @strdict a β γ nMoves nThermMoves metroStep startingPoint points)
    wsave(datadir("Henucleus", "results", results_savename), @strdict a β γ nMoves nThermMoves metroStep startingPoint mean sigma roa)

end