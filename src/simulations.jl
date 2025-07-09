using DrWatson
@quickactivate "VariationalMonteCarlo"

includet(srcdir("MonteCarlo.jl"))
includet(srcdir("1DHarm_functions.jl"))
includet(srcdir("He_functions.jl"))

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

    points_savename = savename("onedharm_points", data, "jld2")
    results_savename = savename("onedharm_results", data, "jld2")

    #Save all points
    #io = open(datadir("onedharm", "generated_points", points_savename), "w") do io
    #    for x in points
    #        println(io, x)
    #    end
    #end
    wsave(datadir("onedharm", "points", points_savename), @strdict alpha nMoves nThermMoves metroStep startingPoint points)
    wsave(datadir("onedharm", "results", results_savename), @strdict alpha nMoves nThermMoves metroStep startingPoint mean sigma roa)


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

end