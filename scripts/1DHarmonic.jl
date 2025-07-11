using Distributed
@everywhere using DrWatson

@everywhere begin
    @quickactivate "VariationalMonteCarlo"
    using Revise: includet
    includet(srcdir("simulations.jl"))
end

#using DrWatson
#@quickactivate "VariationalMonteCarlo"
#includet(srcdir("simulations.jl"))

function onedharm()

    # Choose an interval for the parameter alpha to analyze
    alphas = collect(range(0.9, 1.1, 21))

    # Decide all parameters for MonteCarlo method
    nMoves = 100000
    nThermMoves = 100
    metroStep = 1.0
    startingPoint = 0.0

    #Create the dictionaries for simulations
    d = Dict()
    d[:nMoves] = nMoves
    d[:nThermMoves] = nThermMoves
    d[:metroStep] = metroStep
    d[:startingPoint] = startingPoint
    d[:alpha] = alphas

    #Then create all possible input combinations
    dicts = dict_list(d)

    #Then evaluate for each combination
    pmap(simulation_onedharm, dicts)

end

function onedharm__corr()
    #Choose a parameter from which to generate the points
    starting_alpha = 1.0
    nPoints = 100
    alphastep = 0.005

    # Decide all parameters for MonteCarlo method
    nMoves = 100000
    nThermMoves = 100
    metroStep = 1.0
    startingPoint = 0.0

    #Create the dictionaries for simulations
    d = Dict()
    d[:nMoves] = nMoves
    d[:nThermMoves] = nThermMoves
    d[:metroStep] = metroStep
    d[:startingPoint] = startingPoint
    d[:starting_alpha] = starting_alpha
    d[:nPoints] = nPoints
    d[:alphastep] = alphastep

    dicts = dict_list(d)

    pmap(simulation_onedharm_corr, dicts)
    
end