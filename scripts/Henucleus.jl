using Distributed
@everywhere using DrWatson

@everywhere begin
    @quickactivate "VariationalMonteCarlo"
    using Revise: includet
    includet(srcdir("simulations.jl"))
end

function Henucleus()

    # Choose an interval for the parameter alpha to analyze
    a = collect(range(-0.8, -0.6, 11))
    β = collect(range(1.7, 2.3, 11))
    γ = collect(range(0.07, 0.12, 11))

    # Decide all parameters for MonteCarlo method
    nMoves = 10000
    nThermMoves = 100
    metroStep = 1.0
    startingPoint = [zeros(12)]

    #Create the dictionaries for the simulationss
    d = Dict()
    d = Dict()
    d[:nMoves] = nMoves
    d[:nThermMoves] = nThermMoves
    d[:metroStep] = metroStep
    d[:startingPoint] = startingPoint
    d[:a] = a
    d[:β] = β
    d[:γ] = γ

    #Then create all possible input combinations
    dicts = dict_list(d)

    #Then evaluate for each combination
    pmap(simulation_Henucleus, dicts)

end
