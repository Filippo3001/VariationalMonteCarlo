using DrWatson
@quickactivate "VariationalMonteCarlo"

using DataFrames, Plots, StatsPlots

# Define a function which plots the results of the MonteCarlo simulations for the 1D Harmonic oscillator
# Only consider data given by arguments
function plot_onedharm(nMoves, nThermMoves, metroStep, startingPoint, alphamin, alphamax)
    #Load results of simulations data folder
    df = collect_results!(datadir("onedharm", "results_table.jld2"), datadir("onedharm", "results"))
    # Select only the data of interest
    subset!(df, :nMoves => nm -> nm .== nMoves, :nThermMoves => ntm -> ntm .== nThermMoves,
        :metroStep => ms -> ms .== metroStep, :startingPoint => sp -> sp .== startingPoint, :alpha => a -> a .<= alphamax, :alpha => a -> a .>= alphamin) 

    #Plot the results in a scatter plot and save it
    gr()
    plot = @df df scatter(:alpha, :mean; yerror=:sigma)
    savefig(plot, plotsdir("onedharm", "scatterplot"))

end

function plot_onedharm_corr(nMoves, nThermMoves, metroStep, startingPoint, starting_alpha)
    #Load results of simulations data folder
    df = collect_results!(datadir("onedharm", "corr_results_table.jld2"), datadir("onedharm", "corr_results"))
    # Select only the data of interest
    subset!(df)
end



#nMoves, nThermMoves, metroStep, startingPoint, alphamin, alphamax
#plot_onedharm(100000, 100, 1., 0., 0.8, 1.2)