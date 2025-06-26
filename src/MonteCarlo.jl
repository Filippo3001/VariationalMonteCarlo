module MonteCarlo

# Define a function which, give an importance sampling f(x) generate nMoves points using metropolis algorithm. Points can be either vectors of Float or Floats
function metropolis(impsampling, nMoves, nThermMoves, metroStep, startingPoint)
    acceptedPoints = Vector()
    point = copy(startingPoint)
    f_r = impsampling(point)
    naccepted = 0
    for i in 1:nMoves
        for j in 1:nThermMoves
            trialStep = point + (rand(Float64, size(point)) .- 0.5) .* metroStep
            newf_r = impsampling(trialStep)
            if newf_r > f_r || rand() < newf_r/f_r
                point = trialStep
                f_r =  newf_r
                naccepted += 1
            end
        end
        push!(acceptedPoints, point)
    end
    rateofacceptance = naccepted / (nMoves * nThermMoves)
    return acceptedPoints, rateofacceptance
end
    
end