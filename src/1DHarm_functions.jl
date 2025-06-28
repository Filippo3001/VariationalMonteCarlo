module OneDHarm

using DifferentiationInterface
import ForwardDiff
using LinearAlgebra
# Define the  function necesary to study the 1D Harmonic oscillator
function harmtrialstate(x, alpha)
    (sqrt(alpha) / pi^(0.25)) * exp(- 0.5 * alpha^2 * x^2)

end

function harmimpsampling(alpha::Float64)
    function inner(x)
        harmtrialstate(x, alpha)^2
    end
end

function harmHamiltonian(x, alpha)
    kinPart = - second_derivative(harmtrialstate, AutoForwardDiff(), x, Constant(alpha))
    potPart = x^2* harmtrialstate(x, alpha)
    kinPart + potPart
end

function harmlocalen(alpha::Float64)
    function inner(x)
        harmHamiltonian(x, alpha) / harmtrialstate(x, alpha)
    end
end

function harmlocalenan(alpha)
    function inner(x)
        alpha^2 + x^2 * (1 - alpha^4)
    end
end

function harmresultsan(alpha)
    function inner(x)

    end
end

end