using DrWatson, Test
@quickactivate "VariationalMonteCarlo"
# Test the metropolis step of the MonteCarlo method using 1D probabilities distributions
includet(srcdir("MonteCarlo.jl"))

using Distributions, HypothesisTests

# Normal distribution
#distrib = Normal()
normalpdf(x) = pdf(Normal(), x)

points, roa = MonteCarlo.metropolis(normalpdf, 10000, 100, 1., 0.)

@testset "Normal distribution tests" begin
    @test pvalue(OneSampleTTest(points, mean(Normal()))) > 0.05
    @test pvalue(OneSampleZTest(points, mean(Normal()))) > 0.05
    @test pvalue(ApproximateOneSampleKSTest(points, Normal())) > 0.05
    @test pvalue(OneSampleADTest(points, Normal())) > 0.05
end

# Cauchy distribution
chipdf(x) = pdf(Chi(3), x)

points2, roa = MonteCarlo.metropolis(chipdf, 10000, 100, 1., 0.)

@testset "Chi(ν = 3) distribution tests" begin
    @test pvalue(OneSampleTTest(points2, mean(Chi(3)))) > 0.05
    @test pvalue(OneSampleZTest(points2, mean(Chi(3)))) > 0.05
    @test pvalue(ApproximateOneSampleKSTest(points2, Chi(3))) > 0.05
    @test pvalue(OneSampleADTest(points2, Chi(3))) > 0.05
end

# Gamma distribution
gammapdf(x) = pdf(Gamma(), x)

points3, roa = MonteCarlo.metropolis(gammapdf, 10000, 100, 1., 0.)

@testset "Gamma distribution tests" begin
    @test pvalue(OneSampleTTest(points3, mean(Gamma()))) > 0.05
    @test pvalue(OneSampleZTest(points3, mean(Gamma()))) > 0.05
    @test pvalue(ApproximateOneSampleKSTest(points3, Gamma())) > 0.05
    @test pvalue(OneSampleADTest(points3, Gamma())) > 0.05
end

# Uniform distribution
uniformpdf(x) = pdf(Uniform(-2,2), x)

points4, roa = MonteCarlo.metropolis(uniformpdf, 10000, 100, 1., 0.)

@testset "Uniform distribution tests" begin
    @test pvalue(OneSampleTTest(points4, mean(Uniform(-2,2)))) > 0.05
    @test pvalue(OneSampleZTest(points4, mean(Uniform(-2,2)))) > 0.05
    @test pvalue(ApproximateOneSampleKSTest(points4, Uniform(-2,2))) > 0.05
    @test pvalue(OneSampleADTest(points4, Uniform(-2,2))) > 0.05
end