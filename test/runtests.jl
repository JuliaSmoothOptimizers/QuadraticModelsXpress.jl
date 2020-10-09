using QPSReader, QuadraticModels, QuadraticModelsXpress
using Test

@testset "QuadraticModelsXpress.jl" begin
    qps1 = readqps("QAFIRO.SIF") #lower bounds
    stats1 = xpress(QuadraticModel(qps1))
    @test isapprox(stats1.objective, -1.59078179, atol=1e-2)

    qps2 = readqps("HS21.SIF") # low/upp bounds
    stats2 = xpress(QuadraticModel(qps2))
    @test isapprox(stats2.objective, -9.99599999e1, atol=1e-2)

    qps3 = readqps("HS52.SIF") # free bounds
    stats3 = xpress(QuadraticModel(qps3))
    @test isapprox(stats3.objective, 5.32664756, atol=1e-2)

    qps4 = readqps("AFIRO.SIF") # test lp
    stats4 = xpress(QuadraticModel(qps4))
    @test isapprox(stats4.objective, -4.6475314286e2, atol=1e-2)
end
