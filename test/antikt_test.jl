using AntiKt
using Test

@testset "AntiKt test" begin
    #define variable event, a vector of Pseudojet objets
    include("event.jl")

    #Antikt clustering distance parameter
    R=0.4

    #Minimum pt in GeV for the clusterinh
    ptmin= 5

    #Run the jet clustering on the event
    jets = antikt(event, R, ptmin)

    #Six jets expected in the event
    @test length(jets) == 6
end


