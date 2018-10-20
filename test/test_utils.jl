using Test, QSimulator

@testset "photons and indices" begin
    dims = [2,2]
    @test tensor_product_states(dims) == [[0,0], [0,1], [1,0], [1,1]]
    dims = [2,3]
    @test tensor_product_states(dims) == [[0,0], [0,1], [0,2], [1,0], [1,1], [1,2]]
    dims = [2, 3, 5, 4]
    for (index, photons) in enumerate(tensor_product_states(dims))
      @test photons_to_index(photons, dims) == index
    end
    @test photons_to_index(dims, dims) == nothing
    @test index_to_photons(prod(dims)+1, dims) == nothing
end

@testset "photon and index to state" begin
    dims = [2, 3]
    photons = [1,0]
    index = photons_to_index(photons, dims)
    @test index == 4
    state = Complex128[0.0, 1.0] ⊗ Complex128[1.0, 0.0, 0.0]
    @test index_to_state(index, prod(dims)) == state
    @test photons_to_state(photons, dims) == state
end
