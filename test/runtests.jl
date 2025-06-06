using LineCableModels
using Test

@testset "LineCableModels.jl tests" begin
	@info "Running BaseParams tests..."
	@testset "BaseParams module" begin
		include("baseparams.jl")
	end
	@info "BaseParams tests completed."

	@testset "DataModel module 1/1" begin
		include("datamodel.jl")
	end
	@info "DataModel tests completed."


	@info "All tests completed."
end
