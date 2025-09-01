# This code adds 'using Revise' to your Julia startup file
startup_file = joinpath(homedir(), ".julia", "config", "startup.jl")
mkpath(dirname(startup_file))
open(startup_file, "a") do f
    # The 'try...catch' block is a robust way to handle this in a startup file
    write(f, "\ntry\n    using Revise\ncatch e\n    @warn \"Error initializing Revise\" e\nend\n")
end

println("Revise.jl has been added to your startup file: ", startup_file)