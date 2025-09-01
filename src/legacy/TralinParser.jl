module TralinParser
 
using LinearAlgebra

export parse_tralin_file, extract_tralin_variable, take_complex_list, clean_variable_list

"""
    parse_tralin_file(filename, ord)

Parse a TRALIN file and extract impedance, admittance, and potential coefficient matrices
for multiple frequency samples.
"""
function parse_tralin_file(filename; ord)
    fileLines = readlines(filename)

    # Get all occurrences of "GROUND WIRES ELIMINATED"
    limited_str = "GROUND WIRES ELIMINATED"
    all_idx = findall(row -> occursin(limited_str, row), fileLines)

    # Initialize arrays to store matrices for all frequency samples
    Z_matrices = Vector{Matrix{ComplexF64}}(undef, length(all_idx))
    Y_matrices = Vector{Matrix{ComplexF64}}(undef, length(all_idx))
    P_matrices = Vector{Matrix{ComplexF64}}(undef, length(all_idx))

    # Loop through each frequency block
    for (k, start_idx) in enumerate(all_idx)

        # Slice the file from the current "GROUND WIRES ELIMINATED" position to end
        block_lines = fileLines[start_idx:end]

        # Extract matrices for this frequency sample, ensuring output is ComplexF64
        Z_matrices[k] = Complex{Float64}.(extract_tralin_variable(block_lines, ord, "SERIES IMPEDANCES - (ohms/kilometer)", "SHUNT ADMITTANCES (microsiemens/kilometer)"))
        Y_matrices[k] = Complex{Float64}.(extract_tralin_variable(block_lines, ord, "SHUNT ADMITTANCES (microsiemens/kilometer)", "SERIES ADMITTANCES (siemens.kilometer)"))
        P_matrices[k] = Complex{Float64}.(extract_tralin_variable(block_lines, ord, "POTENTIAL COEFFICIENTS (meghoms.kilometer)", "SERIES IMPEDANCES - (ohms/kilometer)"))
    end

    # Convert lists of matrices into 3D arrays for each matrix type
    Z_stack = reshape(hcat(Z_matrices...), ord, ord, length(Z_matrices))
    Y_stack = reshape(hcat(Y_matrices...), ord, ord, length(Y_matrices))
    P_stack = reshape(hcat(P_matrices...), ord, ord, length(P_matrices))

    Z_stack = Z_stack ./ 1000
    Y_stack = Y_stack .* 1e-6 ./ 1000
    P_stack = P_stack .* 1e6 .* 1000

    return Z_stack, Y_stack, P_stack
end

"""
    extract_tralin_variable(fileLines, order, str_init, str_final)

Extracts matrix data between specified headers in `fileLines`, handling complex formatting.
"""
function extract_tralin_variable(fileLines, order, str_init, str_final)
    # Locate header and footer lines
    variable_init = findfirst(line -> occursin(str_init, line), fileLines)
    variable_final = findfirst(line -> occursin(str_final, line), fileLines)

    if isnothing(variable_init) || isnothing(variable_final)
        println("Could not locate start or end of the block.")
        return zeros(ComplexF64, order, order)
    end

    # Parse the relevant lines into a list of complex numbers
    variable_list_number = []
    for line in fileLines[variable_init+15:variable_final-1]
        numbers = take_complex_list(line)
        if !isempty(numbers)
            push!(variable_list_number, numbers)
        end
    end

    # Process, clean, and arrange data into matrix form
    variable_list_number = clean_variable_list(variable_list_number, order)

    # Initialize matrix and fill, with padding if necessary
    matrix = zeros(ComplexF64, order, order)
    for (i, row) in enumerate(variable_list_number)
        matrix[i, 1:length(row)] = row
    end

    # Make symmetric by filling lower triangle
    matrix += tril(matrix, -1)'

    return matrix
end


"""
    take_complex_list(s)

Parses a string to identify real and complex numbers, with conditional scaling for scientific notation.
"""
function take_complex_list(s)
    numbers = []

    # Match the first real number (decimal, integer, or scientific notation)
    first_real_pattern = r"([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?|\d+)"
    first_real_match = match(first_real_pattern, s)
    if !isnothing(first_real_match)
        real_part_str = strip(first_real_match.match)
        real_value = occursin(r"[Ee]", real_part_str) ? parse(Float64, real_part_str) : parse(Float64, real_part_str) * 1
        push!(numbers, real_value)
    end

    # Match complex numbers (handles scientific notation or regular float, allowing extra whitespace before 'j')
    complex_pattern = r"([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?|\d+)\s*\+\s*j\s*([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?|\d+)"
    for m in eachmatch(complex_pattern, s)
        real_part_str, imag_part_str = m.captures
        real_value = occursin(r"[Ee]", real_part_str) ? parse(Float64, real_part_str) : parse(Float64, real_part_str) * 1
        imag_value = occursin(r"[Ee]", imag_part_str) ? parse(Float64, imag_part_str) : parse(Float64, imag_part_str) * 1
        push!(numbers, Complex(real_value, imag_value))
    end

    return numbers
end


"""
    clean_variable_list(variable_list_number, order)

Cleans and arranges extracted list into a proper matrix format.
"""
function clean_variable_list(data, order)
    # Remove entries that lack values, filter short lists
    filter!(lst -> length(lst) > 1, data)
    
    # Trim row label elements and only keep the actual data
    data = [lst[2:end] for lst in data]

    # Apply padding to each row as needed to align with specified order
    data_padded = [vcat(lst, fill(0.0 + 0.0im, order - length(lst))) for lst in data]

    # Ensure `data_padded` has `order` rows; add extra rows of zeros if required
    if length(data_padded) < order
        for _ in 1:(order - length(data_padded))
            push!(data_padded, fill(0.0 + 0.0im, order))
        end
    end

    return data_padded
end

end # of module
