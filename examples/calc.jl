
# # area_filaments_mm2 = [0.5, 0.75,1,1.5,2.5,4,6,10,16,25,35,50,70,95,120,150,185,240,300,400] # mm²
# # area_filaments = area_filaments_mm2 * 1e-6 # mm²
# # diam_filaments = sqrt.(area_filaments ./ π)*2 # m
# # diam_fix = sqrt(2500/pi)*1e-3*2 # m


# # for (idx, diam) in enumerate(diam_filaments)
# #     for n_layers in 2:1:10000
# #         radius = (n_layers-1)*diam + diam/2
# #         error = radius*2 - diam_fix
        
# #         if abs(error) < 1e-4
# #             println("The filament with an diameter equal to $(diam) m with $n_layers layers is the closest to $diam_fix, with $(radius*2) m of diameter")
# #             break
# #         end
# #     end
# # end


# Input data
# area_filaments_mm2 = [0.5, 0.75, 1, 1.5, 2.5, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 25, 28, 29, 30, 32, 34, 35, 36, 38, 39, 40, 42, 44, 46,48, 50, 52,54,56] # mm²
area_filaments_mm2 = [i for i in range(1,100)] 
# Convert units and calculate diameters
area_filaments_m2 = area_filaments_mm2 .* 1e-6 # m²
diam_filaments_m = @. sqrt(area_filaments_m2 / π) * 2 # m

# Target diameter to match
diam_target_m = 57.8e-3#sqrt(2500 / π) * 1e-3 * 2 # m

# Flag to exit the outer loop once a match is found
match_found = false

# Iterate to find the first match
for (idx, diam) in enumerate(diam_filaments_m)
    for n_layers in 2:10000
        radius = (n_layers - 1) * diam + (diam / 2)
        number_of_strands = 1 + 3 * n_layers * (n_layers - 1) # Total strands for n layers

        total_diam = 2 * radius
        
        if abs(total_diam - diam_target_m) < 1e-3
            # A match was found, print the details
            println("Match Found:")
            println("  - Idx: $(idx)")
            println("  - Filament Area: $(area_filaments_mm2[idx]) mm²")
            println("  - Strand Diameter: $(sqrt(area_filaments_mm2[idx]/(pi))/1000) m")
            println("  - Number of Layers: $(n_layers)")
            println("  - Number of Strands: $(number_of_strands)")
            println("  - Calculated Diameter: $(round(total_diam * 1000, digits=5)) mm")
            println("  - Target Diameter:     $(round(diam_target_m * 1000, digits=5)) mm")
            
        end
    end
end

real_diam = 57.44769
rin = real_diam/2+0.6+1.5+21.3+1.4+0.7+3+2.5
rin =0.04855
println(2*pi*rin/5e-3)



# # --- EXAMPLE USAGE ---
# filepath = "C:\\ATP\\work\\18kV_1000mm2_equivalent_1.lis"
# Z_atp = read_atp_data(filepath, cable_system)
# # display(Ze)
# display(Z_atp)


1+1