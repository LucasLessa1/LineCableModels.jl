using Printf
using LineCableModels
using LineCableModels.ImportExport
using LineCableModels.Commons: f₀
# using Revise
# include("C:/Users/Amauri/OneDrive/Documentos/UnB/Mestrado/LineCableModels.jl/src/commons/consts.jl")

fullfile(filename) = joinpath(@__DIR__, filename); #hide
set_logger!(0); #hide

function mod_prin(len::Int, num::T) where {T<:Number}
    s_general = @sprintf("%g", num)
    y = ""
    k = findfirst('e', s_general)

    if !isnothing(k)
        mantissa_str = s_general[1:k-1]
        exponent = parse(Int, s_general[k+1:end])
        
        dec_point_idx = findfirst('.', mantissa_str)
        if !isnothing(dec_point_idx)
            num_after_dot = length(mantissa_str) - dec_point_idx
            new_exponent = exponent - num_after_dot + 1
            mantissa_no_dot = replace(mantissa_str, "." => "")
            
            if length(mantissa_no_dot) > 5
                mantissa_no_dot = mantissa_no_dot[1:5]
            end
            y = string(mantissa_no_dot, ".e", new_exponent)
        else
            y = s_general
        end
    else
        y = string(num)
    end
    
    if length(y) > len
        y = y[1:len]
    end
    
    return lpad(y, len)
end

# The function signature is now parameterized for type stability.
function generate_atp_file(
    cable_system::LineCableSystem,
    config::NamedTuple, 
    earth_params::EarthModel; 
    base_freq=f₀, 
    file_name::String="cable_lcc.dat"
)::Union{String,Nothing}

    file_name = isabspath(file_name) ? file_name : joinpath(@__DIR__, file_name)
    # We write to an in-memory IOBuffer to minimize string allocations. This is much faster.
    buffer = IOBuffer()

    println(buffer, "BEGIN NEW DATA CASE")
    println(buffer, "CABLE CONSTANTS")
    println(buffer, "CABLE PARAMETERS")

    write(buffer, "BRANCH  ")
    for (i, cable) in enumerate(cable_system.cables)
        write(buffer, "IN_", lpad(i, 3, '0'), "OUT", lpad(i, 3, '0'))
    end
    println(buffer)

    num_cables = length(cable_system.cables)
    num_components = [length(cable.design_data.components) for cable in cable_system.cables]
    vert_comps = [cable.vert for cable in cable_system.cables]# Air=1, Surface=0, Ground=-1
    position = sign(count(x -> x < 0, vert_comps) > length(cable_system.cables)/2 ? -1 : 1)
    model = 1 # Bergeron=1

    println(buffer, "    2   ", mod_prin(2, position), "    ", num_cables, "         ", model, "    1            ", (sum(num_components) - num_cables), "    0    ", config.addCG)
    
    join(buffer, ("    " * string(n) for n in num_components))
    println(buffer)

    for cable in cable_system.cables
        for comp in cable.design_data.components
            write(buffer, mod_prin(10, comp.conductor_group.radius_in))
            write(buffer, mod_prin(10, comp.conductor_group.radius_ext))
        end
        write(buffer, mod_prin(10, cable.design_data.components[end].insulator_group.radius_ext))
        println(buffer)

        let n_written = 0
            for comp in cable.design_data.components
                for val in (comp.conductor_props.rho, comp.conductor_props.mu_r, comp.insulator_group.layers[1].material_props.mu_r, comp.insulator_group.layers[1].material_props.eps_r)
                    write(buffer, mod_prin(10, val))
                    n_written += 10
                    if n_written >= 80
                        println(buffer)
                        n_written = 0
                    end
                end
            end
            if n_written > 0; println(buffer); end
        end

        if config.addCG > 0
            println(buffer, "C Additional conductance and capacitance card")
            let n_written = 0
                for comp in cable.design_data.components
                    ins_group = comp.insulator_group
                    if config.addCG == 1
                        vals = ("", ins_group.shunt_capacitance)
                    elseif config.addCG == 2
                        vals = (ins_group.shunt_conductance, "")
                    else # config.addCG == 3
                        vals = (ins_group.shunt_conductance, ins_group.shunt_capacitance)
                    end

                    for val in vals
                        str_to_write = val isa Number ? mod_prin(10, val) : lpad("", 10)
                        write(buffer, str_to_write)
                        n_written += 10
                        if n_written >= 80
                            println(buffer)
                            n_written = 0
                        end
                    end
                end
                if n_written > 0; println(buffer); end
            end
        end
    end

    let n_written = 0
        for cable in cable_system.cables
            for val in (cable.vert, cable.horz)
                write(buffer, mod_prin(10, val))
                n_written += 10
                if n_written >= 80
                    println(buffer)
                    n_written = 0
                end
            end
        end
        if n_written > 0; println(buffer); end
    end

    println(buffer, "     ", mod_prin(10, earth_params.layers[end].base_rho_g), "     ", mod_prin(10, base_freq), "     ", mod_prin(10, cable_system.line_length), "     ", config.modal_z_type)

    println(buffer, "BLANK CARD ENDING FREQUENCY CARDS")
    println(buffer, "\$PUNCH")
    println(buffer, "BLANK CARD ENDING CABLE CONSTANTS")
    println(buffer, "BEGIN NEW DATA CASE")
    println(buffer, "BLANK CARD")
    
    # Write the entire buffer to the file in a single, efficient operation.
    output_string = String(take!(buffer))
    open(file_name, "w") do file
        write(file, output_string)
    end

    @info "ATP file generated: $file_name"
    
    return file_name
end

using Base.Filesystem

"""
    run_dat_file(dat_file::String; ...)

Executes an ATP simulation using the definitive and robust method for legacy
solvers: it creates a temporary 'startup' file, changes the working directory
to the data file's location, runs the solver, captures the output, and then
cleans up all temporary files while safely returning to the original directory.
"""
function run_dat_file(
    dat_file::String; 
    atp_executable::String = "C:\\ATPDraw\\ATP_solver2023\\tpgigm.exe"
)
    # The required content of the 'startup' file.
    startup_content = """
    1  RHIGH   EPSZNO  EPWARN  EPSTOP  EPSUBA  EPDGEL  EPOMEG    SZPLT   SZBED   TENFLZ
      1.D+10   1.D-8   1.D-3     0.1   100.  1.D-16  1.D-15     10.0    72.0     10.
    2 SIGMAX  TENERG  DEGMIN  DEGMAX  ZNOLIM(1), (2)  STATFR  ZNVREF  XMAXMX    AINCR
         4.0  1.D+20     0.0   360.     1.0     1.5    60.0   1.D-6     2.0     .05
    3 FREQFR  HLETT1  Unused     VHS     VS     VH  TAXISL  VAXISL    FILL1    FILL2
               0.25            8.0     1.0     10.    20.0     8.0     6.      7.0
    4 TOLRCE   FHTAX   FXSUP   FYSUP   FXTIT   FYTIT  VPLOTS  VPLOTL  FACTVI  FTCARR
      8.E-5     0.5     .25     .03    0.10     0.1     1.0     5.0     0.0     1.5
    5 FXNUMV  FXNUMH  FVAXTT  FXVERT  UNIXON  TIMTAC  OVRLAP  FLZERO  EPSILN  FLTINF
         1.5     5.0    -1.5     0.0     0.0     0.0     0.5  1.D-12  1.D-12  1.D+19
    6 XHEADM  YHEADM  HGTHDM  XCASTI  YCASTI  HGTCST  XLEGND  YLEGND  HGTLGN  TSTALL
         2.5    7.95     .55     0.5     7.3     .35     0.5    1.30     .25    -0.0
    7 XALPHA  YALPHA  HGTALF  D4FACT  PEKEXP  EPSLRT  EPSPIV  PLMARK  FACOSC  Unused
         1.5     6.5     .25     0.0     43.  1.E-12  1.E-16     1.0     0.3
    8 NMAUTO  INTINF  KOL132  MUNIT5  MAXZNO  IPRSPY  IPRSUP   LNPIN  MINHAR  MAXHAR
           1 9999999     132       1      50       0       0       6       0      20
    9 NFORS2  NIOMAX   MRGN  LINLIM   MPAGE  MODE28  KPGRID KPEN(1) KPEN(2) KPEN(3)
         30      10       2     100       0       1       3      12      10      11
    10 ..(4)  KOMLEV   NSMTH  MODSCR  KOLALP  MAXFLG  LIMCRD  NOBLAN  MOUSET  NOTPPL
         14      -1      50       2       5       1   30000       0       0       1
    1 NOCOMM  NOHELP  NEWPL4  JDELAY  NOTMAX  NSMPLT  KOLWID  KOLSEP  JCOLU1  KSLOWR
           0       0       0       0       0      50      11       1       0      25
    2 KSYMBL  NOBACK  KOLEXM   LTEK   NCUT1   NCUT2  INCHPX  INCHPY  NODPCX  LCHLIM
         200       1      60       1      13      11       2       2       0       0
    13 NORUN  JTURBO  MAXSYM     IHS  LIMCOL   KLEVL   KEXTR  NOHPGL  NOPOST  NOSM59
           0       5       3       3      79       0       0       0       0       0
    4 LEFTA6  LENREC  LU6VRT   LRLIM  KASEND  LUNDAT  KTRPL4  JORIEN  LIMPNL  LUNTEX
           0       0   32768      75       5       3   -6666       0     200     -11
    5 KINSEN  LISTON  LIMTAC  NOCALC  MFLUSH  L4BYTE  KOMPAR  LIST01   NOGNU  KROSEC
           1       0      25       0    1000       1       0       0       0       0
    6 LUNIT1  LUNIT2  LUNIT3  LUNIT4  LUNIT5  LUNIT6  LUNIT7  LUNIT8  LUNIT9  LUNT10
          21      22       3      -4       1       6       7       8       9      10
    17 KS(1)   KS(2)   KS(3)   KS(4)   KP(1)   KP(2)   KP(3)   KP(4)  KOLROV  NUMHLD
           0       0      12      10       7      14       0       0      18
    8 L4FULL  NOQUOT  JJEATS  NUMBUS   NOTAB  NOPISA   MSCSV  MAXL31  LIM132  MAXMVC
           0       0       0      -1       0       0       0     400       0      80
    18 Name of language font file  ] Window] Root name for SPY @K usage    ]
    blockd51.bin                  junk    inclspy .dat
    9 SSONLY  CHEFLD  TEXNAM  CHVBAR  BRANCH  TXCOPY  USERID  -TRASH  -TERRA  CHRCOM
      PHASOR  E       DUM     |       NAME    COPY    Hannov  ......  TERRA   C {}\$,
    0 DATTYP  LISTYP  PCHTYP  PL4TYP  EFIELD  FMTPL4  PSCTYP  DBGTYP  BINTYP  EXTTYP
      .dat    .lis    .pch    .pl4                      .ps     .dbg    .bin    .ext
    C    After regular STARTUP comes optional VMS-like symbol definitions that are
    C    used for input data file name in response to the opening prompt.
    scott:==c:\\atp\\    { 1st of 2 remote directories
    tsu:==c:\\tsu-huei\\    { 2nd of 2
    \$EOF    { Software end-of-file terminates last of 20 or fewer VMS-like symbols
    """

    println("Executing ATP for: ", basename(dat_file))
    
    dat_dir = dirname(dat_file)
    temp_startup_path = joinpath(dat_dir, "startup")
    base_path_name = splitext(dat_file)[1]

    captured_output = ""

    try
        # --- 1. Create the temporary 'startup' file in the data directory ---
        write(temp_startup_path, startup_content)
        
        # --- 2. Change to the data directory to run the solver ---
        # The `cd(path) do ... end` block automatically returns to the original directory.
        cd(dat_dir) do
            println("Changed directory to: $(pwd())")
            # The command uses the full path to the executable, but just the
            # filename for the .dat file since we are now in its directory.
            atp_cmd = Cmd([atp_executable, basename(dat_file)])
            
            # --- 3. Capture the output ---
            captured_output = readchomp(atp_cmd)
        end # Automatically returns to the original directory here

        println("SUCCESS: ATP execution complete and output captured.")

    catch e
        println("ERROR: Failed during ATP execution.")
        println(e)
    finally
        # --- COMPREHENSIVE CLEANUP ---
        println("Cleaning up all temporary files...")
        
        # 1. Clean up files with known names based on the input .dat file
        base_path_name = splitext(dat_file)[1]
        rm(joinpath(dat_dir, "startup"), force=true)
        rm(base_path_name * ".dbg", force=true)
        rm(base_path_name * ".pch", force=true)
        rm(temp_startup_path, force=true)
        
        # 2. Scan the directory and clean up files with known extensions
        extensions_to_delete = [".tmp", ".bin", ".49"]
        for filename in readdir(dat_dir)
            if any(ext -> endswith(filename, ext), extensions_to_delete)
                rm(joinpath(dat_dir, filename), force=true)
            end
        end
        
        println("Cleanup complete.")
    end
    return captured_output
end

config = (
    addCG=3,        # C=1, G=2, CG=3
    modal_z_type=2 # 2-sqrt(Im{Z}/Im{Y}), 3-Re(sqrt(Z/Y))
)
output_file = fullfile("my_cable.dat")
file = generate_atp_file(cable_system, config, earth_params, file_name=output_file)
lis_content = run_dat_file(file)
Z = read_data(Val(:atp), lis_content, cable_system)
# atp_executable="C:/ATP/tools/runATP.exe"
# atp_cmd = `$atp_executable $file`

# # run(atp_cmd)
# output = readchomp(atp_cmd)
# display(output)
# # println(file)
println("File saved to: ", output_file)
display(Z)
