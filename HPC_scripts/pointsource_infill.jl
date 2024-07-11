first_time = false
ENV["JULIA_CONDAPKG_BACKEND"] = "MicroMamba";
if first_time
    Pkg.add(url="https://github.com/andrew-saydjari/CloudClean.jl")
    Pkg.add("FITSIO")
    Pkg.add(["CondaPkg","PythonCall","PythonPlot"])
    using CondaPkg
    CondaPkg.add("colorcet")
end
 
import Pkg
Pkg.add("BenchmarkTools")
using BenchmarkTools
using CloudClean, FITSIO
using PyPlot
using Plots
using Pickle
using Glob

home_dir = "/xdisk/hamden/hina0830/venv39"

println("ARGS: ", ARGS) # Get the command line argument; first arguement is array job number and second argument is the path to the input file

# Get the first argument
index = parse(Int, ARGS[1]) # This correponds to the array number, which is index number + 1

# Get the second argument
inputfile_path = ARGS[2]


# Function that reads the input file
function read_and_split_file(filename)
    split_lines = []
    open(filename, "r") do file
        for line in eachline(file)
            push!(split_lines, split(line)) # Get 4 integers in the line and add them to the list
        end
    end
    return split_lines
end

# Read the line in the input file and get values (batch#, index#, RA, and DEC)
result = read_and_split_file(inputfile_path)[index]

Batch, Index, RA, DEC = result[1], result[2], result[3], result[4] # Define batch#, index#, RA, and DEC
foldername = "/Batch_" * Batch * "/Index" * Index * "-RA" * RA * "-DEC" * DEC 
path_init = home_dir * "/raw_files_newformat" * foldername
println("path_init: ", path_init)

coord_fn, out_fn, bimage_fn, orig_fn = sort(glob("*_coord.pkl", path_init), rev=false), sort(glob("*_masked.fits", path_init), rev=false), sort(glob("*_mask.fits", path_init), rev=false), sort(glob("*-int_Pinfilled.fits", path_init), rev=false)
println(length(coord_fn))

cd(path_init)
println("path changed", pwd())

# Function to print thread ID
function thread_print(string)
    coord_fn, out_fn, bimage_fn, orig_fn = string
    coords = Pickle.load(coord_fn; proto = 5)
    x_locs = [sub[1] for sub in coords]
    y_locs = [sub[2] for sub in coords]
    println("Tuple set. ", length(x_locs), " stars.")

    # Open files 
    filename_original = out_fn 
    f = FITS(filename_original)
    init_image = read(f[1])

    f = FITS(out_fn) 
    out_image = read(f[1])
    close(f)

    f = FITS(bimage_fn) 
    bimage = read(f[1])
    close(f)

    # Convert image to Boolean matrix
    bimage = [iszero(element) for element in bimage]
    bimage_bool = !=(1).(bimage)

    # Parameters
    Np = 95
    halfNp = (Np-1)÷2
    dv = halfNp
    shiftx = 0
    rlim = 20^2

    ndraw0 = 1
    widx = 600 

    # Run the infilling algorithm 
    star_stats = proc_discrete(x_locs.+1 , y_locs.+1 , out_image, bimage_bool, Np=Np, rlim=Inf, tilex=8, ftype=64, widx=widx, seed=2022, ndraw=ndraw0);
   
    # Open the original FITS file
    f_original = FITS(filename_original, "r")

    # Read the data and header
    data = read(f_original[1])
    mean = star_stats[1]
    draw = star_stats[2][:, :, 1]

    # Update the header to reflect the changes in the data
    header = read_header(f_original[1])
    header["NAXIS1"] = size(data, 2)
    header["NAXIS2"] = size(data, 1)

    # Close the original FITS file
    close(f_original)

    # Create a new FITS file for writing with the modified data and header
    new_filename = replace(filename_original, ".fits" => "_Np$(Np)_widx$(widx)_mean.fits" ) 
    println(new_filename)
    f_modified = FITS(new_filename, "w")
    write(f_modified, mean, header=header)

    # Create a new FITS file for writing with the modified data and header
    new_filename = replace(filename_original, ".fits" => "_Np$(Np)_widx$(widx)_draw.fits" ) 
    println(new_filename)
    f_modified = FITS(new_filename, "w")
    write(f_modified, draw, header=header)

    # Close the new FITS file
    close(f_modified)
    
    println("check point3 (files closed)")
    
end

# Define a function to handle the processing of each file
function process_file(file_tuple)
    try
        thread_print(file_tuple)
    catch ex
        println("Error encountered!")
        println(ex)
        println("Error file:", file_tuple)
    end
end


a = zeros(length(coord_fn))


# Run the algorithm using multi-threading 
Threads.@threads for i = 1:length(coord_fn)
    println("iteration $i on thread $(Threads.threadid())")
    a[i] = Threads.threadid()
    file_tuple = coord_fn[i], out_fn[i], bimage_fn[i], orig_fn[i]
    #@time thread_print(file_tuple)
    @time process_file(file_tuple)
    GC.gc()
    
end

println(a)
