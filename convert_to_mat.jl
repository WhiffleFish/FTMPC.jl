using JLD2, MAT

# Name of file
file_name = "nonrobust" # nominal" # "unit" #"nonrobust" 
mode = "7" #"2" 

# Load the JLD2 file
jld2_data = JLD2.load(file_name * mode * ".jld2")

# Extract the data
# data_to_convert = jld2_data["nominalhist2"]
data_to_convert = jld2_data[file_name * "hist"* mode]

dfx = data_to_convert.x
dft = data_to_convert.t

# Specify the MAT file
mat_file_name = file_name * mode* ".mat"

# Create a MAT file and save the variables
matopen(mat_file_name, "w") do file
    write(file, "dfx_"*file_name*"_"*mode, dfx)
    write(file, "dft_"*file_name*"_"*mode, dft)
end

println("Variables saved to $mat_file_name")