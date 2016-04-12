# Script to sample the bounds file -

# Script to max product -
include("DataFile.jl")
include("Solve.jl")
include("Types.jl")

# Load the data file with the default bounds -
data_dictionary = DataFile(0,0,0);

# Set the max -
data_dictionary["MIN_FLAG"] = false;

# Get the flux model dictionary -
flux_model_dictionary = data_dictionary["FLUX_MODEL_DICTIONARY"]
flux_model_dictionary["PROTEIN_export_deGFP"].flux_obj_coeff = 1.0;

# set the bound -
data_dictionary["R_malS_upper_bound"] = 10;
data_dictionary["R_malS_lower_bound"] = 0;
data_dictionary["M_asn_L_c_exchange_reverse_upper_bound"] = 10;
data_dictionary["M_o2_c_exchange_reverse"] = 10
data_dictionary["translation_deGFP_switch_bound"] = false
production_time = 9;

# Set the plasmid concentration -
gene_copy_number = (10e-6/1e9)*(6.02e23)*10.0;
data_dictionary["deGFP_gene_copies"] = gene_copy_number;

# Load the bounds data -
bounds_array = readdlm("Samples.txt")
(nr,nc) = size(bounds_array)

# Setup performance array -
performance_array = zeros(nr)

# Ok, we need to enumerate through this list, and calculate the performance of the model for each value -
for index in collect(1:nr)

  # Set the bounds array -
  data_dictionary["sample_bounds_array"] = bounds_array[index,:]

  # Run the simulation -
  (OV,FA,DA,UA,EF) = Solve(data_dictionary);

  # calculate the protein concentration -
  # deGFP_concentration = (1e3)*OV*production_time;

  # capture -
  performance_array[index] = (1e3)*OV
  @show index
end

# Write the concentration file -
writedlm("./performance.txt",performance_array)
