# Script to max product -
include("DataFile.jl")
include("Solve.jl")
include("Types.jl")

# Setup the objective -
data_dictionary = DataFile(0,0,0);

# Set the max -
data_dictionary["MIN_FLAG"] = false;

# Get the flux model dictionary -
flux_model_dictionary = data_dictionary["FLUX_MODEL_DICTIONARY"]
flux_model_dictionary["PROTEIN_export_deGFP"].flux_obj_coeff = 1.0;

# Call solve -
(OV,FA,DA,UA,EF) = Solve(data_dictionary);
