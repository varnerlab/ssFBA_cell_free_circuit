# Script to max product -
include("DataFile.jl")
include("Solve.jl")
include("Types.jl")

# PyPlot -
using PyPlot

# Setup the objective -
data_dictionary = DataFile(0,0,0);
number_of_fluxes = data_dictionary["NUMBER_OF_FLUXES"];
number_of_samples = 100;

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
data_dictionary["sample_bounds_array"] = 10.0*ones(20)

# Setup steady-state array -
plasmind_concentration_array = collect(linspace(0,10,100));
production_time = 9;
deGFP_array = zeros(length(plasmind_concentration_array))
flux_array = zeros(number_of_samples,number_of_fluxes);

for (index,plasmid_concentration) in enumerate(plasmind_concentration_array)

  # calculate the gene copy number -
  gene_copy_number = (10e-6/1e9)*(6.02e23)*plasmid_concentration;

  # set this value in data dictionary -
  data_dictionary["deGFP_gene_copies"] = gene_copy_number;

  # Run the simulation -
  (OV,FA,DA,UA,EF) = Solve(data_dictionary);

  # calculate the protein concentration -
  deGFP_concentration = (1e3)*OV*production_time;

  # capture -
  deGFP_array[index] = deGFP_concentration;

  for flux_index in collect(1:number_of_fluxes)
    flux_array[index,flux_index] = FA[flux_index]
  end

end

# plot -
#clf()
data = readdlm("Plasmid-IC-SS-P70.dat");
plot(data[:,1],data[:,2],"ko")
plot(plasmind_concentration_array,deGFP_array,linewidth="2.0")
xlabel("Plasmid concentration [nM]",fontsize="14")
ylabel(L"deGFP concentration [$\mu$M]",fontsize="14")
