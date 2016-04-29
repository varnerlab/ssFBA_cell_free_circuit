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

# Setup bounds dictionary -
bounds_dictionary = Dict()
bounds_dictionary["M_ala_L_c_exchange_reverse"] = 10.0
bounds_dictionary["M_arg_L_c_exchange_reverse"] = 10.0
bounds_dictionary["M_asn_L_c_exchange_reverse"] = 10.0
bounds_dictionary["M_asp_L_c_exchange_reverse"] = 10.0
bounds_dictionary["M_cys_L_c_exchange_reverse"] = 10.0
bounds_dictionary["M_glu_L_c_exchange_reverse"] = 10.0
bounds_dictionary["M_gln_L_c_exchange_reverse"] = 10.0
bounds_dictionary["M_gly_L_c_exchange_reverse"] = 10.0
bounds_dictionary["M_ile_L_c_exchange_reverse"] = 10.0
bounds_dictionary["M_leu_L_c_exchange_reverse"] = 10.0
bounds_dictionary["M_his_L_c_exchange_reverse"] = 10.0
bounds_dictionary["M_lys_L_c_exchange_reverse"] = 10.0
bounds_dictionary["M_met_L_c_exchange_reverse"] = 10.0
bounds_dictionary["M_phe_L_c_exchange_reverse"] = 10.0
bounds_dictionary["M_pro_L_c_exchange_reverse"] = 10.0
bounds_dictionary["M_ser_L_c_exchange_reverse"] = 10.0
bounds_dictionary["M_thr_L_c_exchange_reverse"] = 10.0
bounds_dictionary["M_trp_L_c_exchange_reverse"] = 10.0
bounds_dictionary["M_tyr_L_c_exchange_reverse"] = 10.0
bounds_dictionary["M_val_L_c_exchange_reverse"] = 10.0
data_dictionary["aa_uptake_bounds_dictionary"] = bounds_dictionary;

# Setup the time scale -
TSTART = 0.0;
TSTOP = 16.0;
Ts = 0.1;
TSIM = collect(TSTART:Ts:TSTOP);

# Setup initial conditions -
maltose_array = zeros(length(TSIM)+1)
maltose_array[1] = 20;

deGFP_array = zeros(length(TSIM)+1);
deGFP_array[1] = 0.0;

# Setup flux_array -
number_of_fluxes = data_dictionary["NUMBER_OF_FLUXES"];
flux_array = zeros(length(TSIM),number_of_fluxes)

# Main loop -
for (index,time_value) in enumerate(TSIM)

  # Setup bound -
  # if (maltose_array[index]>0)
  #   R_malS_upper_bound = 10;
  #   R_malS_lower_bound = maltose_array[1]/8;
  #   M_asn_L_c_exchange_reverse_upper_bound = 10;
  #   M_o2_c_exchange_reverse_upper_bound = 10;
  #   translation_deGFP_switch_bound = false
  # else
  #   R_malS_upper_bound = 0;
  #   R_malS_lower_bound = 0;
  #   M_asn_L_c_exchange_reverse_upper_bound = 0.0;
  #   M_o2_c_exchange_reverse_upper_bound = 0.0;
  #   translation_deGFP_switch_bound = true
  # end

  # set the bound -
  data_dictionary["R_malS_upper_bound"] = R_malS_upper_bound;
  data_dictionary["R_malS_lower_bound"] = R_malS_lower_bound;
  data_dictionary["M_asn_L_c_exchange_reverse_upper_bound"] = M_asn_L_c_exchange_reverse_upper_bound;
  data_dictionary["M_o2_c_exchange_reverse"] = M_o2_c_exchange_reverse_upper_bound
  data_dictionary["translation_deGFP_switch_bound"] = translation_deGFP_switch_bound

  # Call solve -
  (OV,FA,DA,UA,EF) = Solve(data_dictionary);

  # Get the uptake flux -
  maltose_uptake_flux = FA[1];

  # Compute next value of maltose -
  maltose_array[index+1] = maltose_array[index]-Ts*maltose_uptake_flux;

  # Compute the next value of deGFP -
  if (index>2)
    deGFP_array[index+1] = deGFP_array[index]+Ts*OV;
  end

  # package the fluxes -
  for flux_index in collect(1:number_of_fluxes)
    flux_array[index,flux_index] = FA[flux_index]
  end

end

# plot -
clf()
data = readdlm("S70-Simple.dat");
plot(data[:,1],data[:,2],"ko")
plot(TSIM,deGFP_array[1:end-1]*1e3,linewidth="2.0")
xlabel("Time [hr]",fontsize="14")
ylabel(L"deGFP concentration [$\mu$M]",fontsize="14")
