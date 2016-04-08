# Simple script to testthe Kwatee LFBA code -

# includes -
include("DataFile.jl");
include("FluxDriver.jl");
include("Utility.jl")

function Solve(data_dictionary)

  # is the data_dictionary empty?
  if (isempty(data_dictionary) == true)

    # Load the data_dictionary -
    data_dictionary = DataFile(0,0,0);
  end

  # Split the stm into internal and external blocks -
  (intracellular_index_array, extracellular_index_array, intracellular_block,extracellular_block) = splitStoichiometricMatrixIntoInternalExternalBlocks(data_dictionary);

  # Build state mapping index -
  state_mapping_index_array = calculateStateMappingIndexArray([intracellular_index_array ; extracellular_index_array]);

  # Get model dictionaries -
  species_model_dictionary = data_dictionary["SPECIES_MODEL_DICTIONARY"];
  dsa = calculateDilutionSelectionArray(species_model_dictionary);
  tau_array = calculateTimeDelayParameterArray(species_model_dictionary);

  # Get system dimension -
  number_of_species = data_dictionary["NUMBER_OF_SPECIES"];
  number_of_fluxes = data_dictionary["NUMBER_OF_FLUXES"];

  # Setup the growth rate -
  specific_growth_rate = 0.0;
  TSTEP = 0.1;
  state_array = calculateInitialConditionArray(species_model_dictionary);
  intracellular_state_array = state_array[intracellular_index_array];
  extracellular_state_array = state_array[extracellular_index_array];
  cellmass = 0.01;

  # Growth coefficients -
  gca = calculateGrowthCoefficientArray(species_model_dictionary);

  # Call the flux driver with a SS = true
  (objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(1, state_array, specific_growth_rate, cellmass, dsa, TSTEP, data_dictionary, steady_state_flag=true);


  # return -
  return (objective_value, flux_array, dual_array, uptake_array, exit_flag);
end


function Solve(TSTART,TSTOP,TSTEP,data_dictionary)

# is the data_dictionary empty?
if (isempty(data_dictionary) == true)

  # Load the data_dictionary -
  data_dictionary = DataFile(TSTART,TSTOP,TSTEP);
end

# Split the stm into internal and external blocks -
(intracellular_index_array, extracellular_index_array, intracellular_block,extracellular_block) = splitStoichiometricMatrixIntoInternalExternalBlocks(data_dictionary);

# Build state mapping index -
state_mapping_index_array = calculateStateMappingIndexArray([intracellular_index_array ; extracellular_index_array]);

# Get model dictionaries -
species_model_dictionary = data_dictionary["SPECIES_MODEL_DICTIONARY"];
dsa = calculateDilutionSelectionArray(species_model_dictionary);
tau_array = calculateTimeDelayParameterArray(species_model_dictionary);

# Get system dimension -
number_of_species = data_dictionary["NUMBER_OF_SPECIES"];
number_of_fluxes = data_dictionary["NUMBER_OF_FLUXES"];

# Setup the growth rate -
specific_growth_rate = 0.0;
state_array = calculateInitialConditionArray(species_model_dictionary);
intracellular_state_array = state_array[intracellular_index_array];
extracellular_state_array = state_array[extracellular_index_array];
cellmass = 0.0;

# Growth coefficients -
gca = calculateGrowthCoefficientArray(species_model_dictionary);

# setup the time scale -
time_step_size = TSTEP;
time_start = TSTART;
time_final = TSTOP;
simulation_time_array = collect(time_start:time_step_size:time_final);

# state cache -
state_cache = zeros(length(simulation_time_array),number_of_species);
flux_cache = zeros(length(simulation_time_array),number_of_fluxes);
cellmass_cache = zeros(length(simulation_time_array));
time_step_index = 1;
for time_step in simulation_time_array

  # Setup modified euler scheme -


  # -- first step -------------------------------------------------------------------------------------------------------------------- #
  # Call the flux driver and calc the half step -
  (objective_value, flux_array_1, uptake_array, exit_flag) = FluxDriver(time_step_index, state_array, specific_growth_rate, cellmass, dsa, time_step_size, data_dictionary, steady_state_flag=false);
  K1_intracellular = time_step_size*(tau_array[intracellular_index_array]).*(intracellular_block*flux_array_1 - specific_growth_rate*(dsa[intracellular_index_array].*state_array[intracellular_index_array]) - gca[intracellular_index_array]*specific_growth_rate);
  K1_extracellular = time_step_size*(tau_array[extracellular_index_array]).*(extracellular_block*flux_array_1*cellmass);
  K1 = [K1_intracellular ; K1_extracellular];
  K1 = K1[state_mapping_index_array];
  # --------------n-------------------------------------------------------------------------------------------------------------------- #
  #
  # -- second step ------------------------------------------------------------------------------------------------------------------- #
  # Call the flux driver and calc the half step -
  (objective_value, flux_array_2, uptake_array, exit_flag) = FluxDriver(time_step_index, state_array+K1, specific_growth_rate, cellmass, dsa, time_step_size, data_dictionary, steady_state_flag=false);
  K2_intracellular = time_step_size*(tau_array[intracellular_index_array]).*(intracellular_block*flux_array_2 - specific_growth_rate*(dsa[intracellular_index_array].*state_array[intracellular_index_array]) - gca[intracellular_index_array]*specific_growth_rate);
  K2_extracellular = time_step_size*(tau_array[extracellular_index_array]).*(extracellular_block*flux_array_2*cellmass);
  # ---------------------------------------------------------------------------------------------------------------------------------- #
  #
  # -- update the intracellular and extracellular state ------------------------------------------------------------------------------ #
  intracellular_state_array = intracellular_state_array + (1/2)*(K1_intracellular+K2_intracellular);
  extracellular_state_array = extracellular_state_array + (1/2)*(K1_extracellular+K2_extracellular);
  flux_array = (1/2)*(flux_array_1+flux_array_2);
  tmp_state_array = [intracellular_state_array ; extracellular_state_array];
  state_array = tmp_state_array[state_mapping_index_array]
  # ---------------------------------------------------------------------------------------------------------------------------------- #
  #
  # -- update cellmass state --------------------------------------------------------------------------------------------------------- #
  cellmass = (1 + time_step_size*specific_growth_rate)*cellmass;
  # ---------------------------------------------------------------------------------------------------------------------------------- #

  # cache the state vector -
  for state_index in collect(1:number_of_species)

    # Check the state -
    if (state_array[state_index]<1e-10)
      state_array[state_index] = 0.0;
    end

    state_cache[time_step_index,state_index] = state_array[state_index];
  end

  # cache the flux -
  for flux_index in collect(1:number_of_fluxes)
    flux_cache[time_step_index,flux_index] = flux_array[flux_index];
  end

  # cache the cellmass -
  cellmass_cache[time_step_index] = cellmass;

  # Almost the last thing ...
  specific_growth_rate = calculateUpdatedSpecificGrowthRate(cellmass, state_array, flux_array, data_dictionary);

  # update my time step index, and go around again ...
  time_step_index+=1;

  @show time_step
end

return (simulation_time_array, state_cache, cellmass_cache, flux_cache)
end

function calculateUpdatedSpecificGrowthRate(cellmass, state_array, flux_vector, data_dictionary)

  species_model_dictionary = data_dictionary["SPECIES_MODEL_DICTIONARY"];
  stoichiometric_matrix = data_dictionary["STOICHIOMETRIC_MATRIX"];
  mugmax = data_dictionary["MAX_SPECIFIC_GROWTH_RATE"];

  # What is the system dimension -
  (number_of_species, number_of_fluxes) = size(stoichiometric_matrix);
  species_term_array = ones(Float64,number_of_species);

  # calculate the uptake_array -
  uptake_array = stoichiometric_matrix*flux_vector;

  # Iterate through my species models, and pull out the growth precursors -
  for (key,species_model::SpeciesModel) in species_model_dictionary

    # Is this species a growth precursor?
    is_biomass_precursor  = species_model.is_biomass_precursor
    if (is_biomass_precursor == true)

      # Get my index -
      species_index = species_model.species_index;

      # my growth coefficient -
      growth_precursor_coefficient = species_model.biomass_precursor_coefficient;

      # Calculate dilution term -
      species_term_array[species_index] = state_array[species_index]^(growth_precursor_coefficient);
    end
  end

  # before we do the calculation - if the prod(species_term_array) = 0, then no growth
  updated_specific_growth_rate = mugmax*prod(species_term_array);

  # return -
  return updated_specific_growth_rate;

end
