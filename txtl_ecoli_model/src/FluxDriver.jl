include("DataFile.jl")
include("Control.jl")
using GLPK

# ----------------------------------------------------------------------------------- #
# Copyright (c) 2016 Varnerlab
# School of Chemical Engineering Purdue University
# W. Lafayette IN 46907 USA

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ----------------------------------------------------------------------------------- #
function FluxDriver(time_index, species_abundance_array, specific_growth_rate, cellmass, dsa, step_size, data_dictionary; steady_state_flag=true)
# ----------------------------------------------------------------------------------- #
# FluxDriver.jl was generated using the Kwatee code generation system.
# FluxDriver: Solves the flux balance analysis problem from TSTART to TSTOP given the model encoded in data_dictionary.
# Username: jeffreyvarner
# Type: NFBA-JULIA
# Version: 1.0
# Generation timestamp: 04-07-2016 18:07:24
#
# Input arguments:
# data_dictionary  - Data dictionary instance (holds species and flux models etc)
#
# Return arguments:
# objective_value	 - value of the objective returned by GLPK
# calculated_flux_array	 - optimal metabolic flux array (number_of_fluxes x 1)
# uptake_array	 - stoichiometrix_matrix*flux_array (number_of_species x 1)
# exit_flag	 - Exit flag returned by GLPK
# ----------------------------------------------------------------------------------- #

# Get the stoichiometric_matrix from data_dictionary -
stoichiometric_matrix = data_dictionary["STOICHIOMETRIC_MATRIX"];
(number_of_species,number_of_fluxes) = size(stoichiometric_matrix);

# Call the control function -
control_array = Control(time, species_abundance_array, data_dictionary);

# Setup the GLPK problem -
lp_problem = GLPK.Prob();
GLPK.set_prob_name(lp_problem, "sample");
GLPK.set_obj_name(lp_problem, "objective")

# Are we doing min -or- max?
min_flag = data_dictionary["MIN_FLAG"];
if min_flag == true
	GLPK.set_obj_dir(lp_problem, GLPK.MIN);
else
	GLPK.set_obj_dir(lp_problem, GLPK.MAX);
end

# Set the number of constraints and fluxes -
GLPK.add_rows(lp_problem, number_of_species);
GLPK.add_cols(lp_problem, number_of_fluxes);

# Setup flux bounds, and objective function -
flux_model_dictionary = data_dictionary["FLUX_MODEL_DICTIONARY"];
for (key,flux_model::FluxModel) in flux_model_dictionary

	# Get the default flux bounds and name -
	flux_index = flux_model.flux_index;
	flux_symbol = flux_model.flux_symbol;
	obj_coeff = flux_model.flux_obj_coeff;

	# Update the bounds for this flux -
	bounds_function = flux_model.flux_bounds_model;
	(flux_lower_bound, flux_upper_bound, flux_constraint_type) = bounds_function(key, flux_model, species_abundance_array, control_array[flux_index]);

	# Set the bounds in GLPK -
	GLPK.set_col_name(lp_problem, flux_index, flux_symbol);
	GLPK.set_col_bnds(lp_problem, flux_index, flux_constraint_type, flux_lower_bound, flux_upper_bound);

	# Set the objective function value in GLPK -
	GLPK.set_obj_coef(lp_problem, flux_index, obj_coeff);
end

# Setup problem constraints for the metabolites -
species_model_dictionary = data_dictionary["SPECIES_MODEL_DICTIONARY"];
for (key,species_model::SpeciesModel) in species_model_dictionary

	# Get data for the GLPK problem from the species model array -
	species_index = species_model.species_index;
	species_symbol = species_model.species_symbol;
	species_constraint_type = species_model.species_constraint_type;
	tau = species_model.species_time_constant;

	# Are we solving for a steady-state flux distribution?
	species_lower_bound = species_model.species_lower_bound;
	species_upper_bound = species_model.species_upper_bound;
	if steady_state_flag == false

		if (species_model.is_species_extracellular == false)
			species_lower_bound = -(1.0/step_size)*(tau - dsa[species_index]*step_size*specific_growth_rate)*species_abundance_array[species_index];
			species_constraint_type = GLPK.LO;
		else
			species_lower_bound = -(1.0/(cellmass*step_size))*(tau*species_abundance_array[species_index])
			species_constraint_type = GLPK.LO;
		end
	end

	# Is this species measured?
	if species_model.is_species_measured == true

		measurement_array = species_model.species_measurement_array;
		measured_value = measurement_array[time_index,2];
		measured_value_std = measurement_array[time_index,3];
		measured_value_upper_bound = measured_value + measured_value_std;
		measured_value_lower_bound = measured_value - measured_value_std;

		# Check: is lower bound negative?
		if (measured_value_lower_bound<0.0)
			measured_value_lower_bound = 0.0;
		end

		species_lower_bound =  (1.0/step_size)*(measured_value_lower_bound - (tau - dsa[species_index]*step_size*specific_growth_rate)*species_abundance_array[species_index]);
		species_upper_bound =  (1.0/step_size)*(measured_value_upper_bound - (tau - dsa[species_index]*step_size*specific_growth_rate)*species_abundance_array[species_index]);
		species_constraint_type = GLPK.DB;
	end

	# Set the species bounds in GLPK -
	GLPK.set_row_name(lp_problem, species_index, species_symbol);
	GLPK.set_row_bnds(lp_problem, species_index, species_constraint_type, species_lower_bound, species_upper_bound);
end



# Setup the stoichiometric array -
counter = 1;
row_index_array = zeros(Int,number_of_species*number_of_fluxes);
col_index_array = zeros(Int,number_of_species*number_of_fluxes);
species_index_vector = collect(1:number_of_species);
flux_index_vector = collect(1:number_of_fluxes);
flat_stoichiometric_array = zeros(Float64,number_of_species*number_of_fluxes);
for species_index in species_index_vector
	for flux_index in flux_index_vector
		row_index_array[counter] = species_index;
		col_index_array[counter] = flux_index;
		flat_stoichiometric_array[counter] = stoichiometric_matrix[species_index,flux_index];
		counter+=1;
	end
end
GLPK.load_matrix(lp_problem, number_of_species*number_of_fluxes, row_index_array, col_index_array, flat_stoichiometric_array);

# Set solver parameters
solver_parameters = GLPK.SimplexParam();
solver_parameters.msg_lev = GLPK.MSG_ERR;
solver_parameters.presolve = GLPK.ON;
GLPK.init_smcp(solver_parameters);

# Call the solver -
exit_flag = GLPK.simplex(lp_problem, solver_parameters);

# Get the objective function value -
objective_value = GLPK.get_obj_val(lp_problem);

# Get the calculated flux values from GLPK -
calculated_flux_array = zeros(Float64,number_of_fluxes);
for flux_index in flux_index_vector
	calculated_flux_array[flux_index] = GLPK.get_col_prim(lp_problem, flux_index);
end

# Get the dual values -
dual_value_array = zeros(Float64,number_of_fluxes);
for flux_index in flux_index_vector
	dual_value_array[flux_index] = GLPK.get_col_dual(lp_problem, flux_index);
end

# Calculate the uptake array -
uptake_array = stoichiometric_matrix*calculated_flux_array;

# Formulate the return tuple -
return (objective_value, calculated_flux_array, dual_value_array, uptake_array, exit_flag);
end
