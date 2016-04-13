include("Types.jl")
include("Bounds.jl")
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

function DataFile(TSTART,TSTOP,Ts)
# ----------------------------------------------------------------------------------- #
# DataFile.jl was generated using the Kwatee code generation system.
# DataFile: Stores model parameters as key - value pairs in a Julia Dict()
# Username: jeffreyvarner
# Type: NFBA-JULIA
# Version: 1.0
# Generation timestamp: 04-07-2016 18:07:24
#
# Input arguments:
# TSTART  - Time start
# TSTOP  - Time stop
# Ts - Time step
#
# Return arguments:
# data_dictionary  - Data dictionary instance (holds model parameters)
# ----------------------------------------------------------------------------------- #

stoichiometric_matrix = float(open(readdlm,"../network/Network.dat"));
(number_of_species,number_of_fluxes) = size(stoichiometric_matrix);

# Generate the species model array -
species_model_dictionary = buildSpeciesModelDictionary();

# Generate the flux model array -
flux_model_dictionary = buildFluxModelDictionary();

# Set the min or max flag (default is min) -
min_flag = true;

# Formulate control parameter array -
control_parameter_array = zeros(0,2);

# Blocked reactions -
blocked_reaction_set = Set{AbstractString}()
push!(blocked_reaction_set,"R_arg_deg1")
push!(blocked_reaction_set,"R_arg_deg2")
push!(blocked_reaction_set,"R_arg_deg3")
push!(blocked_reaction_set,"R_arg_deg4")
push!(blocked_reaction_set,"R_arg_deg5")
push!(blocked_reaction_set,"R_arg_deg6")
push!(blocked_reaction_set,"R_thr_deg1")
push!(blocked_reaction_set,"R_thr_deg2")
push!(blocked_reaction_set,"R_thr_deg3")
push!(blocked_reaction_set,"R_thr_deg4")
push!(blocked_reaction_set,"R_thr_deg5")
push!(blocked_reaction_set,"R_asn_deg")
push!(blocked_reaction_set,"R_gaba_deg1")
push!(blocked_reaction_set,"R_gaba_deg2")
push!(blocked_reaction_set,"R_glu_deg")
push!(blocked_reaction_set,"R_gln_deg")
push!(blocked_reaction_set,"R_lys_deg")
push!(blocked_reaction_set,"R_ala_deg")
push!(blocked_reaction_set,"R_cys_deg")
push!(blocked_reaction_set,"R_trp_deg")
push!(blocked_reaction_set,"R_pro_deg")
push!(blocked_reaction_set,"R_ser_deg")
push!(blocked_reaction_set,"R_PTS")

# setup promoter models -
p70a_model = zeros(6);
p70a_model[1] = 1.0   # Hill parameter
p70a_model[2] = 130   # KD
p70a_model[3] = 0.2   # FMAX
p70a_model[4] = 0.014 # K1
p70a_model[5] = 10    # K2

# Setup the mRNA elongation rate, and global translation
RNAP_copies_per_cell = 1900                           # Copies/cell
RNAP_elongation_rate = 65                             # NT/s (ACS SynBio Garamella 2016)
RIBOSOME_copies_per_cell = 23600                       # Copies/cell (upper bound at 72,000)
RIBOSOME_elongation_rate = 16                         # AA/s (ACS SynBio Garamella 2016)
number_of_cells = 2e8;                                # 1e9 cells/ml

# ---------------------------- DO NOT EDIT BELOW THIS LINE -------------------------- #
data_dictionary = Dict();
data_dictionary["STOICHIOMETRIC_MATRIX"] = stoichiometric_matrix;
data_dictionary["CONTROL_PARAMETER_ARRAY"] = control_parameter_array;
data_dictionary["MIN_FLAG"] = min_flag;
data_dictionary["SPECIES_MODEL_DICTIONARY"] = species_model_dictionary;
data_dictionary["FLUX_MODEL_DICTIONARY"] = flux_model_dictionary;
data_dictionary["NUMBER_OF_SPECIES"] = number_of_species;
data_dictionary["NUMBER_OF_FLUXES"] = number_of_fluxes;

data_dictionary["p70a_model"] = p70a_model;
data_dictionary["RNAP_copies_per_cell"] = RNAP_copies_per_cell;
data_dictionary["RNAP_elongation_rate"] = RNAP_elongation_rate;
data_dictionary["RIBOSOME_copies_per_cell"] = RIBOSOME_copies_per_cell;
data_dictionary["RIBOSOME_elongation_rate"] = RIBOSOME_elongation_rate;
data_dictionary["number_of_cells"] = number_of_cells;
data_dictionary["deGFP_mRNA_length"] = 683;
data_dictionary["deGFP_protein_length"] = 229;
data_dictionary["deGFP_gene_copies"] = 3.125e10;
data_dictionary["mRNA_degradation_constant"] = 2.44;
data_dictionary["polysome_gain"] = 16;
data_dictionary["volume"] = 10e-6;

# Inducer levels -
data_dictionary["p70_level"] = 35;

# DNA saturation coeffcient -
data_dictionary["plasmid_saturation_coefficient"] = 4.0

# Reactions that we block ...
data_dictionary["blocked_reaction_set"] = blocked_reaction_set
# ----------------------------------------------------------------------------------- #
return data_dictionary;
end

# ----------------------------------------------------------------------------------- #
# Helper function: buildSpeciesModelDictionary
# Constructs a dictionary of species models
# Generated using the Kwatee code generation system
#
# Input arguments:
# N/A
#
# Return arguments:
# species_model_dictionary  - Dictionary of SpeciesModels key'd by species symbol
# ----------------------------------------------------------------------------------- #
function buildSpeciesModelDictionary()

# function variables -
species_model_dictionary = Dict{AbstractString,SpeciesModel}();

# species_symbol: 1 M_maltose_c -
M_maltose_c_model = SpeciesModel();
M_maltose_c_model.species_index = 1;
M_maltose_c_model.species_symbol = string("M_maltose_c");
M_maltose_c_model.species_constraint_type = GLPK.FX;
M_maltose_c_model.species_lower_bound = 0.0;
M_maltose_c_model.species_upper_bound = 0.0;
M_maltose_c_model.is_species_measured = false;
M_maltose_c_model.is_biomass_precursor = false;
M_maltose_c_model.biomass_precursor_coefficient = 0.0;
M_maltose_c_model.species_time_constant = 1.0;
M_maltose_c_model.species_initial_condition = 0.0;
M_maltose_c_model.is_species_diluted = true;
M_maltose_c_model.is_species_extracellular = false;
species_model_dictionary["M_maltose_c"] = M_maltose_c_model;
M_maltose_c_model = 0;

# species_symbol: 2 M_h2o_c -
M_h2o_c_model = SpeciesModel();
M_h2o_c_model.species_index = 2;
M_h2o_c_model.species_symbol = string("M_h2o_c");
M_h2o_c_model.species_constraint_type = GLPK.FX;
M_h2o_c_model.species_lower_bound = 0.0;
M_h2o_c_model.species_upper_bound = 0.0;
M_h2o_c_model.is_species_measured = false;
M_h2o_c_model.is_biomass_precursor = false;
M_h2o_c_model.biomass_precursor_coefficient = 0.0;
M_h2o_c_model.species_time_constant = 1.0;
M_h2o_c_model.species_initial_condition = 0.0;
M_h2o_c_model.is_species_diluted = true;
M_h2o_c_model.is_species_extracellular = false;
species_model_dictionary["M_h2o_c"] = M_h2o_c_model;
M_h2o_c_model = 0;

# species_symbol: 3 M_glc_D_c -
M_glc_D_c_model = SpeciesModel();
M_glc_D_c_model.species_index = 3;
M_glc_D_c_model.species_symbol = string("M_glc_D_c");
M_glc_D_c_model.species_constraint_type = GLPK.FX;
M_glc_D_c_model.species_lower_bound = 0.0;
M_glc_D_c_model.species_upper_bound = 0.0;
M_glc_D_c_model.is_species_measured = false;
M_glc_D_c_model.is_biomass_precursor = false;
M_glc_D_c_model.biomass_precursor_coefficient = 0.0;
M_glc_D_c_model.species_time_constant = 1.0;
M_glc_D_c_model.species_initial_condition = 0.0;
M_glc_D_c_model.is_species_diluted = true;
M_glc_D_c_model.is_species_extracellular = false;
species_model_dictionary["M_glc_D_c"] = M_glc_D_c_model;
M_glc_D_c_model = 0;

# species_symbol: 4 M_pep_c -
M_pep_c_model = SpeciesModel();
M_pep_c_model.species_index = 4;
M_pep_c_model.species_symbol = string("M_pep_c");
M_pep_c_model.species_constraint_type = GLPK.FX;
M_pep_c_model.species_lower_bound = 0.0;
M_pep_c_model.species_upper_bound = 0.0;
M_pep_c_model.is_species_measured = false;
M_pep_c_model.is_biomass_precursor = false;
M_pep_c_model.biomass_precursor_coefficient = 0.0;
M_pep_c_model.species_time_constant = 1.0;
M_pep_c_model.species_initial_condition = 0.0;
M_pep_c_model.is_species_diluted = true;
M_pep_c_model.is_species_extracellular = false;
species_model_dictionary["M_pep_c"] = M_pep_c_model;
M_pep_c_model = 0;

# species_symbol: 5 M_g6p_c -
M_g6p_c_model = SpeciesModel();
M_g6p_c_model.species_index = 5;
M_g6p_c_model.species_symbol = string("M_g6p_c");
M_g6p_c_model.species_constraint_type = GLPK.FX;
M_g6p_c_model.species_lower_bound = 0.0;
M_g6p_c_model.species_upper_bound = 0.0;
M_g6p_c_model.is_species_measured = false;
M_g6p_c_model.is_biomass_precursor = false;
M_g6p_c_model.biomass_precursor_coefficient = 0.0;
M_g6p_c_model.species_time_constant = 1.0;
M_g6p_c_model.species_initial_condition = 0.0;
M_g6p_c_model.is_species_diluted = true;
M_g6p_c_model.is_species_extracellular = false;
species_model_dictionary["M_g6p_c"] = M_g6p_c_model;
M_g6p_c_model = 0;

# species_symbol: 6 M_pyr_c -
M_pyr_c_model = SpeciesModel();
M_pyr_c_model.species_index = 6;
M_pyr_c_model.species_symbol = string("M_pyr_c");
M_pyr_c_model.species_constraint_type = GLPK.FX;
M_pyr_c_model.species_lower_bound = 0.0;
M_pyr_c_model.species_upper_bound = 0.0;
M_pyr_c_model.is_species_measured = false;
M_pyr_c_model.is_biomass_precursor = false;
M_pyr_c_model.biomass_precursor_coefficient = 0.0;
M_pyr_c_model.species_time_constant = 1.0;
M_pyr_c_model.species_initial_condition = 0.0;
M_pyr_c_model.is_species_diluted = true;
M_pyr_c_model.is_species_extracellular = false;
species_model_dictionary["M_pyr_c"] = M_pyr_c_model;
M_pyr_c_model = 0;

# species_symbol: 7 M_atp_c -
M_atp_c_model = SpeciesModel();
M_atp_c_model.species_index = 7;
M_atp_c_model.species_symbol = string("M_atp_c");
M_atp_c_model.species_constraint_type = GLPK.FX;
M_atp_c_model.species_lower_bound = 0.0;
M_atp_c_model.species_upper_bound = 0.0;
M_atp_c_model.is_species_measured = false;
M_atp_c_model.is_biomass_precursor = false;
M_atp_c_model.biomass_precursor_coefficient = 0.0;
M_atp_c_model.species_time_constant = 1.0;
M_atp_c_model.species_initial_condition = 0.0;
M_atp_c_model.is_species_diluted = true;
M_atp_c_model.is_species_extracellular = false;
species_model_dictionary["M_atp_c"] = M_atp_c_model;
M_atp_c_model = 0;

# species_symbol: 8 M_adp_c -
M_adp_c_model = SpeciesModel();
M_adp_c_model.species_index = 8;
M_adp_c_model.species_symbol = string("M_adp_c");
M_adp_c_model.species_constraint_type = GLPK.FX;
M_adp_c_model.species_lower_bound = 0.0;
M_adp_c_model.species_upper_bound = 0.0;
M_adp_c_model.is_species_measured = false;
M_adp_c_model.is_biomass_precursor = false;
M_adp_c_model.biomass_precursor_coefficient = 0.0;
M_adp_c_model.species_time_constant = 1.0;
M_adp_c_model.species_initial_condition = 0.0;
M_adp_c_model.is_species_diluted = true;
M_adp_c_model.is_species_extracellular = false;
species_model_dictionary["M_adp_c"] = M_adp_c_model;
M_adp_c_model = 0;

# species_symbol: 9 M_h_c -
M_h_c_model = SpeciesModel();
M_h_c_model.species_index = 9;
M_h_c_model.species_symbol = string("M_h_c");
M_h_c_model.species_constraint_type = GLPK.FX;
M_h_c_model.species_lower_bound = 0.0;
M_h_c_model.species_upper_bound = 0.0;
M_h_c_model.is_species_measured = false;
M_h_c_model.is_biomass_precursor = false;
M_h_c_model.biomass_precursor_coefficient = 0.0;
M_h_c_model.species_time_constant = 1.0;
M_h_c_model.species_initial_condition = 0.0;
M_h_c_model.is_species_diluted = true;
M_h_c_model.is_species_extracellular = false;
species_model_dictionary["M_h_c"] = M_h_c_model;
M_h_c_model = 0;

# species_symbol: 10 M_utp_c -
M_utp_c_model = SpeciesModel();
M_utp_c_model.species_index = 10;
M_utp_c_model.species_symbol = string("M_utp_c");
M_utp_c_model.species_constraint_type = GLPK.FX;
M_utp_c_model.species_lower_bound = 0.0;
M_utp_c_model.species_upper_bound = 0.0;
M_utp_c_model.is_species_measured = false;
M_utp_c_model.is_biomass_precursor = false;
M_utp_c_model.biomass_precursor_coefficient = 0.0;
M_utp_c_model.species_time_constant = 1.0;
M_utp_c_model.species_initial_condition = 0.0;
M_utp_c_model.is_species_diluted = true;
M_utp_c_model.is_species_extracellular = false;
species_model_dictionary["M_utp_c"] = M_utp_c_model;
M_utp_c_model = 0;

# species_symbol: 11 M_udp_c -
M_udp_c_model = SpeciesModel();
M_udp_c_model.species_index = 11;
M_udp_c_model.species_symbol = string("M_udp_c");
M_udp_c_model.species_constraint_type = GLPK.FX;
M_udp_c_model.species_lower_bound = 0.0;
M_udp_c_model.species_upper_bound = 0.0;
M_udp_c_model.is_species_measured = false;
M_udp_c_model.is_biomass_precursor = false;
M_udp_c_model.biomass_precursor_coefficient = 0.0;
M_udp_c_model.species_time_constant = 1.0;
M_udp_c_model.species_initial_condition = 0.0;
M_udp_c_model.is_species_diluted = true;
M_udp_c_model.is_species_extracellular = false;
species_model_dictionary["M_udp_c"] = M_udp_c_model;
M_udp_c_model = 0;

# species_symbol: 12 M_ctp_c -
M_ctp_c_model = SpeciesModel();
M_ctp_c_model.species_index = 12;
M_ctp_c_model.species_symbol = string("M_ctp_c");
M_ctp_c_model.species_constraint_type = GLPK.FX;
M_ctp_c_model.species_lower_bound = 0.0;
M_ctp_c_model.species_upper_bound = 0.0;
M_ctp_c_model.is_species_measured = false;
M_ctp_c_model.is_biomass_precursor = false;
M_ctp_c_model.biomass_precursor_coefficient = 0.0;
M_ctp_c_model.species_time_constant = 1.0;
M_ctp_c_model.species_initial_condition = 0.0;
M_ctp_c_model.is_species_diluted = true;
M_ctp_c_model.is_species_extracellular = false;
species_model_dictionary["M_ctp_c"] = M_ctp_c_model;
M_ctp_c_model = 0;

# species_symbol: 13 M_cdp_c -
M_cdp_c_model = SpeciesModel();
M_cdp_c_model.species_index = 13;
M_cdp_c_model.species_symbol = string("M_cdp_c");
M_cdp_c_model.species_constraint_type = GLPK.FX;
M_cdp_c_model.species_lower_bound = 0.0;
M_cdp_c_model.species_upper_bound = 0.0;
M_cdp_c_model.is_species_measured = false;
M_cdp_c_model.is_biomass_precursor = false;
M_cdp_c_model.biomass_precursor_coefficient = 0.0;
M_cdp_c_model.species_time_constant = 1.0;
M_cdp_c_model.species_initial_condition = 0.0;
M_cdp_c_model.is_species_diluted = true;
M_cdp_c_model.is_species_extracellular = false;
species_model_dictionary["M_cdp_c"] = M_cdp_c_model;
M_cdp_c_model = 0;

# species_symbol: 14 M_gtp_c -
M_gtp_c_model = SpeciesModel();
M_gtp_c_model.species_index = 14;
M_gtp_c_model.species_symbol = string("M_gtp_c");
M_gtp_c_model.species_constraint_type = GLPK.FX;
M_gtp_c_model.species_lower_bound = 0.0;
M_gtp_c_model.species_upper_bound = 0.0;
M_gtp_c_model.is_species_measured = false;
M_gtp_c_model.is_biomass_precursor = false;
M_gtp_c_model.biomass_precursor_coefficient = 0.0;
M_gtp_c_model.species_time_constant = 1.0;
M_gtp_c_model.species_initial_condition = 0.0;
M_gtp_c_model.is_species_diluted = true;
M_gtp_c_model.is_species_extracellular = false;
species_model_dictionary["M_gtp_c"] = M_gtp_c_model;
M_gtp_c_model = 0;

# species_symbol: 15 M_gdp_c -
M_gdp_c_model = SpeciesModel();
M_gdp_c_model.species_index = 15;
M_gdp_c_model.species_symbol = string("M_gdp_c");
M_gdp_c_model.species_constraint_type = GLPK.FX;
M_gdp_c_model.species_lower_bound = 0.0;
M_gdp_c_model.species_upper_bound = 0.0;
M_gdp_c_model.is_species_measured = false;
M_gdp_c_model.is_biomass_precursor = false;
M_gdp_c_model.biomass_precursor_coefficient = 0.0;
M_gdp_c_model.species_time_constant = 1.0;
M_gdp_c_model.species_initial_condition = 0.0;
M_gdp_c_model.is_species_diluted = true;
M_gdp_c_model.is_species_extracellular = false;
species_model_dictionary["M_gdp_c"] = M_gdp_c_model;
M_gdp_c_model = 0;

# species_symbol: 16 M_f6p_c -
M_f6p_c_model = SpeciesModel();
M_f6p_c_model.species_index = 16;
M_f6p_c_model.species_symbol = string("M_f6p_c");
M_f6p_c_model.species_constraint_type = GLPK.FX;
M_f6p_c_model.species_lower_bound = 0.0;
M_f6p_c_model.species_upper_bound = 0.0;
M_f6p_c_model.is_species_measured = false;
M_f6p_c_model.is_biomass_precursor = false;
M_f6p_c_model.biomass_precursor_coefficient = 0.0;
M_f6p_c_model.species_time_constant = 1.0;
M_f6p_c_model.species_initial_condition = 0.0;
M_f6p_c_model.is_species_diluted = true;
M_f6p_c_model.is_species_extracellular = false;
species_model_dictionary["M_f6p_c"] = M_f6p_c_model;
M_f6p_c_model = 0;

# species_symbol: 17 M_fdp_c -
M_fdp_c_model = SpeciesModel();
M_fdp_c_model.species_index = 17;
M_fdp_c_model.species_symbol = string("M_fdp_c");
M_fdp_c_model.species_constraint_type = GLPK.FX;
M_fdp_c_model.species_lower_bound = 0.0;
M_fdp_c_model.species_upper_bound = 0.0;
M_fdp_c_model.is_species_measured = false;
M_fdp_c_model.is_biomass_precursor = false;
M_fdp_c_model.biomass_precursor_coefficient = 0.0;
M_fdp_c_model.species_time_constant = 1.0;
M_fdp_c_model.species_initial_condition = 0.0;
M_fdp_c_model.is_species_diluted = true;
M_fdp_c_model.is_species_extracellular = false;
species_model_dictionary["M_fdp_c"] = M_fdp_c_model;
M_fdp_c_model = 0;

# species_symbol: 18 M_pi_c -
M_pi_c_model = SpeciesModel();
M_pi_c_model.species_index = 18;
M_pi_c_model.species_symbol = string("M_pi_c");
M_pi_c_model.species_constraint_type = GLPK.FX;
M_pi_c_model.species_lower_bound = 0.0;
M_pi_c_model.species_upper_bound = 0.0;
M_pi_c_model.is_species_measured = false;
M_pi_c_model.is_biomass_precursor = false;
M_pi_c_model.biomass_precursor_coefficient = 0.0;
M_pi_c_model.species_time_constant = 1.0;
M_pi_c_model.species_initial_condition = 0.0;
M_pi_c_model.is_species_diluted = true;
M_pi_c_model.is_species_extracellular = false;
species_model_dictionary["M_pi_c"] = M_pi_c_model;
M_pi_c_model = 0;

# species_symbol: 19 M_dhap_c -
M_dhap_c_model = SpeciesModel();
M_dhap_c_model.species_index = 19;
M_dhap_c_model.species_symbol = string("M_dhap_c");
M_dhap_c_model.species_constraint_type = GLPK.FX;
M_dhap_c_model.species_lower_bound = 0.0;
M_dhap_c_model.species_upper_bound = 0.0;
M_dhap_c_model.is_species_measured = false;
M_dhap_c_model.is_biomass_precursor = false;
M_dhap_c_model.biomass_precursor_coefficient = 0.0;
M_dhap_c_model.species_time_constant = 1.0;
M_dhap_c_model.species_initial_condition = 0.0;
M_dhap_c_model.is_species_diluted = true;
M_dhap_c_model.is_species_extracellular = false;
species_model_dictionary["M_dhap_c"] = M_dhap_c_model;
M_dhap_c_model = 0;

# species_symbol: 20 M_g3p_c -
M_g3p_c_model = SpeciesModel();
M_g3p_c_model.species_index = 20;
M_g3p_c_model.species_symbol = string("M_g3p_c");
M_g3p_c_model.species_constraint_type = GLPK.FX;
M_g3p_c_model.species_lower_bound = 0.0;
M_g3p_c_model.species_upper_bound = 0.0;
M_g3p_c_model.is_species_measured = false;
M_g3p_c_model.is_biomass_precursor = false;
M_g3p_c_model.biomass_precursor_coefficient = 0.0;
M_g3p_c_model.species_time_constant = 1.0;
M_g3p_c_model.species_initial_condition = 0.0;
M_g3p_c_model.is_species_diluted = true;
M_g3p_c_model.is_species_extracellular = false;
species_model_dictionary["M_g3p_c"] = M_g3p_c_model;
M_g3p_c_model = 0;

# species_symbol: 21 M_nad_c -
M_nad_c_model = SpeciesModel();
M_nad_c_model.species_index = 21;
M_nad_c_model.species_symbol = string("M_nad_c");
M_nad_c_model.species_constraint_type = GLPK.FX;
M_nad_c_model.species_lower_bound = 0.0;
M_nad_c_model.species_upper_bound = 0.0;
M_nad_c_model.is_species_measured = false;
M_nad_c_model.is_biomass_precursor = false;
M_nad_c_model.biomass_precursor_coefficient = 0.0;
M_nad_c_model.species_time_constant = 1.0;
M_nad_c_model.species_initial_condition = 0.0;
M_nad_c_model.is_species_diluted = true;
M_nad_c_model.is_species_extracellular = false;
species_model_dictionary["M_nad_c"] = M_nad_c_model;
M_nad_c_model = 0;

# species_symbol: 22 M_13dpg_c -
M_13dpg_c_model = SpeciesModel();
M_13dpg_c_model.species_index = 22;
M_13dpg_c_model.species_symbol = string("M_13dpg_c");
M_13dpg_c_model.species_constraint_type = GLPK.FX;
M_13dpg_c_model.species_lower_bound = 0.0;
M_13dpg_c_model.species_upper_bound = 0.0;
M_13dpg_c_model.is_species_measured = false;
M_13dpg_c_model.is_biomass_precursor = false;
M_13dpg_c_model.biomass_precursor_coefficient = 0.0;
M_13dpg_c_model.species_time_constant = 1.0;
M_13dpg_c_model.species_initial_condition = 0.0;
M_13dpg_c_model.is_species_diluted = true;
M_13dpg_c_model.is_species_extracellular = false;
species_model_dictionary["M_13dpg_c"] = M_13dpg_c_model;
M_13dpg_c_model = 0;

# species_symbol: 23 M_nadh_c -
M_nadh_c_model = SpeciesModel();
M_nadh_c_model.species_index = 23;
M_nadh_c_model.species_symbol = string("M_nadh_c");
M_nadh_c_model.species_constraint_type = GLPK.FX;
M_nadh_c_model.species_lower_bound = 0.0;
M_nadh_c_model.species_upper_bound = 0.0;
M_nadh_c_model.is_species_measured = false;
M_nadh_c_model.is_biomass_precursor = false;
M_nadh_c_model.biomass_precursor_coefficient = 0.0;
M_nadh_c_model.species_time_constant = 1.0;
M_nadh_c_model.species_initial_condition = 0.0;
M_nadh_c_model.is_species_diluted = true;
M_nadh_c_model.is_species_extracellular = false;
species_model_dictionary["M_nadh_c"] = M_nadh_c_model;
M_nadh_c_model = 0;

# species_symbol: 24 M_3pg_c -
M_3pg_c_model = SpeciesModel();
M_3pg_c_model.species_index = 24;
M_3pg_c_model.species_symbol = string("M_3pg_c");
M_3pg_c_model.species_constraint_type = GLPK.FX;
M_3pg_c_model.species_lower_bound = 0.0;
M_3pg_c_model.species_upper_bound = 0.0;
M_3pg_c_model.is_species_measured = false;
M_3pg_c_model.is_biomass_precursor = false;
M_3pg_c_model.biomass_precursor_coefficient = 0.0;
M_3pg_c_model.species_time_constant = 1.0;
M_3pg_c_model.species_initial_condition = 0.0;
M_3pg_c_model.is_species_diluted = true;
M_3pg_c_model.is_species_extracellular = false;
species_model_dictionary["M_3pg_c"] = M_3pg_c_model;
M_3pg_c_model = 0;

# species_symbol: 25 M_2pg_c -
M_2pg_c_model = SpeciesModel();
M_2pg_c_model.species_index = 25;
M_2pg_c_model.species_symbol = string("M_2pg_c");
M_2pg_c_model.species_constraint_type = GLPK.FX;
M_2pg_c_model.species_lower_bound = 0.0;
M_2pg_c_model.species_upper_bound = 0.0;
M_2pg_c_model.is_species_measured = false;
M_2pg_c_model.is_biomass_precursor = false;
M_2pg_c_model.biomass_precursor_coefficient = 0.0;
M_2pg_c_model.species_time_constant = 1.0;
M_2pg_c_model.species_initial_condition = 0.0;
M_2pg_c_model.is_species_diluted = true;
M_2pg_c_model.is_species_extracellular = false;
species_model_dictionary["M_2pg_c"] = M_2pg_c_model;
M_2pg_c_model = 0;

# species_symbol: 26 M_oaa_c -
M_oaa_c_model = SpeciesModel();
M_oaa_c_model.species_index = 26;
M_oaa_c_model.species_symbol = string("M_oaa_c");
M_oaa_c_model.species_constraint_type = GLPK.FX;
M_oaa_c_model.species_lower_bound = 0.0;
M_oaa_c_model.species_upper_bound = 0.0;
M_oaa_c_model.is_species_measured = false;
M_oaa_c_model.is_biomass_precursor = false;
M_oaa_c_model.biomass_precursor_coefficient = 0.0;
M_oaa_c_model.species_time_constant = 1.0;
M_oaa_c_model.species_initial_condition = 0.0;
M_oaa_c_model.is_species_diluted = true;
M_oaa_c_model.is_species_extracellular = false;
species_model_dictionary["M_oaa_c"] = M_oaa_c_model;
M_oaa_c_model = 0;

# species_symbol: 27 M_co2_c -
M_co2_c_model = SpeciesModel();
M_co2_c_model.species_index = 27;
M_co2_c_model.species_symbol = string("M_co2_c");
M_co2_c_model.species_constraint_type = GLPK.FX;
M_co2_c_model.species_lower_bound = 0.0;
M_co2_c_model.species_upper_bound = 0.0;
M_co2_c_model.is_species_measured = false;
M_co2_c_model.is_biomass_precursor = false;
M_co2_c_model.biomass_precursor_coefficient = 0.0;
M_co2_c_model.species_time_constant = 1.0;
M_co2_c_model.species_initial_condition = 0.0;
M_co2_c_model.is_species_diluted = true;
M_co2_c_model.is_species_extracellular = false;
species_model_dictionary["M_co2_c"] = M_co2_c_model;
M_co2_c_model = 0;

# species_symbol: 28 M_coa_c -
M_coa_c_model = SpeciesModel();
M_coa_c_model.species_index = 28;
M_coa_c_model.species_symbol = string("M_coa_c");
M_coa_c_model.species_constraint_type = GLPK.FX;
M_coa_c_model.species_lower_bound = 0.0;
M_coa_c_model.species_upper_bound = 0.0;
M_coa_c_model.is_species_measured = false;
M_coa_c_model.is_biomass_precursor = false;
M_coa_c_model.biomass_precursor_coefficient = 0.0;
M_coa_c_model.species_time_constant = 1.0;
M_coa_c_model.species_initial_condition = 0.0;
M_coa_c_model.is_species_diluted = true;
M_coa_c_model.is_species_extracellular = false;
species_model_dictionary["M_coa_c"] = M_coa_c_model;
M_coa_c_model = 0;

# species_symbol: 29 M_accoa_c -
M_accoa_c_model = SpeciesModel();
M_accoa_c_model.species_index = 29;
M_accoa_c_model.species_symbol = string("M_accoa_c");
M_accoa_c_model.species_constraint_type = GLPK.FX;
M_accoa_c_model.species_lower_bound = 0.0;
M_accoa_c_model.species_upper_bound = 0.0;
M_accoa_c_model.is_species_measured = false;
M_accoa_c_model.is_biomass_precursor = false;
M_accoa_c_model.biomass_precursor_coefficient = 0.0;
M_accoa_c_model.species_time_constant = 1.0;
M_accoa_c_model.species_initial_condition = 0.0;
M_accoa_c_model.is_species_diluted = true;
M_accoa_c_model.is_species_extracellular = false;
species_model_dictionary["M_accoa_c"] = M_accoa_c_model;
M_accoa_c_model = 0;

# species_symbol: 30 M_amp_c -
M_amp_c_model = SpeciesModel();
M_amp_c_model.species_index = 30;
M_amp_c_model.species_symbol = string("M_amp_c");
M_amp_c_model.species_constraint_type = GLPK.FX;
M_amp_c_model.species_lower_bound = 0.0;
M_amp_c_model.species_upper_bound = 0.0;
M_amp_c_model.is_species_measured = false;
M_amp_c_model.is_biomass_precursor = false;
M_amp_c_model.biomass_precursor_coefficient = 0.0;
M_amp_c_model.species_time_constant = 1.0;
M_amp_c_model.species_initial_condition = 0.0;
M_amp_c_model.is_species_diluted = true;
M_amp_c_model.is_species_extracellular = false;
species_model_dictionary["M_amp_c"] = M_amp_c_model;
M_amp_c_model = 0;

# species_symbol: 31 M_nadp_c -
M_nadp_c_model = SpeciesModel();
M_nadp_c_model.species_index = 31;
M_nadp_c_model.species_symbol = string("M_nadp_c");
M_nadp_c_model.species_constraint_type = GLPK.FX;
M_nadp_c_model.species_lower_bound = 0.0;
M_nadp_c_model.species_upper_bound = 0.0;
M_nadp_c_model.is_species_measured = false;
M_nadp_c_model.is_biomass_precursor = false;
M_nadp_c_model.biomass_precursor_coefficient = 0.0;
M_nadp_c_model.species_time_constant = 1.0;
M_nadp_c_model.species_initial_condition = 0.0;
M_nadp_c_model.is_species_diluted = true;
M_nadp_c_model.is_species_extracellular = false;
species_model_dictionary["M_nadp_c"] = M_nadp_c_model;
M_nadp_c_model = 0;

# species_symbol: 32 M_6pgl_c -
M_6pgl_c_model = SpeciesModel();
M_6pgl_c_model.species_index = 32;
M_6pgl_c_model.species_symbol = string("M_6pgl_c");
M_6pgl_c_model.species_constraint_type = GLPK.FX;
M_6pgl_c_model.species_lower_bound = 0.0;
M_6pgl_c_model.species_upper_bound = 0.0;
M_6pgl_c_model.is_species_measured = false;
M_6pgl_c_model.is_biomass_precursor = false;
M_6pgl_c_model.biomass_precursor_coefficient = 0.0;
M_6pgl_c_model.species_time_constant = 1.0;
M_6pgl_c_model.species_initial_condition = 0.0;
M_6pgl_c_model.is_species_diluted = true;
M_6pgl_c_model.is_species_extracellular = false;
species_model_dictionary["M_6pgl_c"] = M_6pgl_c_model;
M_6pgl_c_model = 0;

# species_symbol: 33 M_nadph_c -
M_nadph_c_model = SpeciesModel();
M_nadph_c_model.species_index = 33;
M_nadph_c_model.species_symbol = string("M_nadph_c");
M_nadph_c_model.species_constraint_type = GLPK.FX;
M_nadph_c_model.species_lower_bound = 0.0;
M_nadph_c_model.species_upper_bound = 0.0;
M_nadph_c_model.is_species_measured = false;
M_nadph_c_model.is_biomass_precursor = false;
M_nadph_c_model.biomass_precursor_coefficient = 0.0;
M_nadph_c_model.species_time_constant = 1.0;
M_nadph_c_model.species_initial_condition = 0.0;
M_nadph_c_model.is_species_diluted = true;
M_nadph_c_model.is_species_extracellular = false;
species_model_dictionary["M_nadph_c"] = M_nadph_c_model;
M_nadph_c_model = 0;

# species_symbol: 34 M_6pgc_c -
M_6pgc_c_model = SpeciesModel();
M_6pgc_c_model.species_index = 34;
M_6pgc_c_model.species_symbol = string("M_6pgc_c");
M_6pgc_c_model.species_constraint_type = GLPK.FX;
M_6pgc_c_model.species_lower_bound = 0.0;
M_6pgc_c_model.species_upper_bound = 0.0;
M_6pgc_c_model.is_species_measured = false;
M_6pgc_c_model.is_biomass_precursor = false;
M_6pgc_c_model.biomass_precursor_coefficient = 0.0;
M_6pgc_c_model.species_time_constant = 1.0;
M_6pgc_c_model.species_initial_condition = 0.0;
M_6pgc_c_model.is_species_diluted = true;
M_6pgc_c_model.is_species_extracellular = false;
species_model_dictionary["M_6pgc_c"] = M_6pgc_c_model;
M_6pgc_c_model = 0;

# species_symbol: 35 M_ru5p_D_c -
M_ru5p_D_c_model = SpeciesModel();
M_ru5p_D_c_model.species_index = 35;
M_ru5p_D_c_model.species_symbol = string("M_ru5p_D_c");
M_ru5p_D_c_model.species_constraint_type = GLPK.FX;
M_ru5p_D_c_model.species_lower_bound = 0.0;
M_ru5p_D_c_model.species_upper_bound = 0.0;
M_ru5p_D_c_model.is_species_measured = false;
M_ru5p_D_c_model.is_biomass_precursor = false;
M_ru5p_D_c_model.biomass_precursor_coefficient = 0.0;
M_ru5p_D_c_model.species_time_constant = 1.0;
M_ru5p_D_c_model.species_initial_condition = 0.0;
M_ru5p_D_c_model.is_species_diluted = true;
M_ru5p_D_c_model.is_species_extracellular = false;
species_model_dictionary["M_ru5p_D_c"] = M_ru5p_D_c_model;
M_ru5p_D_c_model = 0;

# species_symbol: 36 M_xu5p_D_c -
M_xu5p_D_c_model = SpeciesModel();
M_xu5p_D_c_model.species_index = 36;
M_xu5p_D_c_model.species_symbol = string("M_xu5p_D_c");
M_xu5p_D_c_model.species_constraint_type = GLPK.FX;
M_xu5p_D_c_model.species_lower_bound = 0.0;
M_xu5p_D_c_model.species_upper_bound = 0.0;
M_xu5p_D_c_model.is_species_measured = false;
M_xu5p_D_c_model.is_biomass_precursor = false;
M_xu5p_D_c_model.biomass_precursor_coefficient = 0.0;
M_xu5p_D_c_model.species_time_constant = 1.0;
M_xu5p_D_c_model.species_initial_condition = 0.0;
M_xu5p_D_c_model.is_species_diluted = true;
M_xu5p_D_c_model.is_species_extracellular = false;
species_model_dictionary["M_xu5p_D_c"] = M_xu5p_D_c_model;
M_xu5p_D_c_model = 0;

# species_symbol: 37 M_r5p_c -
M_r5p_c_model = SpeciesModel();
M_r5p_c_model.species_index = 37;
M_r5p_c_model.species_symbol = string("M_r5p_c");
M_r5p_c_model.species_constraint_type = GLPK.FX;
M_r5p_c_model.species_lower_bound = 0.0;
M_r5p_c_model.species_upper_bound = 0.0;
M_r5p_c_model.is_species_measured = false;
M_r5p_c_model.is_biomass_precursor = false;
M_r5p_c_model.biomass_precursor_coefficient = 0.0;
M_r5p_c_model.species_time_constant = 1.0;
M_r5p_c_model.species_initial_condition = 0.0;
M_r5p_c_model.is_species_diluted = true;
M_r5p_c_model.is_species_extracellular = false;
species_model_dictionary["M_r5p_c"] = M_r5p_c_model;
M_r5p_c_model = 0;

# species_symbol: 38 M_s7p_c -
M_s7p_c_model = SpeciesModel();
M_s7p_c_model.species_index = 38;
M_s7p_c_model.species_symbol = string("M_s7p_c");
M_s7p_c_model.species_constraint_type = GLPK.FX;
M_s7p_c_model.species_lower_bound = 0.0;
M_s7p_c_model.species_upper_bound = 0.0;
M_s7p_c_model.is_species_measured = false;
M_s7p_c_model.is_biomass_precursor = false;
M_s7p_c_model.biomass_precursor_coefficient = 0.0;
M_s7p_c_model.species_time_constant = 1.0;
M_s7p_c_model.species_initial_condition = 0.0;
M_s7p_c_model.is_species_diluted = true;
M_s7p_c_model.is_species_extracellular = false;
species_model_dictionary["M_s7p_c"] = M_s7p_c_model;
M_s7p_c_model = 0;

# species_symbol: 39 M_e4p_c -
M_e4p_c_model = SpeciesModel();
M_e4p_c_model.species_index = 39;
M_e4p_c_model.species_symbol = string("M_e4p_c");
M_e4p_c_model.species_constraint_type = GLPK.FX;
M_e4p_c_model.species_lower_bound = 0.0;
M_e4p_c_model.species_upper_bound = 0.0;
M_e4p_c_model.is_species_measured = false;
M_e4p_c_model.is_biomass_precursor = false;
M_e4p_c_model.biomass_precursor_coefficient = 0.0;
M_e4p_c_model.species_time_constant = 1.0;
M_e4p_c_model.species_initial_condition = 0.0;
M_e4p_c_model.is_species_diluted = false;
M_e4p_c_model.is_species_extracellular = true;
species_model_dictionary["M_e4p_c"] = M_e4p_c_model;
M_e4p_c_model = 0;

# species_symbol: 40 M_2ddg6p_c -
M_2ddg6p_c_model = SpeciesModel();
M_2ddg6p_c_model.species_index = 40;
M_2ddg6p_c_model.species_symbol = string("M_2ddg6p_c");
M_2ddg6p_c_model.species_constraint_type = GLPK.FX;
M_2ddg6p_c_model.species_lower_bound = 0.0;
M_2ddg6p_c_model.species_upper_bound = 0.0;
M_2ddg6p_c_model.is_species_measured = false;
M_2ddg6p_c_model.is_biomass_precursor = false;
M_2ddg6p_c_model.biomass_precursor_coefficient = 0.0;
M_2ddg6p_c_model.species_time_constant = 1.0;
M_2ddg6p_c_model.species_initial_condition = 0.0;
M_2ddg6p_c_model.is_species_diluted = true;
M_2ddg6p_c_model.is_species_extracellular = false;
species_model_dictionary["M_2ddg6p_c"] = M_2ddg6p_c_model;
M_2ddg6p_c_model = 0;

# species_symbol: 41 M_cit_c -
M_cit_c_model = SpeciesModel();
M_cit_c_model.species_index = 41;
M_cit_c_model.species_symbol = string("M_cit_c");
M_cit_c_model.species_constraint_type = GLPK.FX;
M_cit_c_model.species_lower_bound = 0.0;
M_cit_c_model.species_upper_bound = 0.0;
M_cit_c_model.is_species_measured = false;
M_cit_c_model.is_biomass_precursor = false;
M_cit_c_model.biomass_precursor_coefficient = 0.0;
M_cit_c_model.species_time_constant = 1.0;
M_cit_c_model.species_initial_condition = 0.0;
M_cit_c_model.is_species_diluted = true;
M_cit_c_model.is_species_extracellular = false;
species_model_dictionary["M_cit_c"] = M_cit_c_model;
M_cit_c_model = 0;

# species_symbol: 42 M_icit_c -
M_icit_c_model = SpeciesModel();
M_icit_c_model.species_index = 42;
M_icit_c_model.species_symbol = string("M_icit_c");
M_icit_c_model.species_constraint_type = GLPK.FX;
M_icit_c_model.species_lower_bound = 0.0;
M_icit_c_model.species_upper_bound = 0.0;
M_icit_c_model.is_species_measured = false;
M_icit_c_model.is_biomass_precursor = false;
M_icit_c_model.biomass_precursor_coefficient = 0.0;
M_icit_c_model.species_time_constant = 1.0;
M_icit_c_model.species_initial_condition = 0.0;
M_icit_c_model.is_species_diluted = true;
M_icit_c_model.is_species_extracellular = false;
species_model_dictionary["M_icit_c"] = M_icit_c_model;
M_icit_c_model = 0;

# species_symbol: 43 M_akg_c -
M_akg_c_model = SpeciesModel();
M_akg_c_model.species_index = 43;
M_akg_c_model.species_symbol = string("M_akg_c");
M_akg_c_model.species_constraint_type = GLPK.FX;
M_akg_c_model.species_lower_bound = 0.0;
M_akg_c_model.species_upper_bound = 0.0;
M_akg_c_model.is_species_measured = false;
M_akg_c_model.is_biomass_precursor = false;
M_akg_c_model.biomass_precursor_coefficient = 0.0;
M_akg_c_model.species_time_constant = 1.0;
M_akg_c_model.species_initial_condition = 0.0;
M_akg_c_model.is_species_diluted = true;
M_akg_c_model.is_species_extracellular = false;
species_model_dictionary["M_akg_c"] = M_akg_c_model;
M_akg_c_model = 0;

# species_symbol: 44 M_succoa_c -
M_succoa_c_model = SpeciesModel();
M_succoa_c_model.species_index = 44;
M_succoa_c_model.species_symbol = string("M_succoa_c");
M_succoa_c_model.species_constraint_type = GLPK.FX;
M_succoa_c_model.species_lower_bound = 0.0;
M_succoa_c_model.species_upper_bound = 0.0;
M_succoa_c_model.is_species_measured = false;
M_succoa_c_model.is_biomass_precursor = false;
M_succoa_c_model.biomass_precursor_coefficient = 0.0;
M_succoa_c_model.species_time_constant = 1.0;
M_succoa_c_model.species_initial_condition = 0.0;
M_succoa_c_model.is_species_diluted = true;
M_succoa_c_model.is_species_extracellular = false;
species_model_dictionary["M_succoa_c"] = M_succoa_c_model;
M_succoa_c_model = 0;

# species_symbol: 45 M_succ_c -
M_succ_c_model = SpeciesModel();
M_succ_c_model.species_index = 45;
M_succ_c_model.species_symbol = string("M_succ_c");
M_succ_c_model.species_constraint_type = GLPK.FX;
M_succ_c_model.species_lower_bound = 0.0;
M_succ_c_model.species_upper_bound = 0.0;
M_succ_c_model.is_species_measured = false;
M_succ_c_model.is_biomass_precursor = false;
M_succ_c_model.biomass_precursor_coefficient = 0.0;
M_succ_c_model.species_time_constant = 1.0;
M_succ_c_model.species_initial_condition = 0.0;
M_succ_c_model.is_species_diluted = true;
M_succ_c_model.is_species_extracellular = false;
species_model_dictionary["M_succ_c"] = M_succ_c_model;
M_succ_c_model = 0;

# species_symbol: 46 M_q8_c -
M_q8_c_model = SpeciesModel();
M_q8_c_model.species_index = 46;
M_q8_c_model.species_symbol = string("M_q8_c");
M_q8_c_model.species_constraint_type = GLPK.FX;
M_q8_c_model.species_lower_bound = 0.0;
M_q8_c_model.species_upper_bound = 0.0;
M_q8_c_model.is_species_measured = false;
M_q8_c_model.is_biomass_precursor = false;
M_q8_c_model.biomass_precursor_coefficient = 0.0;
M_q8_c_model.species_time_constant = 1.0;
M_q8_c_model.species_initial_condition = 0.0;
M_q8_c_model.is_species_diluted = true;
M_q8_c_model.is_species_extracellular = false;
species_model_dictionary["M_q8_c"] = M_q8_c_model;
M_q8_c_model = 0;

# species_symbol: 47 M_fum_c -
M_fum_c_model = SpeciesModel();
M_fum_c_model.species_index = 47;
M_fum_c_model.species_symbol = string("M_fum_c");
M_fum_c_model.species_constraint_type = GLPK.FX;
M_fum_c_model.species_lower_bound = 0.0;
M_fum_c_model.species_upper_bound = 0.0;
M_fum_c_model.is_species_measured = false;
M_fum_c_model.is_biomass_precursor = false;
M_fum_c_model.biomass_precursor_coefficient = 0.0;
M_fum_c_model.species_time_constant = 1.0;
M_fum_c_model.species_initial_condition = 0.0;
M_fum_c_model.is_species_diluted = true;
M_fum_c_model.is_species_extracellular = false;
species_model_dictionary["M_fum_c"] = M_fum_c_model;
M_fum_c_model = 0;

# species_symbol: 48 M_q8h2_c -
M_q8h2_c_model = SpeciesModel();
M_q8h2_c_model.species_index = 48;
M_q8h2_c_model.species_symbol = string("M_q8h2_c");
M_q8h2_c_model.species_constraint_type = GLPK.FX;
M_q8h2_c_model.species_lower_bound = 0.0;
M_q8h2_c_model.species_upper_bound = 0.0;
M_q8h2_c_model.is_species_measured = false;
M_q8h2_c_model.is_biomass_precursor = false;
M_q8h2_c_model.biomass_precursor_coefficient = 0.0;
M_q8h2_c_model.species_time_constant = 1.0;
M_q8h2_c_model.species_initial_condition = 0.0;
M_q8h2_c_model.is_species_diluted = true;
M_q8h2_c_model.is_species_extracellular = false;
species_model_dictionary["M_q8h2_c"] = M_q8h2_c_model;
M_q8h2_c_model = 0;

# species_symbol: 49 M_mql8_c -
M_mql8_c_model = SpeciesModel();
M_mql8_c_model.species_index = 49;
M_mql8_c_model.species_symbol = string("M_mql8_c");
M_mql8_c_model.species_constraint_type = GLPK.FX;
M_mql8_c_model.species_lower_bound = 0.0;
M_mql8_c_model.species_upper_bound = 0.0;
M_mql8_c_model.is_species_measured = false;
M_mql8_c_model.is_biomass_precursor = false;
M_mql8_c_model.biomass_precursor_coefficient = 0.0;
M_mql8_c_model.species_time_constant = 1.0;
M_mql8_c_model.species_initial_condition = 0.0;
M_mql8_c_model.is_species_diluted = true;
M_mql8_c_model.is_species_extracellular = false;
species_model_dictionary["M_mql8_c"] = M_mql8_c_model;
M_mql8_c_model = 0;

# species_symbol: 50 M_mqn8_c -
M_mqn8_c_model = SpeciesModel();
M_mqn8_c_model.species_index = 50;
M_mqn8_c_model.species_symbol = string("M_mqn8_c");
M_mqn8_c_model.species_constraint_type = GLPK.FX;
M_mqn8_c_model.species_lower_bound = 0.0;
M_mqn8_c_model.species_upper_bound = 0.0;
M_mqn8_c_model.is_species_measured = false;
M_mqn8_c_model.is_biomass_precursor = false;
M_mqn8_c_model.biomass_precursor_coefficient = 0.0;
M_mqn8_c_model.species_time_constant = 1.0;
M_mqn8_c_model.species_initial_condition = 0.0;
M_mqn8_c_model.is_species_diluted = true;
M_mqn8_c_model.is_species_extracellular = false;
species_model_dictionary["M_mqn8_c"] = M_mqn8_c_model;
M_mqn8_c_model = 0;

# species_symbol: 51 M_mal_L_c -
M_mal_L_c_model = SpeciesModel();
M_mal_L_c_model.species_index = 51;
M_mal_L_c_model.species_symbol = string("M_mal_L_c");
M_mal_L_c_model.species_constraint_type = GLPK.FX;
M_mal_L_c_model.species_lower_bound = 0.0;
M_mal_L_c_model.species_upper_bound = 0.0;
M_mal_L_c_model.is_species_measured = false;
M_mal_L_c_model.is_biomass_precursor = false;
M_mal_L_c_model.biomass_precursor_coefficient = 0.0;
M_mal_L_c_model.species_time_constant = 1.0;
M_mal_L_c_model.species_initial_condition = 0.0;
M_mal_L_c_model.is_species_diluted = true;
M_mal_L_c_model.is_species_extracellular = false;
species_model_dictionary["M_mal_L_c"] = M_mal_L_c_model;
M_mal_L_c_model = 0;

# species_symbol: 52 M_glx_c -
M_glx_c_model = SpeciesModel();
M_glx_c_model.species_index = 52;
M_glx_c_model.species_symbol = string("M_glx_c");
M_glx_c_model.species_constraint_type = GLPK.FX;
M_glx_c_model.species_lower_bound = 0.0;
M_glx_c_model.species_upper_bound = 0.0;
M_glx_c_model.is_species_measured = false;
M_glx_c_model.is_biomass_precursor = false;
M_glx_c_model.biomass_precursor_coefficient = 0.0;
M_glx_c_model.species_time_constant = 1.0;
M_glx_c_model.species_initial_condition = 0.0;
M_glx_c_model.is_species_diluted = true;
M_glx_c_model.is_species_extracellular = false;
species_model_dictionary["M_glx_c"] = M_glx_c_model;
M_glx_c_model = 0;

# species_symbol: 53 M_actp_c -
M_actp_c_model = SpeciesModel();
M_actp_c_model.species_index = 53;
M_actp_c_model.species_symbol = string("M_actp_c");
M_actp_c_model.species_constraint_type = GLPK.FX;
M_actp_c_model.species_lower_bound = 0.0;
M_actp_c_model.species_upper_bound = 0.0;
M_actp_c_model.is_species_measured = false;
M_actp_c_model.is_biomass_precursor = false;
M_actp_c_model.biomass_precursor_coefficient = 0.0;
M_actp_c_model.species_time_constant = 1.0;
M_actp_c_model.species_initial_condition = 0.0;
M_actp_c_model.is_species_diluted = true;
M_actp_c_model.is_species_extracellular = false;
species_model_dictionary["M_actp_c"] = M_actp_c_model;
M_actp_c_model = 0;

# species_symbol: 54 M_ac_c -
M_ac_c_model = SpeciesModel();
M_ac_c_model.species_index = 54;
M_ac_c_model.species_symbol = string("M_ac_c");
M_ac_c_model.species_constraint_type = GLPK.FX;
M_ac_c_model.species_lower_bound = 0.0;
M_ac_c_model.species_upper_bound = 0.0;
M_ac_c_model.is_species_measured = false;
M_ac_c_model.is_biomass_precursor = false;
M_ac_c_model.biomass_precursor_coefficient = 0.0;
M_ac_c_model.species_time_constant = 1.0;
M_ac_c_model.species_initial_condition = 0.0;
M_ac_c_model.is_species_diluted = true;
M_ac_c_model.is_species_extracellular = false;
species_model_dictionary["M_ac_c"] = M_ac_c_model;
M_ac_c_model = 0;

# species_symbol: 55 M_ppi_c -
M_ppi_c_model = SpeciesModel();
M_ppi_c_model.species_index = 55;
M_ppi_c_model.species_symbol = string("M_ppi_c");
M_ppi_c_model.species_constraint_type = GLPK.FX;
M_ppi_c_model.species_lower_bound = 0.0;
M_ppi_c_model.species_upper_bound = 0.0;
M_ppi_c_model.is_species_measured = false;
M_ppi_c_model.is_biomass_precursor = false;
M_ppi_c_model.biomass_precursor_coefficient = 0.0;
M_ppi_c_model.species_time_constant = 1.0;
M_ppi_c_model.species_initial_condition = 0.0;
M_ppi_c_model.is_species_diluted = true;
M_ppi_c_model.is_species_extracellular = false;
species_model_dictionary["M_ppi_c"] = M_ppi_c_model;
M_ppi_c_model = 0;

# species_symbol: 56 M_etoh_c -
M_etoh_c_model = SpeciesModel();
M_etoh_c_model.species_index = 56;
M_etoh_c_model.species_symbol = string("M_etoh_c");
M_etoh_c_model.species_constraint_type = GLPK.FX;
M_etoh_c_model.species_lower_bound = 0.0;
M_etoh_c_model.species_upper_bound = 0.0;
M_etoh_c_model.is_species_measured = false;
M_etoh_c_model.is_biomass_precursor = false;
M_etoh_c_model.biomass_precursor_coefficient = 0.0;
M_etoh_c_model.species_time_constant = 1.0;
M_etoh_c_model.species_initial_condition = 0.0;
M_etoh_c_model.is_species_diluted = false;
M_etoh_c_model.is_species_extracellular = true;
species_model_dictionary["M_etoh_c"] = M_etoh_c_model;
M_etoh_c_model = 0;

# species_symbol: 57 M_lac_D_c -
M_lac_D_c_model = SpeciesModel();
M_lac_D_c_model.species_index = 57;
M_lac_D_c_model.species_symbol = string("M_lac_D_c");
M_lac_D_c_model.species_constraint_type = GLPK.FX;
M_lac_D_c_model.species_lower_bound = 0.0;
M_lac_D_c_model.species_upper_bound = 0.0;
M_lac_D_c_model.is_species_measured = false;
M_lac_D_c_model.is_biomass_precursor = false;
M_lac_D_c_model.biomass_precursor_coefficient = 0.0;
M_lac_D_c_model.species_time_constant = 1.0;
M_lac_D_c_model.species_initial_condition = 0.0;
M_lac_D_c_model.is_species_diluted = true;
M_lac_D_c_model.is_species_extracellular = false;
species_model_dictionary["M_lac_D_c"] = M_lac_D_c_model;
M_lac_D_c_model = 0;

# species_symbol: 58 M_for_c -
M_for_c_model = SpeciesModel();
M_for_c_model.species_index = 58;
M_for_c_model.species_symbol = string("M_for_c");
M_for_c_model.species_constraint_type = GLPK.FX;
M_for_c_model.species_lower_bound = 0.0;
M_for_c_model.species_upper_bound = 0.0;
M_for_c_model.is_species_measured = false;
M_for_c_model.is_biomass_precursor = false;
M_for_c_model.biomass_precursor_coefficient = 0.0;
M_for_c_model.species_time_constant = 1.0;
M_for_c_model.species_initial_condition = 0.0;
M_for_c_model.is_species_diluted = true;
M_for_c_model.is_species_extracellular = false;
species_model_dictionary["M_for_c"] = M_for_c_model;
M_for_c_model = 0;

# species_symbol: 59 M_o2_c -
M_o2_c_model = SpeciesModel();
M_o2_c_model.species_index = 59;
M_o2_c_model.species_symbol = string("M_o2_c");
M_o2_c_model.species_constraint_type = GLPK.FX;
M_o2_c_model.species_lower_bound = 0.0;
M_o2_c_model.species_upper_bound = 0.0;
M_o2_c_model.is_species_measured = false;
M_o2_c_model.is_biomass_precursor = false;
M_o2_c_model.biomass_precursor_coefficient = 0.0;
M_o2_c_model.species_time_constant = 1.0;
M_o2_c_model.species_initial_condition = 0.0;
M_o2_c_model.is_species_diluted = true;
M_o2_c_model.is_species_extracellular = false;
species_model_dictionary["M_o2_c"] = M_o2_c_model;
M_o2_c_model = 0;

# species_symbol: 60 M_h_e -
M_h_e_model = SpeciesModel();
M_h_e_model.species_index = 60;
M_h_e_model.species_symbol = string("M_h_e");
M_h_e_model.species_constraint_type = GLPK.FX;
M_h_e_model.species_lower_bound = 0.0;
M_h_e_model.species_upper_bound = 0.0;
M_h_e_model.is_species_measured = false;
M_h_e_model.is_biomass_precursor = false;
M_h_e_model.biomass_precursor_coefficient = 0.0;
M_h_e_model.species_time_constant = 1.0;
M_h_e_model.species_initial_condition = 0.0;
M_h_e_model.is_species_diluted = false;
M_h_e_model.is_species_extracellular = true;
species_model_dictionary["M_h_e"] = M_h_e_model;
M_h_e_model = 0;

# species_symbol: 61 M_chor_c -
M_chor_c_model = SpeciesModel();
M_chor_c_model.species_index = 61;
M_chor_c_model.species_symbol = string("M_chor_c");
M_chor_c_model.species_constraint_type = GLPK.FX;
M_chor_c_model.species_lower_bound = 0.0;
M_chor_c_model.species_upper_bound = 0.0;
M_chor_c_model.is_species_measured = false;
M_chor_c_model.is_biomass_precursor = false;
M_chor_c_model.biomass_precursor_coefficient = 0.0;
M_chor_c_model.species_time_constant = 1.0;
M_chor_c_model.species_initial_condition = 0.0;
M_chor_c_model.is_species_diluted = true;
M_chor_c_model.is_species_extracellular = false;
species_model_dictionary["M_chor_c"] = M_chor_c_model;
M_chor_c_model = 0;

# species_symbol: 62 M_gln_L_c -
M_gln_L_c_model = SpeciesModel();
M_gln_L_c_model.species_index = 62;
M_gln_L_c_model.species_symbol = string("M_gln_L_c");
M_gln_L_c_model.species_constraint_type = GLPK.FX;
M_gln_L_c_model.species_lower_bound = 0.0;
M_gln_L_c_model.species_upper_bound = 0.0;
M_gln_L_c_model.is_species_measured = false;
M_gln_L_c_model.is_biomass_precursor = false;
M_gln_L_c_model.biomass_precursor_coefficient = 0.0;
M_gln_L_c_model.species_time_constant = 1.0;
M_gln_L_c_model.species_initial_condition = 0.0;
M_gln_L_c_model.is_species_diluted = true;
M_gln_L_c_model.is_species_extracellular = false;
species_model_dictionary["M_gln_L_c"] = M_gln_L_c_model;
M_gln_L_c_model = 0;

# species_symbol: 63 M_gly_L_c -
M_gly_L_c_model = SpeciesModel();
M_gly_L_c_model.species_index = 63;
M_gly_L_c_model.species_symbol = string("M_gly_L_c");
M_gly_L_c_model.species_constraint_type = GLPK.FX;
M_gly_L_c_model.species_lower_bound = 0.0;
M_gly_L_c_model.species_upper_bound = 0.0;
M_gly_L_c_model.is_species_measured = false;
M_gly_L_c_model.is_biomass_precursor = false;
M_gly_L_c_model.biomass_precursor_coefficient = 0.0;
M_gly_L_c_model.species_time_constant = 1.0;
M_gly_L_c_model.species_initial_condition = 0.0;
M_gly_L_c_model.is_species_diluted = true;
M_gly_L_c_model.is_species_extracellular = false;
species_model_dictionary["M_gly_L_c"] = M_gly_L_c_model;
M_gly_L_c_model = 0;

# species_symbol: 64 M_gar_c -
M_gar_c_model = SpeciesModel();
M_gar_c_model.species_index = 64;
M_gar_c_model.species_symbol = string("M_gar_c");
M_gar_c_model.species_constraint_type = GLPK.FX;
M_gar_c_model.species_lower_bound = 0.0;
M_gar_c_model.species_upper_bound = 0.0;
M_gar_c_model.is_species_measured = false;
M_gar_c_model.is_biomass_precursor = false;
M_gar_c_model.biomass_precursor_coefficient = 0.0;
M_gar_c_model.species_time_constant = 1.0;
M_gar_c_model.species_initial_condition = 0.0;
M_gar_c_model.is_species_diluted = true;
M_gar_c_model.is_species_extracellular = false;
species_model_dictionary["M_gar_c"] = M_gar_c_model;
M_gar_c_model = 0;

# species_symbol: 65 M_glu_L_c -
M_glu_L_c_model = SpeciesModel();
M_glu_L_c_model.species_index = 65;
M_glu_L_c_model.species_symbol = string("M_glu_L_c");
M_glu_L_c_model.species_constraint_type = GLPK.FX;
M_glu_L_c_model.species_lower_bound = 0.0;
M_glu_L_c_model.species_upper_bound = 0.0;
M_glu_L_c_model.is_species_measured = false;
M_glu_L_c_model.is_biomass_precursor = false;
M_glu_L_c_model.biomass_precursor_coefficient = 0.0;
M_glu_L_c_model.species_time_constant = 1.0;
M_glu_L_c_model.species_initial_condition = 0.0;
M_glu_L_c_model.is_species_diluted = true;
M_glu_L_c_model.is_species_extracellular = false;
species_model_dictionary["M_glu_L_c"] = M_glu_L_c_model;
M_glu_L_c_model = 0;

# species_symbol: 66 M_10fthf_c -
M_10fthf_c_model = SpeciesModel();
M_10fthf_c_model.species_index = 66;
M_10fthf_c_model.species_symbol = string("M_10fthf_c");
M_10fthf_c_model.species_constraint_type = GLPK.FX;
M_10fthf_c_model.species_lower_bound = 0.0;
M_10fthf_c_model.species_upper_bound = 0.0;
M_10fthf_c_model.is_species_measured = false;
M_10fthf_c_model.is_biomass_precursor = false;
M_10fthf_c_model.biomass_precursor_coefficient = 0.0;
M_10fthf_c_model.species_time_constant = 1.0;
M_10fthf_c_model.species_initial_condition = 0.0;
M_10fthf_c_model.is_species_diluted = true;
M_10fthf_c_model.is_species_extracellular = false;
species_model_dictionary["M_10fthf_c"] = M_10fthf_c_model;
M_10fthf_c_model = 0;

# species_symbol: 67 M_air_c -
M_air_c_model = SpeciesModel();
M_air_c_model.species_index = 67;
M_air_c_model.species_symbol = string("M_air_c");
M_air_c_model.species_constraint_type = GLPK.FX;
M_air_c_model.species_lower_bound = 0.0;
M_air_c_model.species_upper_bound = 0.0;
M_air_c_model.is_species_measured = false;
M_air_c_model.is_biomass_precursor = false;
M_air_c_model.biomass_precursor_coefficient = 0.0;
M_air_c_model.species_time_constant = 1.0;
M_air_c_model.species_initial_condition = 0.0;
M_air_c_model.is_species_diluted = true;
M_air_c_model.is_species_extracellular = false;
species_model_dictionary["M_air_c"] = M_air_c_model;
M_air_c_model = 0;

# species_symbol: 68 M_thf_c -
M_thf_c_model = SpeciesModel();
M_thf_c_model.species_index = 68;
M_thf_c_model.species_symbol = string("M_thf_c");
M_thf_c_model.species_constraint_type = GLPK.FX;
M_thf_c_model.species_lower_bound = 0.0;
M_thf_c_model.species_upper_bound = 0.0;
M_thf_c_model.is_species_measured = false;
M_thf_c_model.is_biomass_precursor = false;
M_thf_c_model.biomass_precursor_coefficient = 0.0;
M_thf_c_model.species_time_constant = 1.0;
M_thf_c_model.species_initial_condition = 0.0;
M_thf_c_model.is_species_diluted = true;
M_thf_c_model.is_species_extracellular = false;
species_model_dictionary["M_thf_c"] = M_thf_c_model;
M_thf_c_model = 0;

# species_symbol: 69 M_asp_L_c -
M_asp_L_c_model = SpeciesModel();
M_asp_L_c_model.species_index = 69;
M_asp_L_c_model.species_symbol = string("M_asp_L_c");
M_asp_L_c_model.species_constraint_type = GLPK.FX;
M_asp_L_c_model.species_lower_bound = 0.0;
M_asp_L_c_model.species_upper_bound = 0.0;
M_asp_L_c_model.is_species_measured = false;
M_asp_L_c_model.is_biomass_precursor = false;
M_asp_L_c_model.biomass_precursor_coefficient = 0.0;
M_asp_L_c_model.species_time_constant = 1.0;
M_asp_L_c_model.species_initial_condition = 0.0;
M_asp_L_c_model.is_species_diluted = true;
M_asp_L_c_model.is_species_extracellular = false;
species_model_dictionary["M_asp_L_c"] = M_asp_L_c_model;
M_asp_L_c_model = 0;

# species_symbol: 70 M_hco3_c -
M_hco3_c_model = SpeciesModel();
M_hco3_c_model.species_index = 70;
M_hco3_c_model.species_symbol = string("M_hco3_c");
M_hco3_c_model.species_constraint_type = GLPK.FX;
M_hco3_c_model.species_lower_bound = 0.0;
M_hco3_c_model.species_upper_bound = 0.0;
M_hco3_c_model.is_species_measured = false;
M_hco3_c_model.is_biomass_precursor = false;
M_hco3_c_model.biomass_precursor_coefficient = 0.0;
M_hco3_c_model.species_time_constant = 1.0;
M_hco3_c_model.species_initial_condition = 0.0;
M_hco3_c_model.is_species_diluted = true;
M_hco3_c_model.is_species_extracellular = false;
species_model_dictionary["M_hco3_c"] = M_hco3_c_model;
M_hco3_c_model = 0;

# species_symbol: 71 M_aicar_c -
M_aicar_c_model = SpeciesModel();
M_aicar_c_model.species_index = 71;
M_aicar_c_model.species_symbol = string("M_aicar_c");
M_aicar_c_model.species_constraint_type = GLPK.FX;
M_aicar_c_model.species_lower_bound = 0.0;
M_aicar_c_model.species_upper_bound = 0.0;
M_aicar_c_model.is_species_measured = false;
M_aicar_c_model.is_biomass_precursor = false;
M_aicar_c_model.biomass_precursor_coefficient = 0.0;
M_aicar_c_model.species_time_constant = 1.0;
M_aicar_c_model.species_initial_condition = 0.0;
M_aicar_c_model.is_species_diluted = true;
M_aicar_c_model.is_species_extracellular = false;
species_model_dictionary["M_aicar_c"] = M_aicar_c_model;
M_aicar_c_model = 0;

# species_symbol: 72 M_imp_c -
M_imp_c_model = SpeciesModel();
M_imp_c_model.species_index = 72;
M_imp_c_model.species_symbol = string("M_imp_c");
M_imp_c_model.species_constraint_type = GLPK.FX;
M_imp_c_model.species_lower_bound = 0.0;
M_imp_c_model.species_upper_bound = 0.0;
M_imp_c_model.is_species_measured = false;
M_imp_c_model.is_biomass_precursor = false;
M_imp_c_model.biomass_precursor_coefficient = 0.0;
M_imp_c_model.species_time_constant = 1.0;
M_imp_c_model.species_initial_condition = 0.0;
M_imp_c_model.is_species_diluted = true;
M_imp_c_model.is_species_extracellular = false;
species_model_dictionary["M_imp_c"] = M_imp_c_model;
M_imp_c_model = 0;

# species_symbol: 73 M_methf_c -
M_methf_c_model = SpeciesModel();
M_methf_c_model.species_index = 73;
M_methf_c_model.species_symbol = string("M_methf_c");
M_methf_c_model.species_constraint_type = GLPK.FX;
M_methf_c_model.species_lower_bound = 0.0;
M_methf_c_model.species_upper_bound = 0.0;
M_methf_c_model.is_species_measured = false;
M_methf_c_model.is_biomass_precursor = false;
M_methf_c_model.biomass_precursor_coefficient = 0.0;
M_methf_c_model.species_time_constant = 1.0;
M_methf_c_model.species_initial_condition = 0.0;
M_methf_c_model.is_species_diluted = true;
M_methf_c_model.is_species_extracellular = false;
species_model_dictionary["M_methf_c"] = M_methf_c_model;
M_methf_c_model = 0;

# species_symbol: 74 M_mlthf_c -
M_mlthf_c_model = SpeciesModel();
M_mlthf_c_model.species_index = 74;
M_mlthf_c_model.species_symbol = string("M_mlthf_c");
M_mlthf_c_model.species_constraint_type = GLPK.FX;
M_mlthf_c_model.species_lower_bound = 0.0;
M_mlthf_c_model.species_upper_bound = 0.0;
M_mlthf_c_model.is_species_measured = false;
M_mlthf_c_model.is_biomass_precursor = false;
M_mlthf_c_model.biomass_precursor_coefficient = 0.0;
M_mlthf_c_model.species_time_constant = 1.0;
M_mlthf_c_model.species_initial_condition = 0.0;
M_mlthf_c_model.is_species_diluted = true;
M_mlthf_c_model.is_species_extracellular = false;
species_model_dictionary["M_mlthf_c"] = M_mlthf_c_model;
M_mlthf_c_model = 0;

# species_symbol: 75 M_5mthf_c -
M_5mthf_c_model = SpeciesModel();
M_5mthf_c_model.species_index = 75;
M_5mthf_c_model.species_symbol = string("M_5mthf_c");
M_5mthf_c_model.species_constraint_type = GLPK.FX;
M_5mthf_c_model.species_lower_bound = 0.0;
M_5mthf_c_model.species_upper_bound = 0.0;
M_5mthf_c_model.is_species_measured = false;
M_5mthf_c_model.is_biomass_precursor = false;
M_5mthf_c_model.biomass_precursor_coefficient = 0.0;
M_5mthf_c_model.species_time_constant = 1.0;
M_5mthf_c_model.species_initial_condition = 0.0;
M_5mthf_c_model.is_species_diluted = true;
M_5mthf_c_model.is_species_extracellular = false;
species_model_dictionary["M_5mthf_c"] = M_5mthf_c_model;
M_5mthf_c_model = 0;

# species_symbol: 76 M_gmp_c -
M_gmp_c_model = SpeciesModel();
M_gmp_c_model.species_index = 76;
M_gmp_c_model.species_symbol = string("M_gmp_c");
M_gmp_c_model.species_constraint_type = GLPK.FX;
M_gmp_c_model.species_lower_bound = 0.0;
M_gmp_c_model.species_upper_bound = 0.0;
M_gmp_c_model.is_species_measured = false;
M_gmp_c_model.is_biomass_precursor = false;
M_gmp_c_model.biomass_precursor_coefficient = 0.0;
M_gmp_c_model.species_time_constant = 1.0;
M_gmp_c_model.species_initial_condition = 0.0;
M_gmp_c_model.is_species_diluted = true;
M_gmp_c_model.is_species_extracellular = false;
species_model_dictionary["M_gmp_c"] = M_gmp_c_model;
M_gmp_c_model = 0;

# species_symbol: 77 M_ump_c -
M_ump_c_model = SpeciesModel();
M_ump_c_model.species_index = 77;
M_ump_c_model.species_symbol = string("M_ump_c");
M_ump_c_model.species_constraint_type = GLPK.FX;
M_ump_c_model.species_lower_bound = 0.0;
M_ump_c_model.species_upper_bound = 0.0;
M_ump_c_model.is_species_measured = false;
M_ump_c_model.is_biomass_precursor = false;
M_ump_c_model.biomass_precursor_coefficient = 0.0;
M_ump_c_model.species_time_constant = 1.0;
M_ump_c_model.species_initial_condition = 0.0;
M_ump_c_model.is_species_diluted = true;
M_ump_c_model.is_species_extracellular = false;
species_model_dictionary["M_ump_c"] = M_ump_c_model;
M_ump_c_model = 0;

# species_symbol: 78 M_cmp_c -
M_cmp_c_model = SpeciesModel();
M_cmp_c_model.species_index = 78;
M_cmp_c_model.species_symbol = string("M_cmp_c");
M_cmp_c_model.species_constraint_type = GLPK.FX;
M_cmp_c_model.species_lower_bound = 0.0;
M_cmp_c_model.species_upper_bound = 0.0;
M_cmp_c_model.is_species_measured = false;
M_cmp_c_model.is_biomass_precursor = false;
M_cmp_c_model.biomass_precursor_coefficient = 0.0;
M_cmp_c_model.species_time_constant = 1.0;
M_cmp_c_model.species_initial_condition = 0.0;
M_cmp_c_model.is_species_diluted = true;
M_cmp_c_model.is_species_extracellular = false;
species_model_dictionary["M_cmp_c"] = M_cmp_c_model;
M_cmp_c_model = 0;

# species_symbol: 79 M_ala_L_c -
M_ala_L_c_model = SpeciesModel();
M_ala_L_c_model.species_index = 79;
M_ala_L_c_model.species_symbol = string("M_ala_L_c");
M_ala_L_c_model.species_constraint_type = GLPK.FX;
M_ala_L_c_model.species_lower_bound = 0.0;
M_ala_L_c_model.species_upper_bound = 0.0;
M_ala_L_c_model.is_species_measured = false;
M_ala_L_c_model.is_biomass_precursor = false;
M_ala_L_c_model.biomass_precursor_coefficient = 0.0;
M_ala_L_c_model.species_time_constant = 1.0;
M_ala_L_c_model.species_initial_condition = 0.0;
M_ala_L_c_model.is_species_diluted = true;
M_ala_L_c_model.is_species_extracellular = false;
species_model_dictionary["M_ala_L_c"] = M_ala_L_c_model;
M_ala_L_c_model = 0;

# species_symbol: 80 M_arg_L_c -
M_arg_L_c_model = SpeciesModel();
M_arg_L_c_model.species_index = 80;
M_arg_L_c_model.species_symbol = string("M_arg_L_c");
M_arg_L_c_model.species_constraint_type = GLPK.FX;
M_arg_L_c_model.species_lower_bound = 0.0;
M_arg_L_c_model.species_upper_bound = 0.0;
M_arg_L_c_model.is_species_measured = false;
M_arg_L_c_model.is_biomass_precursor = false;
M_arg_L_c_model.biomass_precursor_coefficient = 0.0;
M_arg_L_c_model.species_time_constant = 1.0;
M_arg_L_c_model.species_initial_condition = 0.0;
M_arg_L_c_model.is_species_diluted = true;
M_arg_L_c_model.is_species_extracellular = false;
species_model_dictionary["M_arg_L_c"] = M_arg_L_c_model;
M_arg_L_c_model = 0;

# species_symbol: 81 M_nh4_c -
M_nh4_c_model = SpeciesModel();
M_nh4_c_model.species_index = 81;
M_nh4_c_model.species_symbol = string("M_nh4_c");
M_nh4_c_model.species_constraint_type = GLPK.FX;
M_nh4_c_model.species_lower_bound = 0.0;
M_nh4_c_model.species_upper_bound = 0.0;
M_nh4_c_model.is_species_measured = false;
M_nh4_c_model.is_biomass_precursor = false;
M_nh4_c_model.biomass_precursor_coefficient = 0.0;
M_nh4_c_model.species_time_constant = 1.0;
M_nh4_c_model.species_initial_condition = 0.0;
M_nh4_c_model.is_species_diluted = true;
M_nh4_c_model.is_species_extracellular = false;
species_model_dictionary["M_nh4_c"] = M_nh4_c_model;
M_nh4_c_model = 0;

# species_symbol: 82 M_asn_L_c -
M_asn_L_c_model = SpeciesModel();
M_asn_L_c_model.species_index = 82;
M_asn_L_c_model.species_symbol = string("M_asn_L_c");
M_asn_L_c_model.species_constraint_type = GLPK.FX;
M_asn_L_c_model.species_lower_bound = 0.0;
M_asn_L_c_model.species_upper_bound = 0.0;
M_asn_L_c_model.is_species_measured = false;
M_asn_L_c_model.is_biomass_precursor = false;
M_asn_L_c_model.biomass_precursor_coefficient = 0.0;
M_asn_L_c_model.species_time_constant = 1.0;
M_asn_L_c_model.species_initial_condition = 0.0;
M_asn_L_c_model.is_species_diluted = true;
M_asn_L_c_model.is_species_extracellular = false;
species_model_dictionary["M_asn_L_c"] = M_asn_L_c_model;
M_asn_L_c_model = 0;

# species_symbol: 83 M_ser_L_c -
M_ser_L_c_model = SpeciesModel();
M_ser_L_c_model.species_index = 83;
M_ser_L_c_model.species_symbol = string("M_ser_L_c");
M_ser_L_c_model.species_constraint_type = GLPK.FX;
M_ser_L_c_model.species_lower_bound = 0.0;
M_ser_L_c_model.species_upper_bound = 0.0;
M_ser_L_c_model.is_species_measured = false;
M_ser_L_c_model.is_biomass_precursor = false;
M_ser_L_c_model.biomass_precursor_coefficient = 0.0;
M_ser_L_c_model.species_time_constant = 1.0;
M_ser_L_c_model.species_initial_condition = 0.0;
M_ser_L_c_model.is_species_diluted = true;
M_ser_L_c_model.is_species_extracellular = false;
species_model_dictionary["M_ser_L_c"] = M_ser_L_c_model;
M_ser_L_c_model = 0;

# species_symbol: 84 M_h2s_c -
M_h2s_c_model = SpeciesModel();
M_h2s_c_model.species_index = 84;
M_h2s_c_model.species_symbol = string("M_h2s_c");
M_h2s_c_model.species_constraint_type = GLPK.FX;
M_h2s_c_model.species_lower_bound = 0.0;
M_h2s_c_model.species_upper_bound = 0.0;
M_h2s_c_model.is_species_measured = false;
M_h2s_c_model.is_biomass_precursor = false;
M_h2s_c_model.biomass_precursor_coefficient = 0.0;
M_h2s_c_model.species_time_constant = 1.0;
M_h2s_c_model.species_initial_condition = 0.0;
M_h2s_c_model.is_species_diluted = true;
M_h2s_c_model.is_species_extracellular = false;
species_model_dictionary["M_h2s_c"] = M_h2s_c_model;
M_h2s_c_model = 0;

# species_symbol: 85 M_cys_L_c -
M_cys_L_c_model = SpeciesModel();
M_cys_L_c_model.species_index = 85;
M_cys_L_c_model.species_symbol = string("M_cys_L_c");
M_cys_L_c_model.species_constraint_type = GLPK.FX;
M_cys_L_c_model.species_lower_bound = 0.0;
M_cys_L_c_model.species_upper_bound = 0.0;
M_cys_L_c_model.is_species_measured = false;
M_cys_L_c_model.is_biomass_precursor = false;
M_cys_L_c_model.biomass_precursor_coefficient = 0.0;
M_cys_L_c_model.species_time_constant = 1.0;
M_cys_L_c_model.species_initial_condition = 0.0;
M_cys_L_c_model.is_species_diluted = true;
M_cys_L_c_model.is_species_extracellular = false;
species_model_dictionary["M_cys_L_c"] = M_cys_L_c_model;
M_cys_L_c_model = 0;

# species_symbol: 86 M_his_L_c -
M_his_L_c_model = SpeciesModel();
M_his_L_c_model.species_index = 86;
M_his_L_c_model.species_symbol = string("M_his_L_c");
M_his_L_c_model.species_constraint_type = GLPK.FX;
M_his_L_c_model.species_lower_bound = 0.0;
M_his_L_c_model.species_upper_bound = 0.0;
M_his_L_c_model.is_species_measured = false;
M_his_L_c_model.is_biomass_precursor = false;
M_his_L_c_model.biomass_precursor_coefficient = 0.0;
M_his_L_c_model.species_time_constant = 1.0;
M_his_L_c_model.species_initial_condition = 0.0;
M_his_L_c_model.is_species_diluted = true;
M_his_L_c_model.is_species_extracellular = false;
species_model_dictionary["M_his_L_c"] = M_his_L_c_model;
M_his_L_c_model = 0;

# species_symbol: 87 M_thr_L_c -
M_thr_L_c_model = SpeciesModel();
M_thr_L_c_model.species_index = 87;
M_thr_L_c_model.species_symbol = string("M_thr_L_c");
M_thr_L_c_model.species_constraint_type = GLPK.FX;
M_thr_L_c_model.species_lower_bound = 0.0;
M_thr_L_c_model.species_upper_bound = 0.0;
M_thr_L_c_model.is_species_measured = false;
M_thr_L_c_model.is_biomass_precursor = false;
M_thr_L_c_model.biomass_precursor_coefficient = 0.0;
M_thr_L_c_model.species_time_constant = 1.0;
M_thr_L_c_model.species_initial_condition = 0.0;
M_thr_L_c_model.is_species_diluted = true;
M_thr_L_c_model.is_species_extracellular = false;
species_model_dictionary["M_thr_L_c"] = M_thr_L_c_model;
M_thr_L_c_model = 0;

# species_symbol: 88 M_ile_L_c -
M_ile_L_c_model = SpeciesModel();
M_ile_L_c_model.species_index = 88;
M_ile_L_c_model.species_symbol = string("M_ile_L_c");
M_ile_L_c_model.species_constraint_type = GLPK.FX;
M_ile_L_c_model.species_lower_bound = 0.0;
M_ile_L_c_model.species_upper_bound = 0.0;
M_ile_L_c_model.is_species_measured = false;
M_ile_L_c_model.is_biomass_precursor = false;
M_ile_L_c_model.biomass_precursor_coefficient = 0.0;
M_ile_L_c_model.species_time_constant = 1.0;
M_ile_L_c_model.species_initial_condition = 0.0;
M_ile_L_c_model.is_species_diluted = true;
M_ile_L_c_model.is_species_extracellular = false;
species_model_dictionary["M_ile_L_c"] = M_ile_L_c_model;
M_ile_L_c_model = 0;

# species_symbol: 89 M_leu_L_c -
M_leu_L_c_model = SpeciesModel();
M_leu_L_c_model.species_index = 89;
M_leu_L_c_model.species_symbol = string("M_leu_L_c");
M_leu_L_c_model.species_constraint_type = GLPK.FX;
M_leu_L_c_model.species_lower_bound = 0.0;
M_leu_L_c_model.species_upper_bound = 0.0;
M_leu_L_c_model.is_species_measured = false;
M_leu_L_c_model.is_biomass_precursor = false;
M_leu_L_c_model.biomass_precursor_coefficient = 0.0;
M_leu_L_c_model.species_time_constant = 1.0;
M_leu_L_c_model.species_initial_condition = 0.0;
M_leu_L_c_model.is_species_diluted = true;
M_leu_L_c_model.is_species_extracellular = false;
species_model_dictionary["M_leu_L_c"] = M_leu_L_c_model;
M_leu_L_c_model = 0;

# species_symbol: 90 M_lys_L_c -
M_lys_L_c_model = SpeciesModel();
M_lys_L_c_model.species_index = 90;
M_lys_L_c_model.species_symbol = string("M_lys_L_c");
M_lys_L_c_model.species_constraint_type = GLPK.FX;
M_lys_L_c_model.species_lower_bound = 0.0;
M_lys_L_c_model.species_upper_bound = 0.0;
M_lys_L_c_model.is_species_measured = false;
M_lys_L_c_model.is_biomass_precursor = false;
M_lys_L_c_model.biomass_precursor_coefficient = 0.0;
M_lys_L_c_model.species_time_constant = 1.0;
M_lys_L_c_model.species_initial_condition = 0.0;
M_lys_L_c_model.is_species_diluted = true;
M_lys_L_c_model.is_species_extracellular = false;
species_model_dictionary["M_lys_L_c"] = M_lys_L_c_model;
M_lys_L_c_model = 0;

# species_symbol: 91 M_met_L_c -
M_met_L_c_model = SpeciesModel();
M_met_L_c_model.species_index = 91;
M_met_L_c_model.species_symbol = string("M_met_L_c");
M_met_L_c_model.species_constraint_type = GLPK.FX;
M_met_L_c_model.species_lower_bound = 0.0;
M_met_L_c_model.species_upper_bound = 0.0;
M_met_L_c_model.is_species_measured = false;
M_met_L_c_model.is_biomass_precursor = false;
M_met_L_c_model.biomass_precursor_coefficient = 0.0;
M_met_L_c_model.species_time_constant = 1.0;
M_met_L_c_model.species_initial_condition = 0.0;
M_met_L_c_model.is_species_diluted = true;
M_met_L_c_model.is_species_extracellular = false;
species_model_dictionary["M_met_L_c"] = M_met_L_c_model;
M_met_L_c_model = 0;

# species_symbol: 92 M_phe_L_c -
M_phe_L_c_model = SpeciesModel();
M_phe_L_c_model.species_index = 92;
M_phe_L_c_model.species_symbol = string("M_phe_L_c");
M_phe_L_c_model.species_constraint_type = GLPK.FX;
M_phe_L_c_model.species_lower_bound = 0.0;
M_phe_L_c_model.species_upper_bound = 0.0;
M_phe_L_c_model.is_species_measured = false;
M_phe_L_c_model.is_biomass_precursor = false;
M_phe_L_c_model.biomass_precursor_coefficient = 0.0;
M_phe_L_c_model.species_time_constant = 1.0;
M_phe_L_c_model.species_initial_condition = 0.0;
M_phe_L_c_model.is_species_diluted = true;
M_phe_L_c_model.is_species_extracellular = false;
species_model_dictionary["M_phe_L_c"] = M_phe_L_c_model;
M_phe_L_c_model = 0;

# species_symbol: 93 M_pro_L_c -
M_pro_L_c_model = SpeciesModel();
M_pro_L_c_model.species_index = 93;
M_pro_L_c_model.species_symbol = string("M_pro_L_c");
M_pro_L_c_model.species_constraint_type = GLPK.FX;
M_pro_L_c_model.species_lower_bound = 0.0;
M_pro_L_c_model.species_upper_bound = 0.0;
M_pro_L_c_model.is_species_measured = false;
M_pro_L_c_model.is_biomass_precursor = false;
M_pro_L_c_model.biomass_precursor_coefficient = 0.0;
M_pro_L_c_model.species_time_constant = 1.0;
M_pro_L_c_model.species_initial_condition = 0.0;
M_pro_L_c_model.is_species_diluted = true;
M_pro_L_c_model.is_species_extracellular = false;
species_model_dictionary["M_pro_L_c"] = M_pro_L_c_model;
M_pro_L_c_model = 0;

# species_symbol: 94 M_trp_L_c -
M_trp_L_c_model = SpeciesModel();
M_trp_L_c_model.species_index = 94;
M_trp_L_c_model.species_symbol = string("M_trp_L_c");
M_trp_L_c_model.species_constraint_type = GLPK.FX;
M_trp_L_c_model.species_lower_bound = 0.0;
M_trp_L_c_model.species_upper_bound = 0.0;
M_trp_L_c_model.is_species_measured = false;
M_trp_L_c_model.is_biomass_precursor = false;
M_trp_L_c_model.biomass_precursor_coefficient = 0.0;
M_trp_L_c_model.species_time_constant = 1.0;
M_trp_L_c_model.species_initial_condition = 0.0;
M_trp_L_c_model.is_species_diluted = true;
M_trp_L_c_model.is_species_extracellular = false;
species_model_dictionary["M_trp_L_c"] = M_trp_L_c_model;
M_trp_L_c_model = 0;

# species_symbol: 95 M_tyr_L_c -
M_tyr_L_c_model = SpeciesModel();
M_tyr_L_c_model.species_index = 95;
M_tyr_L_c_model.species_symbol = string("M_tyr_L_c");
M_tyr_L_c_model.species_constraint_type = GLPK.FX;
M_tyr_L_c_model.species_lower_bound = 0.0;
M_tyr_L_c_model.species_upper_bound = 0.0;
M_tyr_L_c_model.is_species_measured = false;
M_tyr_L_c_model.is_biomass_precursor = false;
M_tyr_L_c_model.biomass_precursor_coefficient = 0.0;
M_tyr_L_c_model.species_time_constant = 1.0;
M_tyr_L_c_model.species_initial_condition = 0.0;
M_tyr_L_c_model.is_species_diluted = true;
M_tyr_L_c_model.is_species_extracellular = false;
species_model_dictionary["M_tyr_L_c"] = M_tyr_L_c_model;
M_tyr_L_c_model = 0;

# species_symbol: 96 M_val_L_c -
M_val_L_c_model = SpeciesModel();
M_val_L_c_model.species_index = 96;
M_val_L_c_model.species_symbol = string("M_val_L_c");
M_val_L_c_model.species_constraint_type = GLPK.FX;
M_val_L_c_model.species_lower_bound = 0.0;
M_val_L_c_model.species_upper_bound = 0.0;
M_val_L_c_model.is_species_measured = false;
M_val_L_c_model.is_biomass_precursor = false;
M_val_L_c_model.biomass_precursor_coefficient = 0.0;
M_val_L_c_model.species_time_constant = 1.0;
M_val_L_c_model.species_initial_condition = 0.0;
M_val_L_c_model.is_species_diluted = true;
M_val_L_c_model.is_species_extracellular = false;
species_model_dictionary["M_val_L_c"] = M_val_L_c_model;
M_val_L_c_model = 0;

# species_symbol: 97 M_urea_c -
M_urea_c_model = SpeciesModel();
M_urea_c_model.species_index = 97;
M_urea_c_model.species_symbol = string("M_urea_c");
M_urea_c_model.species_constraint_type = GLPK.FX;
M_urea_c_model.species_lower_bound = 0.0;
M_urea_c_model.species_upper_bound = 0.0;
M_urea_c_model.is_species_measured = false;
M_urea_c_model.is_biomass_precursor = false;
M_urea_c_model.biomass_precursor_coefficient = 0.0;
M_urea_c_model.species_time_constant = 1.0;
M_urea_c_model.species_initial_condition = 0.0;
M_urea_c_model.is_species_diluted = true;
M_urea_c_model.is_species_extracellular = false;
species_model_dictionary["M_urea_c"] = M_urea_c_model;
M_urea_c_model = 0;

# species_symbol: 98 M_h2o2_c -
M_h2o2_c_model = SpeciesModel();
M_h2o2_c_model.species_index = 98;
M_h2o2_c_model.species_symbol = string("M_h2o2_c");
M_h2o2_c_model.species_constraint_type = GLPK.FX;
M_h2o2_c_model.species_lower_bound = 0.0;
M_h2o2_c_model.species_upper_bound = 0.0;
M_h2o2_c_model.is_species_measured = false;
M_h2o2_c_model.is_biomass_precursor = false;
M_h2o2_c_model.biomass_precursor_coefficient = 0.0;
M_h2o2_c_model.species_time_constant = 1.0;
M_h2o2_c_model.species_initial_condition = 0.0;
M_h2o2_c_model.is_species_diluted = true;
M_h2o2_c_model.is_species_extracellular = false;
species_model_dictionary["M_h2o2_c"] = M_h2o2_c_model;
M_h2o2_c_model = 0;

# species_symbol: 99 M_mglx_c -
M_mglx_c_model = SpeciesModel();
M_mglx_c_model.species_index = 99;
M_mglx_c_model.species_symbol = string("M_mglx_c");
M_mglx_c_model.species_constraint_type = GLPK.FX;
M_mglx_c_model.species_lower_bound = 0.0;
M_mglx_c_model.species_upper_bound = 0.0;
M_mglx_c_model.is_species_measured = false;
M_mglx_c_model.is_biomass_precursor = false;
M_mglx_c_model.biomass_precursor_coefficient = 0.0;
M_mglx_c_model.species_time_constant = 1.0;
M_mglx_c_model.species_initial_condition = 0.0;
M_mglx_c_model.is_species_diluted = true;
M_mglx_c_model.is_species_extracellular = false;
species_model_dictionary["M_mglx_c"] = M_mglx_c_model;
M_mglx_c_model = 0;

# species_symbol: 100 M_prop_c -
M_prop_c_model = SpeciesModel();
M_prop_c_model.species_index = 100;
M_prop_c_model.species_symbol = string("M_prop_c");
M_prop_c_model.species_constraint_type = GLPK.FX;
M_prop_c_model.species_lower_bound = 0.0;
M_prop_c_model.species_upper_bound = 0.0;
M_prop_c_model.is_species_measured = false;
M_prop_c_model.is_biomass_precursor = false;
M_prop_c_model.biomass_precursor_coefficient = 0.0;
M_prop_c_model.species_time_constant = 1.0;
M_prop_c_model.species_initial_condition = 0.0;
M_prop_c_model.is_species_diluted = true;
M_prop_c_model.is_species_extracellular = false;
species_model_dictionary["M_prop_c"] = M_prop_c_model;
M_prop_c_model = 0;

# species_symbol: 101 M_indole_c -
M_indole_c_model = SpeciesModel();
M_indole_c_model.species_index = 101;
M_indole_c_model.species_symbol = string("M_indole_c");
M_indole_c_model.species_constraint_type = GLPK.FX;
M_indole_c_model.species_lower_bound = 0.0;
M_indole_c_model.species_upper_bound = 0.0;
M_indole_c_model.is_species_measured = false;
M_indole_c_model.is_biomass_precursor = false;
M_indole_c_model.biomass_precursor_coefficient = 0.0;
M_indole_c_model.species_time_constant = 1.0;
M_indole_c_model.species_initial_condition = 0.0;
M_indole_c_model.is_species_diluted = true;
M_indole_c_model.is_species_extracellular = false;
species_model_dictionary["M_indole_c"] = M_indole_c_model;
M_indole_c_model = 0;

# species_symbol: 102 M_cadav_c -
M_cadav_c_model = SpeciesModel();
M_cadav_c_model.species_index = 102;
M_cadav_c_model.species_symbol = string("M_cadav_c");
M_cadav_c_model.species_constraint_type = GLPK.FX;
M_cadav_c_model.species_lower_bound = 0.0;
M_cadav_c_model.species_upper_bound = 0.0;
M_cadav_c_model.is_species_measured = false;
M_cadav_c_model.is_biomass_precursor = false;
M_cadav_c_model.biomass_precursor_coefficient = 0.0;
M_cadav_c_model.species_time_constant = 1.0;
M_cadav_c_model.species_initial_condition = 0.0;
M_cadav_c_model.is_species_diluted = true;
M_cadav_c_model.is_species_extracellular = false;
species_model_dictionary["M_cadav_c"] = M_cadav_c_model;
M_cadav_c_model = 0;

# species_symbol: 103 M_gaba_c -
M_gaba_c_model = SpeciesModel();
M_gaba_c_model.species_index = 103;
M_gaba_c_model.species_symbol = string("M_gaba_c");
M_gaba_c_model.species_constraint_type = GLPK.FX;
M_gaba_c_model.species_lower_bound = 0.0;
M_gaba_c_model.species_upper_bound = 0.0;
M_gaba_c_model.is_species_measured = false;
M_gaba_c_model.is_biomass_precursor = false;
M_gaba_c_model.biomass_precursor_coefficient = 0.0;
M_gaba_c_model.species_time_constant = 1.0;
M_gaba_c_model.species_initial_condition = 0.0;
M_gaba_c_model.is_species_diluted = true;
M_gaba_c_model.is_species_extracellular = false;
species_model_dictionary["M_gaba_c"] = M_gaba_c_model;
M_gaba_c_model = 0;

# species_symbol: 104 GENE_deGFP -
GENE_deGFP_model = SpeciesModel();
GENE_deGFP_model.species_index = 104;
GENE_deGFP_model.species_symbol = string("GENE_deGFP");
GENE_deGFP_model.species_constraint_type = GLPK.FX;
GENE_deGFP_model.species_lower_bound = 0.0;
GENE_deGFP_model.species_upper_bound = 0.0;
GENE_deGFP_model.is_species_measured = false;
GENE_deGFP_model.is_biomass_precursor = false;
GENE_deGFP_model.biomass_precursor_coefficient = 0.0;
GENE_deGFP_model.species_time_constant = 1.0;
GENE_deGFP_model.species_initial_condition = 0.0;
GENE_deGFP_model.is_species_diluted = true;
GENE_deGFP_model.is_species_extracellular = false;
species_model_dictionary["GENE_deGFP"] = GENE_deGFP_model;
GENE_deGFP_model = 0;

# species_symbol: 105 RNAP -
RNAP_model = SpeciesModel();
RNAP_model.species_index = 105;
RNAP_model.species_symbol = string("RNAP");
RNAP_model.species_constraint_type = GLPK.FX;
RNAP_model.species_lower_bound = 0.0;
RNAP_model.species_upper_bound = 0.0;
RNAP_model.is_species_measured = false;
RNAP_model.is_biomass_precursor = false;
RNAP_model.biomass_precursor_coefficient = 0.0;
RNAP_model.species_time_constant = 1.0;
RNAP_model.species_initial_condition = 0.0;
RNAP_model.is_species_diluted = true;
RNAP_model.is_species_extracellular = false;
species_model_dictionary["RNAP"] = RNAP_model;
RNAP_model = 0;

# species_symbol: 106 OPEN_GENE_deGFP -
OPEN_GENE_deGFP_model = SpeciesModel();
OPEN_GENE_deGFP_model.species_index = 106;
OPEN_GENE_deGFP_model.species_symbol = string("OPEN_GENE_deGFP");
OPEN_GENE_deGFP_model.species_constraint_type = GLPK.FX;
OPEN_GENE_deGFP_model.species_lower_bound = 0.0;
OPEN_GENE_deGFP_model.species_upper_bound = 0.0;
OPEN_GENE_deGFP_model.is_species_measured = false;
OPEN_GENE_deGFP_model.is_biomass_precursor = false;
OPEN_GENE_deGFP_model.biomass_precursor_coefficient = 0.0;
OPEN_GENE_deGFP_model.species_time_constant = 1.0;
OPEN_GENE_deGFP_model.species_initial_condition = 0.0;
OPEN_GENE_deGFP_model.is_species_diluted = true;
OPEN_GENE_deGFP_model.is_species_extracellular = false;
species_model_dictionary["OPEN_GENE_deGFP"] = OPEN_GENE_deGFP_model;
OPEN_GENE_deGFP_model = 0;

# species_symbol: 107 mRNA_deGFP -
mRNA_deGFP_model = SpeciesModel();
mRNA_deGFP_model.species_index = 107;
mRNA_deGFP_model.species_symbol = string("mRNA_deGFP");
mRNA_deGFP_model.species_constraint_type = GLPK.FX;
mRNA_deGFP_model.species_lower_bound = 0.0;
mRNA_deGFP_model.species_upper_bound = 0.0;
mRNA_deGFP_model.is_species_measured = false;
mRNA_deGFP_model.is_biomass_precursor = false;
mRNA_deGFP_model.biomass_precursor_coefficient = 0.0;
mRNA_deGFP_model.species_time_constant = 1.0;
mRNA_deGFP_model.species_initial_condition = 0.0;
mRNA_deGFP_model.is_species_diluted = true;
mRNA_deGFP_model.is_species_extracellular = false;
species_model_dictionary["mRNA_deGFP"] = mRNA_deGFP_model;
mRNA_deGFP_model = 0;

# species_symbol: 108 RIBOSOME -
RIBOSOME_model = SpeciesModel();
RIBOSOME_model.species_index = 108;
RIBOSOME_model.species_symbol = string("RIBOSOME");
RIBOSOME_model.species_constraint_type = GLPK.FX;
RIBOSOME_model.species_lower_bound = 0.0;
RIBOSOME_model.species_upper_bound = 0.0;
RIBOSOME_model.is_species_measured = false;
RIBOSOME_model.is_biomass_precursor = false;
RIBOSOME_model.biomass_precursor_coefficient = 0.0;
RIBOSOME_model.species_time_constant = 1.0;
RIBOSOME_model.species_initial_condition = 0.0;
RIBOSOME_model.is_species_diluted = true;
RIBOSOME_model.is_species_extracellular = false;
species_model_dictionary["RIBOSOME"] = RIBOSOME_model;
RIBOSOME_model = 0;

# species_symbol: 109 RIBOSOME_START_deGFP -
RIBOSOME_START_deGFP_model = SpeciesModel();
RIBOSOME_START_deGFP_model.species_index = 109;
RIBOSOME_START_deGFP_model.species_symbol = string("RIBOSOME_START_deGFP");
RIBOSOME_START_deGFP_model.species_constraint_type = GLPK.FX;
RIBOSOME_START_deGFP_model.species_lower_bound = 0.0;
RIBOSOME_START_deGFP_model.species_upper_bound = 0.0;
RIBOSOME_START_deGFP_model.is_species_measured = false;
RIBOSOME_START_deGFP_model.is_biomass_precursor = false;
RIBOSOME_START_deGFP_model.biomass_precursor_coefficient = 0.0;
RIBOSOME_START_deGFP_model.species_time_constant = 1.0;
RIBOSOME_START_deGFP_model.species_initial_condition = 0.0;
RIBOSOME_START_deGFP_model.is_species_diluted = true;
RIBOSOME_START_deGFP_model.is_species_extracellular = false;
species_model_dictionary["RIBOSOME_START_deGFP"] = RIBOSOME_START_deGFP_model;
RIBOSOME_START_deGFP_model = 0;

# species_symbol: 110 M_ala_L_c_tRNA -
M_ala_L_c_tRNA_model = SpeciesModel();
M_ala_L_c_tRNA_model.species_index = 110;
M_ala_L_c_tRNA_model.species_symbol = string("M_ala_L_c_tRNA");
M_ala_L_c_tRNA_model.species_constraint_type = GLPK.FX;
M_ala_L_c_tRNA_model.species_lower_bound = 0.0;
M_ala_L_c_tRNA_model.species_upper_bound = 0.0;
M_ala_L_c_tRNA_model.is_species_measured = false;
M_ala_L_c_tRNA_model.is_biomass_precursor = false;
M_ala_L_c_tRNA_model.biomass_precursor_coefficient = 0.0;
M_ala_L_c_tRNA_model.species_time_constant = 1.0;
M_ala_L_c_tRNA_model.species_initial_condition = 0.0;
M_ala_L_c_tRNA_model.is_species_diluted = true;
M_ala_L_c_tRNA_model.is_species_extracellular = false;
species_model_dictionary["M_ala_L_c_tRNA"] = M_ala_L_c_tRNA_model;
M_ala_L_c_tRNA_model = 0;

# species_symbol: 111 M_arg_L_c_tRNA -
M_arg_L_c_tRNA_model = SpeciesModel();
M_arg_L_c_tRNA_model.species_index = 111;
M_arg_L_c_tRNA_model.species_symbol = string("M_arg_L_c_tRNA");
M_arg_L_c_tRNA_model.species_constraint_type = GLPK.FX;
M_arg_L_c_tRNA_model.species_lower_bound = 0.0;
M_arg_L_c_tRNA_model.species_upper_bound = 0.0;
M_arg_L_c_tRNA_model.is_species_measured = false;
M_arg_L_c_tRNA_model.is_biomass_precursor = false;
M_arg_L_c_tRNA_model.biomass_precursor_coefficient = 0.0;
M_arg_L_c_tRNA_model.species_time_constant = 1.0;
M_arg_L_c_tRNA_model.species_initial_condition = 0.0;
M_arg_L_c_tRNA_model.is_species_diluted = true;
M_arg_L_c_tRNA_model.is_species_extracellular = false;
species_model_dictionary["M_arg_L_c_tRNA"] = M_arg_L_c_tRNA_model;
M_arg_L_c_tRNA_model = 0;

# species_symbol: 112 M_asn_L_c_tRNA -
M_asn_L_c_tRNA_model = SpeciesModel();
M_asn_L_c_tRNA_model.species_index = 112;
M_asn_L_c_tRNA_model.species_symbol = string("M_asn_L_c_tRNA");
M_asn_L_c_tRNA_model.species_constraint_type = GLPK.FX;
M_asn_L_c_tRNA_model.species_lower_bound = 0.0;
M_asn_L_c_tRNA_model.species_upper_bound = 0.0;
M_asn_L_c_tRNA_model.is_species_measured = false;
M_asn_L_c_tRNA_model.is_biomass_precursor = false;
M_asn_L_c_tRNA_model.biomass_precursor_coefficient = 0.0;
M_asn_L_c_tRNA_model.species_time_constant = 1.0;
M_asn_L_c_tRNA_model.species_initial_condition = 0.0;
M_asn_L_c_tRNA_model.is_species_diluted = true;
M_asn_L_c_tRNA_model.is_species_extracellular = false;
species_model_dictionary["M_asn_L_c_tRNA"] = M_asn_L_c_tRNA_model;
M_asn_L_c_tRNA_model = 0;

# species_symbol: 113 M_asp_L_c_tRNA -
M_asp_L_c_tRNA_model = SpeciesModel();
M_asp_L_c_tRNA_model.species_index = 113;
M_asp_L_c_tRNA_model.species_symbol = string("M_asp_L_c_tRNA");
M_asp_L_c_tRNA_model.species_constraint_type = GLPK.FX;
M_asp_L_c_tRNA_model.species_lower_bound = 0.0;
M_asp_L_c_tRNA_model.species_upper_bound = 0.0;
M_asp_L_c_tRNA_model.is_species_measured = false;
M_asp_L_c_tRNA_model.is_biomass_precursor = false;
M_asp_L_c_tRNA_model.biomass_precursor_coefficient = 0.0;
M_asp_L_c_tRNA_model.species_time_constant = 1.0;
M_asp_L_c_tRNA_model.species_initial_condition = 0.0;
M_asp_L_c_tRNA_model.is_species_diluted = true;
M_asp_L_c_tRNA_model.is_species_extracellular = false;
species_model_dictionary["M_asp_L_c_tRNA"] = M_asp_L_c_tRNA_model;
M_asp_L_c_tRNA_model = 0;

# species_symbol: 114 M_cys_L_c_tRNA -
M_cys_L_c_tRNA_model = SpeciesModel();
M_cys_L_c_tRNA_model.species_index = 114;
M_cys_L_c_tRNA_model.species_symbol = string("M_cys_L_c_tRNA");
M_cys_L_c_tRNA_model.species_constraint_type = GLPK.FX;
M_cys_L_c_tRNA_model.species_lower_bound = 0.0;
M_cys_L_c_tRNA_model.species_upper_bound = 0.0;
M_cys_L_c_tRNA_model.is_species_measured = false;
M_cys_L_c_tRNA_model.is_biomass_precursor = false;
M_cys_L_c_tRNA_model.biomass_precursor_coefficient = 0.0;
M_cys_L_c_tRNA_model.species_time_constant = 1.0;
M_cys_L_c_tRNA_model.species_initial_condition = 0.0;
M_cys_L_c_tRNA_model.is_species_diluted = true;
M_cys_L_c_tRNA_model.is_species_extracellular = false;
species_model_dictionary["M_cys_L_c_tRNA"] = M_cys_L_c_tRNA_model;
M_cys_L_c_tRNA_model = 0;

# species_symbol: 115 M_glu_L_c_tRNA -
M_glu_L_c_tRNA_model = SpeciesModel();
M_glu_L_c_tRNA_model.species_index = 115;
M_glu_L_c_tRNA_model.species_symbol = string("M_glu_L_c_tRNA");
M_glu_L_c_tRNA_model.species_constraint_type = GLPK.FX;
M_glu_L_c_tRNA_model.species_lower_bound = 0.0;
M_glu_L_c_tRNA_model.species_upper_bound = 0.0;
M_glu_L_c_tRNA_model.is_species_measured = false;
M_glu_L_c_tRNA_model.is_biomass_precursor = false;
M_glu_L_c_tRNA_model.biomass_precursor_coefficient = 0.0;
M_glu_L_c_tRNA_model.species_time_constant = 1.0;
M_glu_L_c_tRNA_model.species_initial_condition = 0.0;
M_glu_L_c_tRNA_model.is_species_diluted = true;
M_glu_L_c_tRNA_model.is_species_extracellular = false;
species_model_dictionary["M_glu_L_c_tRNA"] = M_glu_L_c_tRNA_model;
M_glu_L_c_tRNA_model = 0;

# species_symbol: 116 M_gln_L_c_tRNA -
M_gln_L_c_tRNA_model = SpeciesModel();
M_gln_L_c_tRNA_model.species_index = 116;
M_gln_L_c_tRNA_model.species_symbol = string("M_gln_L_c_tRNA");
M_gln_L_c_tRNA_model.species_constraint_type = GLPK.FX;
M_gln_L_c_tRNA_model.species_lower_bound = 0.0;
M_gln_L_c_tRNA_model.species_upper_bound = 0.0;
M_gln_L_c_tRNA_model.is_species_measured = false;
M_gln_L_c_tRNA_model.is_biomass_precursor = false;
M_gln_L_c_tRNA_model.biomass_precursor_coefficient = 0.0;
M_gln_L_c_tRNA_model.species_time_constant = 1.0;
M_gln_L_c_tRNA_model.species_initial_condition = 0.0;
M_gln_L_c_tRNA_model.is_species_diluted = true;
M_gln_L_c_tRNA_model.is_species_extracellular = false;
species_model_dictionary["M_gln_L_c_tRNA"] = M_gln_L_c_tRNA_model;
M_gln_L_c_tRNA_model = 0;

# species_symbol: 117 M_gly_L_c_tRNA -
M_gly_L_c_tRNA_model = SpeciesModel();
M_gly_L_c_tRNA_model.species_index = 117;
M_gly_L_c_tRNA_model.species_symbol = string("M_gly_L_c_tRNA");
M_gly_L_c_tRNA_model.species_constraint_type = GLPK.FX;
M_gly_L_c_tRNA_model.species_lower_bound = 0.0;
M_gly_L_c_tRNA_model.species_upper_bound = 0.0;
M_gly_L_c_tRNA_model.is_species_measured = false;
M_gly_L_c_tRNA_model.is_biomass_precursor = false;
M_gly_L_c_tRNA_model.biomass_precursor_coefficient = 0.0;
M_gly_L_c_tRNA_model.species_time_constant = 1.0;
M_gly_L_c_tRNA_model.species_initial_condition = 0.0;
M_gly_L_c_tRNA_model.is_species_diluted = true;
M_gly_L_c_tRNA_model.is_species_extracellular = false;
species_model_dictionary["M_gly_L_c_tRNA"] = M_gly_L_c_tRNA_model;
M_gly_L_c_tRNA_model = 0;

# species_symbol: 118 M_his_L_c_tRNA -
M_his_L_c_tRNA_model = SpeciesModel();
M_his_L_c_tRNA_model.species_index = 118;
M_his_L_c_tRNA_model.species_symbol = string("M_his_L_c_tRNA");
M_his_L_c_tRNA_model.species_constraint_type = GLPK.FX;
M_his_L_c_tRNA_model.species_lower_bound = 0.0;
M_his_L_c_tRNA_model.species_upper_bound = 0.0;
M_his_L_c_tRNA_model.is_species_measured = false;
M_his_L_c_tRNA_model.is_biomass_precursor = false;
M_his_L_c_tRNA_model.biomass_precursor_coefficient = 0.0;
M_his_L_c_tRNA_model.species_time_constant = 1.0;
M_his_L_c_tRNA_model.species_initial_condition = 0.0;
M_his_L_c_tRNA_model.is_species_diluted = true;
M_his_L_c_tRNA_model.is_species_extracellular = false;
species_model_dictionary["M_his_L_c_tRNA"] = M_his_L_c_tRNA_model;
M_his_L_c_tRNA_model = 0;

# species_symbol: 119 M_ile_L_c_tRNA -
M_ile_L_c_tRNA_model = SpeciesModel();
M_ile_L_c_tRNA_model.species_index = 119;
M_ile_L_c_tRNA_model.species_symbol = string("M_ile_L_c_tRNA");
M_ile_L_c_tRNA_model.species_constraint_type = GLPK.FX;
M_ile_L_c_tRNA_model.species_lower_bound = 0.0;
M_ile_L_c_tRNA_model.species_upper_bound = 0.0;
M_ile_L_c_tRNA_model.is_species_measured = false;
M_ile_L_c_tRNA_model.is_biomass_precursor = false;
M_ile_L_c_tRNA_model.biomass_precursor_coefficient = 0.0;
M_ile_L_c_tRNA_model.species_time_constant = 1.0;
M_ile_L_c_tRNA_model.species_initial_condition = 0.0;
M_ile_L_c_tRNA_model.is_species_diluted = true;
M_ile_L_c_tRNA_model.is_species_extracellular = false;
species_model_dictionary["M_ile_L_c_tRNA"] = M_ile_L_c_tRNA_model;
M_ile_L_c_tRNA_model = 0;

# species_symbol: 120 M_leu_L_c_tRNA -
M_leu_L_c_tRNA_model = SpeciesModel();
M_leu_L_c_tRNA_model.species_index = 120;
M_leu_L_c_tRNA_model.species_symbol = string("M_leu_L_c_tRNA");
M_leu_L_c_tRNA_model.species_constraint_type = GLPK.FX;
M_leu_L_c_tRNA_model.species_lower_bound = 0.0;
M_leu_L_c_tRNA_model.species_upper_bound = 0.0;
M_leu_L_c_tRNA_model.is_species_measured = false;
M_leu_L_c_tRNA_model.is_biomass_precursor = false;
M_leu_L_c_tRNA_model.biomass_precursor_coefficient = 0.0;
M_leu_L_c_tRNA_model.species_time_constant = 1.0;
M_leu_L_c_tRNA_model.species_initial_condition = 0.0;
M_leu_L_c_tRNA_model.is_species_diluted = true;
M_leu_L_c_tRNA_model.is_species_extracellular = false;
species_model_dictionary["M_leu_L_c_tRNA"] = M_leu_L_c_tRNA_model;
M_leu_L_c_tRNA_model = 0;

# species_symbol: 121 M_lys_L_c_tRNA -
M_lys_L_c_tRNA_model = SpeciesModel();
M_lys_L_c_tRNA_model.species_index = 121;
M_lys_L_c_tRNA_model.species_symbol = string("M_lys_L_c_tRNA");
M_lys_L_c_tRNA_model.species_constraint_type = GLPK.FX;
M_lys_L_c_tRNA_model.species_lower_bound = 0.0;
M_lys_L_c_tRNA_model.species_upper_bound = 0.0;
M_lys_L_c_tRNA_model.is_species_measured = false;
M_lys_L_c_tRNA_model.is_biomass_precursor = false;
M_lys_L_c_tRNA_model.biomass_precursor_coefficient = 0.0;
M_lys_L_c_tRNA_model.species_time_constant = 1.0;
M_lys_L_c_tRNA_model.species_initial_condition = 0.0;
M_lys_L_c_tRNA_model.is_species_diluted = true;
M_lys_L_c_tRNA_model.is_species_extracellular = false;
species_model_dictionary["M_lys_L_c_tRNA"] = M_lys_L_c_tRNA_model;
M_lys_L_c_tRNA_model = 0;

# species_symbol: 122 M_met_L_c_tRNA -
M_met_L_c_tRNA_model = SpeciesModel();
M_met_L_c_tRNA_model.species_index = 122;
M_met_L_c_tRNA_model.species_symbol = string("M_met_L_c_tRNA");
M_met_L_c_tRNA_model.species_constraint_type = GLPK.FX;
M_met_L_c_tRNA_model.species_lower_bound = 0.0;
M_met_L_c_tRNA_model.species_upper_bound = 0.0;
M_met_L_c_tRNA_model.is_species_measured = false;
M_met_L_c_tRNA_model.is_biomass_precursor = false;
M_met_L_c_tRNA_model.biomass_precursor_coefficient = 0.0;
M_met_L_c_tRNA_model.species_time_constant = 1.0;
M_met_L_c_tRNA_model.species_initial_condition = 0.0;
M_met_L_c_tRNA_model.is_species_diluted = true;
M_met_L_c_tRNA_model.is_species_extracellular = false;
species_model_dictionary["M_met_L_c_tRNA"] = M_met_L_c_tRNA_model;
M_met_L_c_tRNA_model = 0;

# species_symbol: 123 M_phe_L_c_tRNA -
M_phe_L_c_tRNA_model = SpeciesModel();
M_phe_L_c_tRNA_model.species_index = 123;
M_phe_L_c_tRNA_model.species_symbol = string("M_phe_L_c_tRNA");
M_phe_L_c_tRNA_model.species_constraint_type = GLPK.FX;
M_phe_L_c_tRNA_model.species_lower_bound = 0.0;
M_phe_L_c_tRNA_model.species_upper_bound = 0.0;
M_phe_L_c_tRNA_model.is_species_measured = false;
M_phe_L_c_tRNA_model.is_biomass_precursor = false;
M_phe_L_c_tRNA_model.biomass_precursor_coefficient = 0.0;
M_phe_L_c_tRNA_model.species_time_constant = 1.0;
M_phe_L_c_tRNA_model.species_initial_condition = 0.0;
M_phe_L_c_tRNA_model.is_species_diluted = true;
M_phe_L_c_tRNA_model.is_species_extracellular = false;
species_model_dictionary["M_phe_L_c_tRNA"] = M_phe_L_c_tRNA_model;
M_phe_L_c_tRNA_model = 0;

# species_symbol: 124 M_pro_L_c_tRNA -
M_pro_L_c_tRNA_model = SpeciesModel();
M_pro_L_c_tRNA_model.species_index = 124;
M_pro_L_c_tRNA_model.species_symbol = string("M_pro_L_c_tRNA");
M_pro_L_c_tRNA_model.species_constraint_type = GLPK.FX;
M_pro_L_c_tRNA_model.species_lower_bound = 0.0;
M_pro_L_c_tRNA_model.species_upper_bound = 0.0;
M_pro_L_c_tRNA_model.is_species_measured = false;
M_pro_L_c_tRNA_model.is_biomass_precursor = false;
M_pro_L_c_tRNA_model.biomass_precursor_coefficient = 0.0;
M_pro_L_c_tRNA_model.species_time_constant = 1.0;
M_pro_L_c_tRNA_model.species_initial_condition = 0.0;
M_pro_L_c_tRNA_model.is_species_diluted = true;
M_pro_L_c_tRNA_model.is_species_extracellular = false;
species_model_dictionary["M_pro_L_c_tRNA"] = M_pro_L_c_tRNA_model;
M_pro_L_c_tRNA_model = 0;

# species_symbol: 125 M_ser_L_c_tRNA -
M_ser_L_c_tRNA_model = SpeciesModel();
M_ser_L_c_tRNA_model.species_index = 125;
M_ser_L_c_tRNA_model.species_symbol = string("M_ser_L_c_tRNA");
M_ser_L_c_tRNA_model.species_constraint_type = GLPK.FX;
M_ser_L_c_tRNA_model.species_lower_bound = 0.0;
M_ser_L_c_tRNA_model.species_upper_bound = 0.0;
M_ser_L_c_tRNA_model.is_species_measured = false;
M_ser_L_c_tRNA_model.is_biomass_precursor = false;
M_ser_L_c_tRNA_model.biomass_precursor_coefficient = 0.0;
M_ser_L_c_tRNA_model.species_time_constant = 1.0;
M_ser_L_c_tRNA_model.species_initial_condition = 0.0;
M_ser_L_c_tRNA_model.is_species_diluted = true;
M_ser_L_c_tRNA_model.is_species_extracellular = false;
species_model_dictionary["M_ser_L_c_tRNA"] = M_ser_L_c_tRNA_model;
M_ser_L_c_tRNA_model = 0;

# species_symbol: 126 M_thr_L_c_tRNA -
M_thr_L_c_tRNA_model = SpeciesModel();
M_thr_L_c_tRNA_model.species_index = 126;
M_thr_L_c_tRNA_model.species_symbol = string("M_thr_L_c_tRNA");
M_thr_L_c_tRNA_model.species_constraint_type = GLPK.FX;
M_thr_L_c_tRNA_model.species_lower_bound = 0.0;
M_thr_L_c_tRNA_model.species_upper_bound = 0.0;
M_thr_L_c_tRNA_model.is_species_measured = false;
M_thr_L_c_tRNA_model.is_biomass_precursor = false;
M_thr_L_c_tRNA_model.biomass_precursor_coefficient = 0.0;
M_thr_L_c_tRNA_model.species_time_constant = 1.0;
M_thr_L_c_tRNA_model.species_initial_condition = 0.0;
M_thr_L_c_tRNA_model.is_species_diluted = true;
M_thr_L_c_tRNA_model.is_species_extracellular = false;
species_model_dictionary["M_thr_L_c_tRNA"] = M_thr_L_c_tRNA_model;
M_thr_L_c_tRNA_model = 0;

# species_symbol: 127 M_trp_L_c_tRNA -
M_trp_L_c_tRNA_model = SpeciesModel();
M_trp_L_c_tRNA_model.species_index = 127;
M_trp_L_c_tRNA_model.species_symbol = string("M_trp_L_c_tRNA");
M_trp_L_c_tRNA_model.species_constraint_type = GLPK.FX;
M_trp_L_c_tRNA_model.species_lower_bound = 0.0;
M_trp_L_c_tRNA_model.species_upper_bound = 0.0;
M_trp_L_c_tRNA_model.is_species_measured = false;
M_trp_L_c_tRNA_model.is_biomass_precursor = false;
M_trp_L_c_tRNA_model.biomass_precursor_coefficient = 0.0;
M_trp_L_c_tRNA_model.species_time_constant = 1.0;
M_trp_L_c_tRNA_model.species_initial_condition = 0.0;
M_trp_L_c_tRNA_model.is_species_diluted = true;
M_trp_L_c_tRNA_model.is_species_extracellular = false;
species_model_dictionary["M_trp_L_c_tRNA"] = M_trp_L_c_tRNA_model;
M_trp_L_c_tRNA_model = 0;

# species_symbol: 128 M_tyr_L_c_tRNA -
M_tyr_L_c_tRNA_model = SpeciesModel();
M_tyr_L_c_tRNA_model.species_index = 128;
M_tyr_L_c_tRNA_model.species_symbol = string("M_tyr_L_c_tRNA");
M_tyr_L_c_tRNA_model.species_constraint_type = GLPK.FX;
M_tyr_L_c_tRNA_model.species_lower_bound = 0.0;
M_tyr_L_c_tRNA_model.species_upper_bound = 0.0;
M_tyr_L_c_tRNA_model.is_species_measured = false;
M_tyr_L_c_tRNA_model.is_biomass_precursor = false;
M_tyr_L_c_tRNA_model.biomass_precursor_coefficient = 0.0;
M_tyr_L_c_tRNA_model.species_time_constant = 1.0;
M_tyr_L_c_tRNA_model.species_initial_condition = 0.0;
M_tyr_L_c_tRNA_model.is_species_diluted = true;
M_tyr_L_c_tRNA_model.is_species_extracellular = false;
species_model_dictionary["M_tyr_L_c_tRNA"] = M_tyr_L_c_tRNA_model;
M_tyr_L_c_tRNA_model = 0;

# species_symbol: 129 M_val_L_c_tRNA -
M_val_L_c_tRNA_model = SpeciesModel();
M_val_L_c_tRNA_model.species_index = 129;
M_val_L_c_tRNA_model.species_symbol = string("M_val_L_c_tRNA");
M_val_L_c_tRNA_model.species_constraint_type = GLPK.FX;
M_val_L_c_tRNA_model.species_lower_bound = 0.0;
M_val_L_c_tRNA_model.species_upper_bound = 0.0;
M_val_L_c_tRNA_model.is_species_measured = false;
M_val_L_c_tRNA_model.is_biomass_precursor = false;
M_val_L_c_tRNA_model.biomass_precursor_coefficient = 0.0;
M_val_L_c_tRNA_model.species_time_constant = 1.0;
M_val_L_c_tRNA_model.species_initial_condition = 0.0;
M_val_L_c_tRNA_model.is_species_diluted = true;
M_val_L_c_tRNA_model.is_species_extracellular = false;
species_model_dictionary["M_val_L_c_tRNA"] = M_val_L_c_tRNA_model;
M_val_L_c_tRNA_model = 0;

# species_symbol: 130 PROTEIN_deGFP -
PROTEIN_deGFP_model = SpeciesModel();
PROTEIN_deGFP_model.species_index = 130;
PROTEIN_deGFP_model.species_symbol = string("PROTEIN_deGFP");
PROTEIN_deGFP_model.species_constraint_type = GLPK.FX;
PROTEIN_deGFP_model.species_lower_bound = 0.0;
PROTEIN_deGFP_model.species_upper_bound = 0.0;
PROTEIN_deGFP_model.is_species_measured = false;
PROTEIN_deGFP_model.is_biomass_precursor = false;
PROTEIN_deGFP_model.biomass_precursor_coefficient = 0.0;
PROTEIN_deGFP_model.species_time_constant = 1.0;
PROTEIN_deGFP_model.species_initial_condition = 0.0;
PROTEIN_deGFP_model.is_species_diluted = true;
PROTEIN_deGFP_model.is_species_extracellular = false;
species_model_dictionary["PROTEIN_deGFP"] = PROTEIN_deGFP_model;
PROTEIN_deGFP_model = 0;

# species_symbol: 131 tRNA -
tRNA_model = SpeciesModel();
tRNA_model.species_index = 131;
tRNA_model.species_symbol = string("tRNA");
tRNA_model.species_constraint_type = GLPK.FX;
tRNA_model.species_lower_bound = 0.0;
tRNA_model.species_upper_bound = 0.0;
tRNA_model.is_species_measured = false;
tRNA_model.is_biomass_precursor = false;
tRNA_model.biomass_precursor_coefficient = 0.0;
tRNA_model.species_time_constant = 1.0;
tRNA_model.species_initial_condition = 0.0;
tRNA_model.is_species_diluted = true;
tRNA_model.is_species_extracellular = false;
species_model_dictionary["tRNA"] = tRNA_model;
tRNA_model = 0;

return species_model_dictionary;
end

# ----------------------------------------------------------------------------------- #
# Helper function: buildFluxModelDictionary
# Constructs a dictionary of flux models
# Generated using the Kwatee code generation system
#
# Input arguments:
# N/A
#
# Return arguments:
# flux_model_dictionary  - Dictionary of FluxModels key'd by flux symbol
# ----------------------------------------------------------------------------------- #
function buildFluxModelDictionary()

# function variables -
flux_model_dictionary = Dict{AbstractString,FluxModel}();

# 1 R_malS: M_maltose_c+M_h2o_c -([])-> 2*M_glc_D_c
R_malS_model = FluxModel();
R_malS_model.flux_index = 1
R_malS_model.flux_symbol = "R_malS"
R_malS_model.flux_constraint_type = GLPK.DB;
R_malS_model.flux_lower_bound = 0.0;
R_malS_model.flux_upper_bound = 1.0;
R_malS_model.flux_bounds_model = Bounds;
R_malS_model.flux_gamma_array = vec([1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_malS_model.flux_bound_alpha = 1.0;
R_malS_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_malS"] = R_malS_model;
R_malS_model = 0;

# 2 R_PTS: M_pep_c+M_glc_D_c -([])-> M_g6p_c+M_pyr_c
R_PTS_model = FluxModel();
R_PTS_model.flux_index = 2
R_PTS_model.flux_symbol = "R_PTS"
R_PTS_model.flux_constraint_type = GLPK.DB;
R_PTS_model.flux_lower_bound = 0.0;
R_PTS_model.flux_upper_bound = 1.0;
R_PTS_model.flux_bounds_model = Bounds;
R_PTS_model.flux_gamma_array = vec([0.0 0.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_PTS_model.flux_bound_alpha = 1.0;
R_PTS_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_PTS"] = R_PTS_model;
R_PTS_model = 0;

# 3 R_glk_atp: M_atp_c+M_glc_D_c -([])-> M_adp_c+M_g6p_c+M_h_c
R_glk_atp_model = FluxModel();
R_glk_atp_model.flux_index = 3
R_glk_atp_model.flux_symbol = "R_glk_atp"
R_glk_atp_model.flux_constraint_type = GLPK.DB;
R_glk_atp_model.flux_lower_bound = 0.0;
R_glk_atp_model.flux_upper_bound = 1.0;
R_glk_atp_model.flux_bounds_model = Bounds;
R_glk_atp_model.flux_gamma_array = vec([0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_glk_atp_model.flux_bound_alpha = 1.0;
R_glk_atp_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_glk_atp"] = R_glk_atp_model;
R_glk_atp_model = 0;

# 4 R_glk_utp: M_utp_c+M_glc_D_c -([])-> M_udp_c+M_g6p_c+M_h_c
R_glk_utp_model = FluxModel();
R_glk_utp_model.flux_index = 4
R_glk_utp_model.flux_symbol = "R_glk_utp"
R_glk_utp_model.flux_constraint_type = GLPK.DB;
R_glk_utp_model.flux_lower_bound = 0.0;
R_glk_utp_model.flux_upper_bound = 1.0;
R_glk_utp_model.flux_bounds_model = Bounds;
R_glk_utp_model.flux_gamma_array = vec([0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_glk_utp_model.flux_bound_alpha = 1.0;
R_glk_utp_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_glk_utp"] = R_glk_utp_model;
R_glk_utp_model = 0;

# 5 R_glk_ctp: M_ctp_c+M_glc_D_c -([])-> M_cdp_c+M_g6p_c+M_h_c
R_glk_ctp_model = FluxModel();
R_glk_ctp_model.flux_index = 5
R_glk_ctp_model.flux_symbol = "R_glk_ctp"
R_glk_ctp_model.flux_constraint_type = GLPK.DB;
R_glk_ctp_model.flux_lower_bound = 0.0;
R_glk_ctp_model.flux_upper_bound = 1.0;
R_glk_ctp_model.flux_bounds_model = Bounds;
R_glk_ctp_model.flux_gamma_array = vec([0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_glk_ctp_model.flux_bound_alpha = 1.0;
R_glk_ctp_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_glk_ctp"] = R_glk_ctp_model;
R_glk_ctp_model = 0;

# 6 R_glk_gtp: M_gtp_c+M_glc_D_c -([])-> M_gdp_c+M_g6p_c+M_h_c
R_glk_gtp_model = FluxModel();
R_glk_gtp_model.flux_index = 6
R_glk_gtp_model.flux_symbol = "R_glk_gtp"
R_glk_gtp_model.flux_constraint_type = GLPK.DB;
R_glk_gtp_model.flux_lower_bound = 0.0;
R_glk_gtp_model.flux_upper_bound = 1.0;
R_glk_gtp_model.flux_bounds_model = Bounds;
R_glk_gtp_model.flux_gamma_array = vec([0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_glk_gtp_model.flux_bound_alpha = 1.0;
R_glk_gtp_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_glk_gtp"] = R_glk_gtp_model;
R_glk_gtp_model = 0;

# 7 R_pgi: M_g6p_c -([])-> M_f6p_c
R_pgi_model = FluxModel();
R_pgi_model.flux_index = 7
R_pgi_model.flux_symbol = "R_pgi"
R_pgi_model.flux_constraint_type = GLPK.DB;
R_pgi_model.flux_lower_bound = 0.0;
R_pgi_model.flux_upper_bound = 1.0;
R_pgi_model.flux_bounds_model = Bounds;
R_pgi_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_pgi_model.flux_bound_alpha = 1.0;
R_pgi_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_pgi"] = R_pgi_model;
R_pgi_model = 0;

# 8 -1*(R_pgi: M_g6p_c -([])-> M_f6p_c)
R_pgi_reverse_model = FluxModel();
R_pgi_reverse_model.flux_index = 8
R_pgi_reverse_model.flux_symbol = "R_pgi_reverse"
R_pgi_reverse_model.flux_constraint_type = GLPK.DB;
R_pgi_reverse_model.flux_lower_bound = 0.0;
R_pgi_reverse_model.flux_upper_bound = 1.0;
R_pgi_reverse_model.flux_bounds_model = Bounds;
R_pgi_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_pgi_reverse_model.flux_bound_alpha = 1.0;
R_pgi_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_pgi_reverse"] = R_pgi_reverse_model;
R_pgi_reverse_model = 0;

# 9 R_pfk: M_atp_c+M_f6p_c -([])-> M_adp_c+M_fdp_c+M_h_c
R_pfk_model = FluxModel();
R_pfk_model.flux_index = 9
R_pfk_model.flux_symbol = "R_pfk"
R_pfk_model.flux_constraint_type = GLPK.DB;
R_pfk_model.flux_lower_bound = 0.0;
R_pfk_model.flux_upper_bound = 1.0;
R_pfk_model.flux_bounds_model = Bounds;
R_pfk_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_pfk_model.flux_bound_alpha = 1.0;
R_pfk_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_pfk"] = R_pfk_model;
R_pfk_model = 0;

# 10 R_fdp: M_fdp_c+M_h2o_c -([])-> M_f6p_c+M_pi_c
R_fdp_model = FluxModel();
R_fdp_model.flux_index = 10
R_fdp_model.flux_symbol = "R_fdp"
R_fdp_model.flux_constraint_type = GLPK.DB;
R_fdp_model.flux_lower_bound = 0.0;
R_fdp_model.flux_upper_bound = 1.0;
R_fdp_model.flux_bounds_model = Bounds;
R_fdp_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_fdp_model.flux_bound_alpha = 1.0;
R_fdp_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_fdp"] = R_fdp_model;
R_fdp_model = 0;

# 11 R_fbaA: M_fdp_c -([])-> M_dhap_c+M_g3p_c
R_fbaA_model = FluxModel();
R_fbaA_model.flux_index = 11
R_fbaA_model.flux_symbol = "R_fbaA"
R_fbaA_model.flux_constraint_type = GLPK.DB;
R_fbaA_model.flux_lower_bound = 0.0;
R_fbaA_model.flux_upper_bound = 1.0;
R_fbaA_model.flux_bounds_model = Bounds;
R_fbaA_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_fbaA_model.flux_bound_alpha = 1.0;
R_fbaA_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_fbaA"] = R_fbaA_model;
R_fbaA_model = 0;

# 12 -1*(R_fbaA: M_fdp_c -([])-> M_dhap_c+M_g3p_c)
R_fbaA_reverse_model = FluxModel();
R_fbaA_reverse_model.flux_index = 12
R_fbaA_reverse_model.flux_symbol = "R_fbaA_reverse"
R_fbaA_reverse_model.flux_constraint_type = GLPK.DB;
R_fbaA_reverse_model.flux_lower_bound = 0.0;
R_fbaA_reverse_model.flux_upper_bound = 1.0;
R_fbaA_reverse_model.flux_bounds_model = Bounds;
R_fbaA_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_fbaA_reverse_model.flux_bound_alpha = 1.0;
R_fbaA_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_fbaA_reverse"] = R_fbaA_reverse_model;
R_fbaA_reverse_model = 0;

# 13 R_tpiA: M_dhap_c -([])-> M_g3p_c
R_tpiA_model = FluxModel();
R_tpiA_model.flux_index = 13
R_tpiA_model.flux_symbol = "R_tpiA"
R_tpiA_model.flux_constraint_type = GLPK.DB;
R_tpiA_model.flux_lower_bound = 0.0;
R_tpiA_model.flux_upper_bound = 1.0;
R_tpiA_model.flux_bounds_model = Bounds;
R_tpiA_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_tpiA_model.flux_bound_alpha = 1.0;
R_tpiA_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_tpiA"] = R_tpiA_model;
R_tpiA_model = 0;

# 14 -1*(R_tpiA: M_dhap_c -([])-> M_g3p_c)
R_tpiA_reverse_model = FluxModel();
R_tpiA_reverse_model.flux_index = 14
R_tpiA_reverse_model.flux_symbol = "R_tpiA_reverse"
R_tpiA_reverse_model.flux_constraint_type = GLPK.DB;
R_tpiA_reverse_model.flux_lower_bound = 0.0;
R_tpiA_reverse_model.flux_upper_bound = 1.0;
R_tpiA_reverse_model.flux_bounds_model = Bounds;
R_tpiA_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_tpiA_reverse_model.flux_bound_alpha = 1.0;
R_tpiA_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_tpiA_reverse"] = R_tpiA_reverse_model;
R_tpiA_reverse_model = 0;

# 15 R_gapA: M_g3p_c+M_nad_c+M_pi_c -([])-> M_13dpg_c+M_h_c+M_nadh_c
R_gapA_model = FluxModel();
R_gapA_model.flux_index = 15
R_gapA_model.flux_symbol = "R_gapA"
R_gapA_model.flux_constraint_type = GLPK.DB;
R_gapA_model.flux_lower_bound = 0.0;
R_gapA_model.flux_upper_bound = 1.0;
R_gapA_model.flux_bounds_model = Bounds;
R_gapA_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_gapA_model.flux_bound_alpha = 1.0;
R_gapA_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_gapA"] = R_gapA_model;
R_gapA_model = 0;

# 16 -1*(R_gapA: M_g3p_c+M_nad_c+M_pi_c -([])-> M_13dpg_c+M_h_c+M_nadh_c)
R_gapA_reverse_model = FluxModel();
R_gapA_reverse_model.flux_index = 16
R_gapA_reverse_model.flux_symbol = "R_gapA_reverse"
R_gapA_reverse_model.flux_constraint_type = GLPK.DB;
R_gapA_reverse_model.flux_lower_bound = 0.0;
R_gapA_reverse_model.flux_upper_bound = 1.0;
R_gapA_reverse_model.flux_bounds_model = Bounds;
R_gapA_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_gapA_reverse_model.flux_bound_alpha = 1.0;
R_gapA_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_gapA_reverse"] = R_gapA_reverse_model;
R_gapA_reverse_model = 0;

# 17 R_pgk: M_13dpg_c+M_adp_c -([])-> M_3pg_c+M_atp_c
R_pgk_model = FluxModel();
R_pgk_model.flux_index = 17
R_pgk_model.flux_symbol = "R_pgk"
R_pgk_model.flux_constraint_type = GLPK.DB;
R_pgk_model.flux_lower_bound = 0.0;
R_pgk_model.flux_upper_bound = 1.0;
R_pgk_model.flux_bounds_model = Bounds;
R_pgk_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_pgk_model.flux_bound_alpha = 1.0;
R_pgk_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_pgk"] = R_pgk_model;
R_pgk_model = 0;

# 18 -1*(R_pgk: M_13dpg_c+M_adp_c -([])-> M_3pg_c+M_atp_c)
R_pgk_reverse_model = FluxModel();
R_pgk_reverse_model.flux_index = 18
R_pgk_reverse_model.flux_symbol = "R_pgk_reverse"
R_pgk_reverse_model.flux_constraint_type = GLPK.DB;
R_pgk_reverse_model.flux_lower_bound = 0.0;
R_pgk_reverse_model.flux_upper_bound = 1.0;
R_pgk_reverse_model.flux_bounds_model = Bounds;
R_pgk_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_pgk_reverse_model.flux_bound_alpha = 1.0;
R_pgk_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_pgk_reverse"] = R_pgk_reverse_model;
R_pgk_reverse_model = 0;

# 19 R_gpm: M_3pg_c -([])-> M_2pg_c
R_gpm_model = FluxModel();
R_gpm_model.flux_index = 19
R_gpm_model.flux_symbol = "R_gpm"
R_gpm_model.flux_constraint_type = GLPK.DB;
R_gpm_model.flux_lower_bound = 0.0;
R_gpm_model.flux_upper_bound = 1.0;
R_gpm_model.flux_bounds_model = Bounds;
R_gpm_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_gpm_model.flux_bound_alpha = 1.0;
R_gpm_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_gpm"] = R_gpm_model;
R_gpm_model = 0;

# 20 -1*(R_gpm: M_3pg_c -([])-> M_2pg_c)
R_gpm_reverse_model = FluxModel();
R_gpm_reverse_model.flux_index = 20
R_gpm_reverse_model.flux_symbol = "R_gpm_reverse"
R_gpm_reverse_model.flux_constraint_type = GLPK.DB;
R_gpm_reverse_model.flux_lower_bound = 0.0;
R_gpm_reverse_model.flux_upper_bound = 1.0;
R_gpm_reverse_model.flux_bounds_model = Bounds;
R_gpm_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_gpm_reverse_model.flux_bound_alpha = 1.0;
R_gpm_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_gpm_reverse"] = R_gpm_reverse_model;
R_gpm_reverse_model = 0;

# 21 R_eno: M_2pg_c -([])-> M_h2o_c+M_pep_c
R_eno_model = FluxModel();
R_eno_model.flux_index = 21
R_eno_model.flux_symbol = "R_eno"
R_eno_model.flux_constraint_type = GLPK.DB;
R_eno_model.flux_lower_bound = 0.0;
R_eno_model.flux_upper_bound = 1.0;
R_eno_model.flux_bounds_model = Bounds;
R_eno_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_eno_model.flux_bound_alpha = 1.0;
R_eno_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_eno"] = R_eno_model;
R_eno_model = 0;

# 22 -1*(R_eno: M_2pg_c -([])-> M_h2o_c+M_pep_c)
R_eno_reverse_model = FluxModel();
R_eno_reverse_model.flux_index = 22
R_eno_reverse_model.flux_symbol = "R_eno_reverse"
R_eno_reverse_model.flux_constraint_type = GLPK.DB;
R_eno_reverse_model.flux_lower_bound = 0.0;
R_eno_reverse_model.flux_upper_bound = 1.0;
R_eno_reverse_model.flux_bounds_model = Bounds;
R_eno_reverse_model.flux_gamma_array = vec([0.0 1.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_eno_reverse_model.flux_bound_alpha = 1.0;
R_eno_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_eno_reverse"] = R_eno_reverse_model;
R_eno_reverse_model = 0;

# 23 R_pyk: M_adp_c+M_h_c+M_pep_c -([])-> M_atp_c+M_pyr_c
R_pyk_model = FluxModel();
R_pyk_model.flux_index = 23
R_pyk_model.flux_symbol = "R_pyk"
R_pyk_model.flux_constraint_type = GLPK.DB;
R_pyk_model.flux_lower_bound = 0.0;
R_pyk_model.flux_upper_bound = 1.0;
R_pyk_model.flux_bounds_model = Bounds;
R_pyk_model.flux_gamma_array = vec([0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_pyk_model.flux_bound_alpha = 1.0;
R_pyk_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_pyk"] = R_pyk_model;
R_pyk_model = 0;

# 24 R_pck: M_atp_c+M_oaa_c -([])-> M_adp_c+M_co2_c+M_pep_c
R_pck_model = FluxModel();
R_pck_model.flux_index = 24
R_pck_model.flux_symbol = "R_pck"
R_pck_model.flux_constraint_type = GLPK.DB;
R_pck_model.flux_lower_bound = 0.0;
R_pck_model.flux_upper_bound = 1.0;
R_pck_model.flux_bounds_model = Bounds;
R_pck_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_pck_model.flux_bound_alpha = 1.0;
R_pck_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_pck"] = R_pck_model;
R_pck_model = 0;

# 25 R_ppc: M_co2_c+M_h2o_c+M_pep_c -([])-> M_h_c+M_oaa_c+M_pi_c
R_ppc_model = FluxModel();
R_ppc_model.flux_index = 25
R_ppc_model.flux_symbol = "R_ppc"
R_ppc_model.flux_constraint_type = GLPK.DB;
R_ppc_model.flux_lower_bound = 0.0;
R_ppc_model.flux_upper_bound = 1.0;
R_ppc_model.flux_bounds_model = Bounds;
R_ppc_model.flux_gamma_array = vec([0.0 1.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_ppc_model.flux_bound_alpha = 1.0;
R_ppc_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_ppc"] = R_ppc_model;
R_ppc_model = 0;

# 26 R_pdh: M_coa_c+M_nad_c+M_pyr_c -([])-> M_accoa_c+M_co2_c+M_nadh_c
R_pdh_model = FluxModel();
R_pdh_model.flux_index = 26
R_pdh_model.flux_symbol = "R_pdh"
R_pdh_model.flux_constraint_type = GLPK.DB;
R_pdh_model.flux_lower_bound = 0.0;
R_pdh_model.flux_upper_bound = 1.0;
R_pdh_model.flux_bounds_model = Bounds;
R_pdh_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_pdh_model.flux_bound_alpha = 1.0;
R_pdh_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_pdh"] = R_pdh_model;
R_pdh_model = 0;

# 27 R_pps: M_atp_c+M_h2o_c+M_pyr_c -([])-> M_amp_c+2*M_h_c+M_pep_c+M_pi_c
R_pps_model = FluxModel();
R_pps_model.flux_index = 27
R_pps_model.flux_symbol = "R_pps"
R_pps_model.flux_constraint_type = GLPK.DB;
R_pps_model.flux_lower_bound = 0.0;
R_pps_model.flux_upper_bound = 1.0;
R_pps_model.flux_bounds_model = Bounds;
R_pps_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_pps_model.flux_bound_alpha = 1.0;
R_pps_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_pps"] = R_pps_model;
R_pps_model = 0;

# 28 R_zwf: M_g6p_c+M_nadp_c -([])-> M_6pgl_c+M_h_c+M_nadph_c
R_zwf_model = FluxModel();
R_zwf_model.flux_index = 28
R_zwf_model.flux_symbol = "R_zwf"
R_zwf_model.flux_constraint_type = GLPK.DB;
R_zwf_model.flux_lower_bound = 0.0;
R_zwf_model.flux_upper_bound = 1.0;
R_zwf_model.flux_bounds_model = Bounds;
R_zwf_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_zwf_model.flux_bound_alpha = 1.0;
R_zwf_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_zwf"] = R_zwf_model;
R_zwf_model = 0;

# 29 -1*(R_zwf: M_g6p_c+M_nadp_c -([])-> M_6pgl_c+M_h_c+M_nadph_c)
R_zwf_reverse_model = FluxModel();
R_zwf_reverse_model.flux_index = 29
R_zwf_reverse_model.flux_symbol = "R_zwf_reverse"
R_zwf_reverse_model.flux_constraint_type = GLPK.DB;
R_zwf_reverse_model.flux_lower_bound = 0.0;
R_zwf_reverse_model.flux_upper_bound = 1.0;
R_zwf_reverse_model.flux_bounds_model = Bounds;
R_zwf_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_zwf_reverse_model.flux_bound_alpha = 1.0;
R_zwf_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_zwf_reverse"] = R_zwf_reverse_model;
R_zwf_reverse_model = 0;

# 30 R_pgl: M_6pgl_c+M_h2o_c -([])-> M_6pgc_c+M_h_c
R_pgl_model = FluxModel();
R_pgl_model.flux_index = 30
R_pgl_model.flux_symbol = "R_pgl"
R_pgl_model.flux_constraint_type = GLPK.DB;
R_pgl_model.flux_lower_bound = 0.0;
R_pgl_model.flux_upper_bound = 1.0;
R_pgl_model.flux_bounds_model = Bounds;
R_pgl_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_pgl_model.flux_bound_alpha = 1.0;
R_pgl_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_pgl"] = R_pgl_model;
R_pgl_model = 0;

# 31 R_gnd: M_6pgc_c+M_nadp_c -([])-> M_co2_c+M_nadph_c+M_ru5p_D_c
R_gnd_model = FluxModel();
R_gnd_model.flux_index = 31
R_gnd_model.flux_symbol = "R_gnd"
R_gnd_model.flux_constraint_type = GLPK.DB;
R_gnd_model.flux_lower_bound = 0.0;
R_gnd_model.flux_upper_bound = 1.0;
R_gnd_model.flux_bounds_model = Bounds;
R_gnd_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_gnd_model.flux_bound_alpha = 1.0;
R_gnd_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_gnd"] = R_gnd_model;
R_gnd_model = 0;

# 32 R_rpe: M_ru5p_D_c -([])-> M_xu5p_D_c
R_rpe_model = FluxModel();
R_rpe_model.flux_index = 32
R_rpe_model.flux_symbol = "R_rpe"
R_rpe_model.flux_constraint_type = GLPK.DB;
R_rpe_model.flux_lower_bound = 0.0;
R_rpe_model.flux_upper_bound = 1.0;
R_rpe_model.flux_bounds_model = Bounds;
R_rpe_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_rpe_model.flux_bound_alpha = 1.0;
R_rpe_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_rpe"] = R_rpe_model;
R_rpe_model = 0;

# 33 -1*(R_rpe: M_ru5p_D_c -([])-> M_xu5p_D_c)
R_rpe_reverse_model = FluxModel();
R_rpe_reverse_model.flux_index = 33
R_rpe_reverse_model.flux_symbol = "R_rpe_reverse"
R_rpe_reverse_model.flux_constraint_type = GLPK.DB;
R_rpe_reverse_model.flux_lower_bound = 0.0;
R_rpe_reverse_model.flux_upper_bound = 1.0;
R_rpe_reverse_model.flux_bounds_model = Bounds;
R_rpe_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_rpe_reverse_model.flux_bound_alpha = 1.0;
R_rpe_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_rpe_reverse"] = R_rpe_reverse_model;
R_rpe_reverse_model = 0;

# 34 R_rpi: M_r5p_c -([])-> M_ru5p_D_c
R_rpi_model = FluxModel();
R_rpi_model.flux_index = 34
R_rpi_model.flux_symbol = "R_rpi"
R_rpi_model.flux_constraint_type = GLPK.DB;
R_rpi_model.flux_lower_bound = 0.0;
R_rpi_model.flux_upper_bound = 1.0;
R_rpi_model.flux_bounds_model = Bounds;
R_rpi_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_rpi_model.flux_bound_alpha = 1.0;
R_rpi_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_rpi"] = R_rpi_model;
R_rpi_model = 0;

# 35 -1*(R_rpi: M_r5p_c -([])-> M_ru5p_D_c)
R_rpi_reverse_model = FluxModel();
R_rpi_reverse_model.flux_index = 35
R_rpi_reverse_model.flux_symbol = "R_rpi_reverse"
R_rpi_reverse_model.flux_constraint_type = GLPK.DB;
R_rpi_reverse_model.flux_lower_bound = 0.0;
R_rpi_reverse_model.flux_upper_bound = 1.0;
R_rpi_reverse_model.flux_bounds_model = Bounds;
R_rpi_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_rpi_reverse_model.flux_bound_alpha = 1.0;
R_rpi_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_rpi_reverse"] = R_rpi_reverse_model;
R_rpi_reverse_model = 0;

# 36 R_talAB: M_g3p_c+M_s7p_c -([])-> M_e4p_c+M_f6p_c
R_talAB_model = FluxModel();
R_talAB_model.flux_index = 36
R_talAB_model.flux_symbol = "R_talAB"
R_talAB_model.flux_constraint_type = GLPK.DB;
R_talAB_model.flux_lower_bound = 0.0;
R_talAB_model.flux_upper_bound = 1.0;
R_talAB_model.flux_bounds_model = Bounds;
R_talAB_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_talAB_model.flux_bound_alpha = 1.0;
R_talAB_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_talAB"] = R_talAB_model;
R_talAB_model = 0;

# 37 -1*(R_talAB: M_g3p_c+M_s7p_c -([])-> M_e4p_c+M_f6p_c)
R_talAB_reverse_model = FluxModel();
R_talAB_reverse_model.flux_index = 37
R_talAB_reverse_model.flux_symbol = "R_talAB_reverse"
R_talAB_reverse_model.flux_constraint_type = GLPK.DB;
R_talAB_reverse_model.flux_lower_bound = 0.0;
R_talAB_reverse_model.flux_upper_bound = 1.0;
R_talAB_reverse_model.flux_bounds_model = Bounds;
R_talAB_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_talAB_reverse_model.flux_bound_alpha = 1.0;
R_talAB_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_talAB_reverse"] = R_talAB_reverse_model;
R_talAB_reverse_model = 0;

# 38 R_tkt1: M_r5p_c+M_xu5p_D_c -([])-> M_g3p_c+M_s7p_c
R_tkt1_model = FluxModel();
R_tkt1_model.flux_index = 38
R_tkt1_model.flux_symbol = "R_tkt1"
R_tkt1_model.flux_constraint_type = GLPK.DB;
R_tkt1_model.flux_lower_bound = 0.0;
R_tkt1_model.flux_upper_bound = 1.0;
R_tkt1_model.flux_bounds_model = Bounds;
R_tkt1_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_tkt1_model.flux_bound_alpha = 1.0;
R_tkt1_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_tkt1"] = R_tkt1_model;
R_tkt1_model = 0;

# 39 -1*(R_tkt1: M_r5p_c+M_xu5p_D_c -([])-> M_g3p_c+M_s7p_c)
R_tkt1_reverse_model = FluxModel();
R_tkt1_reverse_model.flux_index = 39
R_tkt1_reverse_model.flux_symbol = "R_tkt1_reverse"
R_tkt1_reverse_model.flux_constraint_type = GLPK.DB;
R_tkt1_reverse_model.flux_lower_bound = 0.0;
R_tkt1_reverse_model.flux_upper_bound = 1.0;
R_tkt1_reverse_model.flux_bounds_model = Bounds;
R_tkt1_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_tkt1_reverse_model.flux_bound_alpha = 1.0;
R_tkt1_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_tkt1_reverse"] = R_tkt1_reverse_model;
R_tkt1_reverse_model = 0;

# 40 R_tkt2: M_e4p_c+M_xu5p_D_c -([])-> M_f6p_c+M_g3p_c
R_tkt2_model = FluxModel();
R_tkt2_model.flux_index = 40
R_tkt2_model.flux_symbol = "R_tkt2"
R_tkt2_model.flux_constraint_type = GLPK.DB;
R_tkt2_model.flux_lower_bound = 0.0;
R_tkt2_model.flux_upper_bound = 1.0;
R_tkt2_model.flux_bounds_model = Bounds;
R_tkt2_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_tkt2_model.flux_bound_alpha = 1.0;
R_tkt2_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_tkt2"] = R_tkt2_model;
R_tkt2_model = 0;

# 41 -1*(R_tkt2: M_e4p_c+M_xu5p_D_c -([])-> M_f6p_c+M_g3p_c)
R_tkt2_reverse_model = FluxModel();
R_tkt2_reverse_model.flux_index = 41
R_tkt2_reverse_model.flux_symbol = "R_tkt2_reverse"
R_tkt2_reverse_model.flux_constraint_type = GLPK.DB;
R_tkt2_reverse_model.flux_lower_bound = 0.0;
R_tkt2_reverse_model.flux_upper_bound = 1.0;
R_tkt2_reverse_model.flux_bounds_model = Bounds;
R_tkt2_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_tkt2_reverse_model.flux_bound_alpha = 1.0;
R_tkt2_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_tkt2_reverse"] = R_tkt2_reverse_model;
R_tkt2_reverse_model = 0;

# 42 R_edd: M_6pgc_c -([])-> M_2ddg6p_c+M_h2o_c
R_edd_model = FluxModel();
R_edd_model.flux_index = 42
R_edd_model.flux_symbol = "R_edd"
R_edd_model.flux_constraint_type = GLPK.DB;
R_edd_model.flux_lower_bound = 0.0;
R_edd_model.flux_upper_bound = 1.0;
R_edd_model.flux_bounds_model = Bounds;
R_edd_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_edd_model.flux_bound_alpha = 1.0;
R_edd_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_edd"] = R_edd_model;
R_edd_model = 0;

# 43 R_eda: M_2ddg6p_c -([])-> M_g3p_c+M_pyr_c
R_eda_model = FluxModel();
R_eda_model.flux_index = 43
R_eda_model.flux_symbol = "R_eda"
R_eda_model.flux_constraint_type = GLPK.DB;
R_eda_model.flux_lower_bound = 0.0;
R_eda_model.flux_upper_bound = 1.0;
R_eda_model.flux_bounds_model = Bounds;
R_eda_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_eda_model.flux_bound_alpha = 1.0;
R_eda_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_eda"] = R_eda_model;
R_eda_model = 0;

# 44 R_gltA: M_accoa_c+M_h2o_c+M_oaa_c -([])-> M_cit_c+M_coa_c+M_h_c
R_gltA_model = FluxModel();
R_gltA_model.flux_index = 44
R_gltA_model.flux_symbol = "R_gltA"
R_gltA_model.flux_constraint_type = GLPK.DB;
R_gltA_model.flux_lower_bound = 0.0;
R_gltA_model.flux_upper_bound = 1.0;
R_gltA_model.flux_bounds_model = Bounds;
R_gltA_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_gltA_model.flux_bound_alpha = 1.0;
R_gltA_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_gltA"] = R_gltA_model;
R_gltA_model = 0;

# 45 R_acn: M_cit_c -([])-> M_icit_c
R_acn_model = FluxModel();
R_acn_model.flux_index = 45
R_acn_model.flux_symbol = "R_acn"
R_acn_model.flux_constraint_type = GLPK.DB;
R_acn_model.flux_lower_bound = 0.0;
R_acn_model.flux_upper_bound = 1.0;
R_acn_model.flux_bounds_model = Bounds;
R_acn_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_acn_model.flux_bound_alpha = 1.0;
R_acn_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_acn"] = R_acn_model;
R_acn_model = 0;

# 46 -1*(R_acn: M_cit_c -([])-> M_icit_c)
R_acn_reverse_model = FluxModel();
R_acn_reverse_model.flux_index = 46
R_acn_reverse_model.flux_symbol = "R_acn_reverse"
R_acn_reverse_model.flux_constraint_type = GLPK.DB;
R_acn_reverse_model.flux_lower_bound = 0.0;
R_acn_reverse_model.flux_upper_bound = 1.0;
R_acn_reverse_model.flux_bounds_model = Bounds;
R_acn_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_acn_reverse_model.flux_bound_alpha = 1.0;
R_acn_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_acn_reverse"] = R_acn_reverse_model;
R_acn_reverse_model = 0;

# 47 R_icd: M_icit_c+M_nadp_c -([])-> M_akg_c+M_co2_c+M_nadph_c
R_icd_model = FluxModel();
R_icd_model.flux_index = 47
R_icd_model.flux_symbol = "R_icd"
R_icd_model.flux_constraint_type = GLPK.DB;
R_icd_model.flux_lower_bound = 0.0;
R_icd_model.flux_upper_bound = 1.0;
R_icd_model.flux_bounds_model = Bounds;
R_icd_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_icd_model.flux_bound_alpha = 1.0;
R_icd_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_icd"] = R_icd_model;
R_icd_model = 0;

# 48 -1*(R_icd: M_icit_c+M_nadp_c -([])-> M_akg_c+M_co2_c+M_nadph_c)
R_icd_reverse_model = FluxModel();
R_icd_reverse_model.flux_index = 48
R_icd_reverse_model.flux_symbol = "R_icd_reverse"
R_icd_reverse_model.flux_constraint_type = GLPK.DB;
R_icd_reverse_model.flux_lower_bound = 0.0;
R_icd_reverse_model.flux_upper_bound = 1.0;
R_icd_reverse_model.flux_bounds_model = Bounds;
R_icd_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_icd_reverse_model.flux_bound_alpha = 1.0;
R_icd_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_icd_reverse"] = R_icd_reverse_model;
R_icd_reverse_model = 0;

# 49 R_sucAB: M_akg_c+M_coa_c+M_nad_c -([])-> M_co2_c+M_nadh_c+M_succoa_c
R_sucAB_model = FluxModel();
R_sucAB_model.flux_index = 49
R_sucAB_model.flux_symbol = "R_sucAB"
R_sucAB_model.flux_constraint_type = GLPK.DB;
R_sucAB_model.flux_lower_bound = 0.0;
R_sucAB_model.flux_upper_bound = 1.0;
R_sucAB_model.flux_bounds_model = Bounds;
R_sucAB_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_sucAB_model.flux_bound_alpha = 1.0;
R_sucAB_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_sucAB"] = R_sucAB_model;
R_sucAB_model = 0;

# 50 R_sucCD: M_atp_c+M_coa_c+M_succ_c -([])-> M_adp_c+M_pi_c+M_succoa_c
R_sucCD_model = FluxModel();
R_sucCD_model.flux_index = 50
R_sucCD_model.flux_symbol = "R_sucCD"
R_sucCD_model.flux_constraint_type = GLPK.DB;
R_sucCD_model.flux_lower_bound = 0.0;
R_sucCD_model.flux_upper_bound = 1.0;
R_sucCD_model.flux_bounds_model = Bounds;
R_sucCD_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_sucCD_model.flux_bound_alpha = 1.0;
R_sucCD_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_sucCD"] = R_sucCD_model;
R_sucCD_model = 0;

# 51 -1*(R_sucCD: M_atp_c+M_coa_c+M_succ_c -([])-> M_adp_c+M_pi_c+M_succoa_c)
R_sucCD_reverse_model = FluxModel();
R_sucCD_reverse_model.flux_index = 51
R_sucCD_reverse_model.flux_symbol = "R_sucCD_reverse"
R_sucCD_reverse_model.flux_constraint_type = GLPK.DB;
R_sucCD_reverse_model.flux_lower_bound = 0.0;
R_sucCD_reverse_model.flux_upper_bound = 1.0;
R_sucCD_reverse_model.flux_bounds_model = Bounds;
R_sucCD_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_sucCD_reverse_model.flux_bound_alpha = 1.0;
R_sucCD_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_sucCD_reverse"] = R_sucCD_reverse_model;
R_sucCD_reverse_model = 0;

# 52 R_sdh: M_q8_c+M_succ_c -([])-> M_fum_c+M_q8h2_c
R_sdh_model = FluxModel();
R_sdh_model.flux_index = 52
R_sdh_model.flux_symbol = "R_sdh"
R_sdh_model.flux_constraint_type = GLPK.DB;
R_sdh_model.flux_lower_bound = 0.0;
R_sdh_model.flux_upper_bound = 1.0;
R_sdh_model.flux_bounds_model = Bounds;
R_sdh_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_sdh_model.flux_bound_alpha = 1.0;
R_sdh_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_sdh"] = R_sdh_model;
R_sdh_model = 0;

# 53 R_frd: M_fum_c+M_mql8_c -([])-> M_mqn8_c+M_succ_c
R_frd_model = FluxModel();
R_frd_model.flux_index = 53
R_frd_model.flux_symbol = "R_frd"
R_frd_model.flux_constraint_type = GLPK.DB;
R_frd_model.flux_lower_bound = 0.0;
R_frd_model.flux_upper_bound = 1.0;
R_frd_model.flux_bounds_model = Bounds;
R_frd_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_frd_model.flux_bound_alpha = 1.0;
R_frd_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_frd"] = R_frd_model;
R_frd_model = 0;

# 54 R_fum: M_fum_c+M_h2o_c -([])-> M_mal_L_c
R_fum_model = FluxModel();
R_fum_model.flux_index = 54
R_fum_model.flux_symbol = "R_fum"
R_fum_model.flux_constraint_type = GLPK.DB;
R_fum_model.flux_lower_bound = 0.0;
R_fum_model.flux_upper_bound = 1.0;
R_fum_model.flux_bounds_model = Bounds;
R_fum_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_fum_model.flux_bound_alpha = 1.0;
R_fum_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_fum"] = R_fum_model;
R_fum_model = 0;

# 55 -1*(R_fum: M_fum_c+M_h2o_c -([])-> M_mal_L_c)
R_fum_reverse_model = FluxModel();
R_fum_reverse_model.flux_index = 55
R_fum_reverse_model.flux_symbol = "R_fum_reverse"
R_fum_reverse_model.flux_constraint_type = GLPK.DB;
R_fum_reverse_model.flux_lower_bound = 0.0;
R_fum_reverse_model.flux_upper_bound = 1.0;
R_fum_reverse_model.flux_bounds_model = Bounds;
R_fum_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_fum_reverse_model.flux_bound_alpha = 1.0;
R_fum_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_fum_reverse"] = R_fum_reverse_model;
R_fum_reverse_model = 0;

# 56 R_mdh: M_mal_L_c+M_nad_c -([])-> M_oaa_c+M_h_c+M_nadh_c
R_mdh_model = FluxModel();
R_mdh_model.flux_index = 56
R_mdh_model.flux_symbol = "R_mdh"
R_mdh_model.flux_constraint_type = GLPK.DB;
R_mdh_model.flux_lower_bound = 0.0;
R_mdh_model.flux_upper_bound = 1.0;
R_mdh_model.flux_bounds_model = Bounds;
R_mdh_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_mdh_model.flux_bound_alpha = 1.0;
R_mdh_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_mdh"] = R_mdh_model;
R_mdh_model = 0;

# 57 -1*(R_mdh: M_mal_L_c+M_nad_c -([])-> M_oaa_c+M_h_c+M_nadh_c)
R_mdh_reverse_model = FluxModel();
R_mdh_reverse_model.flux_index = 57
R_mdh_reverse_model.flux_symbol = "R_mdh_reverse"
R_mdh_reverse_model.flux_constraint_type = GLPK.DB;
R_mdh_reverse_model.flux_lower_bound = 0.0;
R_mdh_reverse_model.flux_upper_bound = 1.0;
R_mdh_reverse_model.flux_bounds_model = Bounds;
R_mdh_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_mdh_reverse_model.flux_bound_alpha = 1.0;
R_mdh_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_mdh_reverse"] = R_mdh_reverse_model;
R_mdh_reverse_model = 0;

# 58 R_aceA: M_icit_c -([])-> M_glx_c+M_succ_c
R_aceA_model = FluxModel();
R_aceA_model.flux_index = 58
R_aceA_model.flux_symbol = "R_aceA"
R_aceA_model.flux_constraint_type = GLPK.DB;
R_aceA_model.flux_lower_bound = 0.0;
R_aceA_model.flux_upper_bound = 1.0;
R_aceA_model.flux_bounds_model = Bounds;
R_aceA_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_aceA_model.flux_bound_alpha = 1.0;
R_aceA_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_aceA"] = R_aceA_model;
R_aceA_model = 0;

# 59 R_aceB: M_accoa_c+M_glx_c+M_h2o_c -([])-> M_coa_c+M_h_c+M_mal_L_c
R_aceB_model = FluxModel();
R_aceB_model.flux_index = 59
R_aceB_model.flux_symbol = "R_aceB"
R_aceB_model.flux_constraint_type = GLPK.DB;
R_aceB_model.flux_lower_bound = 0.0;
R_aceB_model.flux_upper_bound = 1.0;
R_aceB_model.flux_bounds_model = Bounds;
R_aceB_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_aceB_model.flux_bound_alpha = 1.0;
R_aceB_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_aceB"] = R_aceB_model;
R_aceB_model = 0;

# 60 R_maeA: M_mal_L_c+M_nad_c -([])-> M_co2_c+M_nadh_c+M_pyr_c
R_maeA_model = FluxModel();
R_maeA_model.flux_index = 60
R_maeA_model.flux_symbol = "R_maeA"
R_maeA_model.flux_constraint_type = GLPK.DB;
R_maeA_model.flux_lower_bound = 0.0;
R_maeA_model.flux_upper_bound = 1.0;
R_maeA_model.flux_bounds_model = Bounds;
R_maeA_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_maeA_model.flux_bound_alpha = 1.0;
R_maeA_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_maeA"] = R_maeA_model;
R_maeA_model = 0;

# 61 R_maeB: M_mal_L_c+M_nadp_c -([])-> M_co2_c+M_nadph_c+M_pyr_c
R_maeB_model = FluxModel();
R_maeB_model.flux_index = 61
R_maeB_model.flux_symbol = "R_maeB"
R_maeB_model.flux_constraint_type = GLPK.DB;
R_maeB_model.flux_lower_bound = 0.0;
R_maeB_model.flux_upper_bound = 1.0;
R_maeB_model.flux_bounds_model = Bounds;
R_maeB_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_maeB_model.flux_bound_alpha = 1.0;
R_maeB_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_maeB"] = R_maeB_model;
R_maeB_model = 0;

# 62 R_pta: M_accoa_c+M_pi_c -([])-> M_actp_c+M_coa_c
R_pta_model = FluxModel();
R_pta_model.flux_index = 62
R_pta_model.flux_symbol = "R_pta"
R_pta_model.flux_constraint_type = GLPK.DB;
R_pta_model.flux_lower_bound = 0.0;
R_pta_model.flux_upper_bound = 1.0;
R_pta_model.flux_bounds_model = Bounds;
R_pta_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_pta_model.flux_bound_alpha = 1.0;
R_pta_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_pta"] = R_pta_model;
R_pta_model = 0;

# 63 -1*(R_pta: M_accoa_c+M_pi_c -([])-> M_actp_c+M_coa_c)
R_pta_reverse_model = FluxModel();
R_pta_reverse_model.flux_index = 63
R_pta_reverse_model.flux_symbol = "R_pta_reverse"
R_pta_reverse_model.flux_constraint_type = GLPK.DB;
R_pta_reverse_model.flux_lower_bound = 0.0;
R_pta_reverse_model.flux_upper_bound = 1.0;
R_pta_reverse_model.flux_bounds_model = Bounds;
R_pta_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_pta_reverse_model.flux_bound_alpha = 1.0;
R_pta_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_pta_reverse"] = R_pta_reverse_model;
R_pta_reverse_model = 0;

# 64 R_ackA: M_actp_c+M_adp_c -([])-> M_ac_c+M_atp_c
R_ackA_model = FluxModel();
R_ackA_model.flux_index = 64
R_ackA_model.flux_symbol = "R_ackA"
R_ackA_model.flux_constraint_type = GLPK.DB;
R_ackA_model.flux_lower_bound = 0.0;
R_ackA_model.flux_upper_bound = 1.0;
R_ackA_model.flux_bounds_model = Bounds;
R_ackA_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_ackA_model.flux_bound_alpha = 1.0;
R_ackA_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_ackA"] = R_ackA_model;
R_ackA_model = 0;

# 65 -1*(R_ackA: M_actp_c+M_adp_c -([])-> M_ac_c+M_atp_c)
R_ackA_reverse_model = FluxModel();
R_ackA_reverse_model.flux_index = 65
R_ackA_reverse_model.flux_symbol = "R_ackA_reverse"
R_ackA_reverse_model.flux_constraint_type = GLPK.DB;
R_ackA_reverse_model.flux_lower_bound = 0.0;
R_ackA_reverse_model.flux_upper_bound = 1.0;
R_ackA_reverse_model.flux_bounds_model = Bounds;
R_ackA_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_ackA_reverse_model.flux_bound_alpha = 1.0;
R_ackA_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_ackA_reverse"] = R_ackA_reverse_model;
R_ackA_reverse_model = 0;

# 66 R_acs: M_ac_c+M_atp_c+M_coa_c -([])-> M_accoa_c+M_amp_c+M_ppi_c
R_acs_model = FluxModel();
R_acs_model.flux_index = 66
R_acs_model.flux_symbol = "R_acs"
R_acs_model.flux_constraint_type = GLPK.DB;
R_acs_model.flux_lower_bound = 0.0;
R_acs_model.flux_upper_bound = 1.0;
R_acs_model.flux_bounds_model = Bounds;
R_acs_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_acs_model.flux_bound_alpha = 1.0;
R_acs_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_acs"] = R_acs_model;
R_acs_model = 0;

# 67 R_adhE: M_accoa_c+2*M_h_c+2*M_nadh_c -([])-> M_coa_c+M_etoh_c+2*M_nad_c
R_adhE_model = FluxModel();
R_adhE_model.flux_index = 67
R_adhE_model.flux_symbol = "R_adhE"
R_adhE_model.flux_constraint_type = GLPK.DB;
R_adhE_model.flux_lower_bound = 0.0;
R_adhE_model.flux_upper_bound = 1.0;
R_adhE_model.flux_bounds_model = Bounds;
R_adhE_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_adhE_model.flux_bound_alpha = 1.0;
R_adhE_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_adhE"] = R_adhE_model;
R_adhE_model = 0;

# 68 -1*(R_adhE: M_accoa_c+2*M_h_c+2*M_nadh_c -([])-> M_coa_c+M_etoh_c+2*M_nad_c)
R_adhE_reverse_model = FluxModel();
R_adhE_reverse_model.flux_index = 68
R_adhE_reverse_model.flux_symbol = "R_adhE_reverse"
R_adhE_reverse_model.flux_constraint_type = GLPK.DB;
R_adhE_reverse_model.flux_lower_bound = 0.0;
R_adhE_reverse_model.flux_upper_bound = 1.0;
R_adhE_reverse_model.flux_bounds_model = Bounds;
R_adhE_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_adhE_reverse_model.flux_bound_alpha = 1.0;
R_adhE_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_adhE_reverse"] = R_adhE_reverse_model;
R_adhE_reverse_model = 0;

# 69 R_ldh_f: M_pyr_c+M_nadh_c+M_h_c -([])-> M_lac_D_c+M_nad_c
R_ldh_f_model = FluxModel();
R_ldh_f_model.flux_index = 69
R_ldh_f_model.flux_symbol = "R_ldh_f"
R_ldh_f_model.flux_constraint_type = GLPK.DB;
R_ldh_f_model.flux_lower_bound = 0.0;
R_ldh_f_model.flux_upper_bound = 1.0;
R_ldh_f_model.flux_bounds_model = Bounds;
R_ldh_f_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_ldh_f_model.flux_bound_alpha = 1.0;
R_ldh_f_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_ldh_f"] = R_ldh_f_model;
R_ldh_f_model = 0;

# 70 R_ldh_r: M_lac_D_c+M_nad_c -([])-> M_pyr_c+M_nadh_c+M_h_c
R_ldh_r_model = FluxModel();
R_ldh_r_model.flux_index = 70
R_ldh_r_model.flux_symbol = "R_ldh_r"
R_ldh_r_model.flux_constraint_type = GLPK.DB;
R_ldh_r_model.flux_lower_bound = 0.0;
R_ldh_r_model.flux_upper_bound = 1.0;
R_ldh_r_model.flux_bounds_model = Bounds;
R_ldh_r_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_ldh_r_model.flux_bound_alpha = 1.0;
R_ldh_r_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_ldh_r"] = R_ldh_r_model;
R_ldh_r_model = 0;

# 71 R_pflAB: M_coa_c+M_pyr_c -([])-> M_accoa_c+M_for_c
R_pflAB_model = FluxModel();
R_pflAB_model.flux_index = 71
R_pflAB_model.flux_symbol = "R_pflAB"
R_pflAB_model.flux_constraint_type = GLPK.DB;
R_pflAB_model.flux_lower_bound = 0.0;
R_pflAB_model.flux_upper_bound = 1.0;
R_pflAB_model.flux_bounds_model = Bounds;
R_pflAB_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_pflAB_model.flux_bound_alpha = 1.0;
R_pflAB_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_pflAB"] = R_pflAB_model;
R_pflAB_model = 0;

# 72 R_cyd: 2*M_h_c+0.5*M_o2_c+M_q8h2_c -([])-> M_h2o_c+M_q8_c+2*M_h_e
R_cyd_model = FluxModel();
R_cyd_model.flux_index = 72
R_cyd_model.flux_symbol = "R_cyd"
R_cyd_model.flux_constraint_type = GLPK.DB;
R_cyd_model.flux_lower_bound = 0.0;
R_cyd_model.flux_upper_bound = 1.0;
R_cyd_model.flux_bounds_model = Bounds;
R_cyd_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_cyd_model.flux_bound_alpha = 1.0;
R_cyd_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_cyd"] = R_cyd_model;
R_cyd_model = 0;

# 73 R_app: 2*M_h_c+M_mql8_c+0.5*M_o2_c -([])-> M_h2o_c+M_mqn8_c+2*M_h_e
R_app_model = FluxModel();
R_app_model.flux_index = 73
R_app_model.flux_symbol = "R_app"
R_app_model.flux_constraint_type = GLPK.DB;
R_app_model.flux_lower_bound = 0.0;
R_app_model.flux_upper_bound = 1.0;
R_app_model.flux_bounds_model = Bounds;
R_app_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_app_model.flux_bound_alpha = 1.0;
R_app_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_app"] = R_app_model;
R_app_model = 0;

# 74 R_atp: M_adp_c+M_pi_c+4*M_h_e -([])-> M_atp_c+3*M_h_c+M_h2o_c
R_atp_model = FluxModel();
R_atp_model.flux_index = 74
R_atp_model.flux_symbol = "R_atp"
R_atp_model.flux_constraint_type = GLPK.DB;
R_atp_model.flux_lower_bound = 0.0;
R_atp_model.flux_upper_bound = 1.0;
R_atp_model.flux_bounds_model = Bounds;
R_atp_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_atp_model.flux_bound_alpha = 1.0;
R_atp_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_atp"] = R_atp_model;
R_atp_model = 0;

# 75 R_nuo: 3*M_h_c+M_nadh_c+M_q8_c -([])-> M_nad_c+M_q8h2_c+2*M_h_e
R_nuo_model = FluxModel();
R_nuo_model.flux_index = 75
R_nuo_model.flux_symbol = "R_nuo"
R_nuo_model.flux_constraint_type = GLPK.DB;
R_nuo_model.flux_lower_bound = 0.0;
R_nuo_model.flux_upper_bound = 1.0;
R_nuo_model.flux_bounds_model = Bounds;
R_nuo_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_nuo_model.flux_bound_alpha = 1.0;
R_nuo_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_nuo"] = R_nuo_model;
R_nuo_model = 0;

# 76 R_pnt1: M_nad_c+M_nadph_c -([])-> M_nadh_c+M_nadp_c
R_pnt1_model = FluxModel();
R_pnt1_model.flux_index = 76
R_pnt1_model.flux_symbol = "R_pnt1"
R_pnt1_model.flux_constraint_type = GLPK.DB;
R_pnt1_model.flux_lower_bound = 0.0;
R_pnt1_model.flux_upper_bound = 1.0;
R_pnt1_model.flux_bounds_model = Bounds;
R_pnt1_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_pnt1_model.flux_bound_alpha = 1.0;
R_pnt1_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_pnt1"] = R_pnt1_model;
R_pnt1_model = 0;

# 77 R_pnt2: M_nadh_c+M_nadp_c+2*M_h_e -([])-> 2*M_h_c+M_nad_c+M_nadph_c
R_pnt2_model = FluxModel();
R_pnt2_model.flux_index = 77
R_pnt2_model.flux_symbol = "R_pnt2"
R_pnt2_model.flux_constraint_type = GLPK.DB;
R_pnt2_model.flux_lower_bound = 0.0;
R_pnt2_model.flux_upper_bound = 1.0;
R_pnt2_model.flux_bounds_model = Bounds;
R_pnt2_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_pnt2_model.flux_bound_alpha = 1.0;
R_pnt2_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_pnt2"] = R_pnt2_model;
R_pnt2_model = 0;

# 78 R_ndh1: M_h_c+M_nadh_c+M_q8_c -([])-> M_nad_c+M_q8h2_c
R_ndh1_model = FluxModel();
R_ndh1_model.flux_index = 78
R_ndh1_model.flux_symbol = "R_ndh1"
R_ndh1_model.flux_constraint_type = GLPK.DB;
R_ndh1_model.flux_lower_bound = 0.0;
R_ndh1_model.flux_upper_bound = 1.0;
R_ndh1_model.flux_bounds_model = Bounds;
R_ndh1_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_ndh1_model.flux_bound_alpha = 1.0;
R_ndh1_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_ndh1"] = R_ndh1_model;
R_ndh1_model = 0;

# 79 R_ndh2: M_h_c+M_mqn8_c+M_nadh_c -([])-> M_mql8_c+M_nad_c
R_ndh2_model = FluxModel();
R_ndh2_model.flux_index = 79
R_ndh2_model.flux_symbol = "R_ndh2"
R_ndh2_model.flux_constraint_type = GLPK.DB;
R_ndh2_model.flux_lower_bound = 0.0;
R_ndh2_model.flux_upper_bound = 1.0;
R_ndh2_model.flux_bounds_model = Bounds;
R_ndh2_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_ndh2_model.flux_bound_alpha = 1.0;
R_ndh2_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_ndh2"] = R_ndh2_model;
R_ndh2_model = 0;

# 80 R_hack1: M_atp_c+M_h2o_c -([])-> M_adp_c+M_h_c+M_pi_c
R_hack1_model = FluxModel();
R_hack1_model.flux_index = 80
R_hack1_model.flux_symbol = "R_hack1"
R_hack1_model.flux_constraint_type = GLPK.DB;
R_hack1_model.flux_lower_bound = 0.0;
R_hack1_model.flux_upper_bound = 1.0;
R_hack1_model.flux_bounds_model = Bounds;
R_hack1_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_hack1_model.flux_bound_alpha = 1.0;
R_hack1_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_hack1"] = R_hack1_model;
R_hack1_model = 0;

# 81 R_ppk: M_atp_c+M_pi_c -([])-> M_adp_c+M_ppi_c
R_ppk_model = FluxModel();
R_ppk_model.flux_index = 81
R_ppk_model.flux_symbol = "R_ppk"
R_ppk_model.flux_constraint_type = GLPK.DB;
R_ppk_model.flux_lower_bound = 0.0;
R_ppk_model.flux_upper_bound = 1.0;
R_ppk_model.flux_bounds_model = Bounds;
R_ppk_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_ppk_model.flux_bound_alpha = 1.0;
R_ppk_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_ppk"] = R_ppk_model;
R_ppk_model = 0;

# 82 R_ppa: M_ppi_c+M_h2o_c -([])-> 2*M_pi_c+M_h_c
R_ppa_model = FluxModel();
R_ppa_model.flux_index = 82
R_ppa_model.flux_symbol = "R_ppa"
R_ppa_model.flux_constraint_type = GLPK.DB;
R_ppa_model.flux_lower_bound = 0.0;
R_ppa_model.flux_upper_bound = 1.0;
R_ppa_model.flux_bounds_model = Bounds;
R_ppa_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_ppa_model.flux_bound_alpha = 1.0;
R_ppa_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_ppa"] = R_ppa_model;
R_ppa_model = 0;

# 83 R_chor: M_e4p_c+2*M_pep_c+M_nadph_c+M_atp_c -([])-> M_chor_c+M_nadp_c+M_adp_c+4*M_pi_c
R_chor_model = FluxModel();
R_chor_model.flux_index = 83
R_chor_model.flux_symbol = "R_chor"
R_chor_model.flux_constraint_type = GLPK.DB;
R_chor_model.flux_lower_bound = 0.0;
R_chor_model.flux_upper_bound = 1.0;
R_chor_model.flux_bounds_model = Bounds;
R_chor_model.flux_gamma_array = vec([0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_chor_model.flux_bound_alpha = 1.0;
R_chor_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_chor"] = R_chor_model;
R_chor_model = 0;

# 84 R_gar: M_r5p_c+M_gln_L_c+M_gly_L_c+2*M_atp_c+M_h2o_c -([])-> M_gar_c+M_glu_L_c+M_adp_c+M_amp_c+M_pi_c+M_ppi_c+7*M_h_c
R_gar_model = FluxModel();
R_gar_model.flux_index = 84
R_gar_model.flux_symbol = "R_gar"
R_gar_model.flux_constraint_type = GLPK.DB;
R_gar_model.flux_lower_bound = 0.0;
R_gar_model.flux_upper_bound = 1.0;
R_gar_model.flux_bounds_model = Bounds;
R_gar_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_gar_model.flux_bound_alpha = 1.0;
R_gar_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_gar"] = R_gar_model;
R_gar_model = 0;

# 85 R_air: M_gar_c+M_10fthf_c+M_gln_L_c+2*M_atp_c+M_h2o_c -([])-> M_air_c+M_thf_c+M_glu_L_c+2*M_adp_c+2*M_pi_c+3*M_h_c
R_air_model = FluxModel();
R_air_model.flux_index = 85
R_air_model.flux_symbol = "R_air"
R_air_model.flux_constraint_type = GLPK.DB;
R_air_model.flux_lower_bound = 0.0;
R_air_model.flux_upper_bound = 1.0;
R_air_model.flux_bounds_model = Bounds;
R_air_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_air_model.flux_bound_alpha = 1.0;
R_air_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_air"] = R_air_model;
R_air_model = 0;

# 86 R_aicar: M_air_c+M_asp_L_c+2*M_atp_c+M_hco3_c -([])-> M_aicar_c+M_fum_c+2*M_adp_c+2*M_h_c+2*M_pi_c
R_aicar_model = FluxModel();
R_aicar_model.flux_index = 86
R_aicar_model.flux_symbol = "R_aicar"
R_aicar_model.flux_constraint_type = GLPK.DB;
R_aicar_model.flux_lower_bound = 0.0;
R_aicar_model.flux_upper_bound = 1.0;
R_aicar_model.flux_bounds_model = Bounds;
R_aicar_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_aicar_model.flux_bound_alpha = 1.0;
R_aicar_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_aicar"] = R_aicar_model;
R_aicar_model = 0;

# 87 R_imp: M_aicar_c+M_10fthf_c -([])-> M_imp_c+M_thf_c+M_h2o_c
R_imp_model = FluxModel();
R_imp_model.flux_index = 87
R_imp_model.flux_symbol = "R_imp"
R_imp_model.flux_constraint_type = GLPK.DB;
R_imp_model.flux_lower_bound = 0.0;
R_imp_model.flux_upper_bound = 1.0;
R_imp_model.flux_bounds_model = Bounds;
R_imp_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_imp_model.flux_bound_alpha = 1.0;
R_imp_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_imp"] = R_imp_model;
R_imp_model = 0;

# 88 R_mthfc: M_h2o_c+M_methf_c -([])-> M_10fthf_c
R_mthfc_model = FluxModel();
R_mthfc_model.flux_index = 88
R_mthfc_model.flux_symbol = "R_mthfc"
R_mthfc_model.flux_constraint_type = GLPK.DB;
R_mthfc_model.flux_lower_bound = 0.0;
R_mthfc_model.flux_upper_bound = 1.0;
R_mthfc_model.flux_bounds_model = Bounds;
R_mthfc_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_mthfc_model.flux_bound_alpha = 1.0;
R_mthfc_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_mthfc"] = R_mthfc_model;
R_mthfc_model = 0;

# 89 -1*(R_mthfc: M_h2o_c+M_methf_c -([])-> M_10fthf_c)
R_mthfc_reverse_model = FluxModel();
R_mthfc_reverse_model.flux_index = 89
R_mthfc_reverse_model.flux_symbol = "R_mthfc_reverse"
R_mthfc_reverse_model.flux_constraint_type = GLPK.DB;
R_mthfc_reverse_model.flux_lower_bound = 0.0;
R_mthfc_reverse_model.flux_upper_bound = 1.0;
R_mthfc_reverse_model.flux_bounds_model = Bounds;
R_mthfc_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_mthfc_reverse_model.flux_bound_alpha = 1.0;
R_mthfc_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_mthfc_reverse"] = R_mthfc_reverse_model;
R_mthfc_reverse_model = 0;

# 90 R_mthfd: M_mlthf_c+M_nadp_c -([])-> M_h_c+M_methf_c+M_nadph_c
R_mthfd_model = FluxModel();
R_mthfd_model.flux_index = 90
R_mthfd_model.flux_symbol = "R_mthfd"
R_mthfd_model.flux_constraint_type = GLPK.DB;
R_mthfd_model.flux_lower_bound = 0.0;
R_mthfd_model.flux_upper_bound = 1.0;
R_mthfd_model.flux_bounds_model = Bounds;
R_mthfd_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_mthfd_model.flux_bound_alpha = 1.0;
R_mthfd_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_mthfd"] = R_mthfd_model;
R_mthfd_model = 0;

# 91 -1*(R_mthfd: M_mlthf_c+M_nadp_c -([])-> M_h_c+M_methf_c+M_nadph_c)
R_mthfd_reverse_model = FluxModel();
R_mthfd_reverse_model.flux_index = 91
R_mthfd_reverse_model.flux_symbol = "R_mthfd_reverse"
R_mthfd_reverse_model.flux_constraint_type = GLPK.DB;
R_mthfd_reverse_model.flux_lower_bound = 0.0;
R_mthfd_reverse_model.flux_upper_bound = 1.0;
R_mthfd_reverse_model.flux_bounds_model = Bounds;
R_mthfd_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_mthfd_reverse_model.flux_bound_alpha = 1.0;
R_mthfd_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_mthfd_reverse"] = R_mthfd_reverse_model;
R_mthfd_reverse_model = 0;

# 92 R_mthfr2: M_mlthf_c+M_h_c+M_nadh_c -([])-> M_5mthf_c+M_nad_c
R_mthfr2_model = FluxModel();
R_mthfr2_model.flux_index = 92
R_mthfr2_model.flux_symbol = "R_mthfr2"
R_mthfr2_model.flux_constraint_type = GLPK.DB;
R_mthfr2_model.flux_lower_bound = 0.0;
R_mthfr2_model.flux_upper_bound = 1.0;
R_mthfr2_model.flux_bounds_model = Bounds;
R_mthfr2_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_mthfr2_model.flux_bound_alpha = 1.0;
R_mthfr2_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_mthfr2"] = R_mthfr2_model;
R_mthfr2_model = 0;

# 93 R_gmp: M_imp_c+M_atp_c+M_gln_L_c+M_nad_c+2*M_h2o_c -([])-> M_gmp_c+M_amp_c+M_glu_L_c+M_nadh_c+3*M_h_c+M_ppi_c
R_gmp_model = FluxModel();
R_gmp_model.flux_index = 93
R_gmp_model.flux_symbol = "R_gmp"
R_gmp_model.flux_constraint_type = GLPK.DB;
R_gmp_model.flux_lower_bound = 0.0;
R_gmp_model.flux_upper_bound = 1.0;
R_gmp_model.flux_bounds_model = Bounds;
R_gmp_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_gmp_model.flux_bound_alpha = 1.0;
R_gmp_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_gmp"] = R_gmp_model;
R_gmp_model = 0;

# 94 R_gdp: M_gmp_c+M_atp_c -([])-> M_gdp_c+M_adp_c
R_gdp_model = FluxModel();
R_gdp_model.flux_index = 94
R_gdp_model.flux_symbol = "R_gdp"
R_gdp_model.flux_constraint_type = GLPK.DB;
R_gdp_model.flux_lower_bound = 0.0;
R_gdp_model.flux_upper_bound = 1.0;
R_gdp_model.flux_bounds_model = Bounds;
R_gdp_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_gdp_model.flux_bound_alpha = 1.0;
R_gdp_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_gdp"] = R_gdp_model;
R_gdp_model = 0;

# 95 R_gtp: M_gdp_c+M_atp_c -([])-> M_gtp_c+M_adp_c
R_gtp_model = FluxModel();
R_gtp_model.flux_index = 95
R_gtp_model.flux_symbol = "R_gtp"
R_gtp_model.flux_constraint_type = GLPK.DB;
R_gtp_model.flux_lower_bound = 0.0;
R_gtp_model.flux_upper_bound = 1.0;
R_gtp_model.flux_bounds_model = Bounds;
R_gtp_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_gtp_model.flux_bound_alpha = 1.0;
R_gtp_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_gtp"] = R_gtp_model;
R_gtp_model = 0;

# 96 R_amp: M_asp_L_c+M_imp_c+M_gtp_c -([])-> M_amp_c+M_gdp_c+M_pi_c+2*M_h_c+M_fum_c
R_amp_model = FluxModel();
R_amp_model.flux_index = 96
R_amp_model.flux_symbol = "R_amp"
R_amp_model.flux_constraint_type = GLPK.DB;
R_amp_model.flux_lower_bound = 0.0;
R_amp_model.flux_upper_bound = 1.0;
R_amp_model.flux_bounds_model = Bounds;
R_amp_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_amp_model.flux_bound_alpha = 1.0;
R_amp_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_amp"] = R_amp_model;
R_amp_model = 0;

# 97 R_adk: M_amp_c+M_atp_c -([])-> 2*M_adp_c
R_adk_model = FluxModel();
R_adk_model.flux_index = 97
R_adk_model.flux_symbol = "R_adk"
R_adk_model.flux_constraint_type = GLPK.DB;
R_adk_model.flux_lower_bound = 0.0;
R_adk_model.flux_upper_bound = 1.0;
R_adk_model.flux_bounds_model = Bounds;
R_adk_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_adk_model.flux_bound_alpha = 1.0;
R_adk_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_adk"] = R_adk_model;
R_adk_model = 0;

# 98 -1*(R_adk: M_amp_c+M_atp_c -([])-> 2*M_adp_c)
R_adk_reverse_model = FluxModel();
R_adk_reverse_model.flux_index = 98
R_adk_reverse_model.flux_symbol = "R_adk_reverse"
R_adk_reverse_model.flux_constraint_type = GLPK.DB;
R_adk_reverse_model.flux_lower_bound = 0.0;
R_adk_reverse_model.flux_upper_bound = 1.0;
R_adk_reverse_model.flux_bounds_model = Bounds;
R_adk_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_adk_reverse_model.flux_bound_alpha = 1.0;
R_adk_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_adk_reverse"] = R_adk_reverse_model;
R_adk_reverse_model = 0;

# 99 R_ump: M_gln_L_c+M_asp_L_c+M_r5p_c+M_q8_c+3*M_atp_c+M_hco3_c -([])-> M_ump_c+M_glu_L_c+M_q8h2_c+2*M_h_c+2*M_adp_c+M_amp_c+2*M_pi_c+M_ppi_c+M_co2_c
R_ump_model = FluxModel();
R_ump_model.flux_index = 99
R_ump_model.flux_symbol = "R_ump"
R_ump_model.flux_constraint_type = GLPK.DB;
R_ump_model.flux_lower_bound = 0.0;
R_ump_model.flux_upper_bound = 1.0;
R_ump_model.flux_bounds_model = Bounds;
R_ump_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_ump_model.flux_bound_alpha = 1.0;
R_ump_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_ump"] = R_ump_model;
R_ump_model = 0;

# 100 R_udp: M_ump_c+M_atp_c -([])-> M_udp_c+M_adp_c
R_udp_model = FluxModel();
R_udp_model.flux_index = 100
R_udp_model.flux_symbol = "R_udp"
R_udp_model.flux_constraint_type = GLPK.DB;
R_udp_model.flux_lower_bound = 0.0;
R_udp_model.flux_upper_bound = 1.0;
R_udp_model.flux_bounds_model = Bounds;
R_udp_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_udp_model.flux_bound_alpha = 1.0;
R_udp_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_udp"] = R_udp_model;
R_udp_model = 0;

# 101 R_utp: M_udp_c+M_atp_c -([])-> M_utp_c+M_adp_c
R_utp_model = FluxModel();
R_utp_model.flux_index = 101
R_utp_model.flux_symbol = "R_utp"
R_utp_model.flux_constraint_type = GLPK.DB;
R_utp_model.flux_lower_bound = 0.0;
R_utp_model.flux_upper_bound = 1.0;
R_utp_model.flux_bounds_model = Bounds;
R_utp_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_utp_model.flux_bound_alpha = 1.0;
R_utp_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_utp"] = R_utp_model;
R_utp_model = 0;

# 102 R_ctp: M_utp_c+M_gln_L_c+M_atp_c+M_h2o_c -([])-> M_ctp_c+M_glu_L_c+M_adp_c+M_pi_c+3*M_h_c
R_ctp_model = FluxModel();
R_ctp_model.flux_index = 102
R_ctp_model.flux_symbol = "R_ctp"
R_ctp_model.flux_constraint_type = GLPK.DB;
R_ctp_model.flux_lower_bound = 0.0;
R_ctp_model.flux_upper_bound = 1.0;
R_ctp_model.flux_bounds_model = Bounds;
R_ctp_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_ctp_model.flux_bound_alpha = 1.0;
R_ctp_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_ctp"] = R_ctp_model;
R_ctp_model = 0;

# 103 R_cdp: M_ctp_c+M_h2o_c -([])-> M_cdp_c+M_pi_c+M_h_c
R_cdp_model = FluxModel();
R_cdp_model.flux_index = 103
R_cdp_model.flux_symbol = "R_cdp"
R_cdp_model.flux_constraint_type = GLPK.DB;
R_cdp_model.flux_lower_bound = 0.0;
R_cdp_model.flux_upper_bound = 1.0;
R_cdp_model.flux_bounds_model = Bounds;
R_cdp_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_cdp_model.flux_bound_alpha = 1.0;
R_cdp_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_cdp"] = R_cdp_model;
R_cdp_model = 0;

# 104 R_cmp: M_cdp_c+M_h2o_c -([])-> M_cmp_c+M_pi_c+M_h_c
R_cmp_model = FluxModel();
R_cmp_model.flux_index = 104
R_cmp_model.flux_symbol = "R_cmp"
R_cmp_model.flux_constraint_type = GLPK.DB;
R_cmp_model.flux_lower_bound = 0.0;
R_cmp_model.flux_upper_bound = 1.0;
R_cmp_model.flux_bounds_model = Bounds;
R_cmp_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_cmp_model.flux_bound_alpha = 1.0;
R_cmp_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_cmp"] = R_cmp_model;
R_cmp_model = 0;

# 105 R_alaAC: M_pyr_c+M_glu_L_c -([])-> M_ala_L_c+M_akg_c
R_alaAC_model = FluxModel();
R_alaAC_model.flux_index = 105
R_alaAC_model.flux_symbol = "R_alaAC"
R_alaAC_model.flux_constraint_type = GLPK.DB;
R_alaAC_model.flux_lower_bound = 0.0;
R_alaAC_model.flux_upper_bound = 1.0;
R_alaAC_model.flux_bounds_model = Bounds;
R_alaAC_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_alaAC_model.flux_bound_alpha = 1.0;
R_alaAC_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_alaAC"] = R_alaAC_model;
R_alaAC_model = 0;

# 106 R_arg: M_glu_L_c+M_accoa_c+4*M_atp_c+M_nadph_c+3*M_h2o_c+M_gln_L_c+M_asp_L_c+M_co2_c -([])-> M_arg_L_c+M_coa_c+5*M_h_c+3*M_adp_c+3*M_pi_c+M_nadp_c+M_akg_c+M_ac_c+M_amp_c+M_ppi_c+M_fum_c
R_arg_model = FluxModel();
R_arg_model.flux_index = 106
R_arg_model.flux_symbol = "R_arg"
R_arg_model.flux_constraint_type = GLPK.DB;
R_arg_model.flux_lower_bound = 0.0;
R_arg_model.flux_upper_bound = 1.0;
R_arg_model.flux_bounds_model = Bounds;
R_arg_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_arg_model.flux_bound_alpha = 1.0;
R_arg_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_arg"] = R_arg_model;
R_arg_model = 0;

# 107 R_aspA: M_fum_c+M_nh4_c -([])-> M_asp_L_c
R_aspA_model = FluxModel();
R_aspA_model.flux_index = 107
R_aspA_model.flux_symbol = "R_aspA"
R_aspA_model.flux_constraint_type = GLPK.DB;
R_aspA_model.flux_lower_bound = 0.0;
R_aspA_model.flux_upper_bound = 1.0;
R_aspA_model.flux_bounds_model = Bounds;
R_aspA_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_aspA_model.flux_bound_alpha = 1.0;
R_aspA_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_aspA"] = R_aspA_model;
R_aspA_model = 0;

# 108 -1*(R_aspA: M_fum_c+M_nh4_c -([])-> M_asp_L_c)
R_aspA_reverse_model = FluxModel();
R_aspA_reverse_model.flux_index = 108
R_aspA_reverse_model.flux_symbol = "R_aspA_reverse"
R_aspA_reverse_model.flux_constraint_type = GLPK.DB;
R_aspA_reverse_model.flux_lower_bound = 0.0;
R_aspA_reverse_model.flux_upper_bound = 1.0;
R_aspA_reverse_model.flux_bounds_model = Bounds;
R_aspA_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_aspA_reverse_model.flux_bound_alpha = 1.0;
R_aspA_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_aspA_reverse"] = R_aspA_reverse_model;
R_aspA_reverse_model = 0;

# 109 R_aspC: M_glu_L_c+M_oaa_c -([])-> M_asp_L_c+M_akg_c
R_aspC_model = FluxModel();
R_aspC_model.flux_index = 109
R_aspC_model.flux_symbol = "R_aspC"
R_aspC_model.flux_constraint_type = GLPK.DB;
R_aspC_model.flux_lower_bound = 0.0;
R_aspC_model.flux_upper_bound = 1.0;
R_aspC_model.flux_bounds_model = Bounds;
R_aspC_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_aspC_model.flux_bound_alpha = 1.0;
R_aspC_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_aspC"] = R_aspC_model;
R_aspC_model = 0;

# 110 R_asnB: M_asp_L_c+M_gln_L_c+M_h2o_c+M_atp_c -([])-> M_asn_L_c+M_glu_L_c+M_h_c+M_ppi_c+M_amp_c
R_asnB_model = FluxModel();
R_asnB_model.flux_index = 110
R_asnB_model.flux_symbol = "R_asnB"
R_asnB_model.flux_constraint_type = GLPK.DB;
R_asnB_model.flux_lower_bound = 0.0;
R_asnB_model.flux_upper_bound = 1.0;
R_asnB_model.flux_bounds_model = Bounds;
R_asnB_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_asnB_model.flux_bound_alpha = 1.0;
R_asnB_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_asnB"] = R_asnB_model;
R_asnB_model = 0;

# 111 R_asnA: M_asp_L_c+M_atp_c+M_nh4_c -([])-> M_asn_L_c+M_h_c+M_ppi_c+M_amp_c
R_asnA_model = FluxModel();
R_asnA_model.flux_index = 111
R_asnA_model.flux_symbol = "R_asnA"
R_asnA_model.flux_constraint_type = GLPK.DB;
R_asnA_model.flux_lower_bound = 0.0;
R_asnA_model.flux_upper_bound = 1.0;
R_asnA_model.flux_bounds_model = Bounds;
R_asnA_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_asnA_model.flux_bound_alpha = 1.0;
R_asnA_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_asnA"] = R_asnA_model;
R_asnA_model = 0;

# 112 R_cysEMK: M_ser_L_c+M_accoa_c+M_h2s_c -([])-> M_cys_L_c+M_coa_c+M_h_c+M_ac_c
R_cysEMK_model = FluxModel();
R_cysEMK_model.flux_index = 112
R_cysEMK_model.flux_symbol = "R_cysEMK"
R_cysEMK_model.flux_constraint_type = GLPK.DB;
R_cysEMK_model.flux_lower_bound = 0.0;
R_cysEMK_model.flux_upper_bound = 1.0;
R_cysEMK_model.flux_bounds_model = Bounds;
R_cysEMK_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_cysEMK_model.flux_bound_alpha = 1.0;
R_cysEMK_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_cysEMK"] = R_cysEMK_model;
R_cysEMK_model = 0;

# 113 R_gltBD: M_gln_L_c+M_akg_c+M_nadph_c+M_h_c -([])-> 2*M_glu_L_c+M_nadp_c
R_gltBD_model = FluxModel();
R_gltBD_model.flux_index = 113
R_gltBD_model.flux_symbol = "R_gltBD"
R_gltBD_model.flux_constraint_type = GLPK.DB;
R_gltBD_model.flux_lower_bound = 0.0;
R_gltBD_model.flux_upper_bound = 1.0;
R_gltBD_model.flux_bounds_model = Bounds;
R_gltBD_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_gltBD_model.flux_bound_alpha = 1.0;
R_gltBD_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_gltBD"] = R_gltBD_model;
R_gltBD_model = 0;

# 114 R_gdhA: M_akg_c+M_nadph_c+M_nh4_c+M_h_c -([])-> M_glu_L_c+M_h2o_c+M_nadp_c
R_gdhA_model = FluxModel();
R_gdhA_model.flux_index = 114
R_gdhA_model.flux_symbol = "R_gdhA"
R_gdhA_model.flux_constraint_type = GLPK.DB;
R_gdhA_model.flux_lower_bound = 0.0;
R_gdhA_model.flux_upper_bound = 1.0;
R_gdhA_model.flux_bounds_model = Bounds;
R_gdhA_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_gdhA_model.flux_bound_alpha = 1.0;
R_gdhA_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_gdhA"] = R_gdhA_model;
R_gdhA_model = 0;

# 115 R_glnA: M_glu_L_c+M_atp_c+M_nh4_c -([])-> M_gln_L_c+M_h_c+M_adp_c+M_pi_c
R_glnA_model = FluxModel();
R_glnA_model.flux_index = 115
R_glnA_model.flux_symbol = "R_glnA"
R_glnA_model.flux_constraint_type = GLPK.DB;
R_glnA_model.flux_lower_bound = 0.0;
R_glnA_model.flux_upper_bound = 1.0;
R_glnA_model.flux_bounds_model = Bounds;
R_glnA_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_glnA_model.flux_bound_alpha = 1.0;
R_glnA_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_glnA"] = R_glnA_model;
R_glnA_model = 0;

# 116 R_glyA: M_ser_L_c+M_thf_c -([])-> M_gly_L_c+M_h2o_c+M_mlthf_c
R_glyA_model = FluxModel();
R_glyA_model.flux_index = 116
R_glyA_model.flux_symbol = "R_glyA"
R_glyA_model.flux_constraint_type = GLPK.DB;
R_glyA_model.flux_lower_bound = 0.0;
R_glyA_model.flux_upper_bound = 1.0;
R_glyA_model.flux_bounds_model = Bounds;
R_glyA_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_glyA_model.flux_bound_alpha = 1.0;
R_glyA_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_glyA"] = R_glyA_model;
R_glyA_model = 0;

# 117 R_his: M_gln_L_c+M_r5p_c+3*M_atp_c+2*M_nad_c+3*M_h2o_c -([])-> M_his_L_c+M_akg_c+M_aicar_c+2*M_adp_c+2*M_nadh_c+M_pi_c+2*M_ppi_c+6*M_h_c
R_his_model = FluxModel();
R_his_model.flux_index = 117
R_his_model.flux_symbol = "R_his"
R_his_model.flux_constraint_type = GLPK.DB;
R_his_model.flux_lower_bound = 0.0;
R_his_model.flux_upper_bound = 1.0;
R_his_model.flux_bounds_model = Bounds;
R_his_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_his_model.flux_bound_alpha = 1.0;
R_his_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_his"] = R_his_model;
R_his_model = 0;

# 118 R_ile: M_thr_L_c+2*M_h_c+M_pyr_c+M_nadph_c+M_glu_L_c -([])-> M_ile_L_c+M_h2o_c+M_nh4_c+M_co2_c+M_nadp_c+M_akg_c
R_ile_model = FluxModel();
R_ile_model.flux_index = 118
R_ile_model.flux_symbol = "R_ile"
R_ile_model.flux_constraint_type = GLPK.DB;
R_ile_model.flux_lower_bound = 0.0;
R_ile_model.flux_upper_bound = 1.0;
R_ile_model.flux_bounds_model = Bounds;
R_ile_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_ile_model.flux_bound_alpha = 1.0;
R_ile_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_ile"] = R_ile_model;
R_ile_model = 0;

# 119 R_leu: 2*M_pyr_c+M_glu_L_c+M_h_c+M_nad_c+M_nadph_c+M_accoa_c -([])-> M_leu_L_c+2*M_co2_c+M_nadp_c+M_coa_c+M_nadh_c+M_akg_c
R_leu_model = FluxModel();
R_leu_model.flux_index = 119
R_leu_model.flux_symbol = "R_leu"
R_leu_model.flux_constraint_type = GLPK.DB;
R_leu_model.flux_lower_bound = 0.0;
R_leu_model.flux_upper_bound = 1.0;
R_leu_model.flux_bounds_model = Bounds;
R_leu_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_leu_model.flux_bound_alpha = 1.0;
R_leu_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_leu"] = R_leu_model;
R_leu_model = 0;

# 120 R_lys: M_asp_L_c+M_atp_c+2*M_nadph_c+2*M_h_c+M_pyr_c+M_succoa_c+M_glu_L_c -([])-> M_lys_L_c+M_adp_c+M_pi_c+2*M_nadp_c+M_coa_c+M_akg_c+M_succ_c+M_co2_c
R_lys_model = FluxModel();
R_lys_model.flux_index = 120
R_lys_model.flux_symbol = "R_lys"
R_lys_model.flux_constraint_type = GLPK.DB;
R_lys_model.flux_lower_bound = 0.0;
R_lys_model.flux_upper_bound = 1.0;
R_lys_model.flux_bounds_model = Bounds;
R_lys_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_lys_model.flux_bound_alpha = 1.0;
R_lys_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_lys"] = R_lys_model;
R_lys_model = 0;

# 121 R_met: M_asp_L_c+M_cys_L_c+M_succoa_c+M_atp_c+2*M_nadph_c+M_5mthf_c -([])-> M_met_L_c+M_coa_c+M_succ_c+M_adp_c+M_pi_c+2*M_nadp_c+M_thf_c+M_nh4_c+M_pyr_c
R_met_model = FluxModel();
R_met_model.flux_index = 121
R_met_model.flux_symbol = "R_met"
R_met_model.flux_constraint_type = GLPK.DB;
R_met_model.flux_lower_bound = 0.0;
R_met_model.flux_upper_bound = 1.0;
R_met_model.flux_bounds_model = Bounds;
R_met_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_met_model.flux_bound_alpha = 1.0;
R_met_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_met"] = R_met_model;
R_met_model = 0;

# 122 R_phe: M_chor_c+M_h_c+M_glu_L_c -([])-> M_phe_L_c+M_co2_c+M_h2o_c+M_akg_c
R_phe_model = FluxModel();
R_phe_model.flux_index = 122
R_phe_model.flux_symbol = "R_phe"
R_phe_model.flux_constraint_type = GLPK.DB;
R_phe_model.flux_lower_bound = 0.0;
R_phe_model.flux_upper_bound = 1.0;
R_phe_model.flux_bounds_model = Bounds;
R_phe_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_phe_model.flux_bound_alpha = 1.0;
R_phe_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_phe"] = R_phe_model;
R_phe_model = 0;

# 123 R_pro: M_glu_L_c+M_atp_c+2*M_h_c+2*M_nadph_c -([])-> M_pro_L_c+M_adp_c+2*M_nadp_c+M_pi_c+M_h2o_c
R_pro_model = FluxModel();
R_pro_model.flux_index = 123
R_pro_model.flux_symbol = "R_pro"
R_pro_model.flux_constraint_type = GLPK.DB;
R_pro_model.flux_lower_bound = 0.0;
R_pro_model.flux_upper_bound = 1.0;
R_pro_model.flux_bounds_model = Bounds;
R_pro_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_pro_model.flux_bound_alpha = 1.0;
R_pro_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_pro"] = R_pro_model;
R_pro_model = 0;

# 124 R_serABC: M_3pg_c+M_nad_c+M_glu_L_c+M_h2o_c -([])-> M_ser_L_c+M_nadh_c+M_h_c+M_akg_c+M_pi_c
R_serABC_model = FluxModel();
R_serABC_model.flux_index = 124
R_serABC_model.flux_symbol = "R_serABC"
R_serABC_model.flux_constraint_type = GLPK.DB;
R_serABC_model.flux_lower_bound = 0.0;
R_serABC_model.flux_upper_bound = 1.0;
R_serABC_model.flux_bounds_model = Bounds;
R_serABC_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_serABC_model.flux_bound_alpha = 1.0;
R_serABC_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_serABC"] = R_serABC_model;
R_serABC_model = 0;

# 125 R_thr: M_asp_L_c+2*M_atp_c+2*M_nadph_c+M_h_c+M_h2o_c -([])-> M_thr_L_c+2*M_adp_c+2*M_pi_c+2*M_nadp_c
R_thr_model = FluxModel();
R_thr_model.flux_index = 125
R_thr_model.flux_symbol = "R_thr"
R_thr_model.flux_constraint_type = GLPK.DB;
R_thr_model.flux_lower_bound = 0.0;
R_thr_model.flux_upper_bound = 1.0;
R_thr_model.flux_bounds_model = Bounds;
R_thr_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 1.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_thr_model.flux_bound_alpha = 1.0;
R_thr_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_thr"] = R_thr_model;
R_thr_model = 0;

# 126 R_trp: M_chor_c+M_gln_L_c+M_ser_L_c+M_r5p_c+2*M_atp_c -([])-> M_trp_L_c+M_glu_L_c+M_pyr_c+M_ppi_c+2*M_h2o_c+M_co2_c+M_g3p_c+2*M_adp_c+M_h_c
R_trp_model = FluxModel();
R_trp_model.flux_index = 126
R_trp_model.flux_symbol = "R_trp"
R_trp_model.flux_constraint_type = GLPK.DB;
R_trp_model.flux_lower_bound = 0.0;
R_trp_model.flux_upper_bound = 1.0;
R_trp_model.flux_bounds_model = Bounds;
R_trp_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_trp_model.flux_bound_alpha = 1.0;
R_trp_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_trp"] = R_trp_model;
R_trp_model = 0;

# 127 R_tyr: M_chor_c+M_glu_L_c+M_nad_c -([])-> M_tyr_L_c+M_akg_c+M_nadh_c+M_co2_c
R_tyr_model = FluxModel();
R_tyr_model.flux_index = 127
R_tyr_model.flux_symbol = "R_tyr"
R_tyr_model.flux_constraint_type = GLPK.DB;
R_tyr_model.flux_lower_bound = 0.0;
R_tyr_model.flux_upper_bound = 1.0;
R_tyr_model.flux_bounds_model = Bounds;
R_tyr_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_tyr_model.flux_bound_alpha = 1.0;
R_tyr_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_tyr"] = R_tyr_model;
R_tyr_model = 0;

# 128 R_val: 2*M_pyr_c+2*M_h_c+M_nadph_c+M_glu_L_c -([])-> M_val_L_c+M_co2_c+M_nadp_c+M_h2o_c+M_akg_c
R_val_model = FluxModel();
R_val_model.flux_index = 128
R_val_model.flux_symbol = "R_val"
R_val_model.flux_constraint_type = GLPK.DB;
R_val_model.flux_lower_bound = 0.0;
R_val_model.flux_upper_bound = 1.0;
R_val_model.flux_bounds_model = Bounds;
R_val_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_val_model.flux_bound_alpha = 1.0;
R_val_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_val"] = R_val_model;
R_val_model = 0;

# 129 R_arg_deg1: M_arg_L_c+5.0*M_h2o_c+M_atp_c+M_o2_c+2.0*M_nad_c+M_akg_c -([])-> 4.0*M_h_c+M_co2_c+M_urea_c+M_glu_L_c+M_pi_c+M_adp_c+M_nh4_c+M_h2o2_c+2.0*M_nadh_c+M_succ_c
R_arg_deg1_model = FluxModel();
R_arg_deg1_model.flux_index = 129
R_arg_deg1_model.flux_symbol = "R_arg_deg1"
R_arg_deg1_model.flux_constraint_type = GLPK.DB;
R_arg_deg1_model.flux_lower_bound = 0.0;
R_arg_deg1_model.flux_upper_bound = 1.0;
R_arg_deg1_model.flux_bounds_model = Bounds;
R_arg_deg1_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_arg_deg1_model.flux_bound_alpha = 1.0;
R_arg_deg1_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_arg_deg1"] = R_arg_deg1_model;
R_arg_deg1_model = 0;

# 130 R_arg_deg2: M_arg_L_c+5.0*M_h2o_c+M_atp_c+M_o2_c+M_nad_c+M_nadp_c+M_akg_c -([])-> 4.0*M_h_c+M_co2_c+M_urea_c+M_glu_L_c+M_pi_c+M_adp_c+M_nh4_c+M_h2o2_c+M_nadh_c+M_nadph_c+M_succ_c
R_arg_deg2_model = FluxModel();
R_arg_deg2_model.flux_index = 130
R_arg_deg2_model.flux_symbol = "R_arg_deg2"
R_arg_deg2_model.flux_constraint_type = GLPK.DB;
R_arg_deg2_model.flux_lower_bound = 0.0;
R_arg_deg2_model.flux_upper_bound = 1.0;
R_arg_deg2_model.flux_bounds_model = Bounds;
R_arg_deg2_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_arg_deg2_model.flux_bound_alpha = 1.0;
R_arg_deg2_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_arg_deg2"] = R_arg_deg2_model;
R_arg_deg2_model = 0;

# 131 R_arg_deg3: M_arg_L_c+5.0*M_h2o_c+M_atp_c+M_o2_c+2.0*M_nadp_c+M_akg_c -([])-> 4.0*M_h_c+M_co2_c+M_urea_c+M_glu_L_c+M_pi_c+M_adp_c+M_nh4_c+M_h2o2_c+2.0*M_nadph_c+M_succ_c
R_arg_deg3_model = FluxModel();
R_arg_deg3_model.flux_index = 131
R_arg_deg3_model.flux_symbol = "R_arg_deg3"
R_arg_deg3_model.flux_constraint_type = GLPK.DB;
R_arg_deg3_model.flux_lower_bound = 0.0;
R_arg_deg3_model.flux_upper_bound = 1.0;
R_arg_deg3_model.flux_bounds_model = Bounds;
R_arg_deg3_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_arg_deg3_model.flux_bound_alpha = 1.0;
R_arg_deg3_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_arg_deg3"] = R_arg_deg3_model;
R_arg_deg3_model = 0;

# 132 R_arg_deg4: M_arg_L_c+3.0*M_h2o_c+2.0*M_akg_c+2.0*M_nad_c -([])-> 3.0*M_h_c+M_co2_c+M_urea_c+2.0*M_glu_L_c+2.0*M_nadh_c+M_succ_c
R_arg_deg4_model = FluxModel();
R_arg_deg4_model.flux_index = 132
R_arg_deg4_model.flux_symbol = "R_arg_deg4"
R_arg_deg4_model.flux_constraint_type = GLPK.DB;
R_arg_deg4_model.flux_lower_bound = 0.0;
R_arg_deg4_model.flux_upper_bound = 1.0;
R_arg_deg4_model.flux_bounds_model = Bounds;
R_arg_deg4_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_arg_deg4_model.flux_bound_alpha = 1.0;
R_arg_deg4_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_arg_deg4"] = R_arg_deg4_model;
R_arg_deg4_model = 0;

# 133 R_arg_deg5: M_arg_L_c+3.0*M_h2o_c+2.0*M_akg_c+M_nad_c+M_nadp_c -([])-> 3.0*M_h_c+M_co2_c+M_urea_c+2.0*M_glu_L_c+M_nadh_c+M_nadph_c+M_succ_c
R_arg_deg5_model = FluxModel();
R_arg_deg5_model.flux_index = 133
R_arg_deg5_model.flux_symbol = "R_arg_deg5"
R_arg_deg5_model.flux_constraint_type = GLPK.DB;
R_arg_deg5_model.flux_lower_bound = 0.0;
R_arg_deg5_model.flux_upper_bound = 1.0;
R_arg_deg5_model.flux_bounds_model = Bounds;
R_arg_deg5_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_arg_deg5_model.flux_bound_alpha = 1.0;
R_arg_deg5_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_arg_deg5"] = R_arg_deg5_model;
R_arg_deg5_model = 0;

# 134 R_arg_deg6: M_arg_L_c+M_accoa_c+4.0*M_h2o_c+M_akg_c+M_nad_c -([])-> M_coa_c+M_h_c+M_co2_c+2.0*M_nh4_c+2.0*M_glu_L_c+M_nadh_c+M_succ_c
R_arg_deg6_model = FluxModel();
R_arg_deg6_model.flux_index = 134
R_arg_deg6_model.flux_symbol = "R_arg_deg6"
R_arg_deg6_model.flux_constraint_type = GLPK.DB;
R_arg_deg6_model.flux_lower_bound = 0.0;
R_arg_deg6_model.flux_upper_bound = 1.0;
R_arg_deg6_model.flux_bounds_model = Bounds;
R_arg_deg6_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_arg_deg6_model.flux_bound_alpha = 1.0;
R_arg_deg6_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_arg_deg6"] = R_arg_deg6_model;
R_arg_deg6_model = 0;

# 135 R_thr_deg1: M_thr_L_c+M_nad_c+M_coa_c -([])-> M_nadh_c+M_h_c+M_accoa_c+M_gly_L_c
R_thr_deg1_model = FluxModel();
R_thr_deg1_model.flux_index = 135
R_thr_deg1_model.flux_symbol = "R_thr_deg1"
R_thr_deg1_model.flux_constraint_type = GLPK.DB;
R_thr_deg1_model.flux_lower_bound = 0.0;
R_thr_deg1_model.flux_upper_bound = 1.0;
R_thr_deg1_model.flux_bounds_model = Bounds;
R_thr_deg1_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_thr_deg1_model.flux_bound_alpha = 1.0;
R_thr_deg1_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_thr_deg1"] = R_thr_deg1_model;
R_thr_deg1_model = 0;

# 136 R_gly_deg: M_gly_L_c+M_accoa_c+M_h_c+M_o2_c+M_h2o_c -([])-> M_coa_c+M_co2_c+M_h2o2_c+M_nh4_c+M_mglx_c
R_gly_deg_model = FluxModel();
R_gly_deg_model.flux_index = 136
R_gly_deg_model.flux_symbol = "R_gly_deg"
R_gly_deg_model.flux_constraint_type = GLPK.DB;
R_gly_deg_model.flux_lower_bound = 0.0;
R_gly_deg_model.flux_upper_bound = 1.0;
R_gly_deg_model.flux_bounds_model = Bounds;
R_gly_deg_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_gly_deg_model.flux_bound_alpha = 1.0;
R_gly_deg_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_gly_deg"] = R_gly_deg_model;
R_gly_deg_model = 0;

# 137 R_thr_deg2: M_thr_L_c+M_nad_c+M_o2_c+M_h2o_c -([])-> M_nadh_c+M_co2_c+M_h2o2_c+M_nh4_c+M_mglx_c
R_thr_deg2_model = FluxModel();
R_thr_deg2_model.flux_index = 137
R_thr_deg2_model.flux_symbol = "R_thr_deg2"
R_thr_deg2_model.flux_constraint_type = GLPK.DB;
R_thr_deg2_model.flux_lower_bound = 0.0;
R_thr_deg2_model.flux_upper_bound = 1.0;
R_thr_deg2_model.flux_bounds_model = Bounds;
R_thr_deg2_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_thr_deg2_model.flux_bound_alpha = 1.0;
R_thr_deg2_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_thr_deg2"] = R_thr_deg2_model;
R_thr_deg2_model = 0;

# 138 R_thr_deg3: M_thr_L_c+M_coa_c+M_nad_c -([])-> M_gly_L_c+M_accoa_c+M_nadh_c+M_h_c
R_thr_deg3_model = FluxModel();
R_thr_deg3_model.flux_index = 138
R_thr_deg3_model.flux_symbol = "R_thr_deg3"
R_thr_deg3_model.flux_constraint_type = GLPK.DB;
R_thr_deg3_model.flux_lower_bound = 0.0;
R_thr_deg3_model.flux_upper_bound = 1.0;
R_thr_deg3_model.flux_bounds_model = Bounds;
R_thr_deg3_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_thr_deg3_model.flux_bound_alpha = 1.0;
R_thr_deg3_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_thr_deg3"] = R_thr_deg3_model;
R_thr_deg3_model = 0;

# 139 R_thr_deg4: M_thr_L_c+M_pi_c+M_adp_c -([])-> M_h_c+M_h2o_c+M_for_c+M_atp_c+M_prop_c
R_thr_deg4_model = FluxModel();
R_thr_deg4_model.flux_index = 139
R_thr_deg4_model.flux_symbol = "R_thr_deg4"
R_thr_deg4_model.flux_constraint_type = GLPK.DB;
R_thr_deg4_model.flux_lower_bound = 0.0;
R_thr_deg4_model.flux_upper_bound = 1.0;
R_thr_deg4_model.flux_bounds_model = Bounds;
R_thr_deg4_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_thr_deg4_model.flux_bound_alpha = 1.0;
R_thr_deg4_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_thr_deg4"] = R_thr_deg4_model;
R_thr_deg4_model = 0;

# 140 R_thr_deg5: M_thr_L_c+M_h_c+M_pyr_c+M_nadph_c+M_glu_L_c -([])-> 2.0*M_h2o_c+M_co2_c+M_nadp_c+M_akg_c+M_ile_L_c
R_thr_deg5_model = FluxModel();
R_thr_deg5_model.flux_index = 140
R_thr_deg5_model.flux_symbol = "R_thr_deg5"
R_thr_deg5_model.flux_constraint_type = GLPK.DB;
R_thr_deg5_model.flux_lower_bound = 0.0;
R_thr_deg5_model.flux_upper_bound = 1.0;
R_thr_deg5_model.flux_bounds_model = Bounds;
R_thr_deg5_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_thr_deg5_model.flux_bound_alpha = 1.0;
R_thr_deg5_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_thr_deg5"] = R_thr_deg5_model;
R_thr_deg5_model = 0;

# 141 R_mglx_deg: M_mglx_c+M_nadp_c+M_h2o_c -([])-> M_pyr_c+M_nadph_c+M_h_c
R_mglx_deg_model = FluxModel();
R_mglx_deg_model.flux_index = 141
R_mglx_deg_model.flux_symbol = "R_mglx_deg"
R_mglx_deg_model.flux_constraint_type = GLPK.DB;
R_mglx_deg_model.flux_lower_bound = 0.0;
R_mglx_deg_model.flux_upper_bound = 1.0;
R_mglx_deg_model.flux_bounds_model = Bounds;
R_mglx_deg_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_mglx_deg_model.flux_bound_alpha = 1.0;
R_mglx_deg_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_mglx_deg"] = R_mglx_deg_model;
R_mglx_deg_model = 0;

# 142 -1*(R_mglx_deg: M_mglx_c+M_nadp_c+M_h2o_c -([])-> M_pyr_c+M_nadph_c+M_h_c)
R_mglx_deg_reverse_model = FluxModel();
R_mglx_deg_reverse_model.flux_index = 142
R_mglx_deg_reverse_model.flux_symbol = "R_mglx_deg_reverse"
R_mglx_deg_reverse_model.flux_constraint_type = GLPK.DB;
R_mglx_deg_reverse_model.flux_lower_bound = 0.0;
R_mglx_deg_reverse_model.flux_upper_bound = 1.0;
R_mglx_deg_reverse_model.flux_bounds_model = Bounds;
R_mglx_deg_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_mglx_deg_reverse_model.flux_bound_alpha = 1.0;
R_mglx_deg_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_mglx_deg_reverse"] = R_mglx_deg_reverse_model;
R_mglx_deg_reverse_model = 0;

# 143 R_ser_deg: M_ser_L_c -([])-> M_nh4_c+M_pyr_c
R_ser_deg_model = FluxModel();
R_ser_deg_model.flux_index = 143
R_ser_deg_model.flux_symbol = "R_ser_deg"
R_ser_deg_model.flux_constraint_type = GLPK.DB;
R_ser_deg_model.flux_lower_bound = 0.0;
R_ser_deg_model.flux_upper_bound = 1.0;
R_ser_deg_model.flux_bounds_model = Bounds;
R_ser_deg_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_ser_deg_model.flux_bound_alpha = 1.0;
R_ser_deg_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_ser_deg"] = R_ser_deg_model;
R_ser_deg_model = 0;

# 144 R_pro_deg: M_pro_L_c+M_q8_c+2.0*M_h2o_c+M_nad_c -([])-> 2.0*M_h_c+M_q8h2_c+M_nadh_c+M_glu_L_c
R_pro_deg_model = FluxModel();
R_pro_deg_model.flux_index = 144
R_pro_deg_model.flux_symbol = "R_pro_deg"
R_pro_deg_model.flux_constraint_type = GLPK.DB;
R_pro_deg_model.flux_lower_bound = 0.0;
R_pro_deg_model.flux_upper_bound = 1.0;
R_pro_deg_model.flux_bounds_model = Bounds;
R_pro_deg_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_pro_deg_model.flux_bound_alpha = 1.0;
R_pro_deg_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_pro_deg"] = R_pro_deg_model;
R_pro_deg_model = 0;

# 145 R_trp_deg: M_trp_L_c+M_h2o_c -([])-> M_indole_c+M_nh4_c+M_pyr_c
R_trp_deg_model = FluxModel();
R_trp_deg_model.flux_index = 145
R_trp_deg_model.flux_symbol = "R_trp_deg"
R_trp_deg_model.flux_constraint_type = GLPK.DB;
R_trp_deg_model.flux_lower_bound = 0.0;
R_trp_deg_model.flux_upper_bound = 1.0;
R_trp_deg_model.flux_bounds_model = Bounds;
R_trp_deg_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_trp_deg_model.flux_bound_alpha = 1.0;
R_trp_deg_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_trp_deg"] = R_trp_deg_model;
R_trp_deg_model = 0;

# 146 R_cys_deg: M_cys_L_c+M_h2o_c -([])-> M_h2s_c+M_nh4_c+M_pyr_c
R_cys_deg_model = FluxModel();
R_cys_deg_model.flux_index = 146
R_cys_deg_model.flux_symbol = "R_cys_deg"
R_cys_deg_model.flux_constraint_type = GLPK.DB;
R_cys_deg_model.flux_lower_bound = 0.0;
R_cys_deg_model.flux_upper_bound = 1.0;
R_cys_deg_model.flux_bounds_model = Bounds;
R_cys_deg_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_cys_deg_model.flux_bound_alpha = 1.0;
R_cys_deg_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_cys_deg"] = R_cys_deg_model;
R_cys_deg_model = 0;

# 147 R_ala_deg: M_ala_L_c+M_h2o_c+M_q8_c -([])-> M_q8h2_c+M_nh4_c+M_pyr_c
R_ala_deg_model = FluxModel();
R_ala_deg_model.flux_index = 147
R_ala_deg_model.flux_symbol = "R_ala_deg"
R_ala_deg_model.flux_constraint_type = GLPK.DB;
R_ala_deg_model.flux_lower_bound = 0.0;
R_ala_deg_model.flux_upper_bound = 1.0;
R_ala_deg_model.flux_bounds_model = Bounds;
R_ala_deg_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_ala_deg_model.flux_bound_alpha = 1.0;
R_ala_deg_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_ala_deg"] = R_ala_deg_model;
R_ala_deg_model = 0;

# 148 R_lys_deg: M_lys_L_c -([])-> M_co2_c+M_cadav_c
R_lys_deg_model = FluxModel();
R_lys_deg_model.flux_index = 148
R_lys_deg_model.flux_symbol = "R_lys_deg"
R_lys_deg_model.flux_constraint_type = GLPK.DB;
R_lys_deg_model.flux_lower_bound = 0.0;
R_lys_deg_model.flux_upper_bound = 1.0;
R_lys_deg_model.flux_bounds_model = Bounds;
R_lys_deg_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_lys_deg_model.flux_bound_alpha = 1.0;
R_lys_deg_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_lys_deg"] = R_lys_deg_model;
R_lys_deg_model = 0;

# 149 R_gln_deg: M_gln_L_c+M_h2o_c -([])-> M_nh4_c+M_glu_L_c
R_gln_deg_model = FluxModel();
R_gln_deg_model.flux_index = 149
R_gln_deg_model.flux_symbol = "R_gln_deg"
R_gln_deg_model.flux_constraint_type = GLPK.DB;
R_gln_deg_model.flux_lower_bound = 0.0;
R_gln_deg_model.flux_upper_bound = 1.0;
R_gln_deg_model.flux_bounds_model = Bounds;
R_gln_deg_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_gln_deg_model.flux_bound_alpha = 1.0;
R_gln_deg_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_gln_deg"] = R_gln_deg_model;
R_gln_deg_model = 0;

# 150 R_glu_deg: M_glu_L_c+M_h_c -([])-> M_co2_c+M_gaba_c
R_glu_deg_model = FluxModel();
R_glu_deg_model.flux_index = 150
R_glu_deg_model.flux_symbol = "R_glu_deg"
R_glu_deg_model.flux_constraint_type = GLPK.DB;
R_glu_deg_model.flux_lower_bound = 0.0;
R_glu_deg_model.flux_upper_bound = 1.0;
R_glu_deg_model.flux_bounds_model = Bounds;
R_glu_deg_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_glu_deg_model.flux_bound_alpha = 1.0;
R_glu_deg_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_glu_deg"] = R_glu_deg_model;
R_glu_deg_model = 0;

# 151 R_gaba_deg1: M_gaba_c+M_akg_c+M_h2o_c+M_nad_c -([])-> M_succ_c+M_glu_L_c+2*M_h_c+M_nadh_c
R_gaba_deg1_model = FluxModel();
R_gaba_deg1_model.flux_index = 151
R_gaba_deg1_model.flux_symbol = "R_gaba_deg1"
R_gaba_deg1_model.flux_constraint_type = GLPK.DB;
R_gaba_deg1_model.flux_lower_bound = 0.0;
R_gaba_deg1_model.flux_upper_bound = 1.0;
R_gaba_deg1_model.flux_bounds_model = Bounds;
R_gaba_deg1_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_gaba_deg1_model.flux_bound_alpha = 1.0;
R_gaba_deg1_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_gaba_deg1"] = R_gaba_deg1_model;
R_gaba_deg1_model = 0;

# 152 R_gaba_deg2: M_gaba_c+M_akg_c+M_h2o_c+M_nadp_c -([])-> M_succ_c+M_glu_L_c+2*M_h_c+M_nadph_c
R_gaba_deg2_model = FluxModel();
R_gaba_deg2_model.flux_index = 152
R_gaba_deg2_model.flux_symbol = "R_gaba_deg2"
R_gaba_deg2_model.flux_constraint_type = GLPK.DB;
R_gaba_deg2_model.flux_lower_bound = 0.0;
R_gaba_deg2_model.flux_upper_bound = 1.0;
R_gaba_deg2_model.flux_bounds_model = Bounds;
R_gaba_deg2_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_gaba_deg2_model.flux_bound_alpha = 1.0;
R_gaba_deg2_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_gaba_deg2"] = R_gaba_deg2_model;
R_gaba_deg2_model = 0;

# 153 R_asn_deg: M_asn_L_c+M_h2o_c+M_adp_c+M_pi_c -([])-> M_nh4_c+M_asp_L_c+M_atp_c
R_asn_deg_model = FluxModel();
R_asn_deg_model.flux_index = 153
R_asn_deg_model.flux_symbol = "R_asn_deg"
R_asn_deg_model.flux_constraint_type = GLPK.DB;
R_asn_deg_model.flux_lower_bound = 0.0;
R_asn_deg_model.flux_upper_bound = 1.0;
R_asn_deg_model.flux_bounds_model = Bounds;
R_asn_deg_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_asn_deg_model.flux_bound_alpha = 1.0;
R_asn_deg_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_asn_deg"] = R_asn_deg_model;
R_asn_deg_model = 0;

# 154 R_amp_ppi: M_amp_c+M_ppi_c+4*M_h_c -([])-> M_atp_c+M_h2o_c
R_amp_ppi_model = FluxModel();
R_amp_ppi_model.flux_index = 154
R_amp_ppi_model.flux_symbol = "R_amp_ppi"
R_amp_ppi_model.flux_constraint_type = GLPK.DB;
R_amp_ppi_model.flux_lower_bound = 0.0;
R_amp_ppi_model.flux_upper_bound = 1.0;
R_amp_ppi_model.flux_bounds_model = Bounds;
R_amp_ppi_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_amp_ppi_model.flux_bound_alpha = 1.0;
R_amp_ppi_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_amp_ppi"] = R_amp_ppi_model;
R_amp_ppi_model = 0;

# 155 R_amp_pi: M_amp_c+2*M_pi_c+6*M_h_c -([])-> M_atp_c+2*M_h2o_c
R_amp_pi_model = FluxModel();
R_amp_pi_model.flux_index = 155
R_amp_pi_model.flux_symbol = "R_amp_pi"
R_amp_pi_model.flux_constraint_type = GLPK.DB;
R_amp_pi_model.flux_lower_bound = 0.0;
R_amp_pi_model.flux_upper_bound = 1.0;
R_amp_pi_model.flux_bounds_model = Bounds;
R_amp_pi_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_amp_pi_model.flux_bound_alpha = 1.0;
R_amp_pi_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_amp_pi"] = R_amp_pi_model;
R_amp_pi_model = 0;

# 156 R_gmp_ppi: M_gmp_c+M_ppi_c+4*M_h_c -([])-> M_gtp_c+M_h2o_c
R_gmp_ppi_model = FluxModel();
R_gmp_ppi_model.flux_index = 156
R_gmp_ppi_model.flux_symbol = "R_gmp_ppi"
R_gmp_ppi_model.flux_constraint_type = GLPK.DB;
R_gmp_ppi_model.flux_lower_bound = 0.0;
R_gmp_ppi_model.flux_upper_bound = 1.0;
R_gmp_ppi_model.flux_bounds_model = Bounds;
R_gmp_ppi_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_gmp_ppi_model.flux_bound_alpha = 1.0;
R_gmp_ppi_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_gmp_ppi"] = R_gmp_ppi_model;
R_gmp_ppi_model = 0;

# 157 R_gmp_pi: M_gmp_c+2*M_pi_c+6*M_h_c -([])-> M_gtp_c+2*M_h2o_c
R_gmp_pi_model = FluxModel();
R_gmp_pi_model.flux_index = 157
R_gmp_pi_model.flux_symbol = "R_gmp_pi"
R_gmp_pi_model.flux_constraint_type = GLPK.DB;
R_gmp_pi_model.flux_lower_bound = 0.0;
R_gmp_pi_model.flux_upper_bound = 1.0;
R_gmp_pi_model.flux_bounds_model = Bounds;
R_gmp_pi_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_gmp_pi_model.flux_bound_alpha = 1.0;
R_gmp_pi_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_gmp_pi"] = R_gmp_pi_model;
R_gmp_pi_model = 0;

# 158 R_cmp_ppi: M_cmp_c+M_ppi_c+4*M_h_c -([])-> M_ctp_c+M_h2o_c
R_cmp_ppi_model = FluxModel();
R_cmp_ppi_model.flux_index = 158
R_cmp_ppi_model.flux_symbol = "R_cmp_ppi"
R_cmp_ppi_model.flux_constraint_type = GLPK.DB;
R_cmp_ppi_model.flux_lower_bound = 0.0;
R_cmp_ppi_model.flux_upper_bound = 1.0;
R_cmp_ppi_model.flux_bounds_model = Bounds;
R_cmp_ppi_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_cmp_ppi_model.flux_bound_alpha = 1.0;
R_cmp_ppi_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_cmp_ppi"] = R_cmp_ppi_model;
R_cmp_ppi_model = 0;

# 159 R_cmp_pi: M_cmp_c+2*M_pi_c+6*M_h_c -([])-> M_ctp_c+2*M_h2o_c
R_cmp_pi_model = FluxModel();
R_cmp_pi_model.flux_index = 159
R_cmp_pi_model.flux_symbol = "R_cmp_pi"
R_cmp_pi_model.flux_constraint_type = GLPK.DB;
R_cmp_pi_model.flux_lower_bound = 0.0;
R_cmp_pi_model.flux_upper_bound = 1.0;
R_cmp_pi_model.flux_bounds_model = Bounds;
R_cmp_pi_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_cmp_pi_model.flux_bound_alpha = 1.0;
R_cmp_pi_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_cmp_pi"] = R_cmp_pi_model;
R_cmp_pi_model = 0;

# 160 R_ump_ppi: M_ump_c+M_ppi_c+4*M_h_c -([])-> M_utp_c+M_h2o_c
R_ump_ppi_model = FluxModel();
R_ump_ppi_model.flux_index = 160
R_ump_ppi_model.flux_symbol = "R_ump_ppi"
R_ump_ppi_model.flux_constraint_type = GLPK.DB;
R_ump_ppi_model.flux_lower_bound = 0.0;
R_ump_ppi_model.flux_upper_bound = 1.0;
R_ump_ppi_model.flux_bounds_model = Bounds;
R_ump_ppi_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_ump_ppi_model.flux_bound_alpha = 1.0;
R_ump_ppi_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_ump_ppi"] = R_ump_ppi_model;
R_ump_ppi_model = 0;

# 161 R_ump_pi: M_ump_c+2*M_pi_c+6*M_h_c -([])-> M_utp_c+2*M_h2o_c
R_ump_pi_model = FluxModel();
R_ump_pi_model.flux_index = 161
R_ump_pi_model.flux_symbol = "R_ump_pi"
R_ump_pi_model.flux_constraint_type = GLPK.DB;
R_ump_pi_model.flux_lower_bound = 0.0;
R_ump_pi_model.flux_upper_bound = 1.0;
R_ump_pi_model.flux_bounds_model = Bounds;
R_ump_pi_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
R_ump_pi_model.flux_bound_alpha = 1.0;
R_ump_pi_model.flux_obj_coeff = 0.0;
flux_model_dictionary["R_ump_pi"] = R_ump_pi_model;
R_ump_pi_model = 0;

# 162 transcriptional_initiation_deGFP: GENE_deGFP+RNAP -([])-> OPEN_GENE_deGFP
transcriptional_initiation_deGFP_model = FluxModel();
transcriptional_initiation_deGFP_model.flux_index = 162
transcriptional_initiation_deGFP_model.flux_symbol = "transcriptional_initiation_deGFP"
transcriptional_initiation_deGFP_model.flux_constraint_type = GLPK.DB;
transcriptional_initiation_deGFP_model.flux_lower_bound = 0.0;
transcriptional_initiation_deGFP_model.flux_upper_bound = 1.0;
transcriptional_initiation_deGFP_model.flux_bounds_model = Bounds;
transcriptional_initiation_deGFP_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
transcriptional_initiation_deGFP_model.flux_bound_alpha = 1.0;
transcriptional_initiation_deGFP_model.flux_obj_coeff = 0.0;
flux_model_dictionary["transcriptional_initiation_deGFP"] = transcriptional_initiation_deGFP_model;
transcriptional_initiation_deGFP_model = 0;

# 163 transcription_deGFP: OPEN_GENE_deGFP+183*M_gtp_c+231*M_ctp_c+101*M_utp_c+163*M_atp_c -([])-> mRNA_deGFP+GENE_deGFP+RNAP+1356*M_pi_c
transcription_deGFP_model = FluxModel();
transcription_deGFP_model.flux_index = 163
transcription_deGFP_model.flux_symbol = "transcription_deGFP"
transcription_deGFP_model.flux_constraint_type = GLPK.DB;
transcription_deGFP_model.flux_lower_bound = 0.0;
transcription_deGFP_model.flux_upper_bound = 1.0;
transcription_deGFP_model.flux_bounds_model = Bounds;
transcription_deGFP_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
transcription_deGFP_model.flux_bound_alpha = 1.0;
transcription_deGFP_model.flux_obj_coeff = 0.0;
flux_model_dictionary["transcription_deGFP"] = transcription_deGFP_model;
transcription_deGFP_model = 0;

# 164 mRNA_degradation_deGFP: mRNA_deGFP -([])-> 183*M_gmp_c+231*M_cmp_c+101*M_ump_c+163*M_amp_c
mRNA_degradation_deGFP_model = FluxModel();
mRNA_degradation_deGFP_model.flux_index = 164
mRNA_degradation_deGFP_model.flux_symbol = "mRNA_degradation_deGFP"
mRNA_degradation_deGFP_model.flux_constraint_type = GLPK.DB;
mRNA_degradation_deGFP_model.flux_lower_bound = 0.0;
mRNA_degradation_deGFP_model.flux_upper_bound = 1.0;
mRNA_degradation_deGFP_model.flux_bounds_model = Bounds;
mRNA_degradation_deGFP_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
mRNA_degradation_deGFP_model.flux_bound_alpha = 1.0;
mRNA_degradation_deGFP_model.flux_obj_coeff = 0.0;
flux_model_dictionary["mRNA_degradation_deGFP"] = mRNA_degradation_deGFP_model;
mRNA_degradation_deGFP_model = 0;

# 165 translation_initiation_deGFP: mRNA_deGFP+RIBOSOME -([])-> RIBOSOME_START_deGFP
translation_initiation_deGFP_model = FluxModel();
translation_initiation_deGFP_model.flux_index = 165
translation_initiation_deGFP_model.flux_symbol = "translation_initiation_deGFP"
translation_initiation_deGFP_model.flux_constraint_type = GLPK.DB;
translation_initiation_deGFP_model.flux_lower_bound = 0.0;
translation_initiation_deGFP_model.flux_upper_bound = 1.0;
translation_initiation_deGFP_model.flux_bounds_model = Bounds;
translation_initiation_deGFP_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
translation_initiation_deGFP_model.flux_bound_alpha = 1.0;
translation_initiation_deGFP_model.flux_obj_coeff = 0.0;
flux_model_dictionary["translation_initiation_deGFP"] = translation_initiation_deGFP_model;
translation_initiation_deGFP_model = 0;

# 166 translation_deGFP: RIBOSOME_START_deGFP+450*M_gtp_c+8.0*M_ala_L_c_tRNA+6.0*M_arg_L_c_tRNA+13.0*M_asn_L_c_tRNA+17.0*M_asp_L_c_tRNA+2.0*M_cys_L_c_tRNA+14.0*M_glu_L_c_tRNA+8.0*M_gln_L_c_tRNA+20.0*M_gly_L_c_tRNA+9.0*M_his_L_c_tRNA+12.0*M_ile_L_c_tRNA+19.0*M_leu_L_c_tRNA+18.0*M_lys_L_c_tRNA+5.0*M_met_L_c_tRNA+12.0*M_phe_L_c_tRNA+10.0*M_pro_L_c_tRNA+9.0*M_ser_L_c_tRNA+15.0*M_thr_L_c_tRNA+1.0*M_trp_L_c_tRNA+10.0*M_tyr_L_c_tRNA+17.0*M_val_L_c_tRNA -([])-> RIBOSOME+mRNA_deGFP+PROTEIN_deGFP+450*M_gdp_c+450*M_pi_c+225*tRNA
translation_deGFP_model = FluxModel();
translation_deGFP_model.flux_index = 166
translation_deGFP_model.flux_symbol = "translation_deGFP"
translation_deGFP_model.flux_constraint_type = GLPK.DB;
translation_deGFP_model.flux_lower_bound = 0.0;
translation_deGFP_model.flux_upper_bound = 1.0;
translation_deGFP_model.flux_bounds_model = Bounds;
translation_deGFP_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 0.0 0.0]);
translation_deGFP_model.flux_bound_alpha = 1.0;
translation_deGFP_model.flux_obj_coeff = 0.0;
flux_model_dictionary["translation_deGFP"] = translation_deGFP_model;
translation_deGFP_model = 0;

# 167 tRNA_charging_M_ala_L_c_deGFP: 8.0*M_ala_L_c+8.0*M_atp_c+8.0*tRNA -([])-> 8.0*M_ala_L_c_tRNA+8.0*M_amp_c+16.0*M_pi_c
tRNA_charging_M_ala_L_c_deGFP_model = FluxModel();
tRNA_charging_M_ala_L_c_deGFP_model.flux_index = 167
tRNA_charging_M_ala_L_c_deGFP_model.flux_symbol = "tRNA_charging_M_ala_L_c_deGFP"
tRNA_charging_M_ala_L_c_deGFP_model.flux_constraint_type = GLPK.DB;
tRNA_charging_M_ala_L_c_deGFP_model.flux_lower_bound = 0.0;
tRNA_charging_M_ala_L_c_deGFP_model.flux_upper_bound = 1.0;
tRNA_charging_M_ala_L_c_deGFP_model.flux_bounds_model = Bounds;
tRNA_charging_M_ala_L_c_deGFP_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0]);
tRNA_charging_M_ala_L_c_deGFP_model.flux_bound_alpha = 1.0;
tRNA_charging_M_ala_L_c_deGFP_model.flux_obj_coeff = 0.0;
flux_model_dictionary["tRNA_charging_M_ala_L_c_deGFP"] = tRNA_charging_M_ala_L_c_deGFP_model;
tRNA_charging_M_ala_L_c_deGFP_model = 0;

# 168 tRNA_charging_M_arg_L_c_deGFP: 6.0*M_arg_L_c+6.0*M_atp_c+6.0*tRNA -([])-> 6.0*M_arg_L_c_tRNA+6.0*M_amp_c+12.0*M_pi_c
tRNA_charging_M_arg_L_c_deGFP_model = FluxModel();
tRNA_charging_M_arg_L_c_deGFP_model.flux_index = 168
tRNA_charging_M_arg_L_c_deGFP_model.flux_symbol = "tRNA_charging_M_arg_L_c_deGFP"
tRNA_charging_M_arg_L_c_deGFP_model.flux_constraint_type = GLPK.DB;
tRNA_charging_M_arg_L_c_deGFP_model.flux_lower_bound = 0.0;
tRNA_charging_M_arg_L_c_deGFP_model.flux_upper_bound = 1.0;
tRNA_charging_M_arg_L_c_deGFP_model.flux_bounds_model = Bounds;
tRNA_charging_M_arg_L_c_deGFP_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0]);
tRNA_charging_M_arg_L_c_deGFP_model.flux_bound_alpha = 1.0;
tRNA_charging_M_arg_L_c_deGFP_model.flux_obj_coeff = 0.0;
flux_model_dictionary["tRNA_charging_M_arg_L_c_deGFP"] = tRNA_charging_M_arg_L_c_deGFP_model;
tRNA_charging_M_arg_L_c_deGFP_model = 0;

# 169 tRNA_charging_M_asn_L_c_deGFP: 13.0*M_asn_L_c+13.0*M_atp_c+13.0*tRNA -([])-> 13.0*M_asn_L_c_tRNA+13.0*M_amp_c+26.0*M_pi_c
tRNA_charging_M_asn_L_c_deGFP_model = FluxModel();
tRNA_charging_M_asn_L_c_deGFP_model.flux_index = 169
tRNA_charging_M_asn_L_c_deGFP_model.flux_symbol = "tRNA_charging_M_asn_L_c_deGFP"
tRNA_charging_M_asn_L_c_deGFP_model.flux_constraint_type = GLPK.DB;
tRNA_charging_M_asn_L_c_deGFP_model.flux_lower_bound = 0.0;
tRNA_charging_M_asn_L_c_deGFP_model.flux_upper_bound = 1.0;
tRNA_charging_M_asn_L_c_deGFP_model.flux_bounds_model = Bounds;
tRNA_charging_M_asn_L_c_deGFP_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0]);
tRNA_charging_M_asn_L_c_deGFP_model.flux_bound_alpha = 1.0;
tRNA_charging_M_asn_L_c_deGFP_model.flux_obj_coeff = 0.0;
flux_model_dictionary["tRNA_charging_M_asn_L_c_deGFP"] = tRNA_charging_M_asn_L_c_deGFP_model;
tRNA_charging_M_asn_L_c_deGFP_model = 0;

# 170 tRNA_charging_M_asp_L_c_deGFP: 17.0*M_asp_L_c+17.0*M_atp_c+17.0*tRNA -([])-> 17.0*M_asp_L_c_tRNA+17.0*M_amp_c+34.0*M_pi_c
tRNA_charging_M_asp_L_c_deGFP_model = FluxModel();
tRNA_charging_M_asp_L_c_deGFP_model.flux_index = 170
tRNA_charging_M_asp_L_c_deGFP_model.flux_symbol = "tRNA_charging_M_asp_L_c_deGFP"
tRNA_charging_M_asp_L_c_deGFP_model.flux_constraint_type = GLPK.DB;
tRNA_charging_M_asp_L_c_deGFP_model.flux_lower_bound = 0.0;
tRNA_charging_M_asp_L_c_deGFP_model.flux_upper_bound = 1.0;
tRNA_charging_M_asp_L_c_deGFP_model.flux_bounds_model = Bounds;
tRNA_charging_M_asp_L_c_deGFP_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0]);
tRNA_charging_M_asp_L_c_deGFP_model.flux_bound_alpha = 1.0;
tRNA_charging_M_asp_L_c_deGFP_model.flux_obj_coeff = 0.0;
flux_model_dictionary["tRNA_charging_M_asp_L_c_deGFP"] = tRNA_charging_M_asp_L_c_deGFP_model;
tRNA_charging_M_asp_L_c_deGFP_model = 0;

# 171 tRNA_charging_M_cys_L_c_deGFP: 2.0*M_cys_L_c+2.0*M_atp_c+2.0*tRNA -([])-> 2.0*M_cys_L_c_tRNA+2.0*M_amp_c+4.0*M_pi_c
tRNA_charging_M_cys_L_c_deGFP_model = FluxModel();
tRNA_charging_M_cys_L_c_deGFP_model.flux_index = 171
tRNA_charging_M_cys_L_c_deGFP_model.flux_symbol = "tRNA_charging_M_cys_L_c_deGFP"
tRNA_charging_M_cys_L_c_deGFP_model.flux_constraint_type = GLPK.DB;
tRNA_charging_M_cys_L_c_deGFP_model.flux_lower_bound = 0.0;
tRNA_charging_M_cys_L_c_deGFP_model.flux_upper_bound = 1.0;
tRNA_charging_M_cys_L_c_deGFP_model.flux_bounds_model = Bounds;
tRNA_charging_M_cys_L_c_deGFP_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0]);
tRNA_charging_M_cys_L_c_deGFP_model.flux_bound_alpha = 1.0;
tRNA_charging_M_cys_L_c_deGFP_model.flux_obj_coeff = 0.0;
flux_model_dictionary["tRNA_charging_M_cys_L_c_deGFP"] = tRNA_charging_M_cys_L_c_deGFP_model;
tRNA_charging_M_cys_L_c_deGFP_model = 0;

# 172 tRNA_charging_M_glu_L_c_deGFP: 14.0*M_glu_L_c+14.0*M_atp_c+14.0*tRNA -([])-> 14.0*M_glu_L_c_tRNA+14.0*M_amp_c+28.0*M_pi_c
tRNA_charging_M_glu_L_c_deGFP_model = FluxModel();
tRNA_charging_M_glu_L_c_deGFP_model.flux_index = 172
tRNA_charging_M_glu_L_c_deGFP_model.flux_symbol = "tRNA_charging_M_glu_L_c_deGFP"
tRNA_charging_M_glu_L_c_deGFP_model.flux_constraint_type = GLPK.DB;
tRNA_charging_M_glu_L_c_deGFP_model.flux_lower_bound = 0.0;
tRNA_charging_M_glu_L_c_deGFP_model.flux_upper_bound = 1.0;
tRNA_charging_M_glu_L_c_deGFP_model.flux_bounds_model = Bounds;
tRNA_charging_M_glu_L_c_deGFP_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0]);
tRNA_charging_M_glu_L_c_deGFP_model.flux_bound_alpha = 1.0;
tRNA_charging_M_glu_L_c_deGFP_model.flux_obj_coeff = 0.0;
flux_model_dictionary["tRNA_charging_M_glu_L_c_deGFP"] = tRNA_charging_M_glu_L_c_deGFP_model;
tRNA_charging_M_glu_L_c_deGFP_model = 0;

# 173 tRNA_charging_M_gln_L_c_deGFP: 8.0*M_gln_L_c+8.0*M_atp_c+8.0*tRNA -([])-> 8.0*M_gln_L_c_tRNA+8.0*M_amp_c+16.0*M_pi_c
tRNA_charging_M_gln_L_c_deGFP_model = FluxModel();
tRNA_charging_M_gln_L_c_deGFP_model.flux_index = 173
tRNA_charging_M_gln_L_c_deGFP_model.flux_symbol = "tRNA_charging_M_gln_L_c_deGFP"
tRNA_charging_M_gln_L_c_deGFP_model.flux_constraint_type = GLPK.DB;
tRNA_charging_M_gln_L_c_deGFP_model.flux_lower_bound = 0.0;
tRNA_charging_M_gln_L_c_deGFP_model.flux_upper_bound = 1.0;
tRNA_charging_M_gln_L_c_deGFP_model.flux_bounds_model = Bounds;
tRNA_charging_M_gln_L_c_deGFP_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0]);
tRNA_charging_M_gln_L_c_deGFP_model.flux_bound_alpha = 1.0;
tRNA_charging_M_gln_L_c_deGFP_model.flux_obj_coeff = 0.0;
flux_model_dictionary["tRNA_charging_M_gln_L_c_deGFP"] = tRNA_charging_M_gln_L_c_deGFP_model;
tRNA_charging_M_gln_L_c_deGFP_model = 0;

# 174 tRNA_charging_M_gly_L_c_deGFP: 20.0*M_gly_L_c+20.0*M_atp_c+20.0*tRNA -([])-> 20.0*M_gly_L_c_tRNA+20.0*M_amp_c+40.0*M_pi_c
tRNA_charging_M_gly_L_c_deGFP_model = FluxModel();
tRNA_charging_M_gly_L_c_deGFP_model.flux_index = 174
tRNA_charging_M_gly_L_c_deGFP_model.flux_symbol = "tRNA_charging_M_gly_L_c_deGFP"
tRNA_charging_M_gly_L_c_deGFP_model.flux_constraint_type = GLPK.DB;
tRNA_charging_M_gly_L_c_deGFP_model.flux_lower_bound = 0.0;
tRNA_charging_M_gly_L_c_deGFP_model.flux_upper_bound = 1.0;
tRNA_charging_M_gly_L_c_deGFP_model.flux_bounds_model = Bounds;
tRNA_charging_M_gly_L_c_deGFP_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0]);
tRNA_charging_M_gly_L_c_deGFP_model.flux_bound_alpha = 1.0;
tRNA_charging_M_gly_L_c_deGFP_model.flux_obj_coeff = 0.0;
flux_model_dictionary["tRNA_charging_M_gly_L_c_deGFP"] = tRNA_charging_M_gly_L_c_deGFP_model;
tRNA_charging_M_gly_L_c_deGFP_model = 0;

# 175 tRNA_charging_M_his_L_c_deGFP: 9.0*M_his_L_c+9.0*M_atp_c+9.0*tRNA -([])-> 9.0*M_his_L_c_tRNA+9.0*M_amp_c+18.0*M_pi_c
tRNA_charging_M_his_L_c_deGFP_model = FluxModel();
tRNA_charging_M_his_L_c_deGFP_model.flux_index = 175
tRNA_charging_M_his_L_c_deGFP_model.flux_symbol = "tRNA_charging_M_his_L_c_deGFP"
tRNA_charging_M_his_L_c_deGFP_model.flux_constraint_type = GLPK.DB;
tRNA_charging_M_his_L_c_deGFP_model.flux_lower_bound = 0.0;
tRNA_charging_M_his_L_c_deGFP_model.flux_upper_bound = 1.0;
tRNA_charging_M_his_L_c_deGFP_model.flux_bounds_model = Bounds;
tRNA_charging_M_his_L_c_deGFP_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0]);
tRNA_charging_M_his_L_c_deGFP_model.flux_bound_alpha = 1.0;
tRNA_charging_M_his_L_c_deGFP_model.flux_obj_coeff = 0.0;
flux_model_dictionary["tRNA_charging_M_his_L_c_deGFP"] = tRNA_charging_M_his_L_c_deGFP_model;
tRNA_charging_M_his_L_c_deGFP_model = 0;

# 176 tRNA_charging_M_ile_L_c_deGFP: 12.0*M_ile_L_c+12.0*M_atp_c+12.0*tRNA -([])-> 12.0*M_ile_L_c_tRNA+12.0*M_amp_c+24.0*M_pi_c
tRNA_charging_M_ile_L_c_deGFP_model = FluxModel();
tRNA_charging_M_ile_L_c_deGFP_model.flux_index = 176
tRNA_charging_M_ile_L_c_deGFP_model.flux_symbol = "tRNA_charging_M_ile_L_c_deGFP"
tRNA_charging_M_ile_L_c_deGFP_model.flux_constraint_type = GLPK.DB;
tRNA_charging_M_ile_L_c_deGFP_model.flux_lower_bound = 0.0;
tRNA_charging_M_ile_L_c_deGFP_model.flux_upper_bound = 1.0;
tRNA_charging_M_ile_L_c_deGFP_model.flux_bounds_model = Bounds;
tRNA_charging_M_ile_L_c_deGFP_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0]);
tRNA_charging_M_ile_L_c_deGFP_model.flux_bound_alpha = 1.0;
tRNA_charging_M_ile_L_c_deGFP_model.flux_obj_coeff = 0.0;
flux_model_dictionary["tRNA_charging_M_ile_L_c_deGFP"] = tRNA_charging_M_ile_L_c_deGFP_model;
tRNA_charging_M_ile_L_c_deGFP_model = 0;

# 177 tRNA_charging_M_leu_L_c_deGFP: 19.0*M_leu_L_c+19.0*M_atp_c+19.0*tRNA -([])-> 19.0*M_leu_L_c_tRNA+19.0*M_amp_c+38.0*M_pi_c
tRNA_charging_M_leu_L_c_deGFP_model = FluxModel();
tRNA_charging_M_leu_L_c_deGFP_model.flux_index = 177
tRNA_charging_M_leu_L_c_deGFP_model.flux_symbol = "tRNA_charging_M_leu_L_c_deGFP"
tRNA_charging_M_leu_L_c_deGFP_model.flux_constraint_type = GLPK.DB;
tRNA_charging_M_leu_L_c_deGFP_model.flux_lower_bound = 0.0;
tRNA_charging_M_leu_L_c_deGFP_model.flux_upper_bound = 1.0;
tRNA_charging_M_leu_L_c_deGFP_model.flux_bounds_model = Bounds;
tRNA_charging_M_leu_L_c_deGFP_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0]);
tRNA_charging_M_leu_L_c_deGFP_model.flux_bound_alpha = 1.0;
tRNA_charging_M_leu_L_c_deGFP_model.flux_obj_coeff = 0.0;
flux_model_dictionary["tRNA_charging_M_leu_L_c_deGFP"] = tRNA_charging_M_leu_L_c_deGFP_model;
tRNA_charging_M_leu_L_c_deGFP_model = 0;

# 178 tRNA_charging_M_lys_L_c_deGFP: 18.0*M_lys_L_c+18.0*M_atp_c+18.0*tRNA -([])-> 18.0*M_lys_L_c_tRNA+18.0*M_amp_c+36.0*M_pi_c
tRNA_charging_M_lys_L_c_deGFP_model = FluxModel();
tRNA_charging_M_lys_L_c_deGFP_model.flux_index = 178
tRNA_charging_M_lys_L_c_deGFP_model.flux_symbol = "tRNA_charging_M_lys_L_c_deGFP"
tRNA_charging_M_lys_L_c_deGFP_model.flux_constraint_type = GLPK.DB;
tRNA_charging_M_lys_L_c_deGFP_model.flux_lower_bound = 0.0;
tRNA_charging_M_lys_L_c_deGFP_model.flux_upper_bound = 1.0;
tRNA_charging_M_lys_L_c_deGFP_model.flux_bounds_model = Bounds;
tRNA_charging_M_lys_L_c_deGFP_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0]);
tRNA_charging_M_lys_L_c_deGFP_model.flux_bound_alpha = 1.0;
tRNA_charging_M_lys_L_c_deGFP_model.flux_obj_coeff = 0.0;
flux_model_dictionary["tRNA_charging_M_lys_L_c_deGFP"] = tRNA_charging_M_lys_L_c_deGFP_model;
tRNA_charging_M_lys_L_c_deGFP_model = 0;

# 179 tRNA_charging_M_met_L_c_deGFP: 5.0*M_met_L_c+5.0*M_atp_c+5.0*tRNA -([])-> 5.0*M_met_L_c_tRNA+5.0*M_amp_c+10.0*M_pi_c
tRNA_charging_M_met_L_c_deGFP_model = FluxModel();
tRNA_charging_M_met_L_c_deGFP_model.flux_index = 179
tRNA_charging_M_met_L_c_deGFP_model.flux_symbol = "tRNA_charging_M_met_L_c_deGFP"
tRNA_charging_M_met_L_c_deGFP_model.flux_constraint_type = GLPK.DB;
tRNA_charging_M_met_L_c_deGFP_model.flux_lower_bound = 0.0;
tRNA_charging_M_met_L_c_deGFP_model.flux_upper_bound = 1.0;
tRNA_charging_M_met_L_c_deGFP_model.flux_bounds_model = Bounds;
tRNA_charging_M_met_L_c_deGFP_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0]);
tRNA_charging_M_met_L_c_deGFP_model.flux_bound_alpha = 1.0;
tRNA_charging_M_met_L_c_deGFP_model.flux_obj_coeff = 0.0;
flux_model_dictionary["tRNA_charging_M_met_L_c_deGFP"] = tRNA_charging_M_met_L_c_deGFP_model;
tRNA_charging_M_met_L_c_deGFP_model = 0;

# 180 tRNA_charging_M_phe_L_c_deGFP: 12.0*M_phe_L_c+12.0*M_atp_c+12.0*tRNA -([])-> 12.0*M_phe_L_c_tRNA+12.0*M_amp_c+24.0*M_pi_c
tRNA_charging_M_phe_L_c_deGFP_model = FluxModel();
tRNA_charging_M_phe_L_c_deGFP_model.flux_index = 180
tRNA_charging_M_phe_L_c_deGFP_model.flux_symbol = "tRNA_charging_M_phe_L_c_deGFP"
tRNA_charging_M_phe_L_c_deGFP_model.flux_constraint_type = GLPK.DB;
tRNA_charging_M_phe_L_c_deGFP_model.flux_lower_bound = 0.0;
tRNA_charging_M_phe_L_c_deGFP_model.flux_upper_bound = 1.0;
tRNA_charging_M_phe_L_c_deGFP_model.flux_bounds_model = Bounds;
tRNA_charging_M_phe_L_c_deGFP_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0]);
tRNA_charging_M_phe_L_c_deGFP_model.flux_bound_alpha = 1.0;
tRNA_charging_M_phe_L_c_deGFP_model.flux_obj_coeff = 0.0;
flux_model_dictionary["tRNA_charging_M_phe_L_c_deGFP"] = tRNA_charging_M_phe_L_c_deGFP_model;
tRNA_charging_M_phe_L_c_deGFP_model = 0;

# 181 tRNA_charging_M_pro_L_c_deGFP: 10.0*M_pro_L_c+10.0*M_atp_c+10.0*tRNA -([])-> 10.0*M_pro_L_c_tRNA+10.0*M_amp_c+20.0*M_pi_c
tRNA_charging_M_pro_L_c_deGFP_model = FluxModel();
tRNA_charging_M_pro_L_c_deGFP_model.flux_index = 181
tRNA_charging_M_pro_L_c_deGFP_model.flux_symbol = "tRNA_charging_M_pro_L_c_deGFP"
tRNA_charging_M_pro_L_c_deGFP_model.flux_constraint_type = GLPK.DB;
tRNA_charging_M_pro_L_c_deGFP_model.flux_lower_bound = 0.0;
tRNA_charging_M_pro_L_c_deGFP_model.flux_upper_bound = 1.0;
tRNA_charging_M_pro_L_c_deGFP_model.flux_bounds_model = Bounds;
tRNA_charging_M_pro_L_c_deGFP_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0]);
tRNA_charging_M_pro_L_c_deGFP_model.flux_bound_alpha = 1.0;
tRNA_charging_M_pro_L_c_deGFP_model.flux_obj_coeff = 0.0;
flux_model_dictionary["tRNA_charging_M_pro_L_c_deGFP"] = tRNA_charging_M_pro_L_c_deGFP_model;
tRNA_charging_M_pro_L_c_deGFP_model = 0;

# 182 tRNA_charging_M_ser_L_c_deGFP: 9.0*M_ser_L_c+9.0*M_atp_c+9.0*tRNA -([])-> 9.0*M_ser_L_c_tRNA+9.0*M_amp_c+18.0*M_pi_c
tRNA_charging_M_ser_L_c_deGFP_model = FluxModel();
tRNA_charging_M_ser_L_c_deGFP_model.flux_index = 182
tRNA_charging_M_ser_L_c_deGFP_model.flux_symbol = "tRNA_charging_M_ser_L_c_deGFP"
tRNA_charging_M_ser_L_c_deGFP_model.flux_constraint_type = GLPK.DB;
tRNA_charging_M_ser_L_c_deGFP_model.flux_lower_bound = 0.0;
tRNA_charging_M_ser_L_c_deGFP_model.flux_upper_bound = 1.0;
tRNA_charging_M_ser_L_c_deGFP_model.flux_bounds_model = Bounds;
tRNA_charging_M_ser_L_c_deGFP_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0]);
tRNA_charging_M_ser_L_c_deGFP_model.flux_bound_alpha = 1.0;
tRNA_charging_M_ser_L_c_deGFP_model.flux_obj_coeff = 0.0;
flux_model_dictionary["tRNA_charging_M_ser_L_c_deGFP"] = tRNA_charging_M_ser_L_c_deGFP_model;
tRNA_charging_M_ser_L_c_deGFP_model = 0;

# 183 tRNA_charging_M_thr_L_c_deGFP: 15.0*M_thr_L_c+15.0*M_atp_c+15.0*tRNA -([])-> 15.0*M_thr_L_c_tRNA+15.0*M_amp_c+30.0*M_pi_c
tRNA_charging_M_thr_L_c_deGFP_model = FluxModel();
tRNA_charging_M_thr_L_c_deGFP_model.flux_index = 183
tRNA_charging_M_thr_L_c_deGFP_model.flux_symbol = "tRNA_charging_M_thr_L_c_deGFP"
tRNA_charging_M_thr_L_c_deGFP_model.flux_constraint_type = GLPK.DB;
tRNA_charging_M_thr_L_c_deGFP_model.flux_lower_bound = 0.0;
tRNA_charging_M_thr_L_c_deGFP_model.flux_upper_bound = 1.0;
tRNA_charging_M_thr_L_c_deGFP_model.flux_bounds_model = Bounds;
tRNA_charging_M_thr_L_c_deGFP_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0]);
tRNA_charging_M_thr_L_c_deGFP_model.flux_bound_alpha = 1.0;
tRNA_charging_M_thr_L_c_deGFP_model.flux_obj_coeff = 0.0;
flux_model_dictionary["tRNA_charging_M_thr_L_c_deGFP"] = tRNA_charging_M_thr_L_c_deGFP_model;
tRNA_charging_M_thr_L_c_deGFP_model = 0;

# 184 tRNA_charging_M_trp_L_c_deGFP: 1.0*M_trp_L_c+1.0*M_atp_c+1.0*tRNA -([])-> 1.0*M_trp_L_c_tRNA+1.0*M_amp_c+2.0*M_pi_c
tRNA_charging_M_trp_L_c_deGFP_model = FluxModel();
tRNA_charging_M_trp_L_c_deGFP_model.flux_index = 184
tRNA_charging_M_trp_L_c_deGFP_model.flux_symbol = "tRNA_charging_M_trp_L_c_deGFP"
tRNA_charging_M_trp_L_c_deGFP_model.flux_constraint_type = GLPK.DB;
tRNA_charging_M_trp_L_c_deGFP_model.flux_lower_bound = 0.0;
tRNA_charging_M_trp_L_c_deGFP_model.flux_upper_bound = 1.0;
tRNA_charging_M_trp_L_c_deGFP_model.flux_bounds_model = Bounds;
tRNA_charging_M_trp_L_c_deGFP_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0]);
tRNA_charging_M_trp_L_c_deGFP_model.flux_bound_alpha = 1.0;
tRNA_charging_M_trp_L_c_deGFP_model.flux_obj_coeff = 0.0;
flux_model_dictionary["tRNA_charging_M_trp_L_c_deGFP"] = tRNA_charging_M_trp_L_c_deGFP_model;
tRNA_charging_M_trp_L_c_deGFP_model = 0;

# 185 tRNA_charging_M_tyr_L_c_deGFP: 10.0*M_tyr_L_c+10.0*M_atp_c+10.0*tRNA -([])-> 10.0*M_tyr_L_c_tRNA+10.0*M_amp_c+20.0*M_pi_c
tRNA_charging_M_tyr_L_c_deGFP_model = FluxModel();
tRNA_charging_M_tyr_L_c_deGFP_model.flux_index = 185
tRNA_charging_M_tyr_L_c_deGFP_model.flux_symbol = "tRNA_charging_M_tyr_L_c_deGFP"
tRNA_charging_M_tyr_L_c_deGFP_model.flux_constraint_type = GLPK.DB;
tRNA_charging_M_tyr_L_c_deGFP_model.flux_lower_bound = 0.0;
tRNA_charging_M_tyr_L_c_deGFP_model.flux_upper_bound = 1.0;
tRNA_charging_M_tyr_L_c_deGFP_model.flux_bounds_model = Bounds;
tRNA_charging_M_tyr_L_c_deGFP_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0]);
tRNA_charging_M_tyr_L_c_deGFP_model.flux_bound_alpha = 1.0;
tRNA_charging_M_tyr_L_c_deGFP_model.flux_obj_coeff = 0.0;
flux_model_dictionary["tRNA_charging_M_tyr_L_c_deGFP"] = tRNA_charging_M_tyr_L_c_deGFP_model;
tRNA_charging_M_tyr_L_c_deGFP_model = 0;

# 186 tRNA_charging_M_val_L_c_deGFP: 17.0*M_val_L_c+17.0*M_atp_c+17.0*tRNA -([])-> 17.0*M_val_L_c_tRNA+17.0*M_amp_c+34.0*M_pi_c
tRNA_charging_M_val_L_c_deGFP_model = FluxModel();
tRNA_charging_M_val_L_c_deGFP_model.flux_index = 186
tRNA_charging_M_val_L_c_deGFP_model.flux_symbol = "tRNA_charging_M_val_L_c_deGFP"
tRNA_charging_M_val_L_c_deGFP_model.flux_constraint_type = GLPK.DB;
tRNA_charging_M_val_L_c_deGFP_model.flux_lower_bound = 0.0;
tRNA_charging_M_val_L_c_deGFP_model.flux_upper_bound = 1.0;
tRNA_charging_M_val_L_c_deGFP_model.flux_bounds_model = Bounds;
tRNA_charging_M_val_L_c_deGFP_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0]);
tRNA_charging_M_val_L_c_deGFP_model.flux_bound_alpha = 1.0;
tRNA_charging_M_val_L_c_deGFP_model.flux_obj_coeff = 0.0;
flux_model_dictionary["tRNA_charging_M_val_L_c_deGFP"] = tRNA_charging_M_val_L_c_deGFP_model;
tRNA_charging_M_val_L_c_deGFP_model = 0;

# 187 tNRA_exchange: tRNA -([])-> []
tNRA_exchange_model = FluxModel();
tNRA_exchange_model.flux_index = 187
tNRA_exchange_model.flux_symbol = "tNRA_exchange"
tNRA_exchange_model.flux_constraint_type = GLPK.DB;
tNRA_exchange_model.flux_lower_bound = 0.0;
tNRA_exchange_model.flux_upper_bound = 1.0;
tNRA_exchange_model.flux_bounds_model = Bounds;
tNRA_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0]);
tNRA_exchange_model.flux_bound_alpha = 1.0;
tNRA_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["tNRA_exchange"] = tNRA_exchange_model;
tNRA_exchange_model = 0;

# 188 -1*(tNRA_exchange: tRNA -([])-> [])
tNRA_exchange_reverse_model = FluxModel();
tNRA_exchange_reverse_model.flux_index = 188
tNRA_exchange_reverse_model.flux_symbol = "tNRA_exchange_reverse"
tNRA_exchange_reverse_model.flux_constraint_type = GLPK.DB;
tNRA_exchange_reverse_model.flux_lower_bound = 0.0;
tNRA_exchange_reverse_model.flux_upper_bound = 1.0;
tNRA_exchange_reverse_model.flux_bounds_model = Bounds;
tNRA_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
tNRA_exchange_reverse_model.flux_bound_alpha = 1.0;
tNRA_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["tNRA_exchange_reverse"] = tNRA_exchange_reverse_model;
tNRA_exchange_reverse_model = 0;

# 189 PROTEIN_export_deGFP: PROTEIN_deGFP -([])-> []
PROTEIN_export_deGFP_model = FluxModel();
PROTEIN_export_deGFP_model.flux_index = 189
PROTEIN_export_deGFP_model.flux_symbol = "PROTEIN_export_deGFP"
PROTEIN_export_deGFP_model.flux_constraint_type = GLPK.DB;
PROTEIN_export_deGFP_model.flux_lower_bound = 0.0;
PROTEIN_export_deGFP_model.flux_upper_bound = 1.0;
PROTEIN_export_deGFP_model.flux_bounds_model = Bounds;
PROTEIN_export_deGFP_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0]);
PROTEIN_export_deGFP_model.flux_bound_alpha = 1.0;
PROTEIN_export_deGFP_model.flux_obj_coeff = 0.0;
flux_model_dictionary["PROTEIN_export_deGFP"] = PROTEIN_export_deGFP_model;
PROTEIN_export_deGFP_model = 0;

# 190 M_ala_L_c_exchange: M_ala_L_c -([])-> []
M_ala_L_c_exchange_model = FluxModel();
M_ala_L_c_exchange_model.flux_index = 190
M_ala_L_c_exchange_model.flux_symbol = "M_ala_L_c_exchange"
M_ala_L_c_exchange_model.flux_constraint_type = GLPK.DB;
M_ala_L_c_exchange_model.flux_lower_bound = 0.0;
M_ala_L_c_exchange_model.flux_upper_bound = 1.0;
M_ala_L_c_exchange_model.flux_bounds_model = Bounds;
M_ala_L_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_ala_L_c_exchange_model.flux_bound_alpha = 1.0;
M_ala_L_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_ala_L_c_exchange"] = M_ala_L_c_exchange_model;
M_ala_L_c_exchange_model = 0;

# 191 -1*(M_ala_L_c_exchange: M_ala_L_c -([])-> [])
M_ala_L_c_exchange_reverse_model = FluxModel();
M_ala_L_c_exchange_reverse_model.flux_index = 191
M_ala_L_c_exchange_reverse_model.flux_symbol = "M_ala_L_c_exchange_reverse"
M_ala_L_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_ala_L_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_ala_L_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_ala_L_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_ala_L_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_ala_L_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_ala_L_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_ala_L_c_exchange_reverse"] = M_ala_L_c_exchange_reverse_model;
M_ala_L_c_exchange_reverse_model = 0;

# 192 M_arg_L_c_exchange: M_arg_L_c -([])-> []
M_arg_L_c_exchange_model = FluxModel();
M_arg_L_c_exchange_model.flux_index = 192
M_arg_L_c_exchange_model.flux_symbol = "M_arg_L_c_exchange"
M_arg_L_c_exchange_model.flux_constraint_type = GLPK.DB;
M_arg_L_c_exchange_model.flux_lower_bound = 0.0;
M_arg_L_c_exchange_model.flux_upper_bound = 1.0;
M_arg_L_c_exchange_model.flux_bounds_model = Bounds;
M_arg_L_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_arg_L_c_exchange_model.flux_bound_alpha = 1.0;
M_arg_L_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_arg_L_c_exchange"] = M_arg_L_c_exchange_model;
M_arg_L_c_exchange_model = 0;

# 193 -1*(M_arg_L_c_exchange: M_arg_L_c -([])-> [])
M_arg_L_c_exchange_reverse_model = FluxModel();
M_arg_L_c_exchange_reverse_model.flux_index = 193
M_arg_L_c_exchange_reverse_model.flux_symbol = "M_arg_L_c_exchange_reverse"
M_arg_L_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_arg_L_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_arg_L_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_arg_L_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_arg_L_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_arg_L_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_arg_L_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_arg_L_c_exchange_reverse"] = M_arg_L_c_exchange_reverse_model;
M_arg_L_c_exchange_reverse_model = 0;

# 194 M_asn_L_c_exchange: M_asn_L_c -([])-> []
M_asn_L_c_exchange_model = FluxModel();
M_asn_L_c_exchange_model.flux_index = 194
M_asn_L_c_exchange_model.flux_symbol = "M_asn_L_c_exchange"
M_asn_L_c_exchange_model.flux_constraint_type = GLPK.DB;
M_asn_L_c_exchange_model.flux_lower_bound = 0.0;
M_asn_L_c_exchange_model.flux_upper_bound = 1.0;
M_asn_L_c_exchange_model.flux_bounds_model = Bounds;
M_asn_L_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_asn_L_c_exchange_model.flux_bound_alpha = 1.0;
M_asn_L_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_asn_L_c_exchange"] = M_asn_L_c_exchange_model;
M_asn_L_c_exchange_model = 0;

# 195 -1*(M_asn_L_c_exchange: M_asn_L_c -([])-> [])
M_asn_L_c_exchange_reverse_model = FluxModel();
M_asn_L_c_exchange_reverse_model.flux_index = 195
M_asn_L_c_exchange_reverse_model.flux_symbol = "M_asn_L_c_exchange_reverse"
M_asn_L_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_asn_L_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_asn_L_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_asn_L_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_asn_L_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_asn_L_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_asn_L_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_asn_L_c_exchange_reverse"] = M_asn_L_c_exchange_reverse_model;
M_asn_L_c_exchange_reverse_model = 0;

# 196 M_asp_L_c_exchange: M_asp_L_c -([])-> []
M_asp_L_c_exchange_model = FluxModel();
M_asp_L_c_exchange_model.flux_index = 196
M_asp_L_c_exchange_model.flux_symbol = "M_asp_L_c_exchange"
M_asp_L_c_exchange_model.flux_constraint_type = GLPK.DB;
M_asp_L_c_exchange_model.flux_lower_bound = 0.0;
M_asp_L_c_exchange_model.flux_upper_bound = 1.0;
M_asp_L_c_exchange_model.flux_bounds_model = Bounds;
M_asp_L_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_asp_L_c_exchange_model.flux_bound_alpha = 1.0;
M_asp_L_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_asp_L_c_exchange"] = M_asp_L_c_exchange_model;
M_asp_L_c_exchange_model = 0;

# 197 -1*(M_asp_L_c_exchange: M_asp_L_c -([])-> [])
M_asp_L_c_exchange_reverse_model = FluxModel();
M_asp_L_c_exchange_reverse_model.flux_index = 197
M_asp_L_c_exchange_reverse_model.flux_symbol = "M_asp_L_c_exchange_reverse"
M_asp_L_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_asp_L_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_asp_L_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_asp_L_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_asp_L_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_asp_L_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_asp_L_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_asp_L_c_exchange_reverse"] = M_asp_L_c_exchange_reverse_model;
M_asp_L_c_exchange_reverse_model = 0;

# 198 M_cys_L_c_exchange: M_cys_L_c -([])-> []
M_cys_L_c_exchange_model = FluxModel();
M_cys_L_c_exchange_model.flux_index = 198
M_cys_L_c_exchange_model.flux_symbol = "M_cys_L_c_exchange"
M_cys_L_c_exchange_model.flux_constraint_type = GLPK.DB;
M_cys_L_c_exchange_model.flux_lower_bound = 0.0;
M_cys_L_c_exchange_model.flux_upper_bound = 1.0;
M_cys_L_c_exchange_model.flux_bounds_model = Bounds;
M_cys_L_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_cys_L_c_exchange_model.flux_bound_alpha = 1.0;
M_cys_L_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_cys_L_c_exchange"] = M_cys_L_c_exchange_model;
M_cys_L_c_exchange_model = 0;

# 199 -1*(M_cys_L_c_exchange: M_cys_L_c -([])-> [])
M_cys_L_c_exchange_reverse_model = FluxModel();
M_cys_L_c_exchange_reverse_model.flux_index = 199
M_cys_L_c_exchange_reverse_model.flux_symbol = "M_cys_L_c_exchange_reverse"
M_cys_L_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_cys_L_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_cys_L_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_cys_L_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_cys_L_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_cys_L_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_cys_L_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_cys_L_c_exchange_reverse"] = M_cys_L_c_exchange_reverse_model;
M_cys_L_c_exchange_reverse_model = 0;

# 200 M_glu_L_c_exchange: M_glu_L_c -([])-> []
M_glu_L_c_exchange_model = FluxModel();
M_glu_L_c_exchange_model.flux_index = 200
M_glu_L_c_exchange_model.flux_symbol = "M_glu_L_c_exchange"
M_glu_L_c_exchange_model.flux_constraint_type = GLPK.DB;
M_glu_L_c_exchange_model.flux_lower_bound = 0.0;
M_glu_L_c_exchange_model.flux_upper_bound = 1.0;
M_glu_L_c_exchange_model.flux_bounds_model = Bounds;
M_glu_L_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_glu_L_c_exchange_model.flux_bound_alpha = 1.0;
M_glu_L_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_glu_L_c_exchange"] = M_glu_L_c_exchange_model;
M_glu_L_c_exchange_model = 0;

# 201 -1*(M_glu_L_c_exchange: M_glu_L_c -([])-> [])
M_glu_L_c_exchange_reverse_model = FluxModel();
M_glu_L_c_exchange_reverse_model.flux_index = 201
M_glu_L_c_exchange_reverse_model.flux_symbol = "M_glu_L_c_exchange_reverse"
M_glu_L_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_glu_L_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_glu_L_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_glu_L_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_glu_L_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_glu_L_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_glu_L_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_glu_L_c_exchange_reverse"] = M_glu_L_c_exchange_reverse_model;
M_glu_L_c_exchange_reverse_model = 0;

# 202 M_gln_L_c_exchange: M_gln_L_c -([])-> []
M_gln_L_c_exchange_model = FluxModel();
M_gln_L_c_exchange_model.flux_index = 202
M_gln_L_c_exchange_model.flux_symbol = "M_gln_L_c_exchange"
M_gln_L_c_exchange_model.flux_constraint_type = GLPK.DB;
M_gln_L_c_exchange_model.flux_lower_bound = 0.0;
M_gln_L_c_exchange_model.flux_upper_bound = 1.0;
M_gln_L_c_exchange_model.flux_bounds_model = Bounds;
M_gln_L_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_gln_L_c_exchange_model.flux_bound_alpha = 1.0;
M_gln_L_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_gln_L_c_exchange"] = M_gln_L_c_exchange_model;
M_gln_L_c_exchange_model = 0;

# 203 -1*(M_gln_L_c_exchange: M_gln_L_c -([])-> [])
M_gln_L_c_exchange_reverse_model = FluxModel();
M_gln_L_c_exchange_reverse_model.flux_index = 203
M_gln_L_c_exchange_reverse_model.flux_symbol = "M_gln_L_c_exchange_reverse"
M_gln_L_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_gln_L_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_gln_L_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_gln_L_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_gln_L_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_gln_L_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_gln_L_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_gln_L_c_exchange_reverse"] = M_gln_L_c_exchange_reverse_model;
M_gln_L_c_exchange_reverse_model = 0;

# 204 M_gly_L_c_exchange: M_gly_L_c -([])-> []
M_gly_L_c_exchange_model = FluxModel();
M_gly_L_c_exchange_model.flux_index = 204
M_gly_L_c_exchange_model.flux_symbol = "M_gly_L_c_exchange"
M_gly_L_c_exchange_model.flux_constraint_type = GLPK.DB;
M_gly_L_c_exchange_model.flux_lower_bound = 0.0;
M_gly_L_c_exchange_model.flux_upper_bound = 1.0;
M_gly_L_c_exchange_model.flux_bounds_model = Bounds;
M_gly_L_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_gly_L_c_exchange_model.flux_bound_alpha = 1.0;
M_gly_L_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_gly_L_c_exchange"] = M_gly_L_c_exchange_model;
M_gly_L_c_exchange_model = 0;

# 205 -1*(M_gly_L_c_exchange: M_gly_L_c -([])-> [])
M_gly_L_c_exchange_reverse_model = FluxModel();
M_gly_L_c_exchange_reverse_model.flux_index = 205
M_gly_L_c_exchange_reverse_model.flux_symbol = "M_gly_L_c_exchange_reverse"
M_gly_L_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_gly_L_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_gly_L_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_gly_L_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_gly_L_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_gly_L_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_gly_L_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_gly_L_c_exchange_reverse"] = M_gly_L_c_exchange_reverse_model;
M_gly_L_c_exchange_reverse_model = 0;

# 206 M_his_L_c_exchange: M_his_L_c -([])-> []
M_his_L_c_exchange_model = FluxModel();
M_his_L_c_exchange_model.flux_index = 206
M_his_L_c_exchange_model.flux_symbol = "M_his_L_c_exchange"
M_his_L_c_exchange_model.flux_constraint_type = GLPK.DB;
M_his_L_c_exchange_model.flux_lower_bound = 0.0;
M_his_L_c_exchange_model.flux_upper_bound = 1.0;
M_his_L_c_exchange_model.flux_bounds_model = Bounds;
M_his_L_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_his_L_c_exchange_model.flux_bound_alpha = 1.0;
M_his_L_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_his_L_c_exchange"] = M_his_L_c_exchange_model;
M_his_L_c_exchange_model = 0;

# 207 -1*(M_his_L_c_exchange: M_his_L_c -([])-> [])
M_his_L_c_exchange_reverse_model = FluxModel();
M_his_L_c_exchange_reverse_model.flux_index = 207
M_his_L_c_exchange_reverse_model.flux_symbol = "M_his_L_c_exchange_reverse"
M_his_L_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_his_L_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_his_L_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_his_L_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_his_L_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_his_L_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_his_L_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_his_L_c_exchange_reverse"] = M_his_L_c_exchange_reverse_model;
M_his_L_c_exchange_reverse_model = 0;

# 208 M_ile_L_c_exchange: M_ile_L_c -([])-> []
M_ile_L_c_exchange_model = FluxModel();
M_ile_L_c_exchange_model.flux_index = 208
M_ile_L_c_exchange_model.flux_symbol = "M_ile_L_c_exchange"
M_ile_L_c_exchange_model.flux_constraint_type = GLPK.DB;
M_ile_L_c_exchange_model.flux_lower_bound = 0.0;
M_ile_L_c_exchange_model.flux_upper_bound = 1.0;
M_ile_L_c_exchange_model.flux_bounds_model = Bounds;
M_ile_L_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_ile_L_c_exchange_model.flux_bound_alpha = 1.0;
M_ile_L_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_ile_L_c_exchange"] = M_ile_L_c_exchange_model;
M_ile_L_c_exchange_model = 0;

# 209 -1*(M_ile_L_c_exchange: M_ile_L_c -([])-> [])
M_ile_L_c_exchange_reverse_model = FluxModel();
M_ile_L_c_exchange_reverse_model.flux_index = 209
M_ile_L_c_exchange_reverse_model.flux_symbol = "M_ile_L_c_exchange_reverse"
M_ile_L_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_ile_L_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_ile_L_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_ile_L_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_ile_L_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_ile_L_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_ile_L_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_ile_L_c_exchange_reverse"] = M_ile_L_c_exchange_reverse_model;
M_ile_L_c_exchange_reverse_model = 0;

# 210 M_leu_L_c_exchange: M_leu_L_c -([])-> []
M_leu_L_c_exchange_model = FluxModel();
M_leu_L_c_exchange_model.flux_index = 210
M_leu_L_c_exchange_model.flux_symbol = "M_leu_L_c_exchange"
M_leu_L_c_exchange_model.flux_constraint_type = GLPK.DB;
M_leu_L_c_exchange_model.flux_lower_bound = 0.0;
M_leu_L_c_exchange_model.flux_upper_bound = 1.0;
M_leu_L_c_exchange_model.flux_bounds_model = Bounds;
M_leu_L_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_leu_L_c_exchange_model.flux_bound_alpha = 1.0;
M_leu_L_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_leu_L_c_exchange"] = M_leu_L_c_exchange_model;
M_leu_L_c_exchange_model = 0;

# 211 -1*(M_leu_L_c_exchange: M_leu_L_c -([])-> [])
M_leu_L_c_exchange_reverse_model = FluxModel();
M_leu_L_c_exchange_reverse_model.flux_index = 211
M_leu_L_c_exchange_reverse_model.flux_symbol = "M_leu_L_c_exchange_reverse"
M_leu_L_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_leu_L_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_leu_L_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_leu_L_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_leu_L_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_leu_L_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_leu_L_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_leu_L_c_exchange_reverse"] = M_leu_L_c_exchange_reverse_model;
M_leu_L_c_exchange_reverse_model = 0;

# 212 M_lys_L_c_exchange: M_lys_L_c -([])-> []
M_lys_L_c_exchange_model = FluxModel();
M_lys_L_c_exchange_model.flux_index = 212
M_lys_L_c_exchange_model.flux_symbol = "M_lys_L_c_exchange"
M_lys_L_c_exchange_model.flux_constraint_type = GLPK.DB;
M_lys_L_c_exchange_model.flux_lower_bound = 0.0;
M_lys_L_c_exchange_model.flux_upper_bound = 1.0;
M_lys_L_c_exchange_model.flux_bounds_model = Bounds;
M_lys_L_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_lys_L_c_exchange_model.flux_bound_alpha = 1.0;
M_lys_L_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_lys_L_c_exchange"] = M_lys_L_c_exchange_model;
M_lys_L_c_exchange_model = 0;

# 213 -1*(M_lys_L_c_exchange: M_lys_L_c -([])-> [])
M_lys_L_c_exchange_reverse_model = FluxModel();
M_lys_L_c_exchange_reverse_model.flux_index = 213
M_lys_L_c_exchange_reverse_model.flux_symbol = "M_lys_L_c_exchange_reverse"
M_lys_L_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_lys_L_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_lys_L_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_lys_L_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_lys_L_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_lys_L_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_lys_L_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_lys_L_c_exchange_reverse"] = M_lys_L_c_exchange_reverse_model;
M_lys_L_c_exchange_reverse_model = 0;

# 214 M_met_L_c_exchange: M_met_L_c -([])-> []
M_met_L_c_exchange_model = FluxModel();
M_met_L_c_exchange_model.flux_index = 214
M_met_L_c_exchange_model.flux_symbol = "M_met_L_c_exchange"
M_met_L_c_exchange_model.flux_constraint_type = GLPK.DB;
M_met_L_c_exchange_model.flux_lower_bound = 0.0;
M_met_L_c_exchange_model.flux_upper_bound = 1.0;
M_met_L_c_exchange_model.flux_bounds_model = Bounds;
M_met_L_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_met_L_c_exchange_model.flux_bound_alpha = 1.0;
M_met_L_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_met_L_c_exchange"] = M_met_L_c_exchange_model;
M_met_L_c_exchange_model = 0;

# 215 -1*(M_met_L_c_exchange: M_met_L_c -([])-> [])
M_met_L_c_exchange_reverse_model = FluxModel();
M_met_L_c_exchange_reverse_model.flux_index = 215
M_met_L_c_exchange_reverse_model.flux_symbol = "M_met_L_c_exchange_reverse"
M_met_L_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_met_L_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_met_L_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_met_L_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_met_L_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_met_L_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_met_L_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_met_L_c_exchange_reverse"] = M_met_L_c_exchange_reverse_model;
M_met_L_c_exchange_reverse_model = 0;

# 216 M_phe_L_c_exchange: M_phe_L_c -([])-> []
M_phe_L_c_exchange_model = FluxModel();
M_phe_L_c_exchange_model.flux_index = 216
M_phe_L_c_exchange_model.flux_symbol = "M_phe_L_c_exchange"
M_phe_L_c_exchange_model.flux_constraint_type = GLPK.DB;
M_phe_L_c_exchange_model.flux_lower_bound = 0.0;
M_phe_L_c_exchange_model.flux_upper_bound = 1.0;
M_phe_L_c_exchange_model.flux_bounds_model = Bounds;
M_phe_L_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_phe_L_c_exchange_model.flux_bound_alpha = 1.0;
M_phe_L_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_phe_L_c_exchange"] = M_phe_L_c_exchange_model;
M_phe_L_c_exchange_model = 0;

# 217 -1*(M_phe_L_c_exchange: M_phe_L_c -([])-> [])
M_phe_L_c_exchange_reverse_model = FluxModel();
M_phe_L_c_exchange_reverse_model.flux_index = 217
M_phe_L_c_exchange_reverse_model.flux_symbol = "M_phe_L_c_exchange_reverse"
M_phe_L_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_phe_L_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_phe_L_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_phe_L_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_phe_L_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_phe_L_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_phe_L_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_phe_L_c_exchange_reverse"] = M_phe_L_c_exchange_reverse_model;
M_phe_L_c_exchange_reverse_model = 0;

# 218 M_pro_L_c_exchange: M_pro_L_c -([])-> []
M_pro_L_c_exchange_model = FluxModel();
M_pro_L_c_exchange_model.flux_index = 218
M_pro_L_c_exchange_model.flux_symbol = "M_pro_L_c_exchange"
M_pro_L_c_exchange_model.flux_constraint_type = GLPK.DB;
M_pro_L_c_exchange_model.flux_lower_bound = 0.0;
M_pro_L_c_exchange_model.flux_upper_bound = 1.0;
M_pro_L_c_exchange_model.flux_bounds_model = Bounds;
M_pro_L_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_pro_L_c_exchange_model.flux_bound_alpha = 1.0;
M_pro_L_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_pro_L_c_exchange"] = M_pro_L_c_exchange_model;
M_pro_L_c_exchange_model = 0;

# 219 -1*(M_pro_L_c_exchange: M_pro_L_c -([])-> [])
M_pro_L_c_exchange_reverse_model = FluxModel();
M_pro_L_c_exchange_reverse_model.flux_index = 219
M_pro_L_c_exchange_reverse_model.flux_symbol = "M_pro_L_c_exchange_reverse"
M_pro_L_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_pro_L_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_pro_L_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_pro_L_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_pro_L_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_pro_L_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_pro_L_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_pro_L_c_exchange_reverse"] = M_pro_L_c_exchange_reverse_model;
M_pro_L_c_exchange_reverse_model = 0;

# 220 M_ser_L_c_exchange: M_ser_L_c -([])-> []
M_ser_L_c_exchange_model = FluxModel();
M_ser_L_c_exchange_model.flux_index = 220
M_ser_L_c_exchange_model.flux_symbol = "M_ser_L_c_exchange"
M_ser_L_c_exchange_model.flux_constraint_type = GLPK.DB;
M_ser_L_c_exchange_model.flux_lower_bound = 0.0;
M_ser_L_c_exchange_model.flux_upper_bound = 1.0;
M_ser_L_c_exchange_model.flux_bounds_model = Bounds;
M_ser_L_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_ser_L_c_exchange_model.flux_bound_alpha = 1.0;
M_ser_L_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_ser_L_c_exchange"] = M_ser_L_c_exchange_model;
M_ser_L_c_exchange_model = 0;

# 221 -1*(M_ser_L_c_exchange: M_ser_L_c -([])-> [])
M_ser_L_c_exchange_reverse_model = FluxModel();
M_ser_L_c_exchange_reverse_model.flux_index = 221
M_ser_L_c_exchange_reverse_model.flux_symbol = "M_ser_L_c_exchange_reverse"
M_ser_L_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_ser_L_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_ser_L_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_ser_L_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_ser_L_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_ser_L_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_ser_L_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_ser_L_c_exchange_reverse"] = M_ser_L_c_exchange_reverse_model;
M_ser_L_c_exchange_reverse_model = 0;

# 222 M_thr_L_c_exchange: M_thr_L_c -([])-> []
M_thr_L_c_exchange_model = FluxModel();
M_thr_L_c_exchange_model.flux_index = 222
M_thr_L_c_exchange_model.flux_symbol = "M_thr_L_c_exchange"
M_thr_L_c_exchange_model.flux_constraint_type = GLPK.DB;
M_thr_L_c_exchange_model.flux_lower_bound = 0.0;
M_thr_L_c_exchange_model.flux_upper_bound = 1.0;
M_thr_L_c_exchange_model.flux_bounds_model = Bounds;
M_thr_L_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_thr_L_c_exchange_model.flux_bound_alpha = 1.0;
M_thr_L_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_thr_L_c_exchange"] = M_thr_L_c_exchange_model;
M_thr_L_c_exchange_model = 0;

# 223 -1*(M_thr_L_c_exchange: M_thr_L_c -([])-> [])
M_thr_L_c_exchange_reverse_model = FluxModel();
M_thr_L_c_exchange_reverse_model.flux_index = 223
M_thr_L_c_exchange_reverse_model.flux_symbol = "M_thr_L_c_exchange_reverse"
M_thr_L_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_thr_L_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_thr_L_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_thr_L_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_thr_L_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_thr_L_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_thr_L_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_thr_L_c_exchange_reverse"] = M_thr_L_c_exchange_reverse_model;
M_thr_L_c_exchange_reverse_model = 0;

# 224 M_trp_L_c_exchange: M_trp_L_c -([])-> []
M_trp_L_c_exchange_model = FluxModel();
M_trp_L_c_exchange_model.flux_index = 224
M_trp_L_c_exchange_model.flux_symbol = "M_trp_L_c_exchange"
M_trp_L_c_exchange_model.flux_constraint_type = GLPK.DB;
M_trp_L_c_exchange_model.flux_lower_bound = 0.0;
M_trp_L_c_exchange_model.flux_upper_bound = 1.0;
M_trp_L_c_exchange_model.flux_bounds_model = Bounds;
M_trp_L_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_trp_L_c_exchange_model.flux_bound_alpha = 1.0;
M_trp_L_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_trp_L_c_exchange"] = M_trp_L_c_exchange_model;
M_trp_L_c_exchange_model = 0;

# 225 -1*(M_trp_L_c_exchange: M_trp_L_c -([])-> [])
M_trp_L_c_exchange_reverse_model = FluxModel();
M_trp_L_c_exchange_reverse_model.flux_index = 225
M_trp_L_c_exchange_reverse_model.flux_symbol = "M_trp_L_c_exchange_reverse"
M_trp_L_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_trp_L_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_trp_L_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_trp_L_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_trp_L_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_trp_L_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_trp_L_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_trp_L_c_exchange_reverse"] = M_trp_L_c_exchange_reverse_model;
M_trp_L_c_exchange_reverse_model = 0;

# 226 M_tyr_L_c_exchange: M_tyr_L_c -([])-> []
M_tyr_L_c_exchange_model = FluxModel();
M_tyr_L_c_exchange_model.flux_index = 226
M_tyr_L_c_exchange_model.flux_symbol = "M_tyr_L_c_exchange"
M_tyr_L_c_exchange_model.flux_constraint_type = GLPK.DB;
M_tyr_L_c_exchange_model.flux_lower_bound = 0.0;
M_tyr_L_c_exchange_model.flux_upper_bound = 1.0;
M_tyr_L_c_exchange_model.flux_bounds_model = Bounds;
M_tyr_L_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_tyr_L_c_exchange_model.flux_bound_alpha = 1.0;
M_tyr_L_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_tyr_L_c_exchange"] = M_tyr_L_c_exchange_model;
M_tyr_L_c_exchange_model = 0;

# 227 -1*(M_tyr_L_c_exchange: M_tyr_L_c -([])-> [])
M_tyr_L_c_exchange_reverse_model = FluxModel();
M_tyr_L_c_exchange_reverse_model.flux_index = 227
M_tyr_L_c_exchange_reverse_model.flux_symbol = "M_tyr_L_c_exchange_reverse"
M_tyr_L_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_tyr_L_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_tyr_L_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_tyr_L_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_tyr_L_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_tyr_L_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_tyr_L_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_tyr_L_c_exchange_reverse"] = M_tyr_L_c_exchange_reverse_model;
M_tyr_L_c_exchange_reverse_model = 0;

# 228 M_val_L_c_exchange: M_val_L_c -([])-> []
M_val_L_c_exchange_model = FluxModel();
M_val_L_c_exchange_model.flux_index = 228
M_val_L_c_exchange_model.flux_symbol = "M_val_L_c_exchange"
M_val_L_c_exchange_model.flux_constraint_type = GLPK.DB;
M_val_L_c_exchange_model.flux_lower_bound = 0.0;
M_val_L_c_exchange_model.flux_upper_bound = 1.0;
M_val_L_c_exchange_model.flux_bounds_model = Bounds;
M_val_L_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_val_L_c_exchange_model.flux_bound_alpha = 1.0;
M_val_L_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_val_L_c_exchange"] = M_val_L_c_exchange_model;
M_val_L_c_exchange_model = 0;

# 229 -1*(M_val_L_c_exchange: M_val_L_c -([])-> [])
M_val_L_c_exchange_reverse_model = FluxModel();
M_val_L_c_exchange_reverse_model.flux_index = 229
M_val_L_c_exchange_reverse_model.flux_symbol = "M_val_L_c_exchange_reverse"
M_val_L_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_val_L_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_val_L_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_val_L_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_val_L_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_val_L_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_val_L_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_val_L_c_exchange_reverse"] = M_val_L_c_exchange_reverse_model;
M_val_L_c_exchange_reverse_model = 0;

# 230 M_o2_c_exchange: M_o2_c -([])-> []
M_o2_c_exchange_model = FluxModel();
M_o2_c_exchange_model.flux_index = 230
M_o2_c_exchange_model.flux_symbol = "M_o2_c_exchange"
M_o2_c_exchange_model.flux_constraint_type = GLPK.DB;
M_o2_c_exchange_model.flux_lower_bound = 0.0;
M_o2_c_exchange_model.flux_upper_bound = 1.0;
M_o2_c_exchange_model.flux_bounds_model = Bounds;
M_o2_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_o2_c_exchange_model.flux_bound_alpha = 1.0;
M_o2_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_o2_c_exchange"] = M_o2_c_exchange_model;
M_o2_c_exchange_model = 0;

# 231 -1*(M_o2_c_exchange: M_o2_c -([])-> [])
M_o2_c_exchange_reverse_model = FluxModel();
M_o2_c_exchange_reverse_model.flux_index = 231
M_o2_c_exchange_reverse_model.flux_symbol = "M_o2_c_exchange_reverse"
M_o2_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_o2_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_o2_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_o2_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_o2_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_o2_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_o2_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_o2_c_exchange_reverse"] = M_o2_c_exchange_reverse_model;
M_o2_c_exchange_reverse_model = 0;

# 232 M_co2_c_exchange: M_co2_c -([])-> []
M_co2_c_exchange_model = FluxModel();
M_co2_c_exchange_model.flux_index = 232
M_co2_c_exchange_model.flux_symbol = "M_co2_c_exchange"
M_co2_c_exchange_model.flux_constraint_type = GLPK.DB;
M_co2_c_exchange_model.flux_lower_bound = 0.0;
M_co2_c_exchange_model.flux_upper_bound = 1.0;
M_co2_c_exchange_model.flux_bounds_model = Bounds;
M_co2_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_co2_c_exchange_model.flux_bound_alpha = 1.0;
M_co2_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_co2_c_exchange"] = M_co2_c_exchange_model;
M_co2_c_exchange_model = 0;

# 233 -1*(M_co2_c_exchange: M_co2_c -([])-> [])
M_co2_c_exchange_reverse_model = FluxModel();
M_co2_c_exchange_reverse_model.flux_index = 233
M_co2_c_exchange_reverse_model.flux_symbol = "M_co2_c_exchange_reverse"
M_co2_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_co2_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_co2_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_co2_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_co2_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_co2_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_co2_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_co2_c_exchange_reverse"] = M_co2_c_exchange_reverse_model;
M_co2_c_exchange_reverse_model = 0;

# 234 M_h2s_c_exchange: M_h2s_c -([])-> []
M_h2s_c_exchange_model = FluxModel();
M_h2s_c_exchange_model.flux_index = 234
M_h2s_c_exchange_model.flux_symbol = "M_h2s_c_exchange"
M_h2s_c_exchange_model.flux_constraint_type = GLPK.DB;
M_h2s_c_exchange_model.flux_lower_bound = 0.0;
M_h2s_c_exchange_model.flux_upper_bound = 1.0;
M_h2s_c_exchange_model.flux_bounds_model = Bounds;
M_h2s_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_h2s_c_exchange_model.flux_bound_alpha = 1.0;
M_h2s_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_h2s_c_exchange"] = M_h2s_c_exchange_model;
M_h2s_c_exchange_model = 0;

# 235 -1*(M_h2s_c_exchange: M_h2s_c -([])-> [])
M_h2s_c_exchange_reverse_model = FluxModel();
M_h2s_c_exchange_reverse_model.flux_index = 235
M_h2s_c_exchange_reverse_model.flux_symbol = "M_h2s_c_exchange_reverse"
M_h2s_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_h2s_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_h2s_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_h2s_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_h2s_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_h2s_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_h2s_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_h2s_c_exchange_reverse"] = M_h2s_c_exchange_reverse_model;
M_h2s_c_exchange_reverse_model = 0;

# 236 M_h_c_exchange: M_h_c -([])-> []
M_h_c_exchange_model = FluxModel();
M_h_c_exchange_model.flux_index = 236
M_h_c_exchange_model.flux_symbol = "M_h_c_exchange"
M_h_c_exchange_model.flux_constraint_type = GLPK.DB;
M_h_c_exchange_model.flux_lower_bound = 0.0;
M_h_c_exchange_model.flux_upper_bound = 1.0;
M_h_c_exchange_model.flux_bounds_model = Bounds;
M_h_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_h_c_exchange_model.flux_bound_alpha = 1.0;
M_h_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_h_c_exchange"] = M_h_c_exchange_model;
M_h_c_exchange_model = 0;

# 237 -1*(M_h_c_exchange: M_h_c -([])-> [])
M_h_c_exchange_reverse_model = FluxModel();
M_h_c_exchange_reverse_model.flux_index = 237
M_h_c_exchange_reverse_model.flux_symbol = "M_h_c_exchange_reverse"
M_h_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_h_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_h_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_h_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_h_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_h_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_h_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_h_c_exchange_reverse"] = M_h_c_exchange_reverse_model;
M_h_c_exchange_reverse_model = 0;

# 238 M_h2o_c_exchange: M_h2o_c -([])-> []
M_h2o_c_exchange_model = FluxModel();
M_h2o_c_exchange_model.flux_index = 238
M_h2o_c_exchange_model.flux_symbol = "M_h2o_c_exchange"
M_h2o_c_exchange_model.flux_constraint_type = GLPK.DB;
M_h2o_c_exchange_model.flux_lower_bound = 0.0;
M_h2o_c_exchange_model.flux_upper_bound = 1.0;
M_h2o_c_exchange_model.flux_bounds_model = Bounds;
M_h2o_c_exchange_model.flux_gamma_array = vec([0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_h2o_c_exchange_model.flux_bound_alpha = 1.0;
M_h2o_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_h2o_c_exchange"] = M_h2o_c_exchange_model;
M_h2o_c_exchange_model = 0;

# 239 -1*(M_h2o_c_exchange: M_h2o_c -([])-> [])
M_h2o_c_exchange_reverse_model = FluxModel();
M_h2o_c_exchange_reverse_model.flux_index = 239
M_h2o_c_exchange_reverse_model.flux_symbol = "M_h2o_c_exchange_reverse"
M_h2o_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_h2o_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_h2o_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_h2o_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_h2o_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_h2o_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_h2o_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_h2o_c_exchange_reverse"] = M_h2o_c_exchange_reverse_model;
M_h2o_c_exchange_reverse_model = 0;

# 240 M_h_e_exchange: M_h_e -([])-> M_h_c
M_h_e_exchange_model = FluxModel();
M_h_e_exchange_model.flux_index = 240
M_h_e_exchange_model.flux_symbol = "M_h_e_exchange"
M_h_e_exchange_model.flux_constraint_type = GLPK.DB;
M_h_e_exchange_model.flux_lower_bound = 0.0;
M_h_e_exchange_model.flux_upper_bound = 1.0;
M_h_e_exchange_model.flux_bounds_model = Bounds;
M_h_e_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_h_e_exchange_model.flux_bound_alpha = 1.0;
M_h_e_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_h_e_exchange"] = M_h_e_exchange_model;
M_h_e_exchange_model = 0;

# 241 -1*(M_h_e_exchange: M_h_e -([])-> M_h_c)
M_h_e_exchange_reverse_model = FluxModel();
M_h_e_exchange_reverse_model.flux_index = 241
M_h_e_exchange_reverse_model.flux_symbol = "M_h_e_exchange_reverse"
M_h_e_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_h_e_exchange_reverse_model.flux_lower_bound = 0.0;
M_h_e_exchange_reverse_model.flux_upper_bound = 1.0;
M_h_e_exchange_reverse_model.flux_bounds_model = Bounds;
M_h_e_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_h_e_exchange_reverse_model.flux_bound_alpha = 1.0;
M_h_e_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_h_e_exchange_reverse"] = M_h_e_exchange_reverse_model;
M_h_e_exchange_reverse_model = 0;

# 242 M_nh4_c_exchange: M_nh4_c -([])-> []
M_nh4_c_exchange_model = FluxModel();
M_nh4_c_exchange_model.flux_index = 242
M_nh4_c_exchange_model.flux_symbol = "M_nh4_c_exchange"
M_nh4_c_exchange_model.flux_constraint_type = GLPK.DB;
M_nh4_c_exchange_model.flux_lower_bound = 0.0;
M_nh4_c_exchange_model.flux_upper_bound = 1.0;
M_nh4_c_exchange_model.flux_bounds_model = Bounds;
M_nh4_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_nh4_c_exchange_model.flux_bound_alpha = 1.0;
M_nh4_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_nh4_c_exchange"] = M_nh4_c_exchange_model;
M_nh4_c_exchange_model = 0;

# 243 -1*(M_nh4_c_exchange: M_nh4_c -([])-> [])
M_nh4_c_exchange_reverse_model = FluxModel();
M_nh4_c_exchange_reverse_model.flux_index = 243
M_nh4_c_exchange_reverse_model.flux_symbol = "M_nh4_c_exchange_reverse"
M_nh4_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_nh4_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_nh4_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_nh4_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_nh4_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_nh4_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_nh4_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_nh4_c_exchange_reverse"] = M_nh4_c_exchange_reverse_model;
M_nh4_c_exchange_reverse_model = 0;

# 244 M_hco3_c_exchange: M_hco3_c -([])-> []
M_hco3_c_exchange_model = FluxModel();
M_hco3_c_exchange_model.flux_index = 244
M_hco3_c_exchange_model.flux_symbol = "M_hco3_c_exchange"
M_hco3_c_exchange_model.flux_constraint_type = GLPK.DB;
M_hco3_c_exchange_model.flux_lower_bound = 0.0;
M_hco3_c_exchange_model.flux_upper_bound = 1.0;
M_hco3_c_exchange_model.flux_bounds_model = Bounds;
M_hco3_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_hco3_c_exchange_model.flux_bound_alpha = 1.0;
M_hco3_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_hco3_c_exchange"] = M_hco3_c_exchange_model;
M_hco3_c_exchange_model = 0;

# 245 -1*(M_hco3_c_exchange: M_hco3_c -([])-> [])
M_hco3_c_exchange_reverse_model = FluxModel();
M_hco3_c_exchange_reverse_model.flux_index = 245
M_hco3_c_exchange_reverse_model.flux_symbol = "M_hco3_c_exchange_reverse"
M_hco3_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_hco3_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_hco3_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_hco3_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_hco3_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_hco3_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_hco3_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_hco3_c_exchange_reverse"] = M_hco3_c_exchange_reverse_model;
M_hco3_c_exchange_reverse_model = 0;

# 246 M_pi_c_exchange: M_pi_c -([])-> []
M_pi_c_exchange_model = FluxModel();
M_pi_c_exchange_model.flux_index = 246
M_pi_c_exchange_model.flux_symbol = "M_pi_c_exchange"
M_pi_c_exchange_model.flux_constraint_type = GLPK.DB;
M_pi_c_exchange_model.flux_lower_bound = 0.0;
M_pi_c_exchange_model.flux_upper_bound = 1.0;
M_pi_c_exchange_model.flux_bounds_model = Bounds;
M_pi_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_pi_c_exchange_model.flux_bound_alpha = 1.0;
M_pi_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_pi_c_exchange"] = M_pi_c_exchange_model;
M_pi_c_exchange_model = 0;

# 247 -1*(M_pi_c_exchange: M_pi_c -([])-> [])
M_pi_c_exchange_reverse_model = FluxModel();
M_pi_c_exchange_reverse_model.flux_index = 247
M_pi_c_exchange_reverse_model.flux_symbol = "M_pi_c_exchange_reverse"
M_pi_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_pi_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_pi_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_pi_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_pi_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_pi_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_pi_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_pi_c_exchange_reverse"] = M_pi_c_exchange_reverse_model;
M_pi_c_exchange_reverse_model = 0;

# 248 M_maltose_c_exchange: M_maltose_c -([])-> []
M_maltose_c_exchange_model = FluxModel();
M_maltose_c_exchange_model.flux_index = 248
M_maltose_c_exchange_model.flux_symbol = "M_maltose_c_exchange"
M_maltose_c_exchange_model.flux_constraint_type = GLPK.DB;
M_maltose_c_exchange_model.flux_lower_bound = 0.0;
M_maltose_c_exchange_model.flux_upper_bound = 1.0;
M_maltose_c_exchange_model.flux_bounds_model = Bounds;
M_maltose_c_exchange_model.flux_gamma_array = vec([1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_maltose_c_exchange_model.flux_bound_alpha = 1.0;
M_maltose_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_maltose_c_exchange"] = M_maltose_c_exchange_model;
M_maltose_c_exchange_model = 0;

# 249 -1*(M_maltose_c_exchange: M_maltose_c -([])-> [])
M_maltose_c_exchange_reverse_model = FluxModel();
M_maltose_c_exchange_reverse_model.flux_index = 249
M_maltose_c_exchange_reverse_model.flux_symbol = "M_maltose_c_exchange_reverse"
M_maltose_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_maltose_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_maltose_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_maltose_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_maltose_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_maltose_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_maltose_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_maltose_c_exchange_reverse"] = M_maltose_c_exchange_reverse_model;
M_maltose_c_exchange_reverse_model = 0;

# 250 M_glc_D_c_exchange: M_glc_D_c -([])-> []
M_glc_D_c_exchange_model = FluxModel();
M_glc_D_c_exchange_model.flux_index = 250
M_glc_D_c_exchange_model.flux_symbol = "M_glc_D_c_exchange"
M_glc_D_c_exchange_model.flux_constraint_type = GLPK.DB;
M_glc_D_c_exchange_model.flux_lower_bound = 0.0;
M_glc_D_c_exchange_model.flux_upper_bound = 1.0;
M_glc_D_c_exchange_model.flux_bounds_model = Bounds;
M_glc_D_c_exchange_model.flux_gamma_array = vec([0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_glc_D_c_exchange_model.flux_bound_alpha = 1.0;
M_glc_D_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_glc_D_c_exchange"] = M_glc_D_c_exchange_model;
M_glc_D_c_exchange_model = 0;

# 251 -1*(M_glc_D_c_exchange: M_glc_D_c -([])-> [])
M_glc_D_c_exchange_reverse_model = FluxModel();
M_glc_D_c_exchange_reverse_model.flux_index = 251
M_glc_D_c_exchange_reverse_model.flux_symbol = "M_glc_D_c_exchange_reverse"
M_glc_D_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_glc_D_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_glc_D_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_glc_D_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_glc_D_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_glc_D_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_glc_D_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_glc_D_c_exchange_reverse"] = M_glc_D_c_exchange_reverse_model;
M_glc_D_c_exchange_reverse_model = 0;

# 252 M_for_c_exchange: M_for_c -([])-> []
M_for_c_exchange_model = FluxModel();
M_for_c_exchange_model.flux_index = 252
M_for_c_exchange_model.flux_symbol = "M_for_c_exchange"
M_for_c_exchange_model.flux_constraint_type = GLPK.DB;
M_for_c_exchange_model.flux_lower_bound = 0.0;
M_for_c_exchange_model.flux_upper_bound = 1.0;
M_for_c_exchange_model.flux_bounds_model = Bounds;
M_for_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_for_c_exchange_model.flux_bound_alpha = 1.0;
M_for_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_for_c_exchange"] = M_for_c_exchange_model;
M_for_c_exchange_model = 0;

# 253 -1*(M_for_c_exchange: M_for_c -([])-> [])
M_for_c_exchange_reverse_model = FluxModel();
M_for_c_exchange_reverse_model.flux_index = 253
M_for_c_exchange_reverse_model.flux_symbol = "M_for_c_exchange_reverse"
M_for_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_for_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_for_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_for_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_for_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_for_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_for_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_for_c_exchange_reverse"] = M_for_c_exchange_reverse_model;
M_for_c_exchange_reverse_model = 0;

# 254 M_lac_D_c_exchange: M_lac_D_c -([])-> []
M_lac_D_c_exchange_model = FluxModel();
M_lac_D_c_exchange_model.flux_index = 254
M_lac_D_c_exchange_model.flux_symbol = "M_lac_D_c_exchange"
M_lac_D_c_exchange_model.flux_constraint_type = GLPK.DB;
M_lac_D_c_exchange_model.flux_lower_bound = 0.0;
M_lac_D_c_exchange_model.flux_upper_bound = 1.0;
M_lac_D_c_exchange_model.flux_bounds_model = Bounds;
M_lac_D_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_lac_D_c_exchange_model.flux_bound_alpha = 1.0;
M_lac_D_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_lac_D_c_exchange"] = M_lac_D_c_exchange_model;
M_lac_D_c_exchange_model = 0;

# 255 -1*(M_lac_D_c_exchange: M_lac_D_c -([])-> [])
M_lac_D_c_exchange_reverse_model = FluxModel();
M_lac_D_c_exchange_reverse_model.flux_index = 255
M_lac_D_c_exchange_reverse_model.flux_symbol = "M_lac_D_c_exchange_reverse"
M_lac_D_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_lac_D_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_lac_D_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_lac_D_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_lac_D_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_lac_D_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_lac_D_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_lac_D_c_exchange_reverse"] = M_lac_D_c_exchange_reverse_model;
M_lac_D_c_exchange_reverse_model = 0;

# 256 M_ac_c_exchange: M_ac_c -([])-> []
M_ac_c_exchange_model = FluxModel();
M_ac_c_exchange_model.flux_index = 256
M_ac_c_exchange_model.flux_symbol = "M_ac_c_exchange"
M_ac_c_exchange_model.flux_constraint_type = GLPK.DB;
M_ac_c_exchange_model.flux_lower_bound = 0.0;
M_ac_c_exchange_model.flux_upper_bound = 1.0;
M_ac_c_exchange_model.flux_bounds_model = Bounds;
M_ac_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_ac_c_exchange_model.flux_bound_alpha = 1.0;
M_ac_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_ac_c_exchange"] = M_ac_c_exchange_model;
M_ac_c_exchange_model = 0;

# 257 -1*(M_ac_c_exchange: M_ac_c -([])-> [])
M_ac_c_exchange_reverse_model = FluxModel();
M_ac_c_exchange_reverse_model.flux_index = 257
M_ac_c_exchange_reverse_model.flux_symbol = "M_ac_c_exchange_reverse"
M_ac_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_ac_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_ac_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_ac_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_ac_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_ac_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_ac_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_ac_c_exchange_reverse"] = M_ac_c_exchange_reverse_model;
M_ac_c_exchange_reverse_model = 0;

# 258 M_etoh_c_exchange: M_etoh_c -([])-> []
M_etoh_c_exchange_model = FluxModel();
M_etoh_c_exchange_model.flux_index = 258
M_etoh_c_exchange_model.flux_symbol = "M_etoh_c_exchange"
M_etoh_c_exchange_model.flux_constraint_type = GLPK.DB;
M_etoh_c_exchange_model.flux_lower_bound = 0.0;
M_etoh_c_exchange_model.flux_upper_bound = 1.0;
M_etoh_c_exchange_model.flux_bounds_model = Bounds;
M_etoh_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_etoh_c_exchange_model.flux_bound_alpha = 1.0;
M_etoh_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_etoh_c_exchange"] = M_etoh_c_exchange_model;
M_etoh_c_exchange_model = 0;

# 259 -1*(M_etoh_c_exchange: M_etoh_c -([])-> [])
M_etoh_c_exchange_reverse_model = FluxModel();
M_etoh_c_exchange_reverse_model.flux_index = 259
M_etoh_c_exchange_reverse_model.flux_symbol = "M_etoh_c_exchange_reverse"
M_etoh_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_etoh_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_etoh_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_etoh_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_etoh_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_etoh_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_etoh_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_etoh_c_exchange_reverse"] = M_etoh_c_exchange_reverse_model;
M_etoh_c_exchange_reverse_model = 0;

# 260 M_mglx_c_exchange: M_mglx_c -([])-> []
M_mglx_c_exchange_model = FluxModel();
M_mglx_c_exchange_model.flux_index = 260
M_mglx_c_exchange_model.flux_symbol = "M_mglx_c_exchange"
M_mglx_c_exchange_model.flux_constraint_type = GLPK.DB;
M_mglx_c_exchange_model.flux_lower_bound = 0.0;
M_mglx_c_exchange_model.flux_upper_bound = 1.0;
M_mglx_c_exchange_model.flux_bounds_model = Bounds;
M_mglx_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_mglx_c_exchange_model.flux_bound_alpha = 1.0;
M_mglx_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_mglx_c_exchange"] = M_mglx_c_exchange_model;
M_mglx_c_exchange_model = 0;

# 261 -1*(M_mglx_c_exchange: M_mglx_c -([])-> [])
M_mglx_c_exchange_reverse_model = FluxModel();
M_mglx_c_exchange_reverse_model.flux_index = 261
M_mglx_c_exchange_reverse_model.flux_symbol = "M_mglx_c_exchange_reverse"
M_mglx_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_mglx_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_mglx_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_mglx_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_mglx_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_mglx_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_mglx_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_mglx_c_exchange_reverse"] = M_mglx_c_exchange_reverse_model;
M_mglx_c_exchange_reverse_model = 0;

# 262 M_prop_c_exchange: M_prop_c -([])-> []
M_prop_c_exchange_model = FluxModel();
M_prop_c_exchange_model.flux_index = 262
M_prop_c_exchange_model.flux_symbol = "M_prop_c_exchange"
M_prop_c_exchange_model.flux_constraint_type = GLPK.DB;
M_prop_c_exchange_model.flux_lower_bound = 0.0;
M_prop_c_exchange_model.flux_upper_bound = 1.0;
M_prop_c_exchange_model.flux_bounds_model = Bounds;
M_prop_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_prop_c_exchange_model.flux_bound_alpha = 1.0;
M_prop_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_prop_c_exchange"] = M_prop_c_exchange_model;
M_prop_c_exchange_model = 0;

# 263 -1*(M_prop_c_exchange: M_prop_c -([])-> [])
M_prop_c_exchange_reverse_model = FluxModel();
M_prop_c_exchange_reverse_model.flux_index = 263
M_prop_c_exchange_reverse_model.flux_symbol = "M_prop_c_exchange_reverse"
M_prop_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_prop_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_prop_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_prop_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_prop_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_prop_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_prop_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_prop_c_exchange_reverse"] = M_prop_c_exchange_reverse_model;
M_prop_c_exchange_reverse_model = 0;

# 264 M_indole_c_exchange: M_indole_c -([])-> []
M_indole_c_exchange_model = FluxModel();
M_indole_c_exchange_model.flux_index = 264
M_indole_c_exchange_model.flux_symbol = "M_indole_c_exchange"
M_indole_c_exchange_model.flux_constraint_type = GLPK.DB;
M_indole_c_exchange_model.flux_lower_bound = 0.0;
M_indole_c_exchange_model.flux_upper_bound = 1.0;
M_indole_c_exchange_model.flux_bounds_model = Bounds;
M_indole_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_indole_c_exchange_model.flux_bound_alpha = 1.0;
M_indole_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_indole_c_exchange"] = M_indole_c_exchange_model;
M_indole_c_exchange_model = 0;

# 265 -1*(M_indole_c_exchange: M_indole_c -([])-> [])
M_indole_c_exchange_reverse_model = FluxModel();
M_indole_c_exchange_reverse_model.flux_index = 265
M_indole_c_exchange_reverse_model.flux_symbol = "M_indole_c_exchange_reverse"
M_indole_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_indole_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_indole_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_indole_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_indole_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_indole_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_indole_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_indole_c_exchange_reverse"] = M_indole_c_exchange_reverse_model;
M_indole_c_exchange_reverse_model = 0;

# 266 M_h2o2_c_exchange: M_h2o2_c -([])-> []
M_h2o2_c_exchange_model = FluxModel();
M_h2o2_c_exchange_model.flux_index = 266
M_h2o2_c_exchange_model.flux_symbol = "M_h2o2_c_exchange"
M_h2o2_c_exchange_model.flux_constraint_type = GLPK.DB;
M_h2o2_c_exchange_model.flux_lower_bound = 0.0;
M_h2o2_c_exchange_model.flux_upper_bound = 1.0;
M_h2o2_c_exchange_model.flux_bounds_model = Bounds;
M_h2o2_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_h2o2_c_exchange_model.flux_bound_alpha = 1.0;
M_h2o2_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_h2o2_c_exchange"] = M_h2o2_c_exchange_model;
M_h2o2_c_exchange_model = 0;

# 267 -1*(M_h2o2_c_exchange: M_h2o2_c -([])-> [])
M_h2o2_c_exchange_reverse_model = FluxModel();
M_h2o2_c_exchange_reverse_model.flux_index = 267
M_h2o2_c_exchange_reverse_model.flux_symbol = "M_h2o2_c_exchange_reverse"
M_h2o2_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_h2o2_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_h2o2_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_h2o2_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_h2o2_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_h2o2_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_h2o2_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_h2o2_c_exchange_reverse"] = M_h2o2_c_exchange_reverse_model;
M_h2o2_c_exchange_reverse_model = 0;

# 268 M_cadav_c_exchange: M_cadav_c -([])-> []
M_cadav_c_exchange_model = FluxModel();
M_cadav_c_exchange_model.flux_index = 268
M_cadav_c_exchange_model.flux_symbol = "M_cadav_c_exchange"
M_cadav_c_exchange_model.flux_constraint_type = GLPK.DB;
M_cadav_c_exchange_model.flux_lower_bound = 0.0;
M_cadav_c_exchange_model.flux_upper_bound = 1.0;
M_cadav_c_exchange_model.flux_bounds_model = Bounds;
M_cadav_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_cadav_c_exchange_model.flux_bound_alpha = 1.0;
M_cadav_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_cadav_c_exchange"] = M_cadav_c_exchange_model;
M_cadav_c_exchange_model = 0;

# 269 -1*(M_cadav_c_exchange: M_cadav_c -([])-> [])
M_cadav_c_exchange_reverse_model = FluxModel();
M_cadav_c_exchange_reverse_model.flux_index = 269
M_cadav_c_exchange_reverse_model.flux_symbol = "M_cadav_c_exchange_reverse"
M_cadav_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_cadav_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_cadav_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_cadav_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_cadav_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_cadav_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_cadav_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_cadav_c_exchange_reverse"] = M_cadav_c_exchange_reverse_model;
M_cadav_c_exchange_reverse_model = 0;

# 270 M_urea_c_exchange: M_urea_c -([])-> []
M_urea_c_exchange_model = FluxModel();
M_urea_c_exchange_model.flux_index = 270
M_urea_c_exchange_model.flux_symbol = "M_urea_c_exchange"
M_urea_c_exchange_model.flux_constraint_type = GLPK.DB;
M_urea_c_exchange_model.flux_lower_bound = 0.0;
M_urea_c_exchange_model.flux_upper_bound = 1.0;
M_urea_c_exchange_model.flux_bounds_model = Bounds;
M_urea_c_exchange_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_urea_c_exchange_model.flux_bound_alpha = 1.0;
M_urea_c_exchange_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_urea_c_exchange"] = M_urea_c_exchange_model;
M_urea_c_exchange_model = 0;

# 271 -1*(M_urea_c_exchange: M_urea_c -([])-> [])
M_urea_c_exchange_reverse_model = FluxModel();
M_urea_c_exchange_reverse_model.flux_index = 271
M_urea_c_exchange_reverse_model.flux_symbol = "M_urea_c_exchange_reverse"
M_urea_c_exchange_reverse_model.flux_constraint_type = GLPK.DB;
M_urea_c_exchange_reverse_model.flux_lower_bound = 0.0;
M_urea_c_exchange_reverse_model.flux_upper_bound = 1.0;
M_urea_c_exchange_reverse_model.flux_bounds_model = Bounds;
M_urea_c_exchange_reverse_model.flux_gamma_array = vec([0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]);
M_urea_c_exchange_reverse_model.flux_bound_alpha = 1.0;
M_urea_c_exchange_reverse_model.flux_obj_coeff = 0.0;
flux_model_dictionary["M_urea_c_exchange_reverse"] = M_urea_c_exchange_reverse_model;
M_urea_c_exchange_reverse_model = 0;

return flux_model_dictionary;
end
