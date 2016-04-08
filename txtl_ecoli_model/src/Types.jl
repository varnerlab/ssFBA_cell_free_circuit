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

# Define the custom species model type - 
type SpeciesModel

	# Model instance variables - 
	species_index::Int
	species_symbol::AbstractString
	species_lower_bound::Float64
	species_upper_bound::Float64
	is_species_measured::Bool
	species_measurement_array::Array{Float64,2}
	species_constraint_type::Int32
	species_initial_condition::Float64
	is_biomass_precursor::Bool
	biomass_precursor_coefficient::Float64
	species_time_constant::Float64
	is_species_diluted::Bool
	is_species_extracellular::Bool

	# Constructor - 
	function SpeciesModel()
		this = new();
	end
end

# Define the custom flux model type - 
type FluxModel

	# Model instance variables - 
	flux_index::Int
	flux_symbol::AbstractString
	flux_lower_bound::Float64
	flux_upper_bound::Float64
	flux_constraint_type::Int32
	flux_gamma_array::Array{Float64,1}
	flux_bound_alpha::Float64
	flux_bounds_model::Function
	flux_obj_coeff::Float64

	# Constructor - 
	function FluxModel()
		this = new();
	end
end
