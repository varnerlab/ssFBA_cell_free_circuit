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
function Control(t,x,data_dictionary)
# ---------------------------------------------------------------------- #
# Control.jl was generated using the Kwatee code generation system.
# Username: jeffreyvarner
# Type: NFBA-JULIA
# Version: 1.0
# Generation timestamp: 04-07-2016 18:07:24
# 
# Arguments: 
# t  - current time 
# x  - state vector 
# data_dictionary  - Data dictionary instance (holds model parameters) 
# ---------------------------------------------------------------------- #

# Set a default value for the allosteric control variables - 
EPSILON = 1.0e-3;
number_of_fluxes = data_dictionary["NUMBER_OF_FLUXES"];
control_vector = ones(number_of_fluxes);
control_parameter_array = data_dictionary["CONTROL_PARAMETER_ARRAY"];

# Alias the species vector - 
M_maltose_c = x[1];
M_h2o_c = x[2];
M_glc_D_c = x[3];
M_pep_c = x[4];
M_g6p_c = x[5];
M_pyr_c = x[6];
M_atp_c = x[7];
M_adp_c = x[8];
M_h_c = x[9];
M_utp_c = x[10];
M_udp_c = x[11];
M_ctp_c = x[12];
M_cdp_c = x[13];
M_gtp_c = x[14];
M_gdp_c = x[15];
M_f6p_c = x[16];
M_fdp_c = x[17];
M_pi_c = x[18];
M_dhap_c = x[19];
M_g3p_c = x[20];
M_nad_c = x[21];
M_13dpg_c = x[22];
M_nadh_c = x[23];
M_3pg_c = x[24];
M_2pg_c = x[25];
M_oaa_c = x[26];
M_co2_c = x[27];
M_coa_c = x[28];
M_accoa_c = x[29];
M_amp_c = x[30];
M_nadp_c = x[31];
M_6pgl_c = x[32];
M_nadph_c = x[33];
M_6pgc_c = x[34];
M_ru5p_D_c = x[35];
M_xu5p_D_c = x[36];
M_r5p_c = x[37];
M_s7p_c = x[38];
M_e4p_c = x[39];
M_2ddg6p_c = x[40];
M_cit_c = x[41];
M_icit_c = x[42];
M_akg_c = x[43];
M_succoa_c = x[44];
M_succ_c = x[45];
M_q8_c = x[46];
M_fum_c = x[47];
M_q8h2_c = x[48];
M_mql8_c = x[49];
M_mqn8_c = x[50];
M_mal_L_c = x[51];
M_glx_c = x[52];
M_actp_c = x[53];
M_ac_c = x[54];
M_ppi_c = x[55];
M_etoh_c = x[56];
M_lac_D_c = x[57];
M_for_c = x[58];
M_o2_c = x[59];
M_h_e = x[60];
M_chor_c = x[61];
M_gln_L_c = x[62];
M_gly_L_c = x[63];
M_gar_c = x[64];
M_glu_L_c = x[65];
M_10fthf_c = x[66];
M_air_c = x[67];
M_thf_c = x[68];
M_asp_L_c = x[69];
M_hco3_c = x[70];
M_aicar_c = x[71];
M_imp_c = x[72];
M_methf_c = x[73];
M_mlthf_c = x[74];
M_5mthf_c = x[75];
M_gmp_c = x[76];
M_ump_c = x[77];
M_cmp_c = x[78];
M_ala_L_c = x[79];
M_arg_L_c = x[80];
M_nh4_c = x[81];
M_asn_L_c = x[82];
M_ser_L_c = x[83];
M_h2s_c = x[84];
M_cys_L_c = x[85];
M_his_L_c = x[86];
M_thr_L_c = x[87];
M_ile_L_c = x[88];
M_leu_L_c = x[89];
M_lys_L_c = x[90];
M_met_L_c = x[91];
M_phe_L_c = x[92];
M_pro_L_c = x[93];
M_trp_L_c = x[94];
M_tyr_L_c = x[95];
M_val_L_c = x[96];
M_urea_c = x[97];
M_h2o2_c = x[98];
M_mglx_c = x[99];
M_prop_c = x[100];
M_indole_c = x[101];
M_cadav_c = x[102];
M_gaba_c = x[103];
GENE_deGFP = x[104];
RNAP = x[105];
OPEN_GENE_deGFP = x[106];
mRNA_deGFP = x[107];
RIBOSOME = x[108];
RIBOSOME_START_deGFP = x[109];
M_ala_L_c_tRNA = x[110];
M_arg_L_c_tRNA = x[111];
M_asn_L_c_tRNA = x[112];
M_asp_L_c_tRNA = x[113];
M_cys_L_c_tRNA = x[114];
M_glu_L_c_tRNA = x[115];
M_gln_L_c_tRNA = x[116];
M_gly_L_c_tRNA = x[117];
M_his_L_c_tRNA = x[118];
M_ile_L_c_tRNA = x[119];
M_leu_L_c_tRNA = x[120];
M_lys_L_c_tRNA = x[121];
M_met_L_c_tRNA = x[122];
M_phe_L_c_tRNA = x[123];
M_pro_L_c_tRNA = x[124];
M_ser_L_c_tRNA = x[125];
M_thr_L_c_tRNA = x[126];
M_trp_L_c_tRNA = x[127];
M_tyr_L_c_tRNA = x[128];
M_val_L_c_tRNA = x[129];
PROTEIN_deGFP = x[130];
tRNA = x[131];

# Build the control vector - 

# Return the control vector - 
return control_vector;
end
