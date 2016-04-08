# Example driver file for the seq analysis library code
include("SeqAnalysisLib.jl")

# Load/write the transcription reactions for tetR -
gene_sequence_deGFP = load_gene_sequence_from_file("./deGFP.gene")
generate_transcription_reaction_for_gene_sequence("deGFP",gene_sequence_deGFP)

# Load/write the translation reactions for invA -
protein_seq_deGFP = load_protein_sequence_from_file("./deGFP.prot")
generate_translation_reactions_for_protein_sequence("deGFP",protein_seq_deGFP)
