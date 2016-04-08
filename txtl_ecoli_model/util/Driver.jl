# Example driver file for the seq analysis library code
include("SeqAnalysisLib.jl")

# Load/write the transcription reactions for rfp -
gene_sequence_rfp = load_gene_sequence_from_file("./rfp.gene")
generate_transcription_reaction_for_gene_sequence("rfp",gene_sequence_rfp)

# Load/write the translation reactions for rfp -
protein_seq_rfp = load_protein_sequence_from_file("./rfp.prot")
generate_translation_reactions_for_protein_sequence("rfp",protein_seq_rfp)

# Load/write the transcription reactions for sicA -
gene_sequence_sicA = load_gene_sequence_from_file("./sicA.gene")
generate_transcription_reaction_for_gene_sequence("sicA",gene_sequence_sicA)

# Load/write the translation reactions for sicA -
protein_seq_sicA = load_protein_sequence_from_file("./sicA.prot")
generate_translation_reactions_for_protein_sequence("sicA",protein_seq_sicA)

# Load/write the transcription reactions for invA -
gene_sequence_invA = load_gene_sequence_from_file("./invA.gene")
generate_transcription_reaction_for_gene_sequence("invA",gene_sequence_invA)

# Load/write the translation reactions for invA -
protein_seq_invA = load_protein_sequence_from_file("./invA.prot")
generate_translation_reactions_for_protein_sequence("invA",protein_seq_invA)

# Load/write the transcription reactions for tetR -
gene_sequence_tetR = load_gene_sequence_from_file("./tetR.gene")
generate_transcription_reaction_for_gene_sequence("tetR",gene_sequence_tetR)

# Load/write the translation reactions for invA -
protein_seq_tetR = load_protein_sequence_from_file("./tetR.prot")
generate_translation_reactions_for_protein_sequence("tetR",protein_seq_tetR)
