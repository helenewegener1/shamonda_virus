
cd /Users/amaqve/Documents/Amandas_PhD/Courses/molevol_phylogenomics/project_work/shamonda_virus/helene/out

# remove one isolate that does not do well in M segment RevTrans alignment

seqconverter -i L_aligned_revtrans.fasta --remname revtransexclude_L.txt > SML_L_noout_align_rtex.fasta
seqconverter -i S_aligned_revtrans.fasta --remname revtransexclude_S.txt > SML_S_noout_align_rtex.fasta
seqconverter -i M_aligned_revtrans.fasta --remname revtransexclude_M.txt > SML_M_noout_align_rtex.fasta

