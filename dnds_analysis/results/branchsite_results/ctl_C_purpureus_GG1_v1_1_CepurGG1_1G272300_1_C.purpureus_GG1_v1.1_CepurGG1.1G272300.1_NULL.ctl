seqfile  = CesA_transx_alignment_aliview_edit2.phy
treefile = RAxML_bestTree.cesa_alignment_trimmed.output
outfile  = branch-site_null_out.txt

noisy    = 3
verbose  = 1
seqtype  = 1
CodonFreq= 7
clock    = 0

model    = 2
NSsites  = 2
fix_omega= 1          * fix omega2 at 1
omega    = 1          * fixed

cleandata= 0
