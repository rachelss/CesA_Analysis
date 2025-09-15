seqfile  = CesA_transx_alignment_aliview_edit2.phy
treefile = RAxML_bestTree.cesa_alignment_trimmed.output
outfile  = branch-site_alt_out.txt

noisy    = 3
verbose  = 1
seqtype  = 1
CodonFreq= 7
clock    = 0

model    = 2          * branch (foreground/background)
NSsites  = 2          * site class with positive selection allowed (Model A)
fix_omega= 0          * estimate omega2 > 1
omega    = 1.5        * starting value

cleandata= 0
