     seqfile = seqfile_replaced          * gets replaced in the script

    treefile = treefile_replaced

     outfile = outfile_replaced

 

       noisy = 3   * 0,1,2,3,9: how much rubbish on the screen

     verbose = 1   * 1: detailed output, 0: concise output

     seqtype = 1   * 1:codons; 2:AAs; 3:codons-->AAs

     ndata = 1     * 1:one gene alignment

   CodonFreq = 7   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table

   estFreq = 0     * used observed frequencies to calculate fitness/freq pairs

       clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree

       model = 2   * models for codons: 0:one, 1:b, 2:2 or more dN/dS ratios for branches

     NSsites = 0   * 0:one w; 1:NearlyNeutral; 2:PositiveSelection; 3:discrete; 4:freqs; 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;10:3normal

       icode = 0   * 0:standard genetic code; 1:mammalian mt; 2-10:see below

   fix_omega = 0   * 1: omega or omega_1 fixed, 0: estimate              

       omega = 1   * initial or fixed omega, for codons or codon-based AAs**    

   cleandata = 0       * remove sites with ambiguity data (1:yes, 0:no)?


