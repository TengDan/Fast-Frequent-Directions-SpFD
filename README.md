# Fast-Frequent-Direction-SpFD

This repository contains the matlab code for the TPAMI paper 'A Fast Frequent Directions Algorithm for Low Rank Approximation' by Dan Teng and Delin Chu. 
1. run_Protein.m is the main script for running the experiments. After running the script, a .mat file will be created which contains all the results, then the figures are plotted by running plot_code_median.m. 
2. The supporting functions: 
   -- Our fast frequent directions method SpFD is in randFreqDirP.m.
   -- The competing algorithms are: DCT (DisDT.m), SpEmb (Sparse.m), norm sampling (normSamp.m), FD (freqDir.m) and Mina's SFD (MSFD.m). 
3. The example attached is dataset Protein.mat.
   -- The pre-runned results are in final_protein_result.mat (without MSFD), to generate figures, run plot_code_median_withoutMSFD.m. 
