The aim of this study is to examine the acceptance for reconstructing the muon and containing the hadronic system in two cases in Near-Detector Liquid Argon (ND-LAr). In the first
case, where all the modules of the the detector are active, the resulting plot is very similar to the corresponding one in the ND Conceptual Design Report (ND-CDR). In the second case,
where the module in the center of the detector is inactive, the resulting plot shows that the change in relative acceptance is uniform as a function of event kinematics.

How to execute the files:

- source ndcaf_setup,sh
- sh ./build_eff.sh
- sh ./pack.sh

and then, in a new terminal and inside a new empty directory (without source-ing anything), we run:

/PATH/TO/run everything.sh FHC 1000

For the creation of the plots, the program efficiency pe.cpp was used. It can be called as follows:

/PATH/TO/efficiency_pe filelist.txt

While, in this study, it was chosen to examine the cases where the inactive module is set at the center of the detector and at the far left side of its central row, the reader can easily
study any other case. To be more specific, the user can change the coordinates of the inactive module, which results in changing its location.
