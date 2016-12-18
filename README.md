### Instructions for reproducing the numerical results in "Robust nonparametric nearest neighbor random process clustering", 2016, by Michael Tschannen and Helmut Bölcskei

Requirements: matlab, make, pdflatex


#### Comparison of NNPC and TSC (see Section 4 in the paper)

Run "synthDataNNPCvsTSC(0.2)" in Matlab (the argument indicates the noise variance). This will produce a figure showing the clustering error as a function of the observation length M for NNPC and TSC clustering observations stemming from three different generative models. The generative model PSDs will be shown in a separate figure.


### Figures 2 and 3:

In terminal, cd to the folder "supp_material_rpc_ext" and type "make fig23.pdf". 


### Table 1:

1. Download the file "amc_to_matrix.m" from "http://mocap.cs.cmu.edu/tools.php" and move it to "robust_rpc".
2. Download the motion capture files listed below from "http://mocap.cs.cmu.edu" and move them to the corresponding folders in "robust_rpc".
3. Run "testMOCAP(1)" and "testMOCAP(2)" in Matlab to reproduce the first and second row of Table 1, respectively.


File list "MOCAP_S16":

16_08.amc
16_11.amc
16_12.amc
16_13.amc
16_14.amc
16_15.amc
16_16.amc
16_17.amc
16_18.amc
16_19.amc
16_20.amc
16_21.amc
16_22.amc
16_23.amc
16_24.amc
16_25.amc
16_26.amc
16_27.amc
16_28.amc
16_29.amc
16_30.amc
16_31.amc
16_32.amc
16_33.amc
16_34.amc
16_35.amc
16_36.amc
16_37.amc
16_38.amc
16_39.amc
16_40.amc
16_41.amc
16_42.amc
16_43.amc
16_44.amc
16_45.amc
16_46.amc
16_47.amc
16_48.amc
16_49.amc
16_50.amc
16_51.amc
16_52.amc
16_53.amc
16_54.amc
16_55.amc
16_56.amc
16_57.amc
16_58.amc


File list "MOCAP_S35":

35_01.amc
35_02.amc
35_03.amc
35_04.amc
35_05.amc
35_06.amc
35_07.amc
35_08.amc
35_09.amc
35_10.amc
35_11.amc
35_12.amc
35_13.amc
35_14.amc
35_15.amc
35_16.amc
35_17.amc
35_18.amc
35_19.amc
35_20.amc
35_21.amc
35_22.amc
35_23.amc
35_24.amc
35_25.amc
35_26.amc
35_28.amc
35_29.amc
35_30.amc
35_31.amc
35_32.amc
35_33.amc
35_34.amc



#### Table 2 and Figure 4:

1. Download Set A ("Z.zip") and Set E ("S.zip") from http://ntsa.upf.edu/downloads/andrzejak-rg-et-al-2001-indications-nonlinear-deterministic-and-finite-dimensional and unpack the zip files in the folder "robust_rpc".
2. In terminal, cd to the folder "robust_rpc" and type "make fig4.pdf". 


