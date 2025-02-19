This text serves as a guide on how to use the supplementary code to verify the results of the paper.

This folder contains two subfolders, one named "Comparison_under_current_disturbance" and the other 
named "Comparison_under_increasing_current_disturbance". The former subfolder corresponds to the 
first subsection "Comparison on control performance" of section "Result" and the latter subfolder 
corresponds to the second subsection "Thrust comparison experiment". These two subfolders are similar 
to each other in most documents, the only different lies on the docmant "AUV_compare.mlx". Therefore, 
in this article we only introduce one supplementary code in detail, and the other supplementary code 
only explains the operation process.

We take the first subfolder "Comparison_under_current_disturbance" as a example. There are 37 
subdocuments in this subfolder. The most important document which presents the result is 
"AUV_compare.mlx", which is also the only document readers need to change if you want to make some
change to test the effectiveness of the algorithm. So we first introduce this document.

There are too many outputs including constants, matrixs and vectors. But many of them are Intermediate
variables, the important ones are "Error_SMC_state", "Error_lqr_state", "Error_COC_state", "F_SMC_state", 
"F_lqr_state", "F_COC_state", "Error_COC_xr_state", "J_SMC_State", "J_lqr_State", "J_COC_State". These 
variables can be seperated into four groups. 
The first group ("Error_SMC_state", "Error_lqr_state", "Error_COC_state") contributes to Fig. 3, which shows 
the state responses of UUV under the control of LQR, SMC and COC. 
The second group ("Error_COC_xr_state") contributes to Fig. 4, which compares the proportion of 
sub-deviations from thruster noise and current disturbance in the total deviation of UUV under current 
disturbance. 
The third group ("F_SMC_state", "F_lqr_state", "F_COC_state") contributes to Fig. 5, which shows thrust 
changing lines of each thruster. 
The last group ("J_SMC_State", "J_lqr_State", "J_COC_State") provides the corresponding data for Figure 6.

To obtain these four groups of variables, readers only need to open "AUV_compare.mlx" in Matlab version 
2022a and above and set the current folder to "Comparison_under_current_disturbance" folder. Then click 
the RUN button, wait a few minutes, the figures will be shown in each subsections.

It is the same as the document "AUV_compare.mlx" of the subfolder 
"Comparison_under_increasing_current_disturbance". To obtain the results, readers only need to open 
"AUV_compare.mlx" in Matlab version 2022a and above and set the current folder to 
"Comparison_under_increasing_current_disturbance" folder. Then click the RUN button, wait a few minutes, 
the figures will be shown in each subsections.
