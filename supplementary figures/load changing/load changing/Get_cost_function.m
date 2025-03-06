function [J_LQR,J_SMC,J_COC] = Get_cost_function(Error_lqr, Error_SMC, Error_COC, Ui_lqr, Ui_SMC, Ui_COC,Q,R)
    Q = 200 * Q;
    J_LQR = Error_lqr'*Q*Error_lqr + Ui_lqr'*R*Ui_lqr; 
    J_SMC = Error_SMC'*Q*Error_SMC + Ui_SMC'*R*Ui_SMC; 
    J_COC= Error_COC'*Q*Error_COC + Ui_COC'*R*Ui_COC;
end