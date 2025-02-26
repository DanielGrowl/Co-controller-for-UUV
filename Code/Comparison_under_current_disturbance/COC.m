function [Ui_co,J_COC_State,L_State,L_det,K] = COC(Ui_COC,model_information,COC_flag,Ko,P_lqr,State_COC,model_information_COC,L_det_state,Q,R,L_standard,switch_gain,R_SMC)
%     L_standard = 0.01;

    A = model_information{1,1};
    B2 = model_information{2,1};    
    B1 = model_information{3,1};
    C1 = model_information{4,1};                               
    D12 = model_information{5,1};   
    D11 = model_information{6,1}; 
    
    Error = COC_flag{1,1};
    Error_pre = COC_flag{2,1};
    J_COC_preState = COC_flag{3,1};
    L_preState = COC_flag{4,1};

    J_COC=Error_pre'*Q*Error_pre + Ui_COC'*R*Ui_COC + Error'*P_lqr*Error;
    L1 = -2*(R*Ui_COC+B2'*P_lqr*Error);
    L2 = inv(R+B2'*P_lqr*B2)*(P_lqr*B2)';
    L = L1'*L2;
    
    if size(L_preState,1)>11 && size(J_COC_preState,1)>11
        J_COC_State = [J_COC_preState; J_COC];
        J_COC_State(1,:) = [];
        L_State = [L_preState; L];
        L_State(1,:) = [];

        L_det = det(L_State);
        L_det_state = [L_det_state, L_det];
        Flag_xr = 0;
        if (L_det < L_standard)
             Error_xr = 0.*Error;
             Flag_xr = 1;
         else
             Error_xr = pinv(L_State)*J_COC_State;
         end
         for j = 1:size(Error_xr,1)
             if abs(Error_xr(j)) > abs(Error(j))
                    Error_xr(j) = Error(j);
             end
         end
    else
        Flag_xr = 1;
        L_det = 0;
        J_COC_State = [J_COC_preState; J_COC];
        L_State = [L_preState; L];
        Error_xr = Error;
    end
    K = Error_xr./Error;
    Ks_r = diag(Error_xr./Error);
    Error_xo = Error - Error_xr;
    Ks_o = diag(Error_xo./Error);
    if Flag_xr == 1
        Ui_xr = zeros(6,1);
        Error_xo = Error;
    else
    Ui_xr = SMC(State_COC,Error_xr,model_information_COC,switch_gain,R_SMC);
    end
    Ui_xo = -Ko*Error_xo;
    Ui_co = Ui_xr + Ui_xo;
end