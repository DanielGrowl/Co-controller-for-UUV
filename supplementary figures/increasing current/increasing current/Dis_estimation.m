function [Error_estimated, Disturbance_estimated,P] = Dis_estimation(State, State_estimated,L,zeta,Mv_pack,Dv_pack,Cv_pack,gv_pack,J_pack,Ui,s_time,P)
    Velocity = State(7:12);
    Error_estimated = State - State_estimated;

    if Error_estimated > 0
        Gamma = zeta * 1;
    elseif Error_estimated < 0
        Gamma = zeta * -1;
    else
        Gamma = 0;
    end

    Mv = Mv_pack{1,1};
    Mv_inv = inv(Mv);
    Dv = Dv_pack{1,1};
    Cv = Cv_pack{1,1};
    gv = gv_pack{1,1};
    J = J_pack{1,1};
    
    F = Mv_inv * (-gv-Cv*Velocity-Dv*Velocity);
    P_dot = -L * Mv_inv * P - L * (Mv_inv*L*Velocity + F + Mv_inv*Ui);
    P = P + P_dot * s_time;
    Disturbance_estimated = L*Velocity + P;
end