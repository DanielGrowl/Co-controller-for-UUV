function dxdt = AUV_Dynamic_model(t,State,Ui,model_information)
    
    M = model_information{1,1};
    
    C = model_information{3,1};
    
    D = model_information{2,1};
    
    fg = model_information{4,1};
                               
    J = model_information{5,1};
    
    tao_wave = model_information{7,1};
    
    tao_rope = model_information{8,1};

    Wave_extend = [0;0;0;0;0;0];

    V = State(7:12);
    dxdt1 = J*V;
    dxdt2 = inv(M)*(Ui  - C*V - D*V - fg - tao_wave - tao_rope);
    dxdt = [dxdt1; dxdt2];
end