function [tao] = SMC(State,Error,model_information,switch_gain,R_SMC)
    M = model_information{1,1};
    C = model_information{3,1};    
    D = model_information{2,1};
    fg = model_information{4,1};                               
    J = model_information{5,1};   
    tao_wave = model_information{7,1}; 
    tao_rope = model_information{8,1};

    yita =State(1:6);
    yita_error = Error(1:6);
    v = State(7:12);
    R = R_SMC;
    Slide_surface = yita_error + R * v;
    tao = C * v + D * v + fg - M/R * J * v - M/R * switch_gain * sign(Slide_surface);
end