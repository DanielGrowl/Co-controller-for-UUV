function dxdt = AUV_Dynamic_model(t,State,Ui,model_information)
    %[Mv,Dv,Cv,gv,J,cable_mass_pre,tao_wave,tao_rope] = Get_real_model(task,State,Wave,time_flag,F_rope_flag);
    
    % Rigid-Body System Inertia Matrix 刚体系统惯性矩阵
    
    % Hydrodynamic System Inertia Matrix(Added term) 水动力系统惯性矩阵（新增项）
    
    M = model_information{1,1};
    
    % Hydrodynamic Coriolis-Centripetal Matrix(Added term) 流体动力学科里奥利向心矩阵（新增项）
    
    % 1). Largrangian parameterizations 拉格朗日参数化
    
    C = model_information{3,1};
    
    % Linear Part 线性部分
    
    % Nonlinear Part 非线性部分
    
    D = model_information{2,1};
    
    % -- [Resorting Force] -- 回转力矩
    fg = model_information{4,1};
                               
    % Transformation Matrix of linear velocity 线速度变换矩阵
    
    
    % Transformation Matrix of angular velocity 角速度变换矩阵
    
    % Transformation Matrix (Body -> World) 转换矩阵（艇体->大地）
    J = model_information{5,1};
    
    tao_wave = model_information{7,1};
    
    tao_rope = model_information{8,1};

%     Wave_extend = model_information{9,1};
    Wave_extend = [0;0;0;0;0;0];

    V = State(7:12);
%     dxdt1 = J*(V-Wave_extend);  %这个是考虑洋流对于流速的影响的, Find_AB模型里头要配套改
    dxdt1 = J*V;  %不考虑洋流对于流速的影响的, Find_AB模型里头要配套改
    dxdt2 = inv(M)*(Ui  - C*V - D*V - fg - tao_wave - tao_rope);
    dxdt = [dxdt1; dxdt2];
end