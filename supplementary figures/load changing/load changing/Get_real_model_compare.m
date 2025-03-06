function [model_information,m] = Get_real_model_compare(task,State,Wave,time_flag,F_rope_flag)
%相比之前的版本，这一版直接讲模型信息封装到model_information函数里
%一些标志位
%[Tasktime_1,i,M_oc,cable_base_mass,v_cable,g_equivalent,cable_mass_pre,mass_loss]
Tf = time_flag(1);
i=time_flag(2);
M_oc = time_flag(3);
cable_base_mass = time_flag(4);    % 假设光缆基座5kg
v_cable = time_flag(5);
g_equivalent = time_flag(6);
cable_mass_pre = time_flag(7);      % 固定的放缆时间
mass_loss = time_flag(8);       % 

%F_rope_flag = [K;alpha;Sea_density;C_f;d;M_oc]
%高度归一化系数；光缆倾角；海水密度；摩擦力系数；光缆直径；光缆单位质量
K = F_rope_flag(1);
alpha = F_rope_flag(2);
Sea_density = F_rope_flag(3);
C_f = F_rope_flag(4);
d = F_rope_flag(5);
M_oc = F_rope_flag(6);

%  x y z 
 z = State(3);
 ph =State(4);
 th =State(5);
 ps =State(6);
 u  =State(7);
 v =State(8);
 w =State(9);
 p =State(10);
 q =State(11);
 r=State(12);

m = 4825;      % kg 艇体质量

if task == 1    %移动到起始位置
    m = 4825;
elseif task == 2    % 放缆及基座
    if i < (Tf*1/3)
        m = 4825;
    elseif i > (Tf*1/3) && i < (Tf*2/3)
        cable_mass = 1 * (i-fix(Tf*1/3)) * v_cable * M_oc + 1 * v_cable * M_oc * randn;
        m = 4825 - (9.8 - g_equivalent) * ( cable_mass + cable_base_mass );      
        cable_mass_pre = 4825 - m;
    elseif i > (Tf*2/3-1)
        m = 4825 - (9.8 - g_equivalent) * cable_mass_pre - g_equivalent * cable_base_mass;    %1.09是海水质量
    end
elseif task == 3    %下降至基座位置
    cable_mass_unit = cable_mass_pre / Tf;
    m = 4825 - i * g_equivalent * cable_mass_unit - g_equivalent * cable_base_mass;
elseif task == 4    %抓取基座
    if i < (Tf*3/5)
        m = 4825 - g_equivalent * ( cable_mass_pre + cable_base_mass );      
    elseif i > (Tf*3/5-1)
        m = 4825 - g_equivalent * cable_mass_pre + g_equivalent * cable_base_mass;
    end
elseif task == 5    %拖拽缆绳至接驳点，这一块得加拖拽拉力
    m = 4825 - g_equivalent * cable_mass_pre + g_equivalent * cable_base_mass;  
elseif task == 6    %基座安装至接驳点
    if i < (Tf*3/5)
        m = 4825 - g_equivalent * cable_mass_pre + g_equivalent * cable_base_mass;      
    elseif i > (Tf*3/5-1)
        m = 4825 - g_equivalent * cable_mass_pre - g_equivalent * cable_base_mass;
    end
elseif task == 7    % 爬升到拖拽起始位置
    cable_mass = i * v_cable * M_oc;
    m = 4825 - g_equivalent * mass_loss - g_equivalent * cable_base_mass - (9.8-g_equivalent) * cable_mass;
    cable_mass_pre = cable_mass;
elseif task == 8    %开始拖拽
    cable_mass = i * v_cable * M_oc;
    % 开始拖拽,最后一项表示拖拽的时候以alpha角度拉住AUV的绳子质量
    m = 4825 - g_equivalent * mass_loss - g_equivalent * cable_base_mass - g_equivalent * cable_mass + g_equivalent * mass_loss * 3/5 * cos(alpha) ; 
    cable_mass_pre = cable_mass;
elseif task == 9    %释放剩余光缆及基座
    if i < (Tf*4/5)
        cable_mass = i * v_cable * M_oc;
        m = 4825  - g_equivalent * mass_loss - (9.8-g_equivalent) * (cable_mass + cable_base_mass) - g_equivalent * cable_base_mass;
        cable_mass_pre = cable_mass;
    elseif i > (Tf*4/5-1)
        m = 4825 - g_equivalent * mass_loss - (9.8-g_equivalent) * cable_mass_pre - 2 * g_equivalent * cable_base_mass;
    end
elseif task == 10   %下降至基座位置
    cable_mass_unit = cable_mass_pre / Tf;
    m = 4825  - g_equivalent * mass_loss - i * g_equivalent * cable_mass_unit - 2 * g_equivalent * cable_base_mass;
elseif task == 11   %抓取基座
    if i < (Tf*3/5)
        m = 4825 - g_equivalent * mass_loss - g_equivalent * ( cable_mass_pre + 2 * cable_base_mass );      
    elseif i > (Tf*3/5-1)
        m = 4825  - g_equivalent * mass_loss - g_equivalent * ( cable_mass_pre + 2 * cable_base_mass ) + g_equivalent * cable_base_mass;
    end
elseif task == 12   %拖拽光缆至接驳点处
    m = 4825 - g_equivalent * mass_loss - g_equivalent * cable_mass_pre + g_equivalent * cable_base_mass - g_equivalent * cable_base_mass;  
elseif task == 13   %安装基座
    if i < (Tf*3/5)
        m = 4825 - g_equivalent * mass_loss - g_equivalent * cable_mass_pre + g_equivalent * cable_base_mass - g_equivalent * cable_base_mass;      
    elseif i > (Tf*3/5-1)
        m = 4825 - g_equivalent * mass_loss - g_equivalent * cable_mass_pre - 2 * g_equivalent * cable_base_mass;
    end
end

L   = 6.0;    
g   = 9.8;
rho = 1000;  
W = m*g;      % N
B = m*g;   % N
Ix = 947;   % kg*m^2
Iy = 15531;   % kg*m^2
Iz = 16063;   % kg*m^2

rg = [0, 0, 0.04];    % m 这个不是很能理解是个啥
rb = [0, 0, 0]; % m
BG = rg - rb;      % m
xg = BG(1);        % m
yg = BG(2);        % m
zg = BG(3);        % m

r2 = rho*L^2/2;
r3 = rho*L^3/2;
r4 = rho*L^4/2;
r5 = rho*L^5/2;

Xud = -2.996e-3*r3;    % kg
Xqd  = 0.373e-3*r4;    % kg*(m/s^2)/(rad/s^2)=kg*m/rad
Xu = -97.380*100;      % N*s/m
Xw  =  19.645*100;     % N*s/m
Xq = 12.944*1000;      % N*s/rad

Yvd = -45.295e-3*r3;  % kg
Ypd =  0.540e-3*r4;   %kg*m/rad
Yrd =  2.148e-3*r4;   %kg*m/rad
Yv  = -457.996*100;   %N*s/m
Yp  =  16.681*1000;   %N*s/rad 
Yr  =  395.180*1000;  %N*s/rad

Zwd = -58.143e-3*r3;   % kg    
Zqd = -4.323e-3*r4;    %kg*m/rad
Zu  = 47.054*100;      %N*s/m
Zw  = -408.494*100;    %N*s/m   
Zq  = -568.187*1000;   %N*s/rad

Kvd =  0.591e-3*r4;   %kg*m^2/s^2 /(m/s^2) = kg*m
Kpd = -0.068e-3*r5;   %kg*m^2/s^2 /(rad/s^2) = kg*m^2/rad
Kp  = -7.118*10000;   %N*m*s/rad

Mud = 0.385e-3*r4;    %kg*m^2/s^2 /(m/s^2) = kg*m
Mwd = -4.141e-3*r4;   %kg*m^2/s^2 /(m/s^2) = kg*m
Mqd = -3.299e-2*r5;    %kg*m^2/s^2 /(rad/s^2) = kg*m^2/rad 
Mw = -303.345*1000;   %N*s
Mq =-909.683*10000;   %N*m*s/rad

Nvd =  1.952e-3*r4;   %kg*m
Nrd = -2.009e-3*r5;   %kg*m^2/rad
Nv  = 199.295*1000;   %N*s 
Nr  = -618.383*10000; %N*m*s/rad

Xuu = 0; % N*s^2/m^2
Yvv = 0; % N*s^2/m^2
Zww = 0; % N*s^2/m^2
Kpp = 0;  % N*s^2/rad^2
Mqq = 0;  % N*s^2/rad^2
Nrr = 0;  % N*s^2/rad^2

% Rigid-Body System Inertia Matrix 刚体系统惯性矩阵
Mrb = [m,     0,     0,     0,  m*zg,     0;
       0,     m,     0, -m*zg,     0,     0;
       0,     0,     m,     0,     0,     0;
       0, -m*zg,     0,    Ix,     0,     0;
    m*zg,     0,     0,     0,    Iy,     0;
       0,     0,     0,     0,     0,    Iz];  
% Hydrodynamic System Inertia Matrix(Added term) 水动力系统惯性矩阵（新增项）
Ma  = - [Xud,   0,   0,   0, Xqd,   0;
           0, Yvd,   0, Ypd,   0, Yrd;
           0,   0, Zwd,   0, Zqd,   0;
           0, Kvd,   0, Kpd,   0,   0;
         Mud,   0, Mwd,   0, Mqd,   0;
           0, Nvd,   0,   0,   0, Nrd];
M = Ma+Mrb;

% Hydrodynamic Coriolis-Centripetal Matrix(Added term) 流体动力学科里奥利向心矩阵（新增项）
Ca = [  0,      0,      0,       0,  -Zwd*w,  Yvd*v;
        0,      0,      0,   Zwd*w,       0, -Xud*u;
        0,      0,      0,  -Yvd*v,   Xud*u,      0;
        0, -Zwd*w,  Yvd*v,       0,  -Nrd*r,  Mqd*q;
    Zwd*w,      0, -Xud*u,   Nrd*r,       0, -Kpd*p;
   -Yvd*v,  Xud*u,      0,  -Mqd*q,    Kpd*p,     0];
% 1). Largrangian parameterizations 拉格朗日参数化
Crb_L = [ 0,     0,     0,       0,     m*w,    -m*v;
          0,     0,     0,    -m*w,       0,     m*u;
          0,     0,     0,     m*v,    -m*u,       0;
          0,   m*w,  -m*v,       0,    Iz*r,   -Iy*q;
       -m*w,     0,   m*u,   -Iz*r,       0,    Ix*p;
        m*v,  -m*u,     0,    Iy*q,   -Ix*p,       0];
C = Crb_L + Ca;

% Linear Part 线性部分
Dl = - [Xu,   0,  Xw,   0, Xq,  0;
         0,  Yv,   0,  Yp,  0, Yr;
        Zu,   0,  Zw,   0, Zq,  0;
         0,   0,   0,  Kp,  0,  0;
         0,   0,  Mw,   0, Mq,  0;
         0,  Nv,   0,   0,  0, Nr];

% Nonlinear Part 非线性部分
Dnl = - [Xuu*abs(u),           0,           0,           0,           0,           0;
                  0,  Yvv*abs(v),           0,           0,           0,           0;
                  0,           0,  Zww*abs(w),           0,           0,           0;
                  0,           0,           0,  Kpp*abs(p),           0,           0;
                  0,           0,           0,           0,  Mqq*abs(q),           0;
                  0,           0,           0,           0,           0,  Nrr*abs(r)];
D = Dl + Dnl;

% -- [Resorting Force] -- 回转力矩
fg = [ (W - B) * sin(th)           ;
     -(W - B) * cos(th) * sin(ph) ; 
     -(W - B) * cos(th) * cos(ph) ;
      zg * W  * cos(th) * sin(ph) ;
      zg * W  * sin(th)           ; 
                                 0];
                             
% Transformation Matrix of linear velocity 线速度变换矩阵
R = [cos(ps)*cos(th), -sin(ps)*cos(ph)+sin(ph)*sin(th)*cos(ps),  sin(ps)*sin(ph)+sin(th)*cos(ps)*cos(ph);
     sin(ps)*cos(th),  cos(ps)*cos(ph)+sin(ph)*sin(th)*sin(ps), -cos(ps)*sin(ph)+sin(th)*sin(ps)*cos(ph);
            -sin(th),                          sin(ph)*cos(th),                          cos(ph)*cos(th)];

% Transformation Matrix of angular velocity 角速度变换矩阵
T = [1,  sin(ph)*tan(th),  cos(ph)*tan(th);
     0,          cos(ph),         -sin(ph);
     0,  sin(ph)/cos(th),  cos(ph)/cos(th)];

% Transformation Matrix (Body -> World) 转换矩阵（艇体->大地）
J = [         R,  zeros(3);
       zeros(3),         T];

%洋流干扰
%Fd = Cd * 1/2 * rou * v^2 * A
%Cd阻力系数，rou流体密度，v流动速度，A物体前部特征面积
if i < (1)
    Wave = [0;0;0;0;0;0];    %没有考虑洋流的转动力矩影响，之后再加
    Wave_extend = Wave; %由于不加旋转洋流
    tao_wave = 0.38 .* 1/2 .* 1050 .* (inv(J) * Wave).^2 .* [pi/4;1.5;pi/4;0;0;0];
else
    Wave = [Wave;0;0;0];    %没有考虑洋流的转动力矩影响，之后再加
    Wave_extend = Wave; %由于不加旋转洋流
    tao_wave = 0.38 .* 1/2 .* 1050 .* (inv(J) * Wave).^2 .* [pi/4;1.5;pi/4;0;0;0];
end

% tao_wave = zeros(6,1);
%没有考虑前部特征面积在洋流方向的投影面积，之后再改

%光缆耦合干扰
%[1] Yuanhui Wang, Xinqian Bian, Xiaoyun Zhang, and Wenbo Xie, 
% “A study on the influence of cable tension on the movement of cable laying ship,” 
% in OCEANS 2010 MTS/IEEE SEATTLE, Seattle, 
% WA: IEEE, Sep. 2010, pp. 1–8. doi: 10.1109/OCEANS.2010.5664404.
if task == 1
    tao_rope = zeros(6,1);
elseif task == 2
    tao_rope = zeros(6,1);
elseif task == 3
   tao_rope = zeros(6,1);
elseif task == 4
    tao_rope = zeros(6,1);
elseif task == 5
    F_rope = - K * z / sin(alpha) * (1/2 * Sea_density * C_f * pi * d * (u^2 + v^2) + ...
    M_oc * sin(alpha) );
    F_rope_AUV = F_rope .* [0; sin(alpha); cos(alpha); 0; 0; 0];
    tao_rope = F_rope_AUV;
elseif task == 6
    tao_rope = zeros(6,1);
elseif task == 7
    tao_rope = zeros(6,1);
elseif task == 8
    F_rope = - K * z / sin(alpha) * (1/2 * Sea_density * C_f * pi * d * (u^2 + v^2) + ...
    M_oc * sin(alpha) );
    F_rope_AUV = F_rope .* [0; sin(alpha); cos(alpha); 0; 0; 0];
    tao_rope = F_rope_AUV;
elseif task == 9
    tao_rope = zeros(6,1);
elseif task == 10
    tao_rope = zeros(6,1);
elseif task == 11
    tao_rope = zeros(6,1);
elseif task == 12
    F_rope = - K * z / sin(alpha) * (1/2 * Sea_density * C_f * pi * d * (u^2 + v^2) + ...
    M_oc * sin(alpha) );
    F_rope_AUV = F_rope .* [0; sin(alpha); cos(alpha); 0; 0; 0];
    tao_rope = F_rope_AUV;
elseif task == 13
    tao_rope = zeros(6,1);
end

route_k = 0;  %证明比lqr好
% route_k = 0.3;  %证明比lqr好

if i > (Tf*1/2)
    Mv = M + route_k * rot90(M) ;
    Cv = C + route_k * rot90(C) ;
    Dv = D + route_k * rot90(D) ;
    gv = fg;
else
    Mv = M;
    Cv = C;
    Dv = D;
    gv = fg;
end

    model_information{1,1} = Mv;model_information{2,1} = Dv;model_information{3,1} = Cv;model_information{4,1} = gv;
    model_information{5,1} = J;model_information{6,1} = cable_mass_pre;model_information{7,1} = tao_wave;
    model_information{8,1} = tao_rope;model_information{9,1} = Wave_extend;
end