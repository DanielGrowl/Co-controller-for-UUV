function [A,B2,B1,C1,D12,D11] = Find_AB(State,Ui,Wave)
    syms x y z ph th ps u v w p q r
    syms u1 u2 u3 u4 u5 u6 d1 d2 d3 d4 d5 d6
    assume([x y z ph th ps u v w p q r],"real")
    assume([u1 u2 u3 u4 u5 u6 d1 d2 d3 d4 d5 d6],"real")
    
    L   = 6.0;    
    g   = 9.8;
    rho = 1000;   
    m = 4825;      % kg
    W = m*g;      % N
    B = m*g;   % N
    Ix = 947;   % kg*m^2
    Iy = 15531;   % kg*m^2
    Iz = 16063;   % kg*m^2
    
    rg = [0, 0, 0.04];    % m
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
    
    % Rigid-Body System Inertia Matrix
    Mrb = [m,     0,     0,     0,  m*zg,     0;
           0,     m,     0, -m*zg,     0,     0;
           0,     0,     m,     0,     0,     0;
           0, -m*zg,     0,    Ix,     0,     0;
        m*zg,     0,     0,     0,    Iy,     0;
           0,     0,     0,     0,     0,    Iz];  
    % Hydrodynamic System Inertia Matrix(Added term)
    Ma  = - [Xud,   0,   0,   0, Xqd,   0;
               0, Yvd,   0, Ypd,   0, Yrd;
               0,   0, Zwd,   0, Zqd,   0;
               0, Kvd,   0, Kpd,   0,   0;
             Mud,   0, Mwd,   0, Mqd,   0;
               0, Nvd,   0,   0,   0, Nrd];
    M = Ma+Mrb;
    
    % Hydrodynamic Coriolis-Centripetal Matrix(Added term)
    Ca = [  0,      0,      0,       0,  -Zwd*w,  Yvd*v;
            0,      0,      0,   Zwd*w,       0, -Xud*u;
            0,      0,      0,  -Yvd*v,   Xud*u,      0;
            0, -Zwd*w,  Yvd*v,       0,  -Nrd*r,  Mqd*q;
        Zwd*w,      0, -Xud*u,   Nrd*r,       0, -Kpd*p;
       -Yvd*v,  Xud*u,      0,  -Mqd*q,    Kpd*p,     0];
    % 1). Largrangian parameterizations
    Crb_L = [ 0,     0,     0,       0,     m*w,    -m*v;
              0,     0,     0,    -m*w,       0,     m*u;
              0,     0,     0,     m*v,    -m*u,       0;
              0,   m*w,  -m*v,       0,    Iz*r,   -Iy*q;
           -m*w,     0,   m*u,   -Iz*r,       0,    Ix*p;
            m*v,  -m*u,     0,    Iy*q,   -Ix*p,       0];
    C = Crb_L + Ca;
    
    % Linear Part
    Dl = - [Xu,   0,  Xw,   0, Xq,  0;
             0,  Yv,   0,  Yp,  0, Yr;
            Zu,   0,  Zw,   0, Zq,  0;
             0,   0,   0,  Kp,  0,  0;
             0,   0,  Mw,   0, Mq,  0;
             0,  Nv,   0,   0,  0, Nr];
    
    % Nonlinear Part
    Dnl = - [Xuu*abs(u),           0,           0,           0,           0,           0;
                      0,  Yvv*abs(v),           0,           0,           0,           0;
                      0,           0,  Zww*abs(w),           0,           0,           0;
                      0,           0,           0,  Kpp*abs(p),           0,           0;
                      0,           0,           0,           0,  Mqq*abs(q),           0;
                      0,           0,           0,           0,           0,  Nrr*abs(r)];
    D = Dl + Dnl;
    
    % -- [Resorting Force]
    g = [ (W - B) * sin(th)           ;
         -(W - B) * cos(th) * sin(ph) ; 
         -(W - B) * cos(th) * cos(ph) ;
          zg * W  * cos(th) * sin(ph) ;
          zg * W  * sin(th)           ; 
                                     0];
                                 
    % Transformation Matrix of linear velocity
    R = [cos(ps)*cos(th), -sin(ps)*cos(ph)+sin(ph)*sin(th)*cos(ps),  sin(ps)*sin(ph)+sin(th)*cos(ps)*cos(ph);
         sin(ps)*cos(th),  cos(ps)*cos(ph)+sin(ph)*sin(th)*sin(ps), -cos(ps)*sin(ph)+sin(th)*sin(ps)*cos(ph);
                -sin(th),                          sin(ph)*cos(th),                          cos(ph)*cos(th)];
    
    % Transformation Matrix of angular velocity
    T = [1,  sin(ph)*tan(th),  cos(ph)*tan(th);
         0,          cos(ph),         -sin(ph);
         0,  sin(ph)/cos(th),  cos(ph)/cos(th)];
    
    % Transformation Matrix (Body -> World)
    J = [         R,  zeros(3);
           zeros(3),         T];
    
    Observe_m = [ eye(6), zeros(6);
          zeros(6), zeros(6)];
    
    Wave = [Wave;0;0;0];
    tao_wave_var1 = 0.38 .* 1/2 .* 1050 .* [pi/4;1.5;pi/4;0;0;0];
    tao_wave_var2 = inv(J);
    
    Real_model_var = inv(M);
    
    U = [u1;u2;u3;u4;u5;u6];
    du = [d1;d2;d3;d4;d5;d6];
    X = [x; y; z; ph; th; ps; u; v; w;  p;  q;  r];
    V = [u; v; w;  p;  q;  r];
    
    f1 = [J*V ;Real_model_var * (U - tao_wave_var1.^(tao_wave_var2 * du).^2 - C*V - D*V - g)]; 
    f2 = Observe_m * X;
    
    A = jacobian(f1,X);
    B = jacobian(f1,U);
    B1 = jacobian(f1,du);
    C1 = jacobian(f2,X);
    D11 = jacobian(f2,du);
    D12 = jacobian(f2,U);
    
    
    u01 = Ui(1);
    u02 = Ui(2);
    u03 = Ui(3);
    u04 = Ui(4);
    u05 = Ui(5);
    u06 = Ui(6);
    
    d01 = Wave(1);
    d02 = Wave(2);
    d03 = Wave(3);
    d04 = Wave(4);
    d05 = Wave(5);
    d06 = Wave(6);
    
    x01 = State(1);
    x02 = State(2);
    x03 = State(3);
    x04 = State(4);
    x05 = State(5);
    x06 = State(6);
    x07 = State(7);
    x08 = State(8);
    x09 = State(9);
    x10 = State(10);
    x11 = State(11);
    x12 = State(12);
    
    Equilibrium_point = [x01, x02, x03, x04, x05, x06, x07, x08, x09, x10, x11, x12, u01, u02, u03, u04, u05, u06, d01, d02, d03, d04, d05, d06];
    
    A = subs(A,{x y z ph th ps u v w p q r u1 u2 u3 u4 u5 u6 d1 d2 d3 d4 d5 d6},Equilibrium_point);
    B2 = subs(B,{x y z ph th ps u v w p q r u1 u2 u3 u4 u5 u6 d1 d2 d3 d4 d5 d6},Equilibrium_point); %B
    B1 = subs(B1,{x y z ph th ps u v w p q r u1 u2 u3 u4 u5 u6 d1 d2 d3 d4 d5 d6},Equilibrium_point); %Bw
    C1 = subs(C1,{x y z ph th ps u v w p q r u1 u2 u3 u4 u5 u6 d1 d2 d3 d4 d5 d6},Equilibrium_point);
    D11 = subs(D11,{x y z ph th ps u v w p q r u1 u2 u3 u4 u5 u6 d1 d2 d3 d4 d5 d6},Equilibrium_point); %Dw
    D12 = subs(D12,{x y z ph th ps u v w p q r u1 u2 u3 u4 u5 u6 d1 d2 d3 d4 d5 d6},Equilibrium_point); %Du
end