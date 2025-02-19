function [Ko,Po,Kr_optimal,Kr_suboptimal,root1,root2_optimal,root2_suboptimal] = Get_control_gain(A,B2,B1,C1,D12,D11,gamma,Q_lqr,R_lqr)
A = double(A);
B2 = double(B2);
B1 = double(B1);
C1 = double(C1);
D12 = double(D12);
D11 = double(D11);

%LQR
Q = Q_lqr;
R = R_lqr;
[Ko,Po] = lqr(A,B2,Q,R); % P is Riccati Equation Solution
root1 = eig(A+B2*Ko);

%LMI H_infinty
N=length(A);
[X,W,gamma_optimal] = LMID(A,B1,B2,C1,D11,D12,N);
Kr_optimal=W*inv(X);
Kr_suboptimal = LMIC(A,B1,B2,C1,D11,D12,gamma,N);
root2_optimal = eig(A+B2*Kr_optimal);
root2_suboptimal = eig(A+B2*Kr_suboptimal);

% %Raccati H_infinty
% B = [B1 B2];
% D = [D11 D12];
% R=[-eye(6), zeros(6); zeros(6), eye(6)];
% H = [A,B;-Q,-R];
% eig(H);
% X_Raccati=dare(A, B, C'*C, R);
% Kr_Raccati=-B2'*X_Raccati; 
% root2_Raccati = eig(A+B2*Kr_Raccati);

% ctrl_flag = 0;
% obsv_flag = 0;
% ctrlMatrix = ctrb(A, B2);
% if rank(ctrlMatrix) == size(A, 1)
%     ctrl_flag = 1;
% else
%     ctrl_flag = 0;
% end
% obsvMatrix = obsv(A, C1);
% if rank(obsvMatrix) == size(A, 1)
%     obsv_flag = 1;
% else
%     obsv_flag = 0;
% end
% 
% sys = ss(A, B2, C1, D12);
% W1 = makeweight(100, 1.5, 0.5);
% W2 = makeweight(0.1, 1, 10);
% P = augw(sys, W1, W2);
% [Kr_hin, CL, gamma_hin] = hinfsyn(P, 12, 6,gamma);
% 
% sys = ltisys(A, B, C1, D);
% [gamma_hl,Kr_hl] = hinflmi(sys,[6,6],gamma);

end