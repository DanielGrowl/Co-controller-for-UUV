function [Ko,Po,Kr_optimal,Kr_suboptimal,root1,root2_optimal,root2_suboptimal] = Get_control_gain(A,B2,B1,C1,D12,D11,gamma,Q_lqr,R_lqr)
A = double(A);
% A = A + 1 * eye(size(A));
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
[X,W,gamma_optimal] = LMID(A,B1,B2,C1,D11,D12,N);     %这个得到的gamma太小了
Kr_optimal=W*inv(X);      %这个得到的gamma太小了
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

% %能观能控性判断
% ctrl_flag = 0;
% obsv_flag = 0;
% ctrlMatrix = ctrb(A, B2); % 计算可控性矩阵
% if rank(ctrlMatrix) == size(A, 1) % 检查可控性矩阵的秩
%     ctrl_flag = 1;
% else
%     ctrl_flag = 0;
% end
% obsvMatrix = obsv(A, C1); % 计算可观性矩阵
% if rank(obsvMatrix) == size(A, 1) % 检查可观性矩阵的秩
%     obsv_flag = 1;
% else
%     obsv_flag = 0;
% end
% 
% %hinfsyn求解,这个也不行
% % 定义状态空间模型
% sys = ss(A, B2, C1, D12);
% % 创建加权函数
% W1 = makeweight(100, 1.5, 0.5); % 误差加权函数
% W2 = makeweight(0.1, 1, 10);    % 控制输入加权函数
% % 扩展系统模型
% P = augw(sys, W1, W2);
% % 设计H∞控制器
% [Kr_hin, CL, gamma_hin] = hinfsyn(P, 12, 6,gamma);
% 
% %hinflmi求解
% % 定义状态空间模型
% sys = ltisys(A, B, C1, D);
% % 设计H∞控制器
% [gamma_hl,Kr_hl] = hinflmi(sys,[6,6],gamma);

end