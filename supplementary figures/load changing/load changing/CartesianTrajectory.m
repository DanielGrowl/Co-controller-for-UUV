function traj = CartesianTrajectory(Xstart, Xend, Tf, method, s_time)   %点对点轨迹规划函数

    N = Tf/s_time;
    timegap = Tf / (N - 1);
    traj = cell(1, N);
    [Rstart, pstart] = TransToRp(Xstart);
    [Rend, pend] = TransToRp(Xend);
    
    for i = 1: N
        if method == 3  %即三次曲线
            S = CubicTimeScaling(Tf,timegap * (i - 1));
            s = S(1);
        else
            s = QuinticTimeScaling(Tf,timegap * (i - 1));
        end
        traj{i} ...
        = [Rstart * MatrixExp3(MatrixLog3(Rstart' * Rend) * s), ... %这一行的写法没看懂，可能是个累乘
           pstart + s * (pend - pstart); 0, 0, 0, 1];   %为什么要增广？是因为跟踪控制吗？
    end
    
end