function s = QuinticTimeScaling(Tf, t) %这个应该是个还没有写完的五次曲线的规划算法，可以避免速度突变
    s = 10 * (t / Tf) ^ 3 - 15 * (t / Tf) ^ 4 + 6 * (t / Tf) ^ 5;
end