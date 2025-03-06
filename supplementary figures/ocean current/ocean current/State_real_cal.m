function State = State_real_cal(State,Ui,Mv,Dv,Cv,gv,J,tao_wave,tao_rope,s_time)
    velocity = State(7:12);
    location = State(1:6);
    Mv_inv = inv(Mv);
    velocity_dot = Mv_inv * (Cv * velocity + Dv * velocity + gv) + Mv_inv * (Ui - tao_wave - tao_rope);
    velocity = velocity + velocity_dot * s_time;
    location_dot = J * velocity;
    location = location + location_dot * s_time;
    State = [location;velocity];
end