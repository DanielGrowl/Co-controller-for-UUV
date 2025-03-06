function State_estimated = State_estimated_cal(State,Ui,s_time,m,Cv)
    accelerate = Ui/m;
    Velocity_estimated = State(7:12) + accelerate * s_time;
    Location_estimated = State(1:6) + (Cv * Velocity_estimated) * s_time;
    State_estimated = [Location_estimated; Velocity_estimated];
end