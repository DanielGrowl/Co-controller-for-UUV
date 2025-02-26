function Tasktime = Timeallocation(Initial_x, Desired_x, unit_time)
    Euldistance = norm(Desired_x-Initial_x);
    Euldistance = ceil(Euldistance);
    Tasktime = unit_time * Euldistance;
end