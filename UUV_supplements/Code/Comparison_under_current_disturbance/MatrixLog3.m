function so3mat = MatrixLog3(R)
    cosinput = (trace(R) - 1) / 2;
    if cosinput >= 1
        so3mat = zeros(3);
    elseif cosinput <= -1
        if ~NearZero(1 + R(3, 3))
            omg = (1 / sqrt(2 * (1 + R(3, 3)))) ...
                  * [R(1, 3); R(2, 3); 1 + R(3, 3)];
        elseif ~NearZero(1 + R(2, 2))
            omg = (1 / sqrt(2 * (1 + R(2, 2)))) ...
                  * [R(1, 2); 1 + R(2, 2); R(3, 2)];
        else
            omg = (1 / sqrt(2 * (1 + R(1, 1)))) ...
                  * [1 + R(1, 1); R(2, 1); R(3, 1)];
        end
        so3mat = VecToso3(pi * omg);
    else
        theta = acos(cosinput);
        so3mat =  (1 / (2 * sin(theta))) * (R - R');
    end
end