function R = MatrixExp3(so3mat)
    omgtheta = so3ToVec(so3mat);
    if NearZero(norm(omgtheta))
        R = eye(3);
    else
        [omghat, theta] = AxisAng3(omgtheta);
        omgmat = so3mat;
        R = cos(theta) *eye(3) + sin(theta) * omgmat + (1 - cos(theta)) * omghat * omghat';
    end
end