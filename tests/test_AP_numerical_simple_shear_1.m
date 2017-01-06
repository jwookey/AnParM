function [Xs, Ls, ts] = test_AP_numerical_simple_shear_1


    [Xs, Ls, ts] = AP_path_integrate( [0.0 1.0 0.0]', 0.0, 10.0,...
        1.0, 0.1, @AP_velocity_simple_shear);
    
    Ls
    Xs
    ts
end