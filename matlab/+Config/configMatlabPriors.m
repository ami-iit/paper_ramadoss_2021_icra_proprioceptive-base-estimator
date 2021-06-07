function priors = configMatlabPriors()
priors = Estimation.Proprioception.PriorsStdDev;
priors.base_position = [1e-3 1e-3 1e-3];
priors.base_orientation = [deg2rad(0.1) deg2rad(0.1) deg2rad(0.1)];
priors.base_linear_velocity = [1e-2 1e-2 1e-2];
priors.base_angular_velocity = [1e-3 1e-3 1e-3];

% walking
priors_dev.imu_orientation = [deg2rad(10) deg2rad(10) deg2rad(1)]';
priors_dev.imu_position = [0.01 0.01 0.01]';
priors_dev.imu_linear_velocity = [0.5 0.5 0.5]';
priors_dev.left_foot_orientation = [deg2rad(10) deg2rad(10) deg2rad(1)]';
priors_dev.left_foot_position = [0.01 0.01 0.01]';
priors_dev.right_foot_orientation = [deg2rad(10) deg2rad(10) deg2rad(1)]';
priors_dev.right_foot_position = [0.01 0.01 0.01]';
priors_dev.accel_bias = [0.01 0.01 0.01]';
priors_dev.gyro_bias = [0.002 0.002 0.002]';

% coordinator
% priors_dev.imu_orientation = [deg2rad(30) deg2rad(30) deg2rad(30)]';
% priors_dev.imu_position = [0.1 0.1 0.1]';
% priors_dev.imu_linear_velocity = [0.5 0.5 0.5]';
% priors_dev.left_foot_orientation = priors_dev.imu_orientation';
% priors_dev.left_foot_position = priors_dev.imu_position';
% priors_dev.right_foot_orientation = priors_dev.imu_orientation';
% priors_dev.right_foot_position = priors_dev.imu_position';
% priors_dev.accel_bias = [0.01 0.01 0.01]';
% priors_dev.gyro_bias = [0.002 0.002 0.002]';

end

