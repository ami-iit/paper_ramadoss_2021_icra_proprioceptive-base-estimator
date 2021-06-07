function qekfparams = configQEKFParams(dt)
qekfparams = iDynTree.AttitudeQuaternionEKFParameters();

qekfparams.time_step_in_seconds = dt;
qekfparams.use_magnetometer_measurements = false; % turning this to true causes filter to explode

% these values were obtained through Allan variance analysis on waist (xsens imu)
% measurements
% qekfparams.accelerometer_noise_variance = 8e-5;
% qekfparams.gyroscope_noise_variance = 1e-7;
% qekfparams.gyro_bias_noise_variance = 1e-10;
% qekfparams.bias_correlation_time_factor = 0.1;
% qekfparams.magnetometer_noise_variance = 1e-4;

qekfparams.accelerometer_noise_variance = 0.1;
qekfparams.gyroscope_noise_variance = 0.01;
qekfparams.gyro_bias_noise_variance = 0.001;
qekfparams.bias_correlation_time_factor = 0.1;
qekfparams.magnetometer_noise_variance = 1e-4;

qekfparams.initial_orientation_error_variance = deg2rad(10);
qekfparams.initial_ang_vel_error_variance = deg2rad(1);
qekfparams.initial_gyro_bias_error_variance = deg2rad(1);
end

