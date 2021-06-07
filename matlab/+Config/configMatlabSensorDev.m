function sensors_dev = configMatlabSensorDev(nr_joints_est)
sensors_dev = Estimation.Proprioception.SensorsStdDev;

sensors_dev.accel_noise = [0.09 0.09 0.09]'; %acc_std;
sensors_dev.gyro_noise = [0.01 0.01 0.01]'; %gyro_std
sensors_dev.accel_bias_noise = [0.01 0.01 0.01]';
sensors_dev.gyro_bias_noise = [0.001 0.001 0.001]';

sensors_dev.contact_foot_linvel_noise = [9e-3 9e-3 9e-3]'; % fix this
sensors_dev.contact_foot_angvel_noise = [4e-3 4e-3 4e-3]'; % fix this
sensors_dev.swing_foot_linvel_noise = [0.15 0.15 0.15]'; % fix this
sensors_dev.swing_foot_angvel_noise = [0.05 0.05 0.05]'; % fix this
sensors_dev.encoders_noise = deg2rad(0.1)*ones(nr_joints_est, 1);

end

