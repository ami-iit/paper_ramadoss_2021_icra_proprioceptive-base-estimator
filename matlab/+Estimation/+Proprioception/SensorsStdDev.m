classdef SensorsStdDev   
    properties
        accel_noise (3, 1) double;
        gyro_noise (3, 1) double;
        accel_bias_noise  (3, 1) double;
        gyro_bias_noise  (3, 1) double;
        
        contact_foot_linvel_noise (3, 1) double;
        contact_foot_angvel_noise (3, 1) double;
        swing_foot_linvel_noise (3, 1) double;
        swing_foot_angvel_noise (3, 1) double;
        
        base_linvel_noise (3, 1) double;
        base_angvel_noise (3, 1) double;
        encoders_noise (:, 1) double;
        velocity_kinematics_noise (6, 1) double;
        
        model_noise (3, 1) double;               
        landmarks_noise (3, 1) double;
    end
    
    methods
        function obj = SensorsStdDev()
            obj.accel_noise = [0.0; 0.0; 0.0];
            obj.gyro_noise = [0.0; 0.0; 0.0];
            obj.accel_bias_noise = [0.0; 0.0; 0.0];
            obj.gyro_bias_noise = [0.0; 0.0; 0.0];
            obj.contact_foot_linvel_noise = [0.0; 0.0; 0.0];
            obj.swing_foot_linvel_noise = [0.0; 0.0; 0.0];            
            obj.contact_foot_angvel_noise = [0.0; 0.0; 0.0];
            obj.swing_foot_angvel_noise = [0.0; 0.0; 0.0];
            obj.base_linvel_noise = [0.0; 0.0; 0.0];
            obj.base_angvel_noise = [0.0; 0.0; 0.0];
            obj.landmarks_noise = [0; 0; 0];
            obj.model_noise = [0; 0; 0];
            obj.velocity_kinematics_noise = zeros(6, 1);
        end
    end
end

