classdef PriorsStdDev   
    properties
        imu_orientation (3, 1) double;
        imu_position (3, 1) double;
        imu_linear_velocity (3, 1) double;
        
        base_orientation (3, 1) double;
        base_position (3, 1) double;
        base_linear_velocity (3, 1) double;
        base_angular_velocity (3, 1) double;
        
        left_foot_orientation (3, 1) double;
        left_foot_position (3, 1) double;
        right_foot_orientation (3, 1) double;
        right_foot_position (3, 1) double;
        accel_bias (3, 1) double;
        gyro_bias (3, 1) double;
        forward_kinematics (6, 1) double;
    end
    
    methods
        function obj = PriorsStdDev()
            obj.imu_orientation = [0.0; 0.0; 0.0];
            obj.imu_position = [0.0; 0.0; 0.0];            
            obj.imu_linear_velocity = [0.0; 0.0; 0.0];
            
            obj.base_orientation = [0.0; 0.0; 0.0];
            obj.base_position = [0.0; 0.0; 0.0];            
            obj.base_linear_velocity = [0.0; 0.0; 0.0];
            obj.base_angular_velocity = [0.0; 0.0; 0.0];
            
            obj.left_foot_orientation = [0.0; 0.0; 0.0];
            obj.left_foot_position = [0.0; 0.0; 0.0];
            obj.right_foot_orientation = [0.0; 0.0; 0.0];
            obj.right_foot_position = [0.0; 0.0; 0.0];
            obj.accel_bias = [0.0; 0.0; 0.0];            
            obj.gyro_bias = [0.0; 0.0; 0.0];
            obj.forward_kinematics = [0.0; 0.0; 0.0; 0.0; 0.0; 0.0];
        end
    end
end