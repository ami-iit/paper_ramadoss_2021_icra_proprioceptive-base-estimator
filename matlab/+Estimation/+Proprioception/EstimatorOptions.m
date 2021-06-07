classdef EstimatorOptions
   
    properties
        nr_joints_est (1, 1) uint32;
        enable_ekf_update (1, 1) logical;
        enable_bias_estimation (1, 1) logical;
        enable_kinematic_meas (1, 1) logical;
        static_bias_initialization (1, 1) logical;
        debug_mode (1, 1) logical;
        
        enable_landmark_measurements (1, 1) logical;
        enable_static_landmarks (1, 1) logical;
    end
    
    methods
        function obj = EstimatorOptions() 
            obj.nr_joints_est = 0;
            obj.enable_ekf_update = true;
            obj.enable_bias_estimation = false;
            obj.enable_kinematic_meas = true;
            obj.static_bias_initialization = false;
            obj.debug_mode = false;
            
            obj.enable_landmark_measurements = false;
            obj.enable_static_landmarks = false;
        end               
    end
end

