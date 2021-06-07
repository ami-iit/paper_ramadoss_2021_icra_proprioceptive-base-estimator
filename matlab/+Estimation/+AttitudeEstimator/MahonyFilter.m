classdef MahonyFilter
properties (SetAccess = private, GetAccess = private)
        timestamp_id (1, 1) uint64;
        quatmahony iDynTree.AttitudeMahonyFilter;
        quatmahonyparams iDynTree.AttitudeMahonyFilterParameters;        
                                
        est_rpy (3, 1) double;
        gyro_bias (3, 1) double;     
        t_prev (1, 1) double;                
                        
        filter_configured (1, 1) logical;
        filter_initialized (1, 1) logical;                
    end
    
    methods
        function obj = MahonyFilter()
            obj.timestamp_id = 0;
            obj.t_prev = 0.0;
            obj.quatmahony = iDynTree.AttitudeMahonyFilter();
            obj.quatmahonyparams = iDynTree.AttitudeMahonyFilterParameters();
            
            obj.est_rpy = zeros(3, 1);
            obj.gyro_bias = zeros(3, 1);
            obj.filter_configured = false;
            obj.filter_initialized = false;
        end
        
        function obj = setup(obj, mahony_params)
            obj.quatmahonyparams = mahony_params;            
            obj.quatmahony.setParameters(mahony_params);
            obj.filter_configured = true;            
        end
                        
        % initialize
        function obj = initialize(obj, x0)            
            xinit = iDynTree.Vector10();            
            xinit.fromMatlab(x0);
            x0span = iDynTree.DynamicSpan(xinit.data, xinit.size);
            
            obj.quatmahony.setInternalState(x0span);
            obj.filter_initialized = true;
        end
        
        % advance
        function [est_rpy, gyro_bias, obj] = advance(obj, t, alpha, omega, mag)
            obj.timestamp_id = obj.timestamp_id + 1;
            
            if (obj.filter_initialized)
                obj.quatmahony.propagateStates();
            end
            
            if (obj.filter_initialized && obj.t_prev > 0.0)
                acc = iDynTree.Vector3();
                acc.fromMatlab(alpha);
    
                gyro = iDynTree.Vector3();
                gyro.fromMatlab(omega);
    
                magmeas = iDynTree.Vector3();
                magmeas.fromMatlab(mag);
                        
                obj.quatmahony.updateFilterWithMeasurements(acc, gyro, magmeas);                
            end
            
            obj.t_prev = t;

            est_rpy_i = iDynTree.Vector3();
            obj.quatmahony.getOrientationEstimateAsRPY(est_rpy_i);
 
            obj.est_rpy = est_rpy_i.toMatlab();
            
            int_state = iDynTree.Vector10();
            xspan = iDynTree.DynamicSpan(int_state.data(), int_state.size());
            obj.quatmahony.getInternalState(xspan);
            int_state_mat = int_state.toMatlab();
            obj.gyro_bias = int_state_mat(5:7);    
            
            est_rpy = obj.est_rpy;
            gyro_bias = obj.gyro_bias;
        end       
        
    end
                                            
 
end

