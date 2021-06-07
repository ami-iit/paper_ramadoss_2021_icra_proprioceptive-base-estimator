classdef DebugOutputs   
    properties
        access (1, 1) logical;
        
        traceP (1, 1) double;
        condP (1, 1) double;
        P_pred (:, :) double;
        Q (:, :) double;
        F  (:, :) double;
        H (:, :) double;        
        K (:, :) double;
        R (:, :) double;
        y (:, :) double;
        z (:, :) double;
        deltay (:, :) double;
        deltax (:, :) double;
        
        predicted_acc (3, 1) double;
        
        isPsymmetric (1, 1) logical;
        isPpositivedefinite (1, 1) logical;
        isPSPD (1, 1) logical;
    end    
    
    methods
        function obj = DebugOutputs()
            obj.access = false;
            obj.traceP = 0;
            obj.condP = 0;
            obj.P_pred = [];
            obj.Q = [];
            obj.F = [];
            obj.H = [];
            obj.K = [];
            obj.y = [];
            obj.z = [];
            obj.deltay = [];
            obj.deltax = [];
            obj.predicted_acc = zeros(3, 1);
            obj.isPsymmetric = false;
            obj.isPpositivedefinite = false;
            obj.isPSPD = false;
        end
    end
end

