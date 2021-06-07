classdef SchmittTrigger    
    properties (SetAccess = private, GetAccess = private)
        %% State
        currentState (1, 1) logical;
        previousTime (1, 1) double;
        timer (1, 1) double;
        
        %% Params
        stableOFFTime (1, 1) double;
        stableONTime (1, 1) double;
        lowValueThreshold (1, 1) double;
        highValueThreshold (1, 1) double;
        
        %% Input
        rawValue (1, 1) double;
        
        %% Verbose
        verbose (1, 1) logical;
    end
    
    methods
        function obj = SchmittTrigger(stableOFFTime, stableONTime, lowValueThreshold, highValueThreshold)
            arguments              
              stableOFFTime (1, 1) double;
              stableONTime (1, 1) double;
              lowValueThreshold (1, 1) double;
              highValueThreshold (1, 1) double;
            end
            obj = obj.configure(stableOFFTime, stableONTime, lowValueThreshold, highValueThreshold);
            obj = obj.resetDevice();
        end
        
        function obj = configure(obj, stableOFFTime, stableONTime, lowValueThreshold, highValueThreshold)
            arguments
              obj Estimation.ContactHandler.SchmittTrigger;
              stableOFFTime (1, 1) double;
              stableONTime (1, 1) double;
              lowValueThreshold (1, 1) double;
              highValueThreshold (1, 1) double;
            end
            
            obj = obj.setStableOFFTime(stableOFFTime);
            obj = obj.setStableONTime(stableONTime);
            obj = obj.setLowValueThreshold(lowValueThreshold);
            obj = obj.setHighValueThreshold(highValueThreshold);
        end
        
        function obj = resetDevice(obj)
            obj.timer = 0;
            obj.currentState = true;
            obj.rawValue = 0.;
            obj.previousTime = -1;
            obj.verbose = 0;
        end
        
        function obj = updateDevice(obj, rawValue, currentTime)
            arguments
                obj Estimation.ContactHandler.SchmittTrigger;
                rawValue (1, 1) double;
                currentTime (1, 1) double;            
            end
            if (obj.previousTime < 0)
                if (currentTime > 0)
                    obj.previousTime = 0;
                else
                    obj.previousTime = currentTime;
                end
            end
            obj.rawValue = rawValue;        
            if (obj.currentState == false)       
                % Check for transition - if valid over a timeframe, then switch
                if (obj.rawValue >= obj.highValueThreshold)
                    if (obj.timer > obj.stableONTime)
                        % rise to high
                        obj.currentState = true;
                    else
                        % wait for timer
                        obj.timer = obj.timer + (currentTime - obj.previousTime);
                    end
                else
                    % stable low - reset timer
                    obj.timer = 0;
                end
            else
                % check for transition - if valid over a timeframe, then switch
                if (obj.rawValue <= obj.lowValueThreshold)
                    if (obj.timer > obj.stableOFFTime)
                        % fall to low
                        obj.currentState = false;                        
                    else
                        % wait for timer
                        obj.timer = obj.timer + (currentTime - obj.previousTime);
                    end                    
                else
                    % stable high - reset timer
                    obj.timer = 0;
                end
            end
            obj.previousTime = currentTime;
        end
        
        function obj = setStableOFFTime(obj, stableOFFTime)
            arguments
                obj Estimation.ContactHandler.SchmittTrigger;
                stableOFFTime (1,1) double;
            end
           obj.stableOFFTime = stableOFFTime;
        end
        
        function obj = setStableONTime(obj, stableONTime)
            arguments
                obj Estimation.ContactHandler.SchmittTrigger;
                stableONTime (1,1) double;
            end
           obj.stableONTime = stableONTime;
        end
        
        function obj = setLowValueThreshold(obj, lowValueThreshold)
            arguments
                obj Estimation.ContactHandler.SchmittTrigger;
                lowValueThreshold (1,1) double;
            end
           obj.lowValueThreshold = lowValueThreshold;
        end
        
        function obj = setHighValueThreshold(obj, highValueThreshold)
            arguments
                obj Estimation.ContactHandler.SchmittTrigger;
                highValueThreshold (1,1) double;
            end
           obj.highValueThreshold = highValueThreshold;
        end
        
        function obj = setInitialState(obj, state)
            arguments
                obj Estimation.ContactHandler.SchmittTrigger;
                state (1,1) logical;
            end
           obj.currentState = state;
        end
        
        function state = getState(obj)
            state = obj.currentState;
        end
        
        function elapsed_time = getElapsedTime(obj)
            elapsed_time = obj.previousTime;
        end
        
        function obj = setVerbose(obj)
           obj.verbose = true;
        end
        
        function obj = unsetVerbose(obj)
           obj.verbose = false;
        end
    end
end

