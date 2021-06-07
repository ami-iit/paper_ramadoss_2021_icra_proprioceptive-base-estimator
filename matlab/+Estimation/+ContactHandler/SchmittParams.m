classdef SchmittParams
    properties
        stableTimeContactMake (1, 1) double;
        stableTimeContactBreak (1, 1) double;
        contactMakeForceThreshold (1, 1) double;
        contactBreakForceThreshold (1, 1) double;
    end
    
    methods
        function obj = SchmittParams(obj)
            obj.stableTimeContactMake = 0.0;
            obj.stableTimeContactBreak = 0.0;
            obj.contactMakeForceThreshold = 0.0;
            obj.contactBreakForceThreshold = 0.0;
        end
    end
end

