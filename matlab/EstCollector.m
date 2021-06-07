classdef EstCollector < handle 
    properties
        name        
        baseRPY;
        basePos;
        baseLinVel;
        baseAngVel;
        lfPos;
        rfPos;
        lfRPY;
        rfRPY;
        biasAcc;
        biasGyro;
    end
    
    methods
        function obj = EstCollector(name)
            obj.name = name;            
        end
        
        function obj = resizeBuffers(obj, size, Nl, Nr)
            obj.baseRPY = zeros(size, 3);
            obj.basePos = zeros(size, 3);
            obj.baseLinVel = zeros(size, 3);
            obj.baseAngVel = zeros(size, 3);
            obj.lfPos = zeros(size, Nl, 3);
            obj.rfPos = zeros(size, Nr, 3);
            obj.lfRPY = zeros(size, 3);
            obj.rfRPY = zeros(size, 3);
            obj.biasAcc = zeros(size, 3);
            obj.biasGyro = zeros(size, 3);
        end
        
        function obj = updateBasePose(obj, iter, rpy, pos, v, w)
            obj.baseRPY(iter, :) = rpy;
            obj.basePos(iter, :) = pos;
            obj.baseLinVel(iter, :) = v;
            obj.baseAngVel(iter, :) = w;
        end
        
        function obj = updateFootPosition(obj, foot, iter, vId, pos)
            if strcmp(foot, 'left')
                obj.lfPos(iter, vId, :) = pos;
            elseif strcmp(foot, 'right')
                obj.rfPos(iter, vId, :) = pos;
            end
        end
        
        function obj = updateFootRotation(obj, foot, iter, rpy)
            if strcmp(foot, 'left')
                obj.lfRPY(iter, :) = rpy;
            elseif strcmp(foot, 'right')
                obj.rfRPY(iter, :) = rpy;
            end
        end
        
        function obj = updateBias(obj, iter, ba, bg)
            obj.biasAcc(iter, :) = ba;
            obj.biasGyro(iter, :) = bg;
        end
        
        function cropped = startFrom(obj, startIter, endIter)
            cropped = EstCollector(obj.name);
            cropped.baseRPY = obj.baseRPY(startIter:endIter, :);
            cropped.basePos = obj.basePos(startIter:endIter, :);
            cropped.baseLinVel = obj.baseLinVel(startIter:endIter, :);
            cropped.baseAngVel = obj.baseAngVel(startIter:endIter, :);
            cropped.lfPos = obj.lfPos(startIter:endIter, :);
            cropped.rfPos = obj.rfPos(startIter:endIter, :);
            cropped.lfRPY = obj.lfRPY(startIter:endIter, :);
            cropped.rfRPY = obj.rfRPY(startIter:endIter, :);
            cropped.biasAcc = obj.biasAcc(startIter:endIter, :);
            cropped.biasGyro = obj.biasGyro(startIter:endIter, :);
        end
        
        function obj = setFootRotations(obj, lfrpy, rfrpy)
            obj.lfRPY = lfrpy;
            obj.rfRPY = rfrpy;
        end
        
        function obj = setFootPositions(obj, lfPos, rfPos)
            obj.lfPos = lfPos;
            obj.rfPos = rfPos;
        end
        
        function obj = setBias(obj, ba, bg)
            obj.biasAcc = ba;
            obj.biasGyro = bg;
        end
    end
end

