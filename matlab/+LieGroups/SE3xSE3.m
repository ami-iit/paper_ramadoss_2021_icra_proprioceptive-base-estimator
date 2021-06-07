classdef SE3xSE3
    methods (Static)
        function Tout = compose(T1, T2)
            [T11, T12] = LieGroups.SE3xSE3.extractSE3bySE3(T1);
            [T21, T22] = LieGroups.SE3xSE3.extractSE3bySE3(T2);
            
            T1out = LieGroups.SE3.compose(T11, T21);
            T2out = LieGroups.SE3.compose(T12, T22);
            Tout = LieGroups.SE3xSE3.constructSE3bySE3(T1out, T2out);
        end
        
        function Tout = inverse(T)
            [T1, T2] = LieGroups.SE3xSE3.extractSE3bySE3(T);
            
            T1inv = LieGroups.SE3.inverse(T1);
            T2inv = LieGroups.SE3.inverse(T2);

            Tout = LieGroups.SE3xSE3.constructSE3bySE3(T1inv, T2inv);
        end
        
        function Tout = identity()
            Tout = eye(8);
        end
        
        function v = vee(g)
            [g1, g2] = LieGroups.SE3xSE3.extractse3byse3(g);
            v1 = LieGroups.SE3.vee(g1);
            v2 = LieGroups.SE3.vee(g2);
            v = [v1; v2];
        end
        
        function g = hat(v)
            v1 = v(1:6);
            v2 = v(7:12);
            g1 = LieGroups.SE3.hat(v1);
            g2 = LieGroups.SE3.hat(v2);
            g = LieGroups.SE3xSE3.constructse3byse3(g1, g2);
        end
        
        function T = exp(g)
            v = LieGroups.SE3xSE3.vee(g);            
            T = LieGroups.SE3xSE3.exphat(v);            
        end
        
        function T = exphat(v)
            v1 = v(1:6);
            v2 = v(7:12);
            
            T1 = LieGroups.SE3.exphat(v1);
            T2 = LieGroups.SE3.exphat(v2);
            T = LieGroups.SE3xSE3.constructSE3bySE3(T1, T2);
        end
        
        function g = log(T)
            v = LieGroups.SE3xSE3.logvee(T);
            g = LieGroups.SE3xSE3.hat(v);
        end
        
        function v = logvee(T)
            [T1, T2] = LieGroups.SE3xSE3.extractSE3bySE3(T);
            v1 = LieGroups.SE3.logvee(T1);
            v2 = LieGroups.SE3.logvee(T2);         
            
            v = [v1; v2];
        end
        
        function AdT = AdjointMatrix(T)
            [T1, T2] = LieGroups.SE3xSE3.extractSE3bySE3(T);
            X1 = LieGroups.SE3.AdjointMatrix(T1);
            X2 = LieGroups.SE3.AdjointMatrix(T2);
            AdT = blkdiag(X1, X2);
        end
        
        function adT = crossProductMatrix(v)
            v1 = v(1:6);
            v2 = v(7:12);
            
            X1 = LieGroups.SE3.crossProductMatrix(v1);
            X2 = LieGroups.SE3.crossProductMatrix(v2);
            adT = blkdiag(X1, X2);
        end        
        
        function Jl = leftJacobian(v)
            v1 = v(1:6);
            v2 = v(7:12);
            
            Jl1 = LieGroups.SE3.leftJacobian(v1);
            Jl2 = LieGroups.SE3.leftJacobian(v2);
            Jl = blkdiag(Jl1, Jl2);
        end
        
        function Jlinv = leftJacobianInverse(v)
            Jl = LieGroups.SE3xSE3.leftJacobian(v);            
            Jlinv = inv(Jl);
        end
        
        function [T1, T2]  = extractSE3bySE3(X)
            T1 = X(1:4, 1:4);
            T2 = X(5:8, 5:8);
        end
        
        function [g1, g2] = extractse3byse3(g)
            g1 = g(1:4, 1:4);
            g2 = g(5:8, 5:8);
        end
        
        function X = constructSE3bySE3(T1, T2)
            X = blkdiag(T1, T2);
        end
        
        function g = constructse3byse3(g1, g2)
            g = blkdiag(g1, g2);
        end
    end
end

