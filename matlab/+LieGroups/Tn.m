classdef Tn
    methods (Static)
        function Tout = compose(T1, T2)
            assert(size(T1, 2) == size(T2, 2), 'size mismatch');
            t1 = LieGroups.Tn.extractTn(T1);
            t2 = LieGroups.Tn.extractTn(T2);
            Tout = LieGroups.Tn.constructTn(t1+t2);
        end
        
        function Tout = inverse(T)
            t = LieGroups.Tn.extractTn(T);
            Tout = LieGroups.Tn.constructTn(-t);
        end
        
        function Tout = identity(n)
            Tout = eye(n+1);
        end
        
        function v = vee(g)
            v = LieGroups.Tn.extractTn(g);
        end
        
        function g = hat(v)            
            g = LieGroups.Tn.constructtn(v);
        end
        
        function T = exp(g)
            v = LieGroups.Tn.vee(g);
            T = LieGroups.Tn.exphat(v);
        end
        
        function T = exphat(v)           
            T = LieGroups.Tn.constructTn(v);
        end
        
        function g = log(T)
            v = LieGroups.Tn.logvee(T);
            g = LieGroups.Tn.hat(v);
        end
        
        function v = logvee(T)
            v = LieGroups.Tn.extractTn(T);
        end
        
        function AdT = AdjointMatrix(T)
            n = size(T, 2) - 1;
            AdT = eye(n);
        end
        
        function adT = crossProductMatrix(v)
            n = length(v);
            adT = zeros(n);
        end
        
        function Jr = rightJacobian(v)
            n = length(v);
            Jr = eye(n);
        end
        
        function Jl = leftJacobian(v)
            n = length(v);
            Jl = eye(n);
        end
        
        function Jlinv = leftJacobianInverse(v)
            n = length(v);
            Jlinv = -eye(n);
        end
        
        function t  = extractTn(T)
            n = size(T, 2) - 1;            
            t = T(1:n, n+1);
        end
        
        
        function T = constructTn(t)
            n = length(t);
            T = eye(n+1);
            T(1:n, n+1) = t;
        end
        
        function g = constructtn(t)
            n = length(t);
            g = zeros(n+1);
            g(1:n, n+1) = t;
        end
    end
end

