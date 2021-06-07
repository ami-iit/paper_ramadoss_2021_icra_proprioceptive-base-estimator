classdef SE_Nplus2_3
        methods (Static)
        function Tout = compose(T1, T2)
            [R1, p1, v1, d13byn] = LieGroups.SE_Nplus2_3.extractSENplus23(T1);
            [R2, p2, v2, d23byn] = LieGroups.SE_Nplus2_3.extractSENplus23(T2);
            Tout = LieGroups.SE_Nplus2_3.constructSENplus23(R1*R2, ...
                                                            R1*p2 + p1, ...
                                                            R1*v2 + v1, ...
                                                            R1*d23byn + d13byn);
        end
        
        function Tout = inverse(T)
            [R, p, v, d3byn] = LieGroups.SE_Nplus2_3.extractSENplus23(T);
            
            Rinv = R';
            pinv = -R'*p;
            vinv = -R'*v;
            dinv = -R'*d3byn;

            Tout = LieGroups.SE_Nplus2_3.constructSENplus23(Rinv, pinv, vinv, dinv);
        end
        
        function Tout = identity(n)
            Tout = eye(5+n);
        end
        
        function v = vee(g)
            [omega_cross, vlin, alin, vdlin3byn] = LieGroups.SE_Nplus2_3.extractseNplus23(g);
            omega = LieGroups.SO3.unskew(omega_cross);
            vdvert = reshape(d, 3*size(vdlin3byn, 2), 1);
            v = [vlin; omega; alin; vdvert];
        end
        
        function g = hat(v)
            vlin = v(1:3);
            omega = v(4:6);
            alin = v(7:9);
            assert(length(v) > 9, 'Use SE_2_3 instead');
            vdlin = v(10:end);
            assert(mod(length(vdlin), 3) == 0, 'vector size mismatch');
            
            omega_cross = LieGroups.SO3.skew(omega);
            n= length(vdlin)/3;
            assert( n ~= 0, 'Use SE_2_3 instead');
            vdlin3byn = reshape(vdlin, 3, n);
            g = LieGroups.SE_Nplus2_3.constructseNplus23(omega_cross, vlin, alin, vdlin3byn);
        end
        
        function T = exp(g)
            v = LieGroups.SE_Nplus2_3.vee(g);
            T = LieGroups.SE_Nplus2_3.exphat(v);
        end
        
        function T = exphat(v)
            vlin = v(1:3);
            omega = v(4:6);
            alin = v(7:9);
            assert(length(v) > 9, 'Use SE_2_3 instead');
                
            vdlin = v(10:end);
            assert(mod(length(vdlin), 3) == 0, 'vector size mismatch');
            n =  length(vdlin)/3;
            assert( n ~= 0, 'Use SE_2_3 instead');
            vdlin3byn = reshape(vdlin, 3, n);
            
            R = LieGroups.SO3.exphat(omega);
            JlSO3 = LieGroups.SO3.leftJacobian(omega);
            
            p = JlSO3*vlin;
            v = JlSO3*alin;
            d3byn = JlSO3*vdlin3byn;
            T = LieGroups.SE_Nplus2_3.constructSENplus23(R, p, v, d3byn);
        end
        
        function g = log(T)
            v = LieGroups.SE_Nplus2_3.logvee(T);
            g = LieGroups.SE_Nplus2_3.hat(v);
        end
        
        function v = logvee(T)
            [R, p, v, d3byn] = LieGroups.SE_Nplus2_3.extractSENplus23(T);
            
            omega = LieGroups.SO3.logvee(R);
            JlinvSO3 = LieGroups.SO3.leftJacobianInverse(omega);
            rho = JlinvSO3*p;
            alpha = JlinvSO3*v;
            vdrho = JlinvSO3*d3byn;
            vdvert = reshape(vdrho, 3*size(d3byn, 2), 1);

            v = [rho; omega; alpha; vdvert];
        end
        
        function AdT = AdjointMatrix(T)
            [R, p, v, d3byn] = LieGroups.SE_Nplus2_3.extractSENplus23(T);
            n = size(d3byn, 2);
            assert( n ~= 0, 'Use SE_2_3 instead');
            AdT = zeros(9 +3*n);
            p_cross = LieGroups.SO3.skew(p);
            v_cross = LieGroups.SO3.skew(v);
            
            AdT(1:3, 1:3) =   R;
            AdT(1:3, 4:6) =   p_cross*R;
            
            AdT(4:6, 4:6) =   R;

            AdT(7:9, 7:9) =   R;
            AdT(7:9, 4:6) =   v_cross*R;

            for idx = 1:n
                d = d3byn(:, idx);
                d_cross = LieGroups.SO3.skew(d);
                
                begin = 10+3*(idx-1);
                fin = begin+2;
                AdT(begin:fin, begin:fin) =   R;
                AdT(begin:fin, 4:6) =   d_cross*R;
            end
        end
        
        function adT = crossProductMatrix(v)
            vlin = v(1:3);
            omega = v(4:6);
            alin = v(7:9);
            assert(length(v) > 9, 'Use SE_2_3 instead');
            
            n = (length(v)  - 9)/3;
            adT = zeros(length(v));
            
            vlin_cross = LieGroups.SO3.skew(vlin);
            omega_cross = LieGroups.SO3.skew(omega);
            alin_cross = LieGroups.SO3.skew(alin);
            
            adT(1:3, 1:3) =   omega_cross;
            adT(1:3, 4:6) =   vlin_cross;
            
            adT(4:6, 4:6) =   omega_cross;

            adT(7:9, 7:9) =   omega_cross;
            adT(7:9, 4:6) =   alin_cross;
            
            for idx = 1:n
                begin = 10+3*(idx-1);
                fin = begin+2;
                vd = v(begin:fin);
                vd_cross = LieGroups.SO3.skew(vd);
                                
                AdT(begin:fin, begin:fin) =   omega_cross;
                AdT(begin:fin, 4:6) =   vd_cross;
            end            
        end
        
        function Jr = rightJacobian(v)
            Jr = LieGroups.SE_Nplus2_3.leftJacobian(-v);
        end
        
        function Jl = leftJacobian(v)
            omega = v(4:6);
            JlSO3 = LieGroups.SO3.leftJacobian(omega);
            
            Jl = zeros(length(v));
            
            Qv = LieGroups.SE3.leftJacobianLinearSubmatrix(v(1:6));
            Qa = LieGroups.SE3.leftJacobianLinearSubmatrix([v(7:9); omega]);
            
            Jl(1:3, 1:3) =   JlSO3;
            Jl(1:3, 4:6) =   Qv;
            
            Jl(4:6, 4:6) =   JlSO3;

            Jl(7:9, 7:9) =   JlSO3;
            Jl(7:9, 4:6) =   Qa;
            
            assert(length(v) > 9, 'Use SE_2_3 instead');            
            n = (length(v)  - 9)/3;
            
            for idx = 1:n
                begin = 10+3*(idx-1);
                fin = begin+2;
                vd = v(begin:fin);
                Qvd = LieGroups.SE3.leftJacobianLinearSubmatrix([vd; omega]);
                                
                Jl(begin:fin, begin:fin) =   JlSO3;
                Jl(begin:fin, 4:6) =   Qvd;
            end 
        end
        
        function Jlinv = leftJacobianInverse(v)
            Jl = LieGroups.SE_Nplus2_3.leftJacobian(v);            
            Jlinv = inv(Jl);
        end
        
        function [R, p, v, d3byn]  = extractSENplus23(T)      
            n = size(T, 2) - 5;
            assert( n ~= 0, 'Use SE_2_3 instead');
            R = T(1:3, 1:3);
            p = T(1:3, 4);
            v = T(1:3, 5);
            d3byn = T(1:3, 6:6+n-1);
        end
        
        function [omega_cross, vlin, alin, vdlin3byn] = extractseNplus23(g)
            n = size(g, 2) - 5;
            assert( n ~= 0, 'Use SE_2_3 instead');
            omega_cross = g(1:3, 1:3);
            vlin = g(1:3, 4);
            alin = g(1:3, 5);
            vdlin3byn = g(1:3, 6:6+n-1);
        end
        
        function T = constructSENplus23(R, p, v, d3byn)
            n = size(d3byn, 2);
            assert( n ~= 0, 'Use SE_2_3 instead');
            T = eye(5+n);
            T(1:3, 1:3) = R;
            T(1:3, 4)  = p;
            T(1:3, 5)  = v;
            T(1:3, 6:6+n-1) = d3byn;
        end
        
        function g = constructseNplus23(omega_cross, vlin, alin, vdlin3byn)
            n = size(vdlin3byn, 2);
            assert( n ~= 0, 'Use SE_2_3 instead');
            g = zeros(5+n);
            g(1:3, 1:3) = omega_cross;
            g(1:3, 4) = vlin;
            g(1:3, 5) = alin;
            g(1:3, 6:6+n-1) = vdlin3byn;
        end
    end
end

