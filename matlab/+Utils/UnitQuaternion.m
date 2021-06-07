classdef UnitQuaternion 
    methods (Static)
                       
        function qinv = inverse(q)
            arguments                
                q (4, 1) double;
            end
            qinv = zeros(4, 1);
            qinv(1) = q(1);
            qinv(2) = -q(2);
            qinv(3) = -q(3);
            qinv(4) = -q(4);
        end
                                
                        
        function q = expQuaternion(omega)
            arguments
                omega (3, 1) double;
            end            
            q = zeros(4, 1);
            
            a = norm(omega);
            q(1) = cos(a/2);
            if (a > 1e-10)
                q(2:4) = sin(a/2)/a.*omega;
            else
                q(2:4) = omega;
            end            
        end
        
        function omega = logQuaternion(q)
            arguments
                q (4, 1) double;
            end
            omega = zeros(3, 1);            
            c = q(1);
            vec = q(2:4);
            s = norm(vec);
            
            if (s > 1e-10)
                a = 2*atan2(s, c);
                omega = (vec.*(a/s));
            else
                omega = (vec.*2);
            end            
        end
        
        function q = composeQuaternion(q1, q2)
            arguments 
                q1 (4, 1) double;
                q2 (4, 1) double;
            end
                        
            % let's follow a left action of q1 on q2
            I4 = eye(4);
                                    
            q1w = q1(1);
            q1v = q1(2:4);
            
            Lq1 = zeros(4, 4);
            Lq1(1, 2:4) = -q1v';
            Lq1(2:4, 1) = q1v;
            Lq1(2:4, 2:4) = Utils.skew(q1v);
            
            Lq1 = Lq1 + (q1w.*I4);
            
            q = Lq1*q2;            
        end
    end
end

