classdef SE3
    methods (Static)
        function Tout = compose(T1, T2)
            [R1, p1] = LieGroups.SE3.extractSE3(T1);
            [R2, p2] = LieGroups.SE3.extractSE3(T2);
            
            Tout = LieGroups.SE3.constructSE3(R1*R2, R1*p2 + p1);
        end
        
        function Tout = inverse(T)
            [R, p] = LieGroups.SE3.extractSE3(T);
            
            Rinv = R';
            pinv = -R'*p;
            Tout = LieGroups.SE3.constructSE3(Rinv, pinv);
        end
        
        function Tout = identity()
            Tout = eye(4);
        end
        
        function v = vee(g)
            [omega_cross, vlin] = LieGroups.SE3.extractse3(g);
            omega = LieGroups.SO3.unskew(omega_cross);
            v = [vlin; omega];
        end
        
        function g = hat(v)
            vlin = v(1:3);
            omega = v(4:6);
            omega_cross = LieGroups.SO3.skew(omega);
            g = LieGroups.SE3.constructse3(omega_cross, vlin);
        end
        
        function T = exp(g)
            v = LieGroups.SE3.vee(g);
            T = LieGroups.SE3.exphat(v);
        end
        
        function T = exphat(v)
            vlin = v(1:3);
            omega = v(4:6);
            
            R = LieGroups.SO3.exphat(omega);
            JlSO3 = LieGroups.SO3.leftJacobian(omega);
            
            p = JlSO3*vlin;
            T = LieGroups.SE3.constructSE3(R, p);
        end
        
        function g = log(T)
            v = LieGroups.SE3.logvee(T);
            g = LieGroups.SE3.hat(v);
        end
        
        function v = logvee(T)
            [R, p] = LieGroups.SE3.extractSE3(T);
            omega = LieGroups.SO3.logvee(R);
            JlinvSO3 = LieGroups.SO3.leftJacobianInverse(omega);
            rho = JlinvSO3*p;
            
            v = [rho; omega];
        end
        
        function AdT = AdjointMatrix(T)
            [R, p] = LieGroups.SE3.extractSE3(T);
            p_cross = LieGroups.SO3.skew(p);
            AdT = [      R  p_cross*R;
                zeros(3)         R];
        end
        
        function adT = crossProductMatrix(v)
            vlin = v(1:3);
            omega = v(4:6);
            
            vlin_cross = LieGroups.SO3.skew(vlin);
            omega_cross = LieGroups.SO3.skew(omega);
            
            adT = [omega_cross vlin_cross;
                zeros(3)    omega_cross];
        end
        
        function Q = leftJacobianLinearSubmatrix(v)
            vlin = v(1:3);
            
            omega = v(4:6);
            theta = norm(omega);
            
            V = LieGroups.SO3.skew(vlin);
            
            Q = 0.5*V;
            if theta < eps
                return;
            end
            
            W = LieGroups.SO3.skew(omega);
            theta_squared = theta*theta;
            
            if theta_squared < eps
                theta_squared = theta;
                theta_cubed = theta;
                theta_fourth = theta;
                theta_fifth = theta;
            else
                theta_cubed = theta_squared*theta;
                theta_fourth = theta_cubed*theta;
                theta_fifth = theta_fourth*theta;
            end
            
            WV = W*V;
            VW = V*W;
            WVW = W*VW;
            
            a = (theta - sin(theta))/theta_cubed;
            Q = Q + a*(WV + ...
                VW + ...
                WVW);
            
            b = (theta_squared + (2*cos(theta)) - 2)/(2*theta_fourth);
            Q = Q + b*(W*WV + ...
                VW*W - ...
                - 3*WVW);
            
            c = ((2*theta) - (3*sin(theta)) + (theta*cos(theta)))/(2*theta_fifth);
            Q = Q + c*(WVW*W + ...
                W*WVW);
        end
        
        function Jl = rightJacobian(v)
            omega = v(4:6);
            JrSO3 = LieGroups.SO3.rightJacobian(-omega);
            Qr = LieGroups.SE3.leftJacobianLinearSubmatrix(-v);
            
            Jl = [JrSO3       Qr;
                zeros(3)  JrSO3];
        end
        
        function Jl = leftJacobian(v)
            omega = v(4:6);
            JlSO3 = LieGroups.SO3.leftJacobian(omega);
            Q = LieGroups.SE3.leftJacobianLinearSubmatrix(v);
            
            Jl = [JlSO3       Q;
                zeros(3)  JlSO3];
        end
        
        function Jlinv = leftJacobianInverse(v)
            omega = v(4:6);
            JlinvSO3 = LieGroups.SO3.leftJacobianInverse(omega);
            Q = LieGroups.SE3.leftJacobianLinearSubmatrix(v);
            
            Jlinv = [JlinvSO3      -JlinvSO3*Q*JlinvSO3;
                zeros(3)                  JlinvSO3];
        end
        
        function Jl = leftJacobianFromCrossProductMatrix(v)
            phi = norm(v);
            Jl = eye(6);
            if phi < 1e-4
                return
            end
            
            phisquared = phi*phi;
            phi_cubed = phisquared*phi;
            phi_fourth = phi_cubed*phi;
            phi_fifth = phi_fourth*phi;
            
            if phisquared < 1e-6
                phisquared = phi;
                phi_cubed = phi;
                phi_fourth = phi;
                phi_fifth = phi;
            end
            
            sinphi = sin(phi);
            cosphi = cos(phi);
            phisinphi = phi*sinphi;
            phicosphi = phi*cosphi;
            
            a = (4 - phisinphi - 4*cosphi)/(2*phisquared);
            ad = LieGroups.SE3.crossProductMatrix(v);
            Jl = Jl + a*ad;
            
            b = ((4*phi) - (5*sinphi) + phicosphi)/(2*phi_cubed);
            adsquared = ad*ad;
            Jl = Jl + b*adsquared;
            
            c = (2 - phisinphi - (2*cosphi))/(2*phi_fourth);
            adcubed = ad*adsquared;
            Jl = Jl + c*adcubed;
            
            d = ((2*phi) - (3*sinphi) + phicosphi)/(2*phi_fifth);
            adfourth = ad*adcubed;
            Jl = Jl + d*adfourth;            
        end       
        
        function out_12by6 =  dexphat_dv_at_identity()    
            out_12by6 = zeros(12, 6);
            dexp_domega_SO3 = LieGroups.SO3.dexphat_domega_at_identity();
            
            out_12by6(1:9, 4:6) = dexp_domega_SO3;
            out_12by6(10:12, 1:3) = eye(3);
        end
        
        function out_12by12 = dcomposeAB_dA(A, B)
            out_12by12 = kron(B, eye(3));
        end
        
        function out_12by12 = dcomposeAB_dB(A, B)
            [Ra, ~] = LieGroups.SE3.extractSE3(A);
            out_12by12 = kron(eye(4), Ra);
        end
        
        function out_12by6 = dcomposeAexphat_dv_at_identity(A)
            dAI_dI = LieGroups.SE3.dcomposeAB_dB(A, eye(4));
            dexp_dv = LieGroups.SE3.dexphat_dv_at_identity();
            
            out_12by6 = dAI_dI*dexp_dv;
        end
        
        function out_6by12 = dlogveeT_dT(T)
            out_6by12 = zeros(6, 12);
            [R, p] = LieGroups.SE3.extractSE3(T);
            dcomposeLeftJinvSO3Withvec_dR = LieGroups.SO3.dcomposeLeftJinvWithvec_dR(R, p);
            omega = LieGroups.SO3.logvee(R);
            JlinvSO3 = LieGroups.SO3.leftJacobianInverse(omega);
            dlogveeR_dR = LieGroups.SO3.dlogveeR_dR(R);
            
            out_6by12(1:3, 1:9) = dcomposeLeftJinvSO3Withvec_dR;
            out_6by12(1:3, 4:6) = JlinvSO3;            
            out_6by12(4:6, 1:9) = dlogveeR_dR;
        end
        
        function [R, p]  = extractSE3(T)
            R = T(1:3, 1:3);
            p = T(1:3, 4);
        end
        
        function [omega_cross, vlin] = extractse3(g)
            omega_cross = g(1:3, 1:3);
            vlin = g(1:3, 4);
        end
        
        function T = constructSE3(R, p)
            T = eye(4);
            T(1:3, 1:3) = R;
            T(1:3, 4)  = p;
        end
        
        function g = constructse3(omega_cross, vlin)
            g = zeros(4);
            g(1:3, 1:3) = omega_cross;
            g(1:3, 4) = vlin;
        end
        
        function dist = geodesicL2Distance(H1, H2)
            H1inv = LieGroups.SE3.inverse(H1);
            pose_error = LieGroups.SE3.compose(H1inv, H2);
            error_vec = LieGroups.SE3.logvee(pose_error);
            dist = norm(error_vec);
        end
        
        function [meanH, found] = geodesicL2WeightedMeanPose(Harray, weights, tolerance, dtRiemannianDescent, max_iter)            
            assert(length(Harray) == length(weights), 'size mismatch');
            
            dt = dtRiemannianDescent;
            meanH = Harray{1};
            optimal_H_found = false;
            iteration = 1;
            total_weights = sum(weights);
            tic
            while ~optimal_H_found
                perturbation = zeros(6, 1);
                for idx = 1:length(Harray)
                    H_i = Harray{idx};
                    if isempty(weights)
                        w_i = 1;
                    else
                        w_i = weights(idx);
                    end
                    
                    meanHinv = LieGroups.SE3.inverse(meanH);
                    error = LieGroups.SE3.compose(meanHinv, H_i);
                    deltar = LieGroups.SE3.logvee(error);
                    
                    perturbation = perturbation + w_i*deltar;
                end
                perturbation = perturbation*(1/total_weights);
                if norm(perturbation) < tolerance
                    optimal_H_found = true;
                    continue;
                end
                
                % else descend along the curve with stepsize
                perturbation = perturbation*dt;
                projSE3 = LieGroups.SE3.exphat(perturbation);
                meanH = LieGroups.SE3.compose(meanH, projSE3);
                                
                if (iteration > max_iter)
                    found = false;
                    meanH = Harray{1};
                    disp('Iterations maxed, returning first pose in input')
                    return
                end
                
                iteration = iteration + 1;
            end
            found = true;
        end
        
    end
end

