classdef SO3
    methods (Static)
        function Rout = compose(R1, R2)
            Rout = R1*R2;
        end
        
        function Rout = inverse(R1)
            Rout = R1';
        end
        
        function Rout = identity()
            Rout = eye(3);
        end
        
        function omega = vee(omega_cross)
            omega = LieGroups.SO3.unskew(omega_cross);
        end
        
        function omega_cross = hat(omega)
            omega_cross = LieGroups.SO3.skew(omega);
        end
        
        function Rout = exp(omega_cross)
            omega = LieGroups.SO3.unskew(omega_cross);
            Rout = exphat(omega);
        end
        
        function Rout = exphat(omega)
            angle = norm(omega);
            if (abs(angle) < sqrt(eps))
                Rout = eye(3);
                return;
            end
            axis = omega/angle;
            I = eye(3);
            across = LieGroups.SO3.skew(axis);
            first = cos(angle)*I;
            second = sin(angle)*across;
            third = (1 - cos(angle))*(axis*axis');
            Rout =  first + second + third;
        end
        
        function omega_cross = log(R)
            costheta = 0.5*(trace(R) - 1);
            theta = acos(costheta) ;
            
            if theta < sqrt(eps)
                omega_cross = zeros(3, 3);
            else                
                scalar = theta/(2*sin(theta));
                omega_cross = scalar*(R - R');
            end
        end
        
        function omega = logvee(R)
            omega_cross = LieGroups.SO3.log(R);
            omega =  LieGroups.SO3.vee(omega_cross);
        end
        
        function AdR = AdjointMatrix(R)            
            AdR = R;
        end
        
        function adR = crossProductMatrix(omega)            
            adR = LieGroups.SO3.skew(omega);
        end
        
        function Jr = rightJacobian(omega)
            % used to translate the changes in position directions in lie algebra of se(3)
            % to the group element in SE(3), this might be particularly
            % useful also while defining probabilities density functions in
            % Lie algebra and transforming them to the Lie groups
            % this is obtained as a consequence of the BCH formula
            % the first order term in the expansion of
            % log(exphat(omega1) exphat(omega2))
            Jr = LieGroups.SO3.leftJacobian(-omega);
%             theta = norm(omega);
%             if theta < sqrt(eps)
%                 Jr = eye(3);
%                 return
%             end
%             
%             a = sin(theta)/theta;
%             b = (1 - cos(theta))/theta;
%             c = 1 - (sin(theta)/theta);
%             
%             phi = omega/theta;
%             phi_cross = LieGroups.SO3.skew(phi);
%             Jr = a*eye(3) - b*phi_cross + c*(phi*phi');
        end
        
        function Jl = leftJacobian(omega)
            % used to translate the changes in position directions in lie algebra of se(3)
            % to the group element in SE(3), this might be particularly
            % useful also while defining probabilities density functions in
            % Lie algebra and transforming them to the Lie groups
            % this is obtained as a consequence of the BCH formula
            % the first order term in the expansion of
            % log(exphat(omega1) exphat(omega2))
            theta = norm(omega);
            if theta < sqrt(eps)
                Jl = eye(3);
                return
            end
            
            a = sin(theta)/theta;
            b = (1 - cos(theta))/theta;
            c = 1 - (sin(theta)/theta);
            
            phi = omega/theta;
            phi_cross = LieGroups.SO3.skew(phi);
            Jl = a*eye(3) + b*phi_cross + c*(phi*phi');
        end
        
        function Jlinv = leftJacobianInverse(omega)
            theta = norm(omega);
            if theta < sqrt(eps)
                Jlinv = eye(3);
                return
            end
            
            halftheta = theta/2;
            cothalftheta = cos(halftheta)/sin(halftheta);
            
            a = halftheta*cothalftheta;
            b = -halftheta;
            c = 1 - a;
            phi = omega/theta;
            phi_cross = LieGroups.SO3.skew(phi);
           
            Jlinv = a*eye(3) - b*phi_cross + c*(phi*phi');            
        end
        
        function out_9by3 =  dexphat_domega_at_identity()
            e1cross = LieGroups.SO3.skew([1; 0; 0]);
            e2cross = LieGroups.SO3.skew([0; 1; 0]);
            e3cross = LieGroups.SO3.skew([0; 0; 1]);
            
            out_9by3 = [-e1cross; -e2cross; -e3cross];            
        end
        
        function out_3by9 = dlogveeR_dR(R)
            costheta = (trace(R) - 1.0)/2;
            
            if costheta > 0.9999
                a = [0; 0; 0];
                b = 0.5;
            else
                theta = acos(costheta);
                sintheta = sqrt(1 - (costheta*costheta));
                
                alpha = ((theta*costheta) - sintheta)/(4*(sintheta^3));
                a = alpha*LieGroups.SO3.unskew((R-R'));
                b = theta/(2*sintheta);
            end
            
            out_3by9 = zeros(3, 9);
            out_3by9(2, 3) = -b;
            out_3by9(3, 4) = -b;
            out_3by9(1, 8) = -b;
            out_3by9(3, 2) = b;
            out_3by9(1, 6) = b;
            out_3by9(2, 7) = b;
            
            out_3by9(1:3, 1) = a;
            out_3by9(1:3, 5) = a;
            out_3by9(1:3, 9) = a;
        end
        
        function out_3by9 = dcomposeLeftJinvWithvec_dR(R, t)
            costheta = (trace(R) - 1.0)/2;
            
            if costheta > 0.9999
                a = [0; 0; 0];
                b = 0.5;
            else
                theta = acos(costheta);
                sintheta = sqrt(1 - (costheta*costheta));
                
                alpha = ((theta*costheta) - sintheta)/(2*theta*(sintheta^2));
                logR = LieGroups.SO3.log(R);
                a = alpha*logR*t;
                b = theta/(2*sintheta);                
            end
            
            out_3by9 = zeros(3, 9);
            out_3by9(1, 2) = -b*t(2);
            out_3by9(1, 3) = -b*t(3);
            out_3by9(1, 4) = -b*t(3);
            out_3by9(1, 7) = b*t(3);
            
            out_3by9(2, 2) = b*t(1);
            out_3by9(2, 6) = -b*t(3);
            out_3by9(2, 8) = b*t(3);
            
            out_3by9(3, 3) = b*t(1);            
            out_3by9(3, 4) = b*t(1);
            out_3by9(3, 6) = b*t(2);
            out_3by9(3, 7) = -b*t(1);
            out_3by9(3, 8) = b*t(2);
            
            out_3by9(1:3, 1) = a;
            out_3by9(1:3, 5) = a;
            out_3by9(1:3, 9) = a;
        end
        
        function S = skew(w)
            %SKEW skew symmetric matrix
            assert(length(w) == 3, 'skew symmetric matrix for 3d vector only')
            S = [0, -w(3), w(2); ...
                w(3), 0.0, -w(1); ...
                -w(2), w(1), 0];
        end
        
        function w = unskew(S)
            %SKEW skew symmetric matrix
%             assert(isequal(S, -S'), 'input not a skew symmetric matrix')
            w = [S(3, 2); S(1, 3); S(2, 1)];
        end
        
        function R = sampleUniform()
            % fast random rotations by James Arvo
            theta = 2*pi*rand; % rotation about the north pole
            phi = 2*pi*rand; % pick a direction to deflect the pole
            z = rand; % amount of pole deflection
            
            % vector for performing reflection
            v = [cos(phi)*sqrt(z); ...
                sin(phi)*sqrt(z); ...
                sqrt(1 - z)];
            
            % find rotation about z axis and then rotate the z axis to
            % random orientation
            Rotz = [cos(theta) sin(theta) 0; ...
                -sin(theta) cos(theta) 0; ...
                0          0  1];
            
            % householder_matrix describing reflection about a plane containing origin
            H = eye(3) - (2*v*v');
            
            R = -H*Rotz;
        end
        
        function Rout = sampleConcentratedGaussian(rpy_mean_rad, sigma_tangent_space)
            roll = rpy_mean_rad(1);
            pitch = rpy_mean_rad(2);
            yaw = rpy_mean_rad(3);
            rotX = [1 0 0 ; 0 cos(roll) -sin(roll); 0  sin(roll) cos(roll)] ;
            rotY = [cos(pitch) 0 sin(pitch); 0 1 0 ;  -sin(pitch) 0 cos(pitch)] ;
            rotZ = [cos(yaw) -sin(yaw) 0; sin(yaw) cos(yaw) 0; 0 0 1 ] ;
            
            Rmean = rotZ*rotY*rotX;
            
            mu = [0 0 0];
            A = chol(sigma_tangent_space);
            while 1
                epsilon = mu + randn*A;
                if norm(epsilon) < pi
                    break;
                end
            end
            
            Rout = Rmean*LieGroups.SO3.exphat(epsilon);
        end
        
        function dist = geodesicL2Distance(R1, R2)
            R1inv = R1';
            rot_error = R1inv*R2;
            error_vec = LieGroups.SO3.logvee(rot_error);
            dist = norm(error_vec);
        end
        
        function [meanR, found] = geodesicL2WeightedMeanRotation(Rarray, weights, tolerance, dtRiemannianDescent, max_iter)            
            assert(length(Rarray) == length(weights), 'size mismatch');
            
            dt = dtRiemannianDescent;
            meanR = Rarray{1};
            optimal_R_found = false;
            iteration = 1;
            total_weights = sum(weights);
            tic
            while ~optimal_R_found
                perturbation = zeros(3, 1);
                for idx = 1:length(Rarray)
                    R_i = Rarray{idx};
                    if isempty(weights)
                        w_i = 1;
                    else
                        w_i = weights(idx);
                    end
                    
                    meanRinv = meanR';
                    error = meanRinv*R_i;
                    deltar = LieGroups.SO3.logvee(error);
                    
                    perturbation = perturbation + w_i*deltar;
                end
                perturbation = perturbation*(1/total_weights);
                if norm(perturbation) < tolerance
                    optimal_R_found = true;
                    continue;
                end
                
                % else descend along the curve with stepsize
                perturbation = perturbation*dt;
                projSO3 = LieGroups.SO3.exphat(perturbation);
                meanR = meanR*projSO3;
                                
                if (iteration > max_iter)
                    found = false;
                    meanR = Rarray{1};
%                     disp('Iterations maxed, returning first rotation in input')
                    return
                end
                
                iteration = iteration + 1;
            end
            found = true;
        end
        
    end
    
    
    
end

