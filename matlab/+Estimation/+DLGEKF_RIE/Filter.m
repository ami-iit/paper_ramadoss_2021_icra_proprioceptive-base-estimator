classdef Filter < handle
    properties (Access = private)
        timestamp_id (1, 1) uint64;
        
        X (:, :) double;
        P (:, :) double;
        
        t_prev (1, 1) double;
        X_prev (:, :) double;
        P_prev (:, :) double;
        alpha_prev (3, 1) double;
        omega_prev (3, 1) double;
        encoders_prev (:, 1) double;
        contacts_prev (2, 1) logical;
        
        X0 (:, :) double;
        P0 (:, :) double;
        
        Qg (3, 3) double; % gyro variance
        Qa (3, 3) double; % accelerometer variance
        Qcontact (6, 6) double; % contact foot pose variance
        Qswing (6, 6) double; % swing foot pose variance
        Qba (3, 3) double; % gyro bias variance
        Qbg (3, 3) double; % accelerometer bias variance
        
        Renc (:, :) double; % encoder variance
        R0fkin (6, 6) double; % forward kinematic pose variance
        
        model_comp ;
        options Estimation.Proprioception.EstimatorOptions;
        
        debugger Estimation.Proprioception.DebugOutputs;
        
        filter_configured (1, 1) logical;
        filter_initialized (1, 1) logical;
        bias_initialized (1, 1) logical;
    end
    
    properties (Access = private, Constant)
        g = [0; 0; -9.80665];
    end
    
    methods
        function obj = Filter()
            obj.timestamp_id = 0;
            obj.Qa = zeros(3);
            obj.Qg = zeros(3);
            obj.Qcontact = zeros(6);
            obj.Qswing = zeros(6);
            obj.Qba = zeros(3);
            obj.Qbg = zeros(3);
            obj.R0fkin = zeros(6);
            obj.t_prev = 0.0;
            obj.alpha_prev = zeros(3, 1);
            obj.omega_prev = zeros(3, 1);
            obj.contacts_prev = false(2, 1);            
            
            obj.options = Estimation.Proprioception.EstimatorOptions();
            obj.debugger = Estimation.Proprioception.DebugOutputs();
            
            obj.filter_configured = false;
            obj.filter_initialized = false;
            obj.bias_initialized = false;
        end
        
        function obj = setup(obj, prior_dev, sensors_dev, model_comp, options)
            arguments
                obj Estimation.DLGEKF_RIE.Filter;
                prior_dev Estimation.Proprioception.PriorsStdDev;
                sensors_dev Estimation.Proprioception.SensorsStdDev;
                model_comp ;
                options Estimation.Proprioception.EstimatorOptions;
            end
            obj.model_comp = model_comp;
            obj.options = options;
            
            obj.Qa = diag(sensors_dev.accel_noise.^2);
            obj.Qg = diag(sensors_dev.gyro_noise.^2);
            obj.Qcontact(1:3, 1:3) = diag(sensors_dev.contact_foot_linvel_noise.^2);
            obj.Qswing(1:3, 1:3) = diag(sensors_dev.swing_foot_linvel_noise.^2);
            obj.Qcontact(4:6, 4:6) = diag(sensors_dev.contact_foot_angvel_noise.^2);
            obj.Qswing(4:6, 4:6) = diag(sensors_dev.swing_foot_angvel_noise.^2);
            obj.Qba = diag(sensors_dev.accel_bias_noise.^2);
            obj.Qbg = diag(sensors_dev.gyro_bias_noise.^2);
            
            assert(options.nr_joints_est == length(sensors_dev.encoders_noise), 'Mismatch #joints in options');
            obj.encoders_prev = zeros(options.nr_joints_est, 1);
            obj.Renc = diag(sensors_dev.encoders_noise.^2);
            
            obj.R0fkin = diag(prior_dev.forward_kinematics.^2);
            
            obj.P0 = blkdiag(diag(prior_dev.imu_position.^2), ...
                             diag(prior_dev.imu_orientation.^2), ...                
                             diag(prior_dev.imu_linear_velocity.^2), ...
                             diag(prior_dev.left_foot_position.^2), ...
                             diag(prior_dev.left_foot_orientation.^2), ... 
                             diag(prior_dev.right_foot_position.^2), ...
                             diag(prior_dev.right_foot_orientation.^2));            
            if (options.enable_bias_estimation)
               obj.P0 = blkdiag(obj.P0, ...
                                diag(prior_dev.accel_bias.^2), ...
                                diag(prior_dev.gyro_bias.^2));                
            end
            
            obj.P_prev = obj.P0;
            
            obj.filter_configured = true;
            obj.filter_initialized = false;
            obj.bias_initialized = false;
        end % setup
        
        function obj = initialize(obj, X0)
            if (obj.filter_configured)
                if (obj.options.enable_bias_estimation)
                    assert(length(X0) == 20, 'State size mismatch');
                else
                    assert(size(X0, 2) == 13, 'State size mismatch');
                end
                
                obj.X0 = X0;
                obj.X = X0;
                obj.X_prev = X0;
                
                obj.P = obj.P0;
                
                obj.filter_initialized = true;
            end
        end
        
        function [X, P, BaseLinkState, DebugOut, obj] = advance(obj, t, omega, alpha, encoders, contacts)
            % TODO intialize bias static pose
            % obj.initializeBias(alpha, omega, obj.x0)
            obj.timestamp_id = obj.timestamp_id + 1;
            obj.X_prev = obj.X;
            obj.P_prev = obj.P;

            if (obj.filter_initialized && obj.t_prev > 0.0)
                dt = t - obj.t_prev;
                obj = obj.predictState(obj.omega_prev, obj.alpha_prev, obj.contacts_prev, dt);
                
                if obj.options.enable_ekf_update
                    if obj.options.enable_kinematic_meas
                        obj = obj.updateKinematics(encoders, contacts, dt);
                    end
                    
                    if obj.options.debug_mode
                        obj.debugger.access = true;
                    end
                end
            end
            
            % store latest values
            obj.omega_prev = omega;
            obj.alpha_prev = alpha;
            obj.encoders_prev = encoders;
            obj.contacts_prev = contacts;
            obj.t_prev = t;
            
            % outputs
            X = obj.X;
            P = obj.P;
            
            % extract base link state
            [R, p, v, ~, ~, ~, ~, ~, bg] = Estimation.DLGEKF_RIE.State.extract(obj.X);
            q = Utils.rot2quat(R);
            [BaseLinkState.q, BaseLinkState.I_p, BaseLinkState.I_pdot, BaseLinkState.I_omega] = obj.model_comp.getBaseStateFromIMUState(q, p, v, omega-bg);
            
            DebugOut = obj.debugger;
        end
    end
    
    methods (Access = private)
        function obj = predictState(obj, omega, alpha, contact, dt)
            [R, ~, v, ~, ~, ~, ~, ba, bg] = Estimation.DLGEKF_RIE.State.extract(obj.X);
            
            %- bias corrected IMU measurements
            alpha_unbiased = alpha - ba;
            omega_unbiased = omega - bg;
            acc = alpha_unbiased + (R'*obj.g);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%-- left parametrized velocity
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Omegap = R'*(v*dt) + 0.5*acc*(dt^2);
            OmegaR = omega_unbiased*dt;
            Omegav = acc*dt;
            Omega = [Omegap; OmegaR; Omegav; zeros(12, 1)];  % default
            
            if obj.options.enable_bias_estimation
                Omega = [Omega; zeros(6, 1)];
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%
            %%--propagate the mean
            %%%%%%%%%%%%%%%%%%%%%%%
            XPrior = obj.X;
            Omega_lifted = Estimation.DLGEKF_RIE.State.exphat(Omega);   
            obj.X = Estimation.DLGEKF_RIE.State.compose(XPrior, Omega_lifted);
                       
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%-- Left parametrized velocity Jacobian slantFk
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            I3dt = eye(3)*dt;
            gCrossdt = R'*Utils.skew(obj.g)*dt; 
            dOmega1_der = 0.5*gCrossdt*dt;
            dOmega1_dev = R'*dt;
            dOmega1_deba = -0.5*I3dt*dt;

            dOmega2_debg = -I3dt;

            dOmega3_der = gCrossdt;
            dOmega3_deba = -I3dt;

            slantFkbase = [zeros(3)    dOmega1_der dOmega1_dev; ...
                           zeros(3)       zeros(3)    zeros(3); ...
                           zeros(3)    dOmega3_der    zeros(3)];
            
            slantFk = [  slantFkbase    zeros(9, 12); ...
                        zeros(12, 9)   zeros(12, 12)];
            
            if obj.options.enable_bias_estimation
                slantFkbasebias = [     dOmega1_deba     zeros(3); ...
                                            zeros(3) dOmega2_debg; ...
                                        dOmega3_deba     zeros(3); ...
                                        zeros(12, 3) zeros(12, 3)];
                
                slantFk = [    slantFk   slantFkbasebias; ...
                          zeros(6, 21)      zeros(6, 6)];
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%-- Map the changes in the group space onto the velocity space
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            AdjX = Estimation.DLGEKF_RIE.State.AdjointMatrix(XPrior);
            leftJacobianOmega = Estimation.DLGEKF_RIE.State.leftJacobian(Omega);
            dimP = size(AdjX, 1);
            Fk = eye(dimP) + (AdjX*leftJacobianOmega*slantFk);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%-- Discrete system noise covariance
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Lk = blkdiag(-0.5*dt*dt*eye(3), -eye(3)*dt, -eye(3)*dt, ...
                         eye(3)*dt, eye(3)*dt, ...
                         eye(3)*dt, -eye(3)*dt);
            Qc = blkdiag(zeros(3), obj.Qg, obj.Qa,...
                         contact(1)*obj.Qcontact + (1-contact(1))*obj.Qswing, ...
                         contact(2)*obj.Qcontact + (1-contact(2))*obj.Qswing);
            if obj.options.enable_bias_estimation
                Lk = blkdiag(Lk, eye(3)*dt, eye(3)*dt);
                Qc = blkdiag(Qc, obj.Qba, obj.Qbg);
            end
            
            Qk = Lk*Qc*Lk';
            Q = AdjX*leftJacobianOmega*Qk*AdjX'*leftJacobianOmega';

            obj.P = Fk*obj.P*Fk'+ Q;
            obj.P = (obj.P + obj.P')./2;
            
            if obj.options.debug_mode
                obj.debugger.predicted_acc = acc;
                obj.debugger.F = Fk;
                obj.debugger.P_pred = obj.P;
                obj.debugger.Q = Qk;
            end
        end % predictState
        
        function obj = updateState(obj, deltaY, H, R)
            S = H*obj.P*H' + R;
            K = (obj.P*H')/S;
            deltaX = K*deltaY;
            
            deltaX_lifted = Estimation.DLGEKF_RIE.State.exphat(deltaX);

            obj.X = Estimation.DLGEKF_RIE.State.compose(deltaX_lifted, obj.X);
            
            I = eye(size(obj.P));
            leftJacobiandeltaX = Estimation.DLGEKF_RIE.State.leftJacobian(deltaX);
            obj.P = leftJacobiandeltaX*(I - (K*H))*obj.P*leftJacobiandeltaX';
            obj.P = (obj.P + obj.P')./2;
            
            if obj.options.debug_mode
                eigval = eig(obj.P);
                condn = max(eigval)/min(eigval);
                
                tol = length(eigval)*eps(max(eigval));
                isposdef = all(eigval > tol);
                
                obj.debugger.traceP = trace(obj.P);
                obj.debugger.condP = condn;
                obj.debugger.isPsymmetric = issymmetric(obj.P);
                obj.debugger.isPpositivedefinite = isposdef;
                obj.debugger.isPSPD = (issymmetric(obj.P) && isposdef);
                
                obj.debugger.H = H;
                obj.debugger.R = R;
                obj.debugger.K = K;
                obj.debugger.deltay = deltaY;
                obj.debugger.deltax = deltaX;
            end
            
        end % updateState
        
        function obj = updateKinematics(obj, encoders, contacts, dt)
            [A_R_IMU, p, ~, A_R_LF, pl, A_R_RF, pr, ~, ~] = Estimation.DLGEKF_RIE.State.extract(obj.X);
            A_H_IMU = LieGroups.SE3.constructSE3(A_R_IMU, p);
            IMU_H_A =   LieGroups.SE3.inverse(A_H_IMU);
            A_H_LF = LieGroups.SE3.constructSE3(A_R_LF, pl);
            A_H_RF = LieGroups.SE3.constructSE3(A_R_RF, pr);
            IMU_H_LF = LieGroups.SE3.compose(IMU_H_A, A_H_LF);
            IMU_H_RF = LieGroups.SE3.compose(IMU_H_A, A_H_RF);
            
            % Jacobian we get here allows us to express the feet velocities
            % in the feet frames, no need for additional transformations
            [y_LF, LF_J_IMULF] = obj.model_comp.relativeKinIMU_to_foot_contactExplicit(encoders, 'left');
            [y_RF, RF_J_IMURF] = obj.model_comp.relativeKinIMU_to_foot_contactExplicit(encoders, 'right');
            
            if (contacts(1) && contacts(2))
                %%%%%%%%%%%%%%%%%%%%%%%%
                %%-- double support
                %%%%%%%%%%%%%%%%%%%%%%%%
                % measurement model                
                h_of_x = LieGroups.SE3xSE3.constructSE3bySE3(IMU_H_LF, IMU_H_RF);
                
                % observation
                y = LieGroups.SE3xSE3.constructSE3bySE3(y_LF, y_RF);
                
                % innovation
                hinv = LieGroups.SE3xSE3.inverse(h_of_x);
                double_composite_pose_error =  LieGroups.SE3xSE3.compose(hinv, y);
                deltaY = LieGroups.SE3xSE3.logvee(double_composite_pose_error);

                Rk = blkdiag(LF_J_IMULF*obj.Renc*LF_J_IMULF', ...
                             RF_J_IMURF*obj.Renc*RF_J_IMURF');
                Rk = Rk/dt;
                
                % measurement model Jacobian
                ZlT_S_of_dl = (A_R_LF')*Utils.skew(pl);
                ZrT_S_of_dr = (A_R_RF')*Utils.skew(pr);
                H_lf = [-A_R_LF'    ZlT_S_of_dl   zeros(3)  A_R_LF'  -ZlT_S_of_dl  zeros(3) zeros(3); ...
                         zeros(3)      -A_R_LF'   zeros(3) zeros(3)       A_R_LF'  zeros(3) zeros(3)];

                H_rf = [-A_R_RF'    ZrT_S_of_dr   zeros(3)   zeros(3)  zeros(3)   A_R_RF'  -ZrT_S_of_dr; ...
                        zeros(3)       -A_R_RF'   zeros(3)   zeros(3)  zeros(3) zeros(3)        A_R_RF'];
                
                if obj.options.enable_bias_estimation
                    H_lf = [H_lf zeros(6)];
                    H_rf = [H_rf zeros(6)];
                end
                
                H = [H_lf; H_rf];
                obj = obj.updateState(deltaY, H, Rk);
            elseif (contacts(1))
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%-- single support left foot
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % measurement model                
                h_of_x = IMU_H_LF;
                
                % observation
                y = y_LF;
                
                % innovation
                hinv = LieGroups.SE3.inverse(h_of_x);
                pose_error =  LieGroups.SE3.compose(hinv, y);
                deltaY = LieGroups.SE3.logvee(pose_error);
                
                Rk = blkdiag(LF_J_IMULF*obj.Renc*LF_J_IMULF');
                Rk = Rk/dt;
                                
                % measurement model Jacobian               
                ZlT_S_of_dl = (A_R_LF')*Utils.skew(pl);                
                H_lf = [-A_R_LF'    ZlT_S_of_dl   zeros(3)  A_R_LF'  -ZlT_S_of_dl  zeros(3) zeros(3); ...
                         zeros(3)      -A_R_LF'   zeros(3) zeros(3)       A_R_LF'  zeros(3) zeros(3)];                

                if obj.options.enable_bias_estimation
                    H_lf = [H_lf zeros(6)];
                end
                
                H = H_lf;
                obj = obj.updateState(deltaY, H, Rk);
            elseif (contacts(2))
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%-- single support right foot
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % measurement model
                h_of_x = IMU_H_RF;
                
                % observation
                y = y_RF;
                
                % innovation
                hinv = LieGroups.SE3.inverse(h_of_x);
                pose_error =  LieGroups.SE3.compose(hinv, y);
                deltaY = LieGroups.SE3.logvee(pose_error);
                
                Rk = blkdiag(RF_J_IMURF*obj.Renc*RF_J_IMURF');
                Rk = Rk/dt;
                
                % measurement model Jacobian                                
                ZrT_S_of_dr = (A_R_RF')*Utils.skew(pr);
                H_rf = [-A_R_RF'    ZrT_S_of_dr   zeros(3)   zeros(3)  zeros(3)   A_R_RF'  -ZrT_S_of_dr; ...
                        zeros(3)        -A_R_RF'   zeros(3)   zeros(3)  zeros(3) zeros(3)      A_R_RF'];

                if obj.options.enable_bias_estimation
                    H_rf = [H_rf zeros(6)];
                end
                
                H = H_rf;
                obj = obj.updateState(deltaY, H, Rk);
            end
            
            if obj.options.debug_mode
                if sum(contacts) > 0
                    obj.debugger.y = y;
                    obj.debugger.z = h_of_x;                    
                end
            end
            
        end % updateKinematics                        
    end
    
end


