classdef Filter < handle
    properties (Access = private)
        timestamp_id (1, 1) uint64;
        
        x (:, 1) double;
        P (:, :) double;
        
        t_prev (1, 1) double;
        x_prev (:, 1) double;
        P_prev (:, :) double;
        alpha_prev (3, 1) double;
        omega_prev (3, 1) double;
        encoders_prev (:, 1) double;
        contacts_prev (2, 1) logical;
        
        x0 (:, 1) double;
        P0 (:, :) double;
        
        ba0 (3, 1) double;
        bg0 (3, 1) double;
        
        %linearization points
        omegaLin (3, 1) double;
        acc2Lin (3, 1) double;
        plfLin (3, 1) double;
        prfLin (3, 1) double;
        thetarfLin (3, 1) double;
        thetalfLin (3, 1) double;
        
        
        Qg (3, 3) double; % gyro variance
        Qa (3, 3) double; % accelerometer variance
        Qcontact (6, 6) double; % contact foot pose variance
        Qswing (6, 6) double; % swing foot pose variance
        Qba (3, 3) double; % gyro bias variance
        Qbg (3, 3) double; % accelerometer bias variance
        
        Renc (:, :) double; % encoder variance
        R0fkin (6, 6) double; % forward kinematic pose variance
        
        model_comp;
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
            
            %linearization points
            obj.omegaLin = zeros(3, 1);
            obj.acc2Lin = zeros(3, 1);
            obj.plfLin = zeros(3, 1);
            obj.prfLin = zeros(3, 1);
            obj.thetarfLin = zeros(3, 1);
            obj.thetalfLin = zeros(3, 1);
            
            obj.ba0 = zeros(3, 1);
            obj.bg0 = zeros(3, 1);
            
            obj.options = Estimation.Proprioception.EstimatorOptions();
            obj.debugger = Estimation.Proprioception.DebugOutputs();
            
            obj.filter_configured = false;
            obj.filter_initialized = false;
            obj.bias_initialized = false;
        end
        
        function obj = setup(obj, prior_dev, sensors_dev, model_comp, options)
            arguments
                obj Estimation.RotellaEstimator.Filter;
                prior_dev Estimation.Proprioception.PriorsStdDev;
                sensors_dev Estimation.Proprioception.SensorsStdDev;
                model_comp ;
                options Estimation.Proprioception.EstimatorOptions;
            end
            obj.model_comp = model_comp;
            obj.options = options;
            
            obj.Qa = diag(sensors_dev.accel_noise.^2);
            obj.Qg = diag(sensors_dev.gyro_noise.^2);
            obj.Qcontact(1:3, 1:3) = diag(sensors_dev.contact_foot_angvel_noise.^2);
            obj.Qswing(1:3, 1:3) = diag(sensors_dev.swing_foot_angvel_noise.^2);
            obj.Qcontact(4:6, 4:6) = diag(sensors_dev.contact_foot_linvel_noise.^2);
            obj.Qswing(4:6, 4:6) = diag(sensors_dev.swing_foot_linvel_noise.^2);
            obj.Qba = diag(sensors_dev.accel_bias_noise.^2);
            obj.Qbg = diag(sensors_dev.gyro_bias_noise.^2);
            
            assert(options.nr_joints_est == length(sensors_dev.encoders_noise), 'Mismatch #joints in options');
            obj.encoders_prev = zeros(options.nr_joints_est, 1);
            obj.Renc = diag(sensors_dev.encoders_noise.^2);
            
            obj.R0fkin = diag(prior_dev.forward_kinematics.^2);
            obj.P0 = Estimation.RotellaEstimator.State.constructVar(prior_dev, options.enable_bias_estimation);
            
            obj.P_prev = obj.P0;
            
            obj.filter_configured = true;
            obj.filter_initialized = false;
            obj.bias_initialized = false;
        end % setup
        
        function obj = initialize(obj, x0, ba0, bg0)
            if (obj.filter_configured)
                if (obj.options.enable_bias_estimation)
                    assert(length(x0) == 30, 'State size mismatch');
                    state_type = 'nonlinear';
                    ba_range = Estimation.RotellaEstimator.State.getRange('bias_acc', state_type, obj.options.enable_bias_estimation);
                    bg_range = Estimation.RotellaEstimator.State.getRange('bias_gyro', state_type, obj.options.enable_bias_estimation);
                    obj.ba0 = x0(ba_range);
                    obj.bg0 = x0(bg_range);
                else
                    assert(length(x0) == 24, 'State size mismatch');
                end
                
                obj.x0 = x0;
                obj.x = x0;
                obj.x_prev = x0;
                
                obj.ba0 = ba0;
                obj.bg0 = bg0;
                
                obj.P = obj.P0;
                
                obj.filter_initialized = true;
            end
        end
        
        function [x, P, BaseLinkState, DebugOut, obj] = advance(obj, t, omega, alpha, encoders, contacts)
            % TODO intialize bias static pose
            % obj.initializeBias(alpha, omega, obj.x0)
            obj.timestamp_id = obj.timestamp_id + 1;
            obj.x_prev = obj.x;
            obj.P_prev = obj.P;
            
            if (obj.filter_initialized && obj.t_prev > 0.0)
                dt = t - obj.t_prev;
                obj = obj.predictState(obj.omega_prev, obj.alpha_prev, obj.encoders_prev, obj.contacts_prev, dt);

                if obj.options.enable_ekf_update
                    if obj.options.enable_kinematic_meas || obj.options.always_enable_kinematic_meas
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
            x = obj.x;
            P = obj.P;
            
            % extract base link state
            [q, p, v, ~, ~, ~, ~, ~, bg] = Estimation.RotellaEstimator.State.extract(obj.x, obj.ba0, obj.bg0);
            [BaseLinkState.q, BaseLinkState.I_p, BaseLinkState.I_pdot, BaseLinkState.I_omega] = obj.model_comp.getBaseStateFromIMUState(q, p, v, omega-bg);
            
            DebugOut = obj.debugger;
        end
    end
    
    methods (Access = private)
        function obj = predictState(obj, omega, alpha, encoders, contact, dt)
            [q, p, v, qlf, plf, qrf, prf, ba, bg] = Estimation.RotellaEstimator.State.extract(obj.x, obj.ba0, obj.bg0);
            R = Utils.quat2rot(q);
            
            %- bias corrected IMU measurements
            alpha_unbiased = alpha - ba;
            omega_unbiased = omega - bg;
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            %%-- nonlinear dynamics
            %%%%%%%%%%%%%%%%%%%%%%%%

            %- base pose dynamics
            quat_correction = Utils.UnitQuaternion.expQuaternion(omega_unbiased*dt);
            q_pred = Utils.UnitQuaternion.composeQuaternion(q, quat_correction);
            if (abs(norm(q_pred) - 0.0) < 10e-3)
                disp('[WARNING][RotellaEstimator::predictState] Malformed quaternion, returning previous state');
                return;
            end
            
            if (abs(norm(q_pred) - 1.0) > 10e-3)
                normalize(q_pred);
            end
            R_pred = Utils.quat2rot(q_pred);
            
            %-- taking the predicted base acceleration at the midpoint of two time
            %-- instants for integrating to get position and velocity
            %-- this might help avoid any errors due to poor rotation estiamte
            %-- used to rotate g into inertial frame
            %             acc = (0.5*(R*alpha_unbiased)) + obj.g;
            %             acc = acc + (0.5*(R_pred*alpha_unbiased));
            %-- otherwise comment the above two lines and uncomment the below line,
            acc = R*alpha_unbiased + obj.g;
            
            v_pred = v + (dt * acc);
            p_pred = p + (dt * v) + ((dt^2)*0.5*acc);
            
            %- foot pose dynamics                        
            qlf_pred = qlf;
            plf_pred = plf;
            
            
            qrf_pred = qrf;
            prf_pred = prf;
            
            %- bias dynamics
            ba_pred = ba;
            bg_pred = bg;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %%-- Linearized dynamics
            %%%%%%%%%%%%%%%%%%%%%%%%%

            obj.omegaLin = omega_unbiased;
            obj.acc2Lin = alpha_unbiased;
            
            
            %- Choose linearization point based on obsevability
            %- constraint: compute the linearization point
            %- from the previous state and current state
            %- this will keep the unobservable directions unobservable
            %- and avoid it from spuriously becoming observable due to
            %- discretization and linearization errors and due to effects of
            %- noise. This will help avoid inconsistencies in the filter.
            %- comment the lines below, if observability constraint need
            %- not be activated
            [q_prev, p_prev, v_prev, ~, plf_prev, ~, prf_prev, ~, ~] = Estimation.RotellaEstimator.State.extract(obj.x_prev, obj.ba0, obj.bg0);
            
            obj.plfLin = plf_prev;
            obj.prfLin = prf_prev;
            
            halfdtsquared = 0.5*(dt^2);
            dtvminus = dt*v_prev;
            q_prev_inv = Utils.UnitQuaternion.inverse(q_prev);
            q_err = Utils.UnitQuaternion.composeQuaternion(q_prev_inv, q_pred);
            
            % observability constrained linearization
            if obj.timestamp_id > 2
                obj.omegaLin = Utils.UnitQuaternion.logQuaternion(q_err);
                obj.acc2Lin = R'*( ( (v_pred - v_prev)./dt) - obj.g );
                
                if (contact(1))
                    obj.plfLin = plf_pred;
                end
                
                if (contact(2))
                    obj.prfLin = prf_pred;                
               end                                              
            end
                        
            % base
            omegaLinCross = Utils.skew(obj.omegaLin);
            acc2LinCross = Utils.skew(obj.acc2Lin);
            
            %- build Fc
            % serialization:  phi p  v  phi_lf p_lf phi_rf p_rf ba bg
            Fc_base_base = [-omegaLinCross   zeros(3)  zeros(3); ...
                zeros(3)     zeros(3)    eye(3); ...
                -R*acc2LinCross   zeros(3)  zeros(3)];
            
            Fc_base_wo_bias = [Fc_base_base zeros(9, 12)];
            Fc_foot_wo_bias = zeros(12, 21);
            
            if obj.options.enable_bias_estimation
                Fc_bias_base = [zeros(3)     -eye(3); ...
                    zeros(3)   zeros(3); ...
                    -R   zeros(3)];
                
                Fc = [Fc_base_wo_bias Fc_bias_base; ...
                    Fc_foot_wo_bias  zeros(12,6); ...
                    zeros(6, 21)     zeros(6)];
            else
                Fc = [Fc_base_wo_bias; ...
                    Fc_foot_wo_bias];
            end
            
            % obtain Fk using zero-order hold on Fc
            linearized_state_dim = size(Fc, 1);
            Fk = eye(linearized_state_dim) + (Fc*dt); % first order Taylor's expansion of Fk = exp(Fc dt)
            
            % build Qk from Lc Fk and Qc
            Lc = zeros(21, 18);
            Lc(1:3, 1:3) = -eye(3);
            Lc(7:21, 4:18) = blkdiag(-R, eye(3), R, eye(3), R);
            if obj.options.enable_bias_estimation
                Lc = blkdiag(Lc, eye(3), eye(3));
            end
            
            Qc = blkdiag(obj.Qg, obj.Qa,...
                contact(1)*obj.Qcontact(1:3, 1:3) + (1-contact(1))*obj.Qswing(1:3, 1:3), ...
                contact(1)*(R*obj.Qcontact(4:6, 4:6)*R') + (1-contact(1))*(R*obj.Qswing(4:6, 4:6)*R'), ...
                contact(1)*obj.Qcontact(1:3, 1:3) + (1-contact(1))*obj.Qswing(1:3, 1:3), ...
                contact(2)*(R*obj.Qcontact(4:6, 4:6)*R') + (1-contact(2))*(R*obj.Qswing(4:6, 4:6))*R');
            if obj.options.enable_bias_estimation
                Qc = blkdiag(Qc, obj.Qba, obj.Qbg);
            end

            Qk = (Fk*Lc*Qc*Lc'*Fk')*dt;
            
            if obj.options.enable_bias_estimation
                obj.x = [q_pred; p_pred; v_pred; ...
                    qlf_pred; plf_pred; ...
                    qrf_pred; prf_pred; ...
                    ba_pred; bg_pred];
            else
                obj.x = [q_pred; p_pred; v_pred; ...
                    qlf_pred; plf_pred; ...
                    qrf_pred; prf_pred];
            end
            obj.P = Fk*obj.P*Fk' + Qk;
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
            [q, p, v, ql, pl, qr, pr, ba, bg] = Estimation.RotellaEstimator.State.extract(obj.x, obj.ba0, obj.bg0);
            [dq, dp, dv, dql, dpl, dqr, dpr, dba, dbg] = Estimation.RotellaEstimator.State.extractLinearizedState(deltaX);
            
            q_upd = Utils.UnitQuaternion.composeQuaternion(q, Utils.UnitQuaternion.expQuaternion(dq));
            p_upd = p + dp;
            v_upd = v + dv;
            ql_upd = Utils.UnitQuaternion.composeQuaternion(ql, Utils.UnitQuaternion.expQuaternion(dql));
            pl_upd = pl + dpl;
            qr_upd = Utils.UnitQuaternion.composeQuaternion(qr, Utils.UnitQuaternion.expQuaternion(dqr));
            pr_upd = pr + dpr;
            
            if (obj.options.enable_bias_estimation)
                ba_upd = ba + dba;
                bg_upd = bg + dbg;
                obj.x = Estimation.RotellaEstimator.State.constructwithbias(q_upd, p_upd, v_upd, ...
                    ql_upd, pl_upd, ...
                    qr_upd, pr_upd, ...
                    ba_upd, bg_upd);
            else
                obj.x = Estimation.RotellaEstimator.State.constructwithoutbias(q_upd, p_upd, v_upd, ...
                    ql_upd, pl_upd, ...
                    qr_upd, pr_upd);
            end
            
            I = eye(size(obj.P));
            
            obj.P = (I - (K*H))*obj.P;
            obj.P = (obj.P + obj.P')./2;
            
            if obj.options.debug_mode
                [U, S, V] = svd(obj.P);
                eigval = diag(S);                
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
            [q, p, ~, ql, pl, qr, pr, ~, ~] = Estimation.RotellaEstimator.State.extract(obj.x, obj.ba0, obj.bg0);
            A_R_IMU = Utils.quat2rot(q);
            A_R_LF = Utils.quat2rot(ql);
            A_R_RF = Utils.quat2rot(qr);
            
            [IMU_q_LF, IMU_p_LF, LF_J_IMULF] = obj.model_comp.relativeKinIMU_to_foot_contactLeftTriv(encoders, 'left');
            [IMU_q_RF, IMU_p_RF, RF_J_IMURF] = obj.model_comp.relativeKinIMU_to_foot_contactLeftTriv(encoders, 'right');
            IMU_R_LF = Utils.quat2rot(IMU_q_LF);
            IMU_R_RF = Utils.quat2rot(IMU_q_RF);

            if (contacts(1) && contacts(2))
                %%%%%%%%%%%%%%%%%%%%%%%%
                %%-- double support
                %%%%%%%%%%%%%%%%%%%%%%%%
                H =  [          -A_R_LF'*A_R_IMU    zeros(3)  zeros(3)   eye(3)  zeros(3) zeros(3) zeros(3); ...
                    Utils.skew(A_R_IMU'*(obj.plfLin - p))   -A_R_IMU'   zeros(3) zeros(3) A_R_IMU' zeros(3) zeros(3); ...
                                -A_R_RF'*A_R_IMU   zeros(3)  zeros(3) zeros(3) zeros(3)   eye(3) zeros(3); ...
                    Utils.skew(A_R_IMU'*(obj.prfLin - p))   -A_R_IMU'  zeros(3) zeros(3) zeros(3) zeros(3) A_R_IMU'];                             

                if (obj.options.enable_bias_estimation)
                    H = [H zeros(12, 6)];
                end
                
                Rk = blkdiag(LF_J_IMULF(4:6, :)*obj.Renc*LF_J_IMULF(4:6, :)', ...
                             LF_J_IMULF(1:3, :)*obj.Renc*LF_J_IMULF(1:3, :)', ...
                             RF_J_IMURF(4:6, :)*obj.Renc*RF_J_IMURF(4:6, :)', ...
                             RF_J_IMURF(1:3, :)*obj.Renc*RF_J_IMURF(1:3, :)');
                Rk = Rk/dt;
                
                zqlf = Utils.UnitQuaternion.composeQuaternion(Utils.UnitQuaternion.inverse(q), ql);
                zqlf_inv = Utils.UnitQuaternion.inverse(zqlf);
                zplf = A_R_IMU'*(pl - p);
                zqrf = Utils.UnitQuaternion.composeQuaternion(Utils.UnitQuaternion.inverse(q), qr);
                zqrf_inv = Utils.UnitQuaternion.inverse(zqrf);
                zprf = A_R_IMU'*(pr - p);
                
                residualqlf = Utils.UnitQuaternion.composeQuaternion(zqlf_inv, IMU_q_LF);
                residualqrf = Utils.UnitQuaternion.composeQuaternion(zqrf_inv, IMU_q_RF);
                
                deltaY = [Utils.UnitQuaternion.logQuaternion(residualqlf); ...
                    IMU_p_LF - zplf; ...
                    Utils.UnitQuaternion.logQuaternion(residualqrf); ...
                    IMU_p_RF - zprf];
                
                obj = obj.updateState(deltaY, H, Rk);
            elseif (contacts(1))
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%-- single support left foot
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                H = [          -A_R_LF'*A_R_IMU    zeros(3)   zeros(3)   eye(3)  zeros(3) zeros(3) zeros(3); ...
                    Utils.skew(A_R_IMU'*(obj.plfLin - p))   -A_R_IMU'   zeros(3) zeros(3) A_R_IMU' zeros(3) zeros(3)];
                
                if (obj.options.enable_bias_estimation)
                    H = [H zeros(6, 6)];
                end
                
                Rk = blkdiag(LF_J_IMULF(4:6, :)*obj.Renc*LF_J_IMULF(4:6, :)', ...
                             LF_J_IMULF(1:3, :)*obj.Renc*LF_J_IMULF(1:3, :)');
                Rk = Rk/dt;
                
                zqlf = Utils.UnitQuaternion.composeQuaternion(Utils.UnitQuaternion.inverse(q), ql);
                zqlf_inv = Utils.UnitQuaternion.inverse(zqlf);
                zplf = A_R_IMU'*(pl - p);
                
                residualqlf = Utils.UnitQuaternion.composeQuaternion(zqlf_inv, IMU_q_LF);
                
                deltaY = [Utils.UnitQuaternion.logQuaternion(residualqlf); ...
                    IMU_p_LF - zplf];
                obj = obj.updateState(deltaY, H, Rk);
            elseif (contacts(2))
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%-- single support right foot
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                H = [          -A_R_RF'*A_R_IMU   zeros(3)  zeros(3) zeros(3) zeros(3)   eye(3) zeros(3); ...
                    Utils.skew(A_R_IMU'*(obj.prfLin - p))   -A_R_IMU'  zeros(3) zeros(3) zeros(3) zeros(3) A_R_IMU'];
                
                if (obj.options.enable_bias_estimation)
                    H = [H zeros(6, 6)];
                end
                
                Rk = blkdiag(RF_J_IMURF(4:6, :)*obj.Renc*RF_J_IMURF(4:6, :)', ...
                             RF_J_IMURF(1:3, :)*obj.Renc*RF_J_IMURF(1:3, :)');
                Rk = Rk/dt;
                
                zqrf = Utils.UnitQuaternion.composeQuaternion(Utils.UnitQuaternion.inverse(q), qr);
                zqrf_inv = Utils.UnitQuaternion.inverse(zqrf);
                zprf = A_R_IMU'*(pr - p);
                
                residualqrf = Utils.UnitQuaternion.composeQuaternion(zqrf_inv, IMU_q_RF);
                
                deltaY = [Utils.UnitQuaternion.logQuaternion(residualqrf); ...
                    IMU_p_RF - zprf];
                
                obj = obj.updateState(deltaY, H, Rk);
            end
            
            if obj.options.debug_mode
                if (contacts(1) && contacts(2))
                    obj.debugger.y = [IMU_q_LF; IMU_p_LF; IMU_q_RF; IMU_p_RF];
                    obj.debugger.z = [zqlf; zplf; zqrf; zprf];
                elseif contacts(1)
                    obj.debugger.y = [IMU_q_LF; IMU_p_LF];
                    obj.debugger.z = [zqlf; zplf];
                elseif contacts(2)
                    obj.debugger.y = [IMU_q_RF; IMU_p_RF];
                    obj.debugger.z = [zqrf; zprf];
                end
            end
            
        end % updateKinematics
        
    end
    
end


