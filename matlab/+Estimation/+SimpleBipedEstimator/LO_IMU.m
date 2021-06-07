classdef LO_IMU < handle
    properties (SetAccess = private, GetAccess = private)
        timestamp_id (1, 1) uint64;
        LO (1, 1) Estimation.LeggedOdom.Filter;
        IMU;
        
        x (:, 1) double;
        x_prev (:, 1) double;
        
        imu_acc_prev (3, 1) double;
        imu_v_prev (3, 1) double;
        imu_p_prev (3, 1) double;
        
        t_prev (1, 1) double;
        
        model_comp;
        
        filter_configured (1, 1) logical;
        filter_initialized (1, 1) logical;
        
        w_R_wimu (3, 3) double;
        base_pose (4, 4) double;
        
        fuse_position (1, 1) logical;
        imu_aligned (1, 1) logical;
        attest_type (1, :) char;
    end
    
    properties (Access = private, Constant)
        g = [0; 0; -9.80665];
    end
    
    methods
        function obj = LO_IMU(attest_type)
            obj.timestamp_id = 0;
            
            obj.LO = Estimation.LeggedOdom.Filter();
            
            switch attest_type
                case 'qekf'
                    obj.IMU =  Estimation.AttitudeEstimator.QEKF();
                case 'mahony'
                    obj.IMU =  Estimation.AttitudeEstimator.MahonyFilter();
                otherwise
                    disp('invalid choice of attitude estimator. constructor failed');
            end
            obj.attest_type = attest_type;
            
            obj.imu_acc_prev = zeros(3, 1);
            obj.imu_v_prev = zeros(3, 1);
            obj.imu_p_prev = zeros(3, 1);
            
            obj.w_R_wimu = eye(3);
            obj.base_pose = eye(4);
            
            obj.fuse_position = false;
            
            obj.x = zeros(27, 1);
            obj.x_prev = zeros(27, 1);
            
            obj.filter_configured = false;
            obj.filter_initialized = false;
            obj.imu_aligned = false;
        end
        
        function obj = setup(obj, model_comp, primary_foot, flat_floor, switching_pattern, attestparams, fuse_position)
            obj.LO = obj.LO.setup(model_comp, primary_foot, flat_floor, switching_pattern);
            
            switch class(attestparams)
                case 'iDynTree.AttitudeMahonyFilterParams'
                    assert(strcmp(obj.attest_type, 'mahony'), 'Improper attitude estimator params');
                case 'iDynTree.AttitudeQuaternionEKFParameters'
                    assert(strcmp(obj.attest_type, 'qekf'), 'Improper attitude estimator params');
            end
            obj.IMU = obj.IMU.setup(attestparams);
            
            obj.fuse_position = fuse_position;
            
            obj.model_comp = model_comp;
            obj.filter_configured = true;
        end
        
        % initialize
        function obj = initialize(obj, ref_frame_for_world, ref_H_w0, initial_joint_pos, attestx0)
            obj.LO = obj.LO.initialize(ref_frame_for_world, ref_H_w0, initial_joint_pos);
            obj.IMU = obj.IMU.initialize(attestx0);
            obj.filter_initialized = true;
        end
        
        % advance
        function [x, BaseLinkState, rpyBaseIMUinLOWorld, rpyBaseLO, obj] = advance(obj, t, alpha, omega, mag, encoders, encoder_speeds, contacts, lf_wrench, rf_wrench)
            obj.timestamp_id = obj.timestamp_id + 1;
            obj.x_prev = obj.x;
            
            if (obj.filter_initialized && obj.t_prev > 0.0)
                [xLO, ~, ~, obj.LO] = obj.LO.advance(t, encoders, encoder_speeds,...
                    contacts, lf_wrench, rf_wrench);
                obj.x = xLO;
                rpyBaseLO = Utils.rot2rpy(Utils.quat2rot(xLO(1:4)));
                
                [rpyBaseIMU, ~, obj.IMU] = obj.IMU.advance(t, alpha, omega, mag);
                
                % wait for IMU to initialize
                if (~obj.imu_aligned)
                    obj = obj.alignIMU(xLO, rpyBaseIMU);
                    obj.x_prev = obj.x;
                end
                
                dt = t - obj.t_prev;
                if (obj.imu_aligned)
                    % fuse
                    [rpyBaseIMUinLOWorld, obj] = obj.fusePoses(xLO, rpyBaseIMU, alpha, dt);
                    obj = obj.fuseVelocities(omega, encoders, encoder_speeds);
                else
                    rpyBaseIMUinLOWorld = zeros(3, 1);
                end
            else
                rpyBaseLO = zeros(3, 1);
                rpyBaseIMUinLOWorld = zeros(3, 1);
            end
            
            % store latest values
            obj.t_prev = t;
            
            % outputs
            x = obj.x;
                                    
            % extract base link state
            [q, p, v, w_omega_imu, ~, ~, ~, ~] = obj.extract(obj.x);
            R = Utils.quat2rot(q);
            omega_imu = R*w_omega_imu;
            [BaseLinkState.q, BaseLinkState.I_p, BaseLinkState.I_pdot, BaseLinkState.I_omega] = obj.model_comp.getBaseStateFromIMUState(q, p, v, omega_imu);
            
            Rb = Utils.quat2rot(BaseLinkState.q);
            obj.base_pose = [Rb BaseLinkState.I_p; zeros(1, 3) 1];
        end
        
        function [q, p, v, omega, qlf, plf, qrf, prf] = extract(obj, x)
            [q, p, v, omega, qlf, plf, qrf, prf] = obj.LO.extract(x);
        end
    end
    
    methods (Access = private)
        function obj = alignIMU(obj, xLO, rpyBaseIMU)
            wimu_R_imu0 = Utils.rpy2rot(rpyBaseIMU(1), rpyBaseIMU(2), rpyBaseIMU(3));
            [qimuLO, ~, ~, ~, ~, ~, ~, ~] = obj.extract(xLO);
            w_R_imu0 = Utils.quat2rot(qimuLO);
            obj.w_R_wimu = w_R_imu0*(wimu_R_imu0');

            obj.imu_aligned = true;
        end
        
        function [rpyBaseIMUinLOWorld, obj] = fusePoses(obj, xLO, rpyBaseIMU, alpha, dt)
            % fuse rotations
            [qimuLO, pimuLO, ~, ~, ~, ~, ~, ~] = obj.extract(xLO);
            w_R_imuLO = Utils.quat2rot(qimuLO);
            
            % rotate rotation from LO to IMU world and replace yaw
            wimu_R_imuLO = (obj.w_R_wimu')*w_R_imuLO;
            rpyimuLO = Utils.rot2rpy(wimu_R_imuLO);
            wimu_R_imu = Utils.rpy2rot(rpyBaseIMU(1), rpyBaseIMU(2), rpyimuLO(3));
            % rotate back to LO world frame (which is the estimator world)
            w_R_imuIMU = obj.w_R_wimu*wimu_R_imu;
            rpyBaseIMUinLOWorld = Utils.rot2rpy(w_R_imuIMU);
            
            weights = [0.5 0.5];
            Rarray = {w_R_imuLO, w_R_imuIMU};
            tol = 1e-3;
            max_iter = 1000;
            
            Rfused = LieGroups.SO3.geodesicL2WeightedMeanRotation(Rarray, weights, tol, dt, max_iter);
            
            qFused = Utils.rot2quat(Rfused);
            obj.x(1:4) = qFused;
            
            % fuse positions
            if obj.fuse_position
                p_prev = obj.x_prev(5:7);
                v_prev = obj.x_prev(8:10);
                acc = Rfused*alpha + obj.g;
                pimuIMU = p_prev + (v_prev*dt) + (0.5*acc*dt*dt);
                
                weight_imu_pos = 0.5;
                pFused = ((1 - weight_imu_pos)*pimuLO) + (weight_imu_pos*pimuIMU);
                obj.x(5:7) = pFused;
            end
        end
        
        function obj = fuseVelocities(obj, omega, encoders, encoder_speeds)
            switch obj.LO.primary_foot
                case 'left'
                    feet_frame = obj.model_comp.LFVertexIds(1);
                case 'right'
                    feet_frame = obj.model_comp.RFVertexIds(1);
                otherwise
                    feet_frame = '';
            end
            
            if isempty(feet_frame)
                JF = [];
                vF = [];
            else                
                JF = obj.model_comp.getFrameFreeFloatingJacobian(feet_frame, obj.base_pose, encoders, encoder_speeds);
                vF = zeros(6, 1);
            end
            
            JIMU = obj.model_comp.getFrameFreeFloatingJacobian(obj.model_comp.base_link_imu, obj.base_pose, encoders, encoder_speeds);
            W_R_IMU = Utils.quat2rot(obj.x(1:4));
            W_omega = W_R_IMU*omega;
            
            nr_joints = length(encoders);
            B = [zeros(nr_joints, 6) eye(nr_joints)];
            
            y = [vF; W_omega; encoder_speeds];
            A = [JF; JIMU(4:6, :); B];
            
            W = eye(length(y));
            reg = 1e-6*eye(6+nr_joints);
            
            qdot_optimum = (A'*W*A + reg)\A'*W*y;
            v_b = qdot_optimum(1:6);

%             Jc_base_inv = inv(JF(1:6, 1:6));
%             Jc_enc = JF(1:6, 7:end);
%             
%             v_b = -Jc_base_inv*Jc_enc*encoder_speeds;
            
            % output expressed as IMU velocity - transform from base to IMU
            % (right trivialized)
            imu_p_b = obj.model_comp.kindyn.getRelativeTransform(obj.model_comp.base_link, obj.model_comp.base_link_imu).inverse().getPosition().toMatlab();
            imu_p_b_dyn = iDynTree.Position();
            imu_p_b_dyn.fromMatlab(imu_p_b);
            w_R_imu_dyn = iDynTree.Rotation();
            w_R_imu_dyn.fromMatlab(W_R_IMU)
            X = iDynTree.Transform(iDynTree.Rotation.Identity(), w_R_imu_dyn*imu_p_b_dyn).asAdjointTransform().toMatlab();
            
            v_imu = X*v_b;
            v_imu(4:6) = W_omega;
            obj.x(8:13) = v_imu;
        end
        
    end
end



