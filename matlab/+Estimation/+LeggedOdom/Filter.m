classdef Filter < handle
    properties (SetAccess = private, GetAccess = private)
        timestamp_id (1, 1) uint64;
        leggedodom (1, 1) iDynTree.SimpleLeggedOdometry;
        
        x (:, 1) double;
        x_prev (:, 1) double;
        
        prev_fixed_frame;        
        fixed_frame;
        
        w_H_b (4, 4) double;
        w_H_b_prev (4, 4) double;
        w_H_imu (4, 4) double;
        w_H_lf (4, 4) double;
        w_H_rf (4, 4) double;
        v_imu (6, 1) double;
        v_b (6, 1) double;
        
        t_prev (1, 1) double;

        encoders_prev (:, 1) double;
        contacts_prev (2, 1) logical; 
        lf_wrench_prev (6, 1) double;
        rf_wrench_prev (6, 1) double;
        
        model_comp ;
        joint_list;
                
        filter_configured (1, 1) logical;
        filter_initialized (1, 1) logical;
        
        switching_pattern (1, :) char;
        flat_floor    (1, 1) logical;
        
        % ros attributes
        joint_publisher;
        lfw_publisher;
        rfw_publisher;
        tftree;
    end
    
    properties (SetAccess = private, GetAccess = public)
        primary_foot;
        soleheight_from_floor = 0.0108;
        to_ros = false;
    end
    
    methods
        function obj = Filter()
            obj.timestamp_id = 0;
            
            obj.leggedodom = iDynTree.SimpleLeggedOdometry();
                       
            obj.prev_fixed_frame = 'unknown';
            obj.primary_foot = 'unknown';            
            obj.fixed_frame = 'unknown';    
            
            obj.w_H_b = eye(4);
            obj.w_H_b_prev = eye(4);
            obj.w_H_imu = eye(4);
            obj.w_H_lf = eye(4);
            obj.w_H_rf = eye(4);
            obj.v_imu = zeros(6, 1);
            obj.v_b = zeros(6, 1);
            
            obj.contacts_prev = false(2, 1);
            obj.lf_wrench_prev = zeros(6, 1);
            obj.rf_wrench_prev = zeros(6, 1);
            obj.switching_pattern = 'alternate';
            obj.flat_floor = false;
            
            
            if (obj.to_ros)
                obj.tftree = rostf;    
                obj.joint_publisher = rospublisher('/joint_states', 'sensor_msgs/JointState');
                obj.lfw_publisher = rospublisher('/lf_wrench', 'geometry_msgs/WrenchStamped');
                obj.rfw_publisher = rospublisher('/rf_wrench', 'geometry_msgs/WrenchStamped');
            end  
        end
        
        function obj = setup(obj, model_comp, primary_foot, flat_floor, switching_pattern)
            switch switching_pattern
                case {'alternate'}
                otherwise
                    obj.filter_configured = false;
                    disp('Specified switching pattern not available, please configure with proper switching pattern');
                    return
            end
            obj.switching_pattern = switching_pattern;                        
            
            obj.model_comp = model_comp;
            obj.primary_foot = primary_foot;
            switch primary_foot
                case 'right'
                    obj.fixed_frame = obj.model_comp.RFVertexIds(1);
                case 'left'
                    obj.fixed_frame = obj.model_comp.LFVertexIds(1);
                otherwise
                    obj.fixed_frame = 'unknown';
            end
            
            obj.flat_floor = flat_floor;
                        
            model = model_comp.kindyn.model.copy();            
            ok = obj.leggedodom.setModel(model);
            if ~ok
                disp('Could not set model');
                return
            end
            
            nr_joints = model.getNrOfJoints();
            obj.joint_list = cell(nr_joints, 1);
            for idx = 1:nr_joints
                obj.joint_list{idx} = model.getJointName(idx-1);
            end
            
            obj.filter_configured = true;
            % set feet polygon
        end
                        
        % initialize
        function obj = initialize(obj, ref_frame_for_world, ref_H_w0, initial_joint_pos)            
            Hrot = iDynTree.Rotation();
            Hrot.fromMatlab(ref_H_w0(1:3, 1:3));
            Hpos = iDynTree.Position();
            Hpos.fromMatlab(ref_H_w0(1:3, 4));
            
            Hdyn = iDynTree.Transform(Hrot, Hpos);
            obj.leggedodom.init(ref_frame_for_world, Hdyn);            
            
            jposDyn = iDynTree.JointPosDoubleArray();
            jposDyn.resize(length(initial_joint_pos));
            jposDyn.fromMatlab(initial_joint_pos);
                
            obj.leggedodom.updateKinematics(jposDyn);
            
            if (obj.flat_floor)
                w_H_fixed = obj.leggedodom.getWorldFrameTransform(obj.leggedodom.model.getFrameIndex(obj.fixed_frame));
                w_p_fixed = w_H_fixed.getPosition();
                w_p_fixed.setVal(2, obj.soleheight_from_floor); % height of sole from floor
                w_H_fixed.setPosition(w_p_fixed);
                obj.leggedodom.changeFixedFrame(obj.fixed_frame, w_H_fixed);
            end
            
            obj.w_H_b = obj.leggedodom.getWorldLinkTransform(obj.model_comp.base_link_idx).asHomogeneousTransform().toMatlab();                        
            obj.w_H_imu = obj.leggedodom.getWorldFrameTransform(obj.model_comp.base_link_imu_idx).asHomogeneousTransform().toMatlab();                        
            obj.w_H_lf = obj.leggedodom.getWorldFrameTransform(obj.model_comp.LFVertexIds(1)).asHomogeneousTransform().toMatlab();
            obj.w_H_rf = obj.leggedodom.getWorldFrameTransform(obj.model_comp.RFVertexIds(1)).asHomogeneousTransform().toMatlab();
            
            obj.x = obj.constructStateVec();
                                  
            obj.filter_initialized = true;
        end
        
        % advance
        function [x, BaseLinkState, fixed_frame, obj] = advance(obj, t, encoders, encoder_speeds, contacts, lf_wrench, rf_wrench)
            obj.timestamp_id = obj.timestamp_id + 1;
            obj.x_prev = obj.x;            
            
            if (obj.filter_initialized && obj.t_prev > 0.0)
                dt = t - obj.t_prev;
                obj = obj.getPrimaryFoot(contacts, lf_wrench, rf_wrench);
                obj = obj.updateLeggedOdometryWithContacts(encoders);
                obj = obj.updateBaseVelocity(encoders, encoder_speeds, contacts);                
            end
                        
            % store latest values
            obj.w_H_b_prev = obj.w_H_b;
            obj.lf_wrench_prev = lf_wrench;
            obj.rf_wrench_prev = rf_wrench;
            obj.encoders_prev = encoders;
            obj.contacts_prev = contacts;
            obj.t_prev = t;
            
            % outputs
            x = obj.x;            
            fixed_frame = obj.fixed_frame;
            
            % extract base link state
            [q, p, v, omega, ~, ~, ~, ~] = obj.extract(obj.x);
            R = Utils.quat2rot(q);
            omega_imu = R*omega;
            [BaseLinkState.q, BaseLinkState.I_p, ~, ~] = obj.model_comp.getBaseStateFromIMUState(q, p, v, omega_imu);            
            BaseLinkState.I_pdot = obj.v_b(1:3);
            BaseLinkState.I_omega = obj.v_b(4:6);
            
            obj = obj.publishToROS(encoders, encoder_speeds, lf_wrench, rf_wrench);
        end
        
        function [q, p, v, omega, qlf, plf, qrf, prf] = extract(obj, x)
            q = x(1:4);
            p = x(5:7);
            v = x(8:10);
            omega = x(11:13);
            qlf = x(14:17);
            plf = x(18:20);
            qrf = x(21:24);
            prf = x(25:27);            
        end 
    end
                                            
    methods (Access = private)
        function obj = getPrimaryFoot(obj, contacts, lf_wrench, rf_wrench)            
            switch obj.switching_pattern
                case 'alternate'
                    obj = obj.inferAlternatePattern(contacts);
                case 'cop'
                    obj = obj.inferCOPDistance(contacts, lf_wrench, rf_wrench);
            end
        end
                
        function obj = updateLeggedOdometryWithContacts(obj, encoders)            
            if (strcmp(obj.primary_foot, 'left'))
                obj.fixed_frame = obj.model_comp.LFVertexIds(1);
            elseif (strcmp(obj.primary_foot, 'right'))
                obj.fixed_frame = obj.model_comp.RFVertexIds(1);
            elseif (strcmp(obj.primary_foot, 'unknown'))
                obj.prev_fixed_frame = obj.fixed_frame;
                return;
            end
                  
            if (strcmp(obj.prev_fixed_frame, obj.fixed_frame) ~= 1)
                if (obj.flat_floor)
                    w_H_fixed = obj.leggedodom.getWorldFrameTransform(obj.leggedodom.model.getFrameIndex(obj.fixed_frame));
                    w_p_fixed = w_H_fixed.getPosition();
                    w_p_fixed.setVal(2, obj.soleheight_from_floor); % height of sole from floor
                    w_H_fixed.setPosition(w_p_fixed);
                    ret = obj.leggedodom.changeFixedFrame(obj.fixed_frame, w_H_fixed);
                else                    
                    ret = obj.leggedodom.changeFixedFrame(obj.fixed_frame);
                end
                
                if (~ret)                    
                    obj.prev_fixed_frame = obj.fixed_frame;
                    return;
                else
%                     disp(['Changed fixed frame to ', obj.fixed_frame]);
                end
            end
            
            jposDyn = iDynTree.JointPosDoubleArray();
            jposDyn.resize(length(encoders));
            jposDyn.zero();
            jposDyn.fromMatlab(encoders);
                
            obj.leggedodom.updateKinematics(jposDyn); 
            
            obj.w_H_b = obj.leggedodom.getWorldLinkTransform(obj.model_comp.base_link_idx).asHomogeneousTransform().toMatlab();                        
            obj.w_H_imu = obj.leggedodom.getWorldFrameTransform(obj.model_comp.base_link_imu_idx).asHomogeneousTransform().toMatlab();                        
            obj.w_H_lf = obj.leggedodom.getWorldFrameTransform(obj.model_comp.LFVertexIds(1)).asHomogeneousTransform().toMatlab();
            obj.w_H_rf = obj.leggedodom.getWorldFrameTransform(obj.model_comp.RFVertexIds(1)).asHomogeneousTransform().toMatlab();
                        
            obj.x = obj.constructStateVec();
            obj.prev_fixed_frame = obj.fixed_frame;
        end        
        
        % compute velocity
        function obj = updateBaseVelocity(obj, encoders, encoder_speeds, contacts)            
            % low pass filter joint velocities?
            Jc = obj.model_comp.getFrameFreeFloatingJacobian(obj.fixed_frame, obj.w_H_b_prev, encoders, encoder_speeds);
            
            Jc_base_inv = inv(Jc(1:6, 1:6));
            Jc_enc = Jc(1:6, 7:end);
            
            obj.v_b = -Jc_base_inv*Jc_enc*encoder_speeds;
          
            % output expressed as IMU velocity - transform from base to IMU
            % (right trivialized)
            imu_p_b = obj.model_comp.kindyn.getRelativeTransform(obj.model_comp.base_link, obj.model_comp.base_link_imu).inverse().getPosition().toMatlab();         
            imu_p_b_dyn = iDynTree.Position();
            imu_p_b_dyn.fromMatlab(imu_p_b);
            w_R_imu = iDynTree.Rotation();
            w_R_imu.fromMatlab(obj.w_H_imu(1:3, 1:3))
            X = iDynTree.Transform(iDynTree.Rotation.Identity(), w_R_imu*imu_p_b_dyn).asAdjointTransform().toMatlab();

            obj.v_imu = X*obj.v_b;
            obj.x = obj.constructStateVec();
        end
                        
        function x = constructStateVec(obj)
            q = Utils.rot2quat(obj.w_H_imu(1:3, 1:3));
            p = obj.w_H_imu(1:3, 4);
            v = obj.v_imu(1:3);
            omega = obj.v_imu(4:6);
            qlf = Utils.rot2quat(obj.w_H_lf(1:3, 1:3));
            plf = obj.w_H_lf(1:3, 4);
            qrf = Utils.rot2quat(obj.w_H_rf(1:3, 1:3));
            prf = obj.w_H_rf(1:3, 4);
            
            x = [q; p; v; omega; qlf; plf; qrf; prf];
        end                       
        
        function obj = inferAlternatePattern(obj, contacts)
            lf_transition = Estimation.ContactHandler.ContactTransition.contactTransitionMode(obj.contacts_prev(1), contacts(1));
            rf_transition = Estimation.ContactHandler.ContactTransition.contactTransitionMode(obj.contacts_prev(2), contacts(2));
            
            switch obj.primary_foot
                case 'left'
                    if (rf_transition == Estimation.ContactHandler.ContactTransition.CONTACT_MAKE &&...
                            contacts(1) == true)
                        obj.primary_foot = 'right';
                    elseif (lf_transition == Estimation.ContactHandler.ContactTransition.STABLE_ONCONTACT || ...
                            rf_transition == Estimation.ContactHandler.ContactTransition.STABLE_ONCONTACT)
                        obj.primary_foot = 'left';
                    elseif (lf_transition == Estimation.ContactHandler.ContactTransition.STABLE_OFFCONTACT)
                        obj.primary_foot = 'right';
                        if (rf_transition == Estimation.ContactHandler.ContactTransition.STABLE_OFFCONTACT)
                            obj.primary_foot = 'unknown';
                        end
                    end
                case 'right'
                    if (lf_transition == Estimation.ContactHandler.ContactTransition.CONTACT_MAKE &&...
                            contacts(2) == true)
                        obj.primary_foot = 'left';
                    elseif (rf_transition == Estimation.ContactHandler.ContactTransition.STABLE_ONCONTACT || ...
                            lf_transition == Estimation.ContactHandler.ContactTransition.STABLE_ONCONTACT)
                        obj.primary_foot = 'right';
                    elseif (rf_transition == Estimation.ContactHandler.ContactTransition.STABLE_OFFCONTACT)
                        obj.primary_foot = 'left';
                        if (lf_transition == Estimation.ContactHandler.ContactTransition.STABLE_OFFCONTACT)
                            obj.primary_foot = 'unknown';
                        end
                    end
                case 'unknown'
                    if (rf_transition == Estimation.ContactHandler.ContactTransition.CONTACT_MAKE ||...
                            rf_transition == Estimation.ContactHandler.ContactTransition.STABLE_ONCONTACT)
                        obj.primary_foot = 'right';
                    else if (lf_transition == Estimation.ContactHandler.ContactTransition.CONTACT_MAKE ||...
                                lf_transition == Estimation.ContactHandler.ContactTransition.STABLE_ONCONTACT)
                            obj.primary_foot = 'left';
                        end
                    end
            end
        end
        
        function obj = inferCOPDistance(obj, contacts, lf_wrench, rf_wrench)
            % check single support
            % send primary foot
            
            % check double support
            % apply product of truncated Gaussians along x and y
            % compute COP - probability of COP in the foot polygon
            % send highest probability foot as primary foot
        end
        
        function obj = publishToROS(obj, encoders, encoder_speeds, lf_wrench, rf_wrench)
            if obj.to_ros                                
                RosUtils.toRosTF(obj.tftree, obj.w_H_b, 'world', obj.model_comp.base_link_frame);
                RosUtils.toRosTF(obj.tftree, obj.w_H_imu, 'world', obj.model_comp.imu_frame);
                RosUtils.toRosTF(obj.tftree, obj.w_H_lf, 'world', obj.model_comp.l_foot_contact_frame);
                RosUtils.toRosTF(obj.tftree, obj.w_H_rf, 'world', obj.model_comp.r_foot_contact_frame);                
                
                % publish joint states
                effort = zeros(size(encoders));
                lf = zeros(6,1); lf(3) = lf_wrench(3);
                rf = zeros(6,1); rf(3) = rf_wrench(3);
                RosUtils.publishJointStates(obj.joint_publisher, obj.joint_list, encoders, encoder_speeds, effort);
                RosUtils.publishWrench(obj.lfw_publisher, lf, obj.model_comp.l_foot_contact_frame);
                RosUtils.publishWrench(obj.rfw_publisher, rf, obj.model_comp.r_foot_contact_frame);
            end                        
        end
        
    end
end

