classdef ModelComputations < handle
    properties
        base_link (1, :) char;
        base_link_imu (1, :) char;
        l_foot_link (1, :) char;
        r_foot_foot_link (1, :) char;
        
        base_link_idx (1, 1) int64;
        base_link_imu_idx (1, 1) int64;
        l_foot_link_idx (1, 1) int64;
        r_foot_link_idx (1, 1) int64;
        
        RFVertexIds,
        LFVertexIds,
        
        kindyn (1, 1) iDynTree.KinDynComputations;
        nr_joints (1, 1) uint32;
    end
    
    methods
        function obj = ModelComputations(kinDyn, base_link_frame, base_link_imu_frame, l_foot_contact_frame, r_foot_contact_frame, LFVertexIDs, RFVertexIDs)
            obj.kindyn = kinDyn;
            obj.nr_joints = kinDyn.getNrOfDegreesOfFreedom();
            assert(obj.kindyn.model.getFrameIndex(base_link_frame) ~= iDynTree.FRAME_INVALID_INDEX, 'base link frame invalid');
            assert(obj.kindyn.model.getFrameIndex(base_link_imu_frame) ~= iDynTree.FRAME_INVALID_INDEX, 'base link IMU frame invalid');
            assert(obj.kindyn.model.getFrameIndex(l_foot_contact_frame) ~= iDynTree.FRAME_INVALID_INDEX, 'left foot contact frame invalid');
            assert(obj.kindyn.model.getFrameIndex(r_foot_contact_frame) ~= iDynTree.FRAME_INVALID_INDEX, 'right foot contact frame invalid');
            
            obj.base_link_imu = base_link_imu_frame;
            obj.base_link = base_link_frame;
            obj.l_foot_link = l_foot_contact_frame;
            obj.r_foot_foot_link = r_foot_contact_frame;
            
            obj.base_link_idx = obj.kindyn.model.getFrameIndex(base_link_frame);
            obj.base_link_imu_idx = obj.kindyn.model.getFrameIndex(base_link_imu_frame);
            obj.l_foot_link_idx = obj.kindyn.model.getFrameIndex(l_foot_contact_frame);
            obj.r_foot_link_idx = obj.kindyn.model.getFrameIndex(r_foot_contact_frame);
            for idx = 1:length(LFVertexIDs)
                assert(obj.kindyn.model.isValidFrameIndex(LFVertexIDs(idx)), 'left foot contact vertex invalid');
            end
            obj.LFVertexIds = LFVertexIDs;
            for idx = 1:length(RFVertexIDs)
                assert(obj.kindyn.model.isValidFrameIndex(RFVertexIDs(idx)), 'right foot contact vertex invalid');
            end
            obj.RFVertexIds = RFVertexIDs;
        end
        
        function setRobotState(obj, base_pose, baseTwistMixed, encoders, encoderSpeeds)
            jointPos = iDynTree.JointPosDoubleArray(obj.kindyn.model);
            jointPos.fromMatlab(encoders);
            
            R = base_pose(1:3, 1:3);
            t = base_pose(1:3, 4);
            
            baseR = iDynTree.Rotation();
            baseR.fromMatlab(R);
            
            baset = iDynTree.Position();
            baset.fromMatlab(t);
            
            grav = iDynTree.Vector3();
            grav.zero();
            grav.setVal(2, -9.8);
            
            jointVel = iDynTree.VectorDynSize();
            jointVel.resize(jointPos.size());
            jointVel.fromMatlab(encoderSpeeds);
            
            base_vel = iDynTree.Twist();
            base_vel.fromMatlab(baseTwistMixed);
            
            obj.kindyn.setRobotState(iDynTree.Transform(baseR, baset), jointPos, base_vel, jointVel, grav);
        end
        
        function [A_quat_b, A_p_b, A_v_B, A_omega_B] = getBaseStateFromIMUState(obj, q, p, v, imu_omega_Aimu)
            imuRot = iDynTree.Rotation();
            imuRot.fromMatlab(Utils.quat2rot(q));
            imuPos = iDynTree.Position();
            imuPos.fromMatlab(p);
            
            base_H_imu = obj.kindyn.getRelativeTransform(obj.base_link, obj.base_link_imu);
            
            A_H_imu = iDynTree.Transform(imuRot, imuPos);
            A_H_b = A_H_imu*(base_H_imu.inverse());
            A_quat_b = A_H_b.getRotation().asQuaternion().toMatlab();
            A_p_b = A_H_b.getPosition().toMatlab();
            
            v_imu = iDynTree.Twist();
            angvelimu = imuRot.toMatlab()*imu_omega_Aimu;
            v_imu.fromMatlab([v; angvelimu]);
            
            X = iDynTree.Transform(iDynTree.Rotation.Identity(), (A_H_b.getRotation()*base_H_imu.getPosition()));
            v_base = X.asAdjointTransform().toMatlab()*v_imu.toMatlab();
            
            A_v_B = v_base(1:3);
            A_omega_B = v_base(4:6);
        end
        
        function [quat, p, J_imuf] = relativeKinIMU_to_foot_contact(obj, encoders, foot)
            jointPos = iDynTree.JointPosDoubleArray(length(encoders));
            jointPos.fromMatlab(encoders);
            obj.kindyn.setJointPos(jointPos);
            
            if (strcmp(foot, 'left') == true)
                foot_frame_idx = obj.LFVertexIds(1);
            elseif (strcmp(foot, 'right') == true)
                foot_frame_idx = obj.RFVertexIds(1);
            end
            
            imu_H_foot = obj.kindyn.getRelativeTransform(obj.base_link_imu_idx, foot_frame_idx);
            p = imu_H_foot.getPosition().toMatlab();
            quat = imu_H_foot.getRotation().asQuaternion().toMatlab();
            
            J_imuf_dyn = iDynTree.MatrixDynSize();
            obj.kindyn.getRelativeJacobian(obj.base_link_imu_idx, foot_frame_idx, J_imuf_dyn);
            J_imuf = J_imuf_dyn.toMatlab();
        end
        
        function [quat, p, J_imuf] = relativeKinIMU_to_foot_contactLeftTriv(obj, encoders, foot)
            jointPos = iDynTree.JointPosDoubleArray(length(encoders));
            jointPos.fromMatlab(encoders);
            obj.kindyn.setJointPos(jointPos);
            
            if (strcmp(foot, 'left') == true)
                foot_frame_idx = obj.LFVertexIds(1);
            elseif (strcmp(foot, 'right') == true)
                foot_frame_idx = obj.RFVertexIds(1);
            end
            
            imu_H_foot = obj.kindyn.getRelativeTransform(obj.base_link_imu_idx, foot_frame_idx);
            p = imu_H_foot.getPosition().toMatlab();
            quat = imu_H_foot.getRotation().asQuaternion().toMatlab();
            
            J_imuf_dyn = iDynTree.MatrixDynSize();
            obj.kindyn.getRelativeJacobianExplicit(obj.base_link_imu_idx, foot_frame_idx, foot_frame_idx, foot_frame_idx, J_imuf_dyn);
            J_imuf = J_imuf_dyn.toMatlab();
        end
        
        function [imu_H_foot, J_imuf] = relativeKinIMU_to_foot_contactExplicit(obj, encoders, foot)
            jointPos = iDynTree.JointPosDoubleArray(length(encoders));
            jointPos.fromMatlab(encoders);
            obj.kindyn.setJointPos(jointPos);
            
            if (strcmp(foot, 'left') == true)
                foot_frame_idx = obj.LFVertexIds(1);
            elseif (strcmp(foot, 'right') == true)
                foot_frame_idx = obj.RFVertexIds(1);
            end
            
            imu_H_foot = obj.kindyn.getRelativeTransform(obj.base_link_imu_idx, foot_frame_idx).asHomogeneousTransform().toMatlab();
            
            J_imuf_dyn = iDynTree.MatrixDynSize();
            obj.kindyn.getRelativeJacobianExplicit(obj.base_link_imu_idx, foot_frame_idx, foot_frame_idx, foot_frame_idx, J_imuf_dyn);
            J_imuf = J_imuf_dyn.toMatlab();
        end
        
        function Jc = getFrameFreeFloatingJacobian(obj, frame, base_pose, encoders, encoder_speeds)
            jointPos = iDynTree.JointPosDoubleArray(length(encoders));
            jointPos.fromMatlab(encoders);
            
            R = base_pose(1:3, 1:3);
            t = base_pose(1:3, 4);
            
            baseR = iDynTree.Rotation();
            baseR.fromMatlab(R);
            
            baset = iDynTree.Position();
            baset.fromMatlab(t);
            
            grav = iDynTree.Vector3();
            grav.zero();
            grav.setVal(2, -9.8);
            
            jointVel = iDynTree.VectorDynSize();
            jointVel.resize(jointPos.size());
            jointVel.zero();
%             jointVel.fromMatlab(encoder_speeds);
            
            base_vel = iDynTree.Twist();
            base_vel.zero();
            
            obj.kindyn.setRobotState(iDynTree.Transform(baseR, baset), jointPos, base_vel, jointVel, grav);
%             obj.kindyn.setJointPos(jointPos);
            
            Jcdyn = iDynTree.MatrixDynSize();
            Jcdyn.resize(6, 6+length(encoders))
            obj.kindyn.getFrameFreeFloatingJacobian(frame, Jcdyn);
            Jc = Jcdyn.toMatlab();
        end
        
    end
    
end


