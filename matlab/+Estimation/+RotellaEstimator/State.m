classdef State
    methods (Static)
        function [quat, p, v, quatl, pl, quatr, pr, ba, bg] = extract(x, ba0, bg0)
            assert((length(x) == 30 || length(x) == 24), 'State size mismatch');
            
            if (length(x) == 30)
                ba = x(25:27);
                bg = x(28:30);
            else
                if (~isempty(ba0))
                    ba = ba0;
                else
                    ba = [0;0;0];
                end
                if (~isempty(bg0))
                    bg = bg0;
                else
                    bg = [0;0;0];
                end
            end
            
            quat = x(1:4);
            p = x(5:7);
            v = x(8:10);
            quatl = x(11:14);
            pl = x(15:17);
            quatr = x(18:21);
            pr = x(22:24);
        end
        
        function [dq, dp, dv, dql, dpl, dqr, dpr, dba, dbg] = extractLinearizedState(dx)
            assert((length(dx) == 27 || length(dx) == 21), 'State size mismatch');
            
            if (length(dx) == 327)
                dba = dx(22:24);
                dbg = dx(25:27);
            else
                dba = [0;0;0];
                dbg = [0;0;0];
            end
            
            dq = dx(1:3);
            dp = dx(4:6);
            dv = dx(7:9);
            dql = dx(10:12);
            dpl = dx(13:15);
            dqr = dx(16:18);
            dpr = dx(19:21);
        end
        
        function x = constructwithbias(quat, p, v, quatl, pl, quatr, pr, ba, bg)
            x = zeros(30, 1);
            x = [quat;p;v;quatl;pl;quatr;pr; ba; bg];
        end
        
        function x = constructwithoutbias(quat, p, v, quatl, pl, quatr, pr)
            x = zeros(24, 1);
            x = [quat;p;v;quatl;pl;quatr;pr];
        end
        
        function P = constructVar(prior_dev, estimate_bias)
            P = blkdiag(diag(prior_dev.imu_orientation.^2), ...
                diag(prior_dev.imu_position.^2), ...
                diag(prior_dev.imu_linear_velocity.^2), ...
                diag(prior_dev.left_foot_orientation.^2), ...
                diag(prior_dev.left_foot_position.^2), ...
                diag(prior_dev.right_foot_orientation.^2), ...
                diag(prior_dev.right_foot_position.^2));
            
            if (estimate_bias)
                P = blkdiag(P, ...
                    diag(prior_dev.accel_bias.^2), ...
                    diag(prior_dev.gyro_bias.^2));
            end
        end
        
        function Psub = extractSubBlockVar(P, state)
            type = 'linearized';
            estimate_bias = true;
            range = Estimation.RotellaEstimator.State.getRange(state, type, estimate_bias);
            Psub = P(range, range);
        end
        
        function range = getRange(state, type, estimate_bias)
            if (~estimate_bias)
                if ( (strcmp(state, 'acc_bias')  == 1) || (strcmp(state, 'gyro_bias')  == 1))
                    range = [-3: -1];
                    error('Accessing bias states while estiamte_bias option disabled');
                end
            end
            
            if (strcmp(type, 'nonlinear'))
                switch state
                    case 'base_rot'
                        range = [1:4];
                    case 'base_pos'
                        range = [5:7];
                    case 'base_vel'
                        range = [8:10];
                    case 'lf_rot'
                        range = [11:14];
                    case 'lf_pos'
                        range = [15:17];
                    case 'rf_rot'
                        range = [18:21];
                    case 'rf_pos'
                        range = [22:24];
                end
                
                if (estimate_bias)
                    switch state
                        case 'bias_acc'
                            range = [25:27];
                        case 'bias_gyro'
                            range = [28:30];
                    end
                end
                
            elseif (strcmp(type, 'linearized'))
                switch state
                    case 'base_rot'
                        range = [1:3];
                    case 'base_pos'
                        range = [4:6];
                    case 'base_vel'
                        range = [7:9];
                    case 'lf_rot'
                        range = [10:12];
                    case 'lf_pos'
                        range = [13:15];
                    case 'rf_rot'
                        range = [16:18];
                    case 'rf_pos'
                        range = [19:21];
                end
                
                if (estimate_bias)
                    switch state
                        case 'bias_acc'
                            range = [22:24];
                        case 'bias_gyro'
                            range = [25:27];
                    end
                end
            end
        end % end function getRange
        
        
        function [Pbrot, Ppos, Pvel, Plfpos, Plfrot, Prfpos, Prfrot, Pba, Pbg] = extractStateVarSubBlockEvolutions(Ptraj, estimate_bias)
            Pbrot = zeros(length(Ptraj), 3, 3);
            Ppos = zeros(length(Ptraj), 3, 3);
            Pvel = zeros(length(Ptraj), 3, 3);
            Plfpos = zeros(length(Ptraj), 3, 3);
            Plfrot = zeros(length(Ptraj), 3, 3);
            Prfpos = zeros(length(Ptraj), 3, 3);
            Prfrot = zeros(length(Ptraj), 3, 3);
            Pba = zeros(length(Ptraj), 3, 3);
            Pbg = zeros(length(Ptraj), 3, 3);
            
            for iter_idx = 1: length(Ptraj)
                Pbrot(iter_idx, :, :) = Estimation.RotellaEstimator.State.extractSubBlockVar(squeeze(Ptraj{iter_idx}), 'base_rot');
                Ppos(iter_idx, :, :) = Estimation.RotellaEstimator.State.extractSubBlockVar(squeeze(Ptraj{iter_idx}), 'base_pos');
                Pvel(iter_idx, :, :) = Estimation.RotellaEstimator.State.extractSubBlockVar(squeeze(Ptraj{iter_idx}), 'base_vel');
                Plfpos(iter_idx, :, :) = Estimation.RotellaEstimator.State.extractSubBlockVar(squeeze(Ptraj{iter_idx}), 'lf_pos');
                Plfrot(iter_idx, :, :) = Estimation.RotellaEstimator.State.extractSubBlockVar(squeeze(Ptraj{iter_idx}), 'lf_rot');
                Prfpos(iter_idx, :, :) = Estimation.RotellaEstimator.State.extractSubBlockVar(squeeze(Ptraj{iter_idx}), 'rf_pos');
                Prfrot(iter_idx, :, :) = Estimation.RotellaEstimator.State.extractSubBlockVar(squeeze(Ptraj{iter_idx}), 'rf_rot');
                if (estimate_bias)
                    Pba(iter_idx, :, :) = Estimation.RotellaEstimator.State.extractSubBlockVar(squeeze(Ptraj{iter_idx}), 'bias_acc');
                    Pbg(iter_idx, :, :) = Estimation.RotellaEstimator.State.extractSubBlockVar(squeeze(Ptraj{iter_idx}), 'bias_gyro');
                end
            end
        end
        
    end
end

