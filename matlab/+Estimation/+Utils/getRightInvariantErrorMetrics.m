function [Base, LF, RF] = getRightInvariantErrorMetrics(matfile, fixedIntervalIterationsForRelativePoseError, noANEES, alignyaw)
%ATErot RMSE absolute trajectory rotation error
%ATEpos RMSE absolute trajectory position error
%ATEvel RMSE absolute trajectory velocity error
%ANEESrot Average normalized estimation error squared for rotation
%ANEESpos Average normalized estimation error squared for position
%ANEESvel Average normalized estimation error squared for velocity
%RPEpos Relative pose error translation
% the error state vectors filling up the covariance matrix are expected to
% be right trivialized
exp = load(matfile);
b_R_imu = exp.b_H_imu(1:3, 1:3);
%% base
if ~noANEES
[Base.rotError, Base.posError, Base.velError, Base.ATErot, Base.ATEpos, Base.ATEvel, ...
 Base.NEESrot, Base.NEESpos, Base.NEESvel, Base.ANEESrot, Base.ANEESpos, Base.ANEESvel, Base.RPErot, Base.RPEpos] = ...
errors(exp.data.simBasePos, exp.data.simBaseRot, exp.data.simBaseLinVel, ...
       exp.EstimatorBaseState.I_p', exp.EstimatorBaseState.rpy', exp.EstimatorBaseState.I_pdot', ...
       exp.EstimatorOut.Ppos, exp.EstimatorOut.Pbrot, exp.EstimatorOut.Pvel, b_R_imu,...
       fixedIntervalIterationsForRelativePoseError, alignyaw);
else
    [Base.rotError, Base.posError, Base.velError, Base.ATErot, Base.ATEpos, Base.ATEvel, ...
  Base.RPErot, Base.RPEpos] = ...
errorsATERPE(exp.data.simBasePos, exp.data.simBaseRot, exp.data.simBaseLinVel, ...
       exp.EstimatorBaseState.I_p', exp.EstimatorBaseState.rpy', exp.EstimatorBaseState.I_pdot', ...
       b_R_imu,...
       fixedIntervalIterationsForRelativePoseError, alignyaw);
end
%% lf
if ~noANEES
[LF.rotError, LF.posError, ~, LF.ATErot, LF.ATEpos, ~, LF.NEESrot, LF.NEESpos, LF.NEESvel, LF.ANEESrot, LF.ANEESpos, ~, LF.RPErot, LF.RPEpos] = ...
errors(exp.simlsolePosTraj, exp.simlsoleRotTraj, zeros(size(exp.simlsolePosTraj)), ...
       exp.EstimatorOut.I_p_LF', exp.EstimatorOut.lfrollpitchyaw', zeros(size(exp.simlsolePosTraj)), ...
       exp.EstimatorOut.Plfpos, exp.EstimatorOut.Plfrot, zeros(size(exp.EstimatorOut.Plfpos)), eye(3), ...
       fixedIntervalIterationsForRelativePoseError, alignyaw);
else
    [LF.rotError, LF.posError, ~, LF.ATErot, LF.ATEpos, ~,  LF.RPErot, LF.RPEpos] = ...
errorsATERPE(exp.simlsolePosTraj, exp.simlsoleRotTraj, zeros(size(exp.simlsolePosTraj)), ...
       exp.EstimatorOut.I_p_LF', exp.EstimatorOut.lfrollpitchyaw', zeros(size(exp.simlsolePosTraj)), ...
       eye(3), ...
       fixedIntervalIterationsForRelativePoseError, alignyaw);
end
%% rf
if ~noANEES
[RF.rotError, RF.posError, ~, RF.ATErot, RF.ATEpos, ~, RF.NEESrot, RF.NEESpos, RF.NEESvel, RF.ANEESrot, RF.ANEESpos, ~, RF.RPErot, RF.RPEpos] = ...
errors(exp.simrsolePosTraj, exp.simrsoleRotTraj, zeros(size(exp.simrsolePosTraj)), ...
       exp.EstimatorOut.I_p_RF', exp.EstimatorOut.rfrollpitchyaw', zeros(size(exp.simrsolePosTraj)), ...
       exp.EstimatorOut.Prfpos, exp.EstimatorOut.Prfrot, zeros(size(exp.EstimatorOut.Prfpos)), eye(3), ...
       fixedIntervalIterationsForRelativePoseError, alignyaw);
else
    [RF.rotError, RF.posError, ~, RF.ATErot, RF.ATEpos, ~, RF.RPErot, RF.RPEpos] = ...
errorsATERPE(exp.simrsolePosTraj, exp.simrsoleRotTraj, zeros(size(exp.simrsolePosTraj)), ...
       exp.EstimatorOut.I_p_RF', exp.EstimatorOut.rfrollpitchyaw', zeros(size(exp.simrsolePosTraj)), ...
       eye(3), ...
       fixedIntervalIterationsForRelativePoseError, alignyaw);
end


    function [rotError, posError, velError, ATErot, ATEpos, ATEvel, NEESrot, NEESpos, NEESvel, ANEESrot, ANEESpos, ANEESvel, RPErot, RPEpos] = errors(pos, rpy, vel, ... 
                                                                                             poshat, rpyhat, velhat, ...
                                                                                             Ppos, Prot, Pvel, R,....
                                                                                             fixedIntervalIterationsForRelativePoseError, alignyaw)
        n = length(pos);
        
        rotError = zeros(n, 3);
        posError = zeros(n, 3);
        velError = zeros(n, 3);
        ATErot = 0;
        ATEpos = 0;
        ATEvel = 0;
        NEESrot = zeros(n, 1);
        NEESpos = zeros(n, 1);
        NEESvel = zeros(n, 1);
        ANEESrot = 0;
        ANEESpos = 0;
        ANEESvel = 0;
        RPEpos = 0;
        RPErot = 0;
        for idx = 1:n
            iplusdelta = fixedIntervalIterationsForRelativePoseError;
            R_i = Core.rpy2rot(rpy(idx, 1), rpy(idx, 2), rpy(idx, 3));
            p_i = pos(idx, :)';
            v_i = vel(idx, :)';
            
            if ~align_yaw
                R_i_hat = Core.rpy2rot(rpyhat(idx, 1), rpyhat(idx, 2), rpyhat(idx, 3));
            else
                R_i_hat = Core.rpy2rot(rpyhat(idx, 1), rpyhat(idx, 2), rpy(idx, 3));
            end
            p_i_hat = poshat(idx, :)';
            v_i_hat = velhat(idx, :)';
            
            % right invariant error ate
            deltaR = R_i*R_i_hat';
            deltaPhi = LieGroups.SO3.logvee(deltaR);            
            deltap = p_i - deltaR*p_i_hat;
            deltav = v_i - deltaR*v_i_hat;
            
            rotError(idx, :) = deltaPhi';
            posError(idx, :) = deltap';
            velError(idx, :) = deltav';
            
            ATErot = ATErot + (deltaPhi(1)^2 + deltaPhi(2)^2 + deltaPhi(3)^2);
            ATEpos = ATEpos + (deltap(1)^2 + deltap(2)^2 + deltap(3)^2);
            ATEvel = ATEvel + (deltav(1)^2 + deltav(2)^2 + deltav(3)^2);
                                   
            % for ANEES - rotate covariance from IMU frame to Base frame
            Ppos_i = R*squeeze(Ppos(idx, :, :))*R';
            Prot_i = R*squeeze(Prot(idx, :, :))*R';
            Pvel_i = R*squeeze(Pvel(idx, :, :))*R';
            
            NEESrot(idx) = (deltaPhi'*inv(Prot_i)*deltaPhi);
            NEESpos(idx) = (deltap'*inv(Ppos_i)*deltap);
            NEESvel(idx) = (deltav'*inv(Pvel_i)*deltav);
            ANEESrot = ANEESrot + NEESrot(idx);
            ANEESpos = ANEESpos + NEESpos(idx);
            ANEESvel = ANEESvel + NEESvel(idx);
            
            % convert after saving
            NEESrot(idx) = rad2deg(NEESrot(idx));            
            
            % for Relative Pose Error
            if (iplusdelta) < n                                
                R_iplusdelta = Core.rpy2rot(rpy(iplusdelta, 1), rpy(iplusdelta, 2), rpy(iplusdelta, 3));
                p_iplusdelta = pos(iplusdelta, :)';                             
                
                if ~align_yaw
                    R_iplusdelta_hat = Core.rpy2rot(rpyhat(iplusdelta, 1), rpyhat(iplusdelta, 2), rpyhat(iplusdelta, 3));
                else
                    R_iplusdelta_hat = Core.rpy2rot(rpyhat(iplusdelta, 1), rpyhat(iplusdelta, 2), rpy(iplusdelta, 3));
                end
                p_iplusdelta_hat = poshat(iplusdelta, :)';               
                
                H_i = LieGroups.SE3.constructSE3(R_i, p_i);                
                H_i_hat = LieGroups.SE3.constructSE3(R_i_hat, p_i_hat);                
                H_i_inv = LieGroups.SE3.inverse(H_i);
                H_i_hat_inv = LieGroups.SE3.inverse(H_i_hat);
                
                H_iplusdelta = LieGroups.SE3.constructSE3(R_iplusdelta, p_iplusdelta);
                H_iplusdelta_inv = LieGroups.SE3.inverse(H_iplusdelta);
                H_iplusdelta_hat = LieGroups.SE3.constructSE3(R_iplusdelta_hat, p_iplusdelta_hat);
                H_iplusdelta_hat_inv = LieGroups.SE3.inverse(H_iplusdelta_hat);
                
                % relative pose
                RP = LieGroups.SE3.compose(H_i, H_iplusdelta_inv);
                RPhat = LieGroups.SE3.compose(H_i_hat, H_iplusdelta_hat_inv);
                RPhatinv = LieGroups.SE3.inverse(RPhat);
                % right invariant relative pose error
                E_i = LieGroups.SE3.compose(RP, RPhatinv);
                [Ei_R, Ei_p] = LieGroups.SE3.extractSE3(E_i);
                Ei_phi = LieGroups.SO3.logvee(Ei_R);
                RPErot = RPErot + (Ei_phi(1)^2 + Ei_phi(2)^2 + Ei_phi(3)^2);
                RPEpos = RPEpos + (Ei_p(1)^2 + Ei_p(2)^2 + Ei_p(3)^2);                
            end
        end
        
        ATErot = rad2deg(sqrt(ATErot/n));
        ATEpos = sqrt(ATEpos/n);
        ATEvel = sqrt(ATEvel/n);
        ANEESrot = rad2deg(ANEESrot/n);
        ANEESpos = ANEESpos/n;
        ANEESvel = ANEESvel/n;
        RPErot = rad2deg(sqrt(RPErot/n));
        RPEpos = sqrt(RPEpos/n);        
    end


function [rotError, posError, velError, ATErot, ATEpos, ATEvel, RPErot, RPEpos] = errorsATERPE(pos, rpy, vel, ... 
                                                                                             poshat, rpyhat, velhat, ...
                                                                                             R,....
                                                                                             fixedIntervalIterationsForRelativePoseError, align_yaw)
        n = length(pos);
        
        rotError = zeros(n, 3);
        posError = zeros(n, 3);
        velError = zeros(n, 3);
        ATErot = 0;
        ATEpos = 0;
        ATEvel = 0;
        RPEpos = 0;
        RPErot = 0;
        for idx = 1:n
            iplusdelta = fixedIntervalIterationsForRelativePoseError;
            R_i = Core.rpy2rot(rpy(idx, 1), rpy(idx, 2), rpy(idx, 3));
            p_i = pos(idx, :)';
            v_i = vel(idx, :)';
            
            if ~align_yaw
                R_i_hat = Core.rpy2rot(rpyhat(idx, 1), rpyhat(idx, 2), rpyhat(idx, 3));
            else
                R_i_hat = Core.rpy2rot(rpyhat(idx, 1), rpyhat(idx, 2), rpy(idx, 3));
            end
            p_i_hat = poshat(idx, :)';
            v_i_hat = velhat(idx, :)';
            
            % right invariant error ate
            deltaR = R_i*R_i_hat';
            deltaPhi = LieGroups.SO3.logvee(deltaR);            
            deltap = p_i - deltaR*p_i_hat;
            deltav = v_i - deltaR*v_i_hat;
            
            rotError(idx, :) = deltaPhi';
            posError(idx, :) = deltap';
            velError(idx, :) = deltav';
            
            ATErot = ATErot + (deltaPhi(1)^2 + deltaPhi(2)^2 + deltaPhi(3)^2);
            ATEpos = ATEpos + (deltap(1)^2 + deltap(2)^2 + deltap(3)^2);
            ATEvel = ATEvel + (deltav(1)^2 + deltav(2)^2 + deltav(3)^2);                                   
            
            % for Relative Pose Error
            if (iplusdelta) < n                                
                R_iplusdelta = Core.rpy2rot(rpy(iplusdelta, 1), rpy(iplusdelta, 2), rpy(iplusdelta, 3));
                p_iplusdelta = pos(iplusdelta, :)';                             
            
                if ~align_yaw
                    R_iplusdelta_hat = Core.rpy2rot(rpyhat(iplusdelta, 1), rpyhat(iplusdelta, 2), rpyhat(iplusdelta, 3));
                else
                    R_iplusdelta_hat = Core.rpy2rot(rpyhat(iplusdelta, 1), rpyhat(iplusdelta, 2), rpy(iplusdelta, 3));
                end
                p_iplusdelta_hat = poshat(iplusdelta, :)';               
                
                H_i = LieGroups.SE3.constructSE3(R_i, p_i);                
                H_i_hat = LieGroups.SE3.constructSE3(R_i_hat, p_i_hat);                
                H_i_inv = LieGroups.SE3.inverse(H_i);
                H_i_hat_inv = LieGroups.SE3.inverse(H_i_hat);
                
                H_iplusdelta = LieGroups.SE3.constructSE3(R_iplusdelta, p_iplusdelta);
                H_iplusdelta_inv = LieGroups.SE3.inverse(H_iplusdelta);
                H_iplusdelta_hat = LieGroups.SE3.constructSE3(R_iplusdelta_hat, p_iplusdelta_hat);
                H_iplusdelta_hat_inv = LieGroups.SE3.inverse(H_iplusdelta_hat);
                
                % relative pose
                RP = LieGroups.SE3.compose(H_i, H_iplusdelta_inv);
                RPhat = LieGroups.SE3.compose(H_i_hat, H_iplusdelta_hat_inv);
                RPhatinv = LieGroups.SE3.inverse(RPhat);
                % right invariant relative pose error
                E_i = LieGroups.SE3.compose(RP, RPhatinv);
                [Ei_R, Ei_p] = LieGroups.SE3.extractSE3(E_i);
                Ei_phi = LieGroups.SO3.logvee(Ei_R);
                RPErot = RPErot + (Ei_phi(1)^2 + Ei_phi(2)^2 + Ei_phi(3)^2);
                RPEpos = RPEpos + (Ei_p(1)^2 + Ei_p(2)^2 + Ei_p(3)^2);                
            end
        end
        
        ATErot = rad2deg(sqrt(ATErot/n));
        ATEpos = sqrt(ATEpos/n);
        ATEvel = sqrt(ATEvel/n);
        RPErot = rad2deg(sqrt(RPErot/n));
        RPEpos = sqrt(RPEpos/n);        
    end



end

