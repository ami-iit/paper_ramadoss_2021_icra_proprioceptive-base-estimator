function [rotError, posError, velError, ATErot, ATEpos, ATEvel, RPErot, RPEpos] = getRightInvariantErrorMetrics(pos, rpy, vel, ...
    poshat, rpyhat, velhat, ...
    fixedIntervalIterationsForRelativePoseError,  align_yaw)
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
    R_i = Utils.rpy2rot(rpy(idx, 1), rpy(idx, 2), rpy(idx, 3));
    p_i = pos(idx, :)';
    v_i = vel(idx, :)';
    
    if ~align_yaw
        R_i_hat = Utils.rpy2rot(rpyhat(idx, 1), rpyhat(idx, 2), rpyhat(idx, 3));
    else
        R_i_hat = Utils.rpy2rot(rpyhat(idx, 1), rpyhat(idx, 2), rpy(idx, 3));
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
        R_iplusdelta = Utils.rpy2rot(rpy(iplusdelta, 1), rpy(iplusdelta, 2), rpy(iplusdelta, 3));
        p_iplusdelta = pos(iplusdelta, :)';
        
        if ~align_yaw
            R_iplusdelta_hat = Utils.rpy2rot(rpyhat(iplusdelta, 1), rpyhat(iplusdelta, 2), rpyhat(iplusdelta, 3));
        else
            R_iplusdelta_hat = Utils.rpy2rot(rpyhat(iplusdelta, 1), rpyhat(iplusdelta, 2), rpy(iplusdelta, 3));
        end
        p_iplusdelta_hat = poshat(iplusdelta, :)';
        
        H_i = LieGroups.SE3.constructSE3(R_i, p_i);
        H_i_hat = LieGroups.SE3.constructSE3(R_i_hat, p_i_hat);
        H_i_inv = LieGroups.SE3.inverse(H_i);
        H_i_hat_inv = LieGroups.SE3.inverse(H_i_hat);
        
        H_iplusdelta = LieGroups.SE3.constructSE3(R_iplusdelta, p_iplusdelta);
        H_iplusdelta_hat = LieGroups.SE3.constructSE3(R_iplusdelta_hat, p_iplusdelta_hat);
        
        % relative pose
        RP = LieGroups.SE3.compose(H_i_inv, H_iplusdelta);
        RPhat = LieGroups.SE3.compose(H_i_hat_inv, H_iplusdelta_hat);
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



