classdef State
    methods (Static)
        % This state uses a different serialization than the one
        % in LieGroups.SE_Nplus2_3 in both states and velocity vector
        % This change in serialization mainly affects the Adjoint
        % matrix order (provides a lower triangular matrix)
        % We will use this state in the filter to be coherent with the original
        % implementation in
        % https://github.com/UMich-BipedLab/Contact-Aided-Invariant-EKF
        % Note that this state does not contain additional landmarks
        function [R, v, p, prf, plf, bg, ba, lm] = extract(X, theta)
            R = X(1:3, 1:3);
            v = X(1:3, 4);
            p = X(1:3, 5);
            prf = X(1:3, 6);
            plf = X(1:3, 7);
            
            if size(X, 2) > 7
                lm = X(:, 8:end);
            else
                lm = [];
            end
            
            bg = theta(1:3);
            ba = theta(4:6);
        end
        
        function  [X, theta] = construct(R, v, p, prf, plf, bg, ba, lm)
            X =  eye(7 + size(lm, 2));
            X(1:7, 1:7) = [         R    v   p   prf   plf; ...
                zeros(1, 3)   1   0     0     0; ...
                zeros(1, 3)   0   1     0     0; ...
                zeros(1, 3)   0   0     1     0; ...
                zeros(1, 3)   0   0     0     1];
            
            if ~isempty(lm)
                X(1:3, 8:end) = lm;
            end
            
            theta = [bg; ba];
        end
        
        function AdjX = Adjoint(X)
            dummy_theta = zeros(6, 1);
            [R, v, p, prf, plf, ~, ~, lm] = Estimation.InvEKF.State.extract(X, dummy_theta);
            n = size(lm, 2);
            AdjX = zeros(15 +3*n);
            v_cross = Utils.skew(v);
            p_cross = Utils.skew(p);
            prf_cross = Utils.skew(prf);
            plf_cross = Utils.skew(plf);
            
            AdjX(1:15, 1:15) =   [         R  zeros(3) zeros(3) zeros(3) zeros(3); ...
                v_cross*R        R  zeros(3) zeros(3) zeros(3); ...
                p_cross*R  zeros(3)       R  zeros(3) zeros(3); ...
                prf_cross*R  zeros(3) zeros(3)       R  zeros(3); ...
                plf_cross*R  zeros(3) zeros(3) zeros(3)       R ];
            
            if ~isempty(lm)                
                for idx = 1:n
                    d = lm(:, idx);
                    d_cross = Utils.skew(d);
                    
                    begin = 16+3*(idx-1);
                    fin = begin+2;
                    AdT(begin:fin, begin:fin) =   R;
                    AdT(begin:fin, 1:3) =   d_cross*R;
                end
            end
        end
        
        
        function [X, theta] = exphat(v_X, v_theta)
            % use closed form solution instead of expm()
            [omega, alin, vlin, vrf, vlf, vlm] = Estimation.InvEKF.State.splitVector(v_X);            
            n = length(vlm/3);
            vlm3byn = reshape(vlm, 3, n);

            JlSO3 = LieGroups.SO3.leftJacobian(omega);
            
            R = LieGroups.SO3.exphat(omega);
            v = JlSO3*alin;
            p = JlSO3*vlin;
            prf = JlSO3*vrf;
            plf = JlSO3*vlf;
            lm = JlSO3*vlm3byn;
            
            bg = v_theta(1:3);
            ba = v_theta(4:6);
                       
            [X, theta] = Estimation.InvEKF.State.construct(R, v, p, prf, plf, bg, ba, lm);
        end
        
        function [omega, alin, vlin, vrf, vlf, vlm] =  splitVector(v)
            omega = v(1:3);
            alin = v(4:6);
            vlin = v(7:9);
            vrf = v(10:12);
            vlf = v(13:15);
            
            vlm = [];
            if length(v) > 15
                assert( mod( length(v(16:end)), 3) == 0, 'vector size mismatch');
                vlm = v(16:end);
            end
        end
        
        
        function [Pbrot, Pvel, Ppos, Prfpos, Plfpos, Plm, Pbg, Pba] = extractStateVarSubBlockEvolutions(Ptraj, estimate_bias)
            Pbrot = zeros(length(Ptraj), 3, 3);
            Ppos = zeros(length(Ptraj), 3, 3);
            Pvel = zeros(length(Ptraj), 3, 3);
            Plfpos = zeros(length(Ptraj), 3, 3);
            Prfpos = zeros(length(Ptraj), 3, 3);
            
            % {TODO} handle landmarks and bias variances
            Plm = [];
            Pbg = zeros(length(Ptraj), 3, 3);
            Pba = zeros(length(Ptraj), 3, 3);
            
            for iter_idx = 1: length(Ptraj)
                P = squeeze(Ptraj{iter_idx});
                Pbrot(iter_idx, :, :) = P(1:3, 1:3);
                Pvel(iter_idx, :, :) = P(4:6, 4:6);
                Ppos(iter_idx, :, :) = P(7:9, 7:9);
                Prfpos(iter_idx, :, :) = P(10:12, 10:12);
                Plfpos(iter_idx, :, :) = P(13:15, 13:15);
                
                lm_offset  = size(P, 2) - 15 - 6;
                Pbg(iter_idx, :, :) = P(16+lm_offset:18+lm_offset, 16+lm_offset:18+lm_offset);
                Pba(iter_idx, :, :) = P(19+lm_offset:21+lm_offset, 19+lm_offset:21+lm_offset);
            end
        end
              
    end
end

