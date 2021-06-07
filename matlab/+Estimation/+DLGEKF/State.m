classdef State
    methods (Static)
        function Tout = compose(X1, X2)
            assert(size(X1, 2) == size(X2, 2), 'size mismatch');
            
            [Rb1, pb1, vb1, Rlf1, plf1, Rrf1, prf1, ba1, bg1] = Estimation.DLGEKF.State.extract(X1);
            [Rb2, pb2, vb2, Rlf2, plf2, Rrf2, prf2, ba2, bg2] = Estimation.DLGEKF.State.extract(X2);
            
            Rbout = Rb1*Rb2;
            pbout = Rb1*pb2 + pb1;
            vbout = Rb1*vb2 + vb1;
            
            Rlfout = Rlf1*Rlf2;
            plfout = Rlf1*plf2 + plf1;
            
            Rrfout = Rrf1*Rrf2;
            prfout = Rrf1*prf2 + prf1;
            
            baout = ba1 + ba2;
            bgout = bg1 + bg2;
            
            estimate_bias = Estimation.DLGEKF.State.is_estimate_bias_enabled(X1);
            Tout = Estimation.DLGEKF.State.construct(Rbout, pbout, vbout, ...
                Rlfout, plfout, ...
                Rrfout, prfout, ...
                baout, bgout, estimate_bias);
        end
        
        function Tout = inverse(X)
            [Rb, pb, vb, Rlf, plf, Rrf, prf, ba, bg] = Estimation.DLGEKF.State.extract(X);
            
            Rbout = Rb';
            pbout = -Rb'*pb;
            vbout = -Rb'*vb;
            
            Rlfout = Rlf';
            plfout = -Rlf'*plf;
            
            Rrfout = Rrf';
            prfout = -Rrf'*prf;
            
            baout = -ba;
            bgout = -bg;
            
            estimate_bias = Estimation.DLGEKF.State.is_estimate_bias_enabled(X);
            Tout = Estimation.DLGEKF.State.construct(Rbout, pbout, vbout, ...
                Rlfout, plfout, ...
                Rrfout, prfout, ...
                baout, bgout, estimate_bias);
        end
        
        function Tout = identity(estimate_bias)
            if ~estimate_bias
                Tout = eye(13);
            else
                Tout = eye(20);
            end
        end
        
        function v = vee(g)
            [gbase, glf, grf, gbiases] = Estimation.DLGEKF.State.extractLieAlgebra(g);
            vbase = LieGroups.SE_2_3.vee(gbase);
            vlf = LieGroups.SE3.vee(glf);
            vrf = LieGroups.SE3.vee(grf);
            
            v = [vbase; vlf; vrf];
            if Estimation.DLGEKF.State.is_estimate_bias_enabled(g)
                vbias = LieGroups.Tn.vee(gbiases);
                v = [v; vbias];
            end
        end
        
        function g = hat(v)
            estimate_bias = Estimation.DLGEKF.State.is_estimate_bias_enabled(v);
            [vbase, vlf, vrf, vbias] = Estimation.DLGEKF.State.splitVector(v, estimate_bias);
            gbase = LieGroups.SE_2_3.hat(vbase);
            glf = LieGroups.SE3.hat(vlf);
            grf = LieGroups.SE3.hat(vrf);
                        
            gbiases = [];
            if estimate_bias
                gbiases = LieGroups.Tn.hat(vbias);
            end
            g = Estimation.DLGEKF.State.constructLieAlgebra(gbase, glf, grf, gbiases, estimate_bias);
        end
        
        function T = exp(g)
            v = Estimation.DLGEKF.State.vee(g);
            T = Estimation.DLGEKF.State.exphat(v);
        end
        
        function T = exphat(v)
            estimate_bias = Estimation.DLGEKF.State.is_estimate_bias_enabled(v);
            [vbase, vlf, vrf, vbias] = Estimation.DLGEKF.State.splitVector(v, estimate_bias);
            
            Xbase = LieGroups.SE_2_3.exphat(vbase);
            [Rb, pb, vb] = LieGroups.SE_2_3.extractSE23(Xbase);

            Xlf = LieGroups.SE3.exphat(vlf);
            [Rlf, plf] = LieGroups.SE3.extractSE3(Xlf);
            
            Xrf = LieGroups.SE3.exphat(vrf);
            [Rrf, prf] = LieGroups.SE3.extractSE3(Xrf);
            
            ba = [];
            bg = [];
            if estimate_bias
                biases = LieGroups.Tn.exphat(vbias);
                ba = biases(1:3, 7);
                bg = biases(4:6, 7);
                if (size(ba, 1) == 1)
                    ba = ba';
                end
                if (size(bg, 1) == 1)
                    bg = bg';
                end
            end
            
            T = Estimation.DLGEKF.State.construct(Rb, pb, vb, ...
                Rlf, plf, ...
                Rrf, prf, ...
                ba, bg, estimate_bias);
        end
        
        function g = log(X)
            v = Estimation.DLGEKF.State.logvee(X);
            g = Estimation.DLGEKF.State.hat(v);
        end
        
        function v = logvee(X)
            estimate_bias = Estimation.DLGEKF.State.is_estimate_bias_enabled(X);
            [Xbase, Xlf, Xrf, Xbias] = Estimation.DLGEKF.State.extractByParts(X);
            vbase = LieGroups.SE_2_3.logvee(Xbase);
            vlf = LieGroups.SE3.logvee(Xlf);
            vrf = LieGroups.SE3.logvee(Xrf);
            
            v = [vbase; vlf; vrf];
            if estimate_bias
                vbias = LieGroups.Tn.logvee(Xbias);
                v = [v; vbias];
            end
        end
        
        function AdT = AdjointMatrix(X)
            estimate_bias = Estimation.DLGEKF.State.is_estimate_bias_enabled(X);
            [Xbase, Xlf, Xrf, Xbias] = Estimation.DLGEKF.State.extractByParts(X);
            
            AdXbase = LieGroups.SE_2_3.AdjointMatrix(Xbase);
            AdXlf = LieGroups.SE3.AdjointMatrix(Xlf);
            AdXrf = LieGroups.SE3.AdjointMatrix(Xrf);
            
            AdT = blkdiag(AdXbase, AdXlf, AdXrf);
            if estimate_bias
                AdXbias = LieGroups.Tn.AdjointMatrix(Xbias);
                AdT = blkdiag(AdT, AdXbias);
            end            
        end
        
        function adT = crossProductMatrix(v)
            estimate_bias = Estimation.DLGEKF.State.is_estimate_bias_enabled(v);
            [vbase, vlf, vrf, vbias] = Estimation.DLGEKF.State.splitVector(v, estimate_bias);
            
            adXbase = LieGroups.SE_2_3.crossProductMatrix(vbase);
            adXlf = LieGroups.SE3.crossProductMatrix(vlf);
            adXrf = LieGroups.SE3.crossProductMatrix(vrf);
            
            adT = blkdiag(adXbase, adXlf, adXrf);
            if estimate_bias
                adXbias = LieGroups.Tn.crossProductMatrix(vbias);
                adT = blkdiag(adT, adXbias);
            end
        end
        
        function Jr = rightJacobian(v)
            estimate_bias = Estimation.DLGEKF.State.is_estimate_bias_enabled(v);
            [vbase, vlf, vrf, vbias] = Estimation.DLGEKF.State.splitVector(v, estimate_bias);
            
            Jrbase = LieGroups.SE_2_3.rightJacobian(vbase);
            Jrlf = LieGroups.SE3.rightJacobian(vlf);
            Jrrf = LieGroups.SE3.rightJacobian(vrf);
            
            Jr = blkdiag(Jrbase, Jrlf, Jrrf);
            if estimate_bias
                Jrbias = LieGroups.Tn.rightJacobian(vbias);
                Jr = blkdiag(Jr, Jrbias);
            end
        end
        
        function Jl = leftJacobian(v)
            estimate_bias = Estimation.DLGEKF.State.is_estimate_bias_enabled(v);
            [vbase, vlf, vrf, vbias] = Estimation.DLGEKF.State.splitVector(v, estimate_bias);
            
            Jlbase = LieGroups.SE_2_3.leftJacobian(vbase);
            Jllf = LieGroups.SE3.leftJacobian(vlf);
            Jlrf = LieGroups.SE3.leftJacobian(vrf);
            
            Jl = blkdiag(Jlbase, Jllf, Jlrf);
            if estimate_bias
                Jlbias = LieGroups.Tn.leftJacobian(vbias);
                Jl = blkdiag(Jl, Jlbias);
            end
        end
        
        function Jlinv = leftJacobianInverse(v)
            estimate_bias = Estimation.DLGEKF.State.is_estimate_bias_enabled(v);
            [vbase, vlf, vrf, vbias] = Estimation.DLGEKF.State.splitVector(v, estimate_bias);
            
            Jlinvbase = LieGroups.SE_2_3.leftJacobianInverse(vbase);
            Jlinvlf = LieGroups.SE3.leftJacobianInverse(vlf);
            Jlinvrf = LieGroups.SE3.leftJacobianInverse(vrf);
            
            Jlinv = blkdiag(Jlinvbase, Jlinvlf, Jlinvrf);
            if estimate_bias
                Jlinvbias = LieGroups.Tn.leftJacobianInverse(vbias);
                Jlinv = blkdiag(Jlinv, Jlinvbias);
            end
        end
        
        function [Rb, pb, vb, Rlf, plf, Rrf, prf, ba, bg]  = extract(X)
            estimate_bias = Estimation.DLGEKF.State.is_estimate_bias_enabled(X);
            [baseSE23, lfootSE3, rfootSE3, biasesTn] = Estimation.DLGEKF.State.extractByParts(X);

            [Rb, pb, vb] = LieGroups.SE_2_3.extractSE23(baseSE23);
            [Rlf, plf] = LieGroups.SE3.extractSE3(lfootSE3);
            [Rrf, prf] = LieGroups.SE3.extractSE3(rfootSE3);
            
            ba = [0; 0; 0];
            bg = [0; 0; 0];
            if estimate_bias
                biases = LieGroups.Tn.extractTn(biasesTn);
                ba = biases(1:3);
                bg = biases(4:6);
            end
        end
        
        function [Xbase, Xlf, Xrf, Xbiases] = extractByParts(X)
            estimate_bias = Estimation.DLGEKF.State.is_estimate_bias_enabled(X);
            Xbase = X(1:5, 1:5);
            Xlf = X(6:9, 6:9);
            Xrf = X(10:13, 10:13);
            Xbiases = [];
            if estimate_bias
                Xbiases = X(14:20, 14:20);
            end
        end
                
        function [gbase, glf, grf, gbiases] = extractLieAlgebra(g)
            estimate_bias = Estimation.DLGEKF.State.is_estimate_bias_enabled(g);
            gbase = g(1:5, 1:5);
            glf = g(6:9, 6:9);
            grf = g(10:13, 10:13);
            
            gbiases = zeros(7);
            if estimate_bias
                gbiases = g(14:20, 14:20);
            end
        end
        
        function X = construct(Rb, pb, vb, Rlf, plf, Rrf, prf, ba, bg, estimate_bias)
            baseSE23 = LieGroups.SE_2_3.constructSE23(Rb, pb, vb);
            lfootSE3 = LieGroups.SE3.constructSE3(Rlf, plf);
            rfootSE3 = LieGroups.SE3.constructSE3(Rrf, prf);
            biases = LieGroups.Tn.constructTn([ba; bg]);
            
            if estimate_bias
                X = blkdiag(baseSE23, lfootSE3, rfootSE3, biases);
            else
                X = blkdiag(baseSE23, lfootSE3, rfootSE3);
            end
        end
        
        function g = constructLieAlgebra(gbase, glf, grf, gbiases, estimate_bias)
            if estimate_bias
                g = blkdiag(gbase, glf, grf, gbiases);
            else
                g = blkdiag(gbase, glf, grf);
            end
        end
        
        function [vbase, vlf, vrf, vbias] = splitVector(v, estimate_bias)
            vbase = v(1:9);
            vlf = v(10:15);
            vrf = v(16:21);
            
            vbias = [];
            if estimate_bias
                vbias = v(22:27);                            
            end
        end
        
        function estimate_bias = is_estimate_bias_enabled(X)
            estimate_bias = false;
            if isvector(X)
                if length(X) == 27
                    estimate_bias = true;
                end
            else
                if size(X, 2) == 20
                    estimate_bias = true;
                end
            end
        end
        
        function [Pbrot, Ppos, Pvel, Plfpos, Plfrot, Prfpos, Prfrot, Pba, Pbg] = extractStateVarSubBlockEvolutions(Ptraj, estimate_bias)
            Pbrot = zeros(length(Ptraj), 3, 3);
            Ppos = zeros(length(Ptraj), 3, 3);
            Pvel = zeros(length(Ptraj), 3, 3);
            Plfrot = zeros(length(Ptraj), 3, 3);
            Plfpos = zeros(length(Ptraj), 3, 3);
            Prfrot = zeros(length(Ptraj), 3, 3);
            Prfpos = zeros(length(Ptraj), 3, 3);
            
            Pbg = zeros(length(Ptraj), 3, 3);
            Pba = zeros(length(Ptraj), 3, 3);
            
            for iter_idx = 1: length(Ptraj)
                P = squeeze(Ptraj{iter_idx});
                Ppos(iter_idx, :, :) = P(1:3, 1:3);
                Pbrot(iter_idx, :, :) = P(4:6, 4:6);
                Pvel(iter_idx, :, :) = P(7:9, 7:9);                
                
                Plfpos(iter_idx, :, :) = P(10:12, 10:12);
                Plfrot(iter_idx, :, :) = P(13:15, 13:15);
                Prfpos(iter_idx, :, :) = P(16:18, 16:18);
                Prfrot(iter_idx, :, :) = P(19:21, 19:21);
                                
                if estimate_bias
                    Pba(iter_idx, :, :) = P(22:24, 22:24);
                    Pbg(iter_idx, :, :) = P(25:27, 25:27);
                end
            end
        end
        
    end
end

