estimate_bias = false;
X = Estimation.DLGEKF_RIE.State.identity(estimate_bias);
perturbation_vec = Estimation.DLGEKF_RIE.State.logvee(X);
n = length(perturbation_vec);


%% state identity, epsilon 0.01
epsilon = 0.01;
X = Estimation.DLGEKF_RIE.State.identity(estimate_bias);
disp('=================================');
disp(['Comparison Numerical Jacobian with State as Identity and perturbation epsilon = ', num2str(epsilon)])
[Jnum, Janalytical] = checkJacobians(X, n, epsilon)
assert((norm(Jnum - Janalytical) < 1e-7), 'Numerical and analytical Jacobians mismatch - X identity');
disp('=================================');

%% randomize the state,  epsilon 0.01
epsilon = 0.0001;
X_rand = getRandomGroupElement(estimate_bias);
disp('=================================');
disp(['Comparison Numerical Jacobian with State as a random group element and perturbation epsilon = ', num2str(epsilon)])
[Jnum_rand, Janalytical_rand] = checkJacobians(X_rand, n, epsilon)
assert((norm(Jnum_rand - Janalytical_rand) < 1e-7), 'Numerical and analytical Jacobians mismatch - X rand');
disp('=================================');
%% randomize epsilon
epsilon_rand = Core.getRandomDouble(1e-7, 0.1);
X_rand1 = getRandomGroupElement(estimate_bias);
disp('=================================');
disp(['Comparison Numerical Jacobian with State as a random group element and perturbation epsilon = ', num2str(epsilon_rand)])
[Jnum_epsrand, Janalytical_epsrand] = checkJacobians(X_rand1, n, epsilon_rand)
assert((norm(Jnum_epsrand - Janalytical_epsrand) < 1e-7), 'Numerical and analytical Jacobians mismatch - epsilon rand');
disp('=================================');
%%
function [Jnum, Janalytical] = checkJacobians(X, n, epsilon)
Jnum = zeros(6, n);
for i=1:n
    perturbation_vec = zeros(n, 1);
    perturbation_vec(i) = perturbation_vec(i) + epsilon;
    exp_pert = Estimation.DLGEKF_RIE.State.exphat(perturbation_vec);
    X_pert = Estimation.DLGEKF_RIE.State.compose(exp_pert, X);
    hinv = hinvOfXSS(X);
    h = hofXSS(X_pert);
    error = LieGroups.SE3.compose(hinv, h);
    Jnum(:, i) = LieGroups.SE3.logvee(error)/epsilon;
end
Janalytical = calcHSS(X);
end

function h = hofXSS(X)
   [Rb, pb, ~, Rlf, plf, ~, ~, ~, ~]  = Estimation.DLGEKF_RIE.State.extract(X);
   h = [Rb'*Rlf Rb'*(plf - pb); zeros(1, 3) 1];
end

function hinv = hinvOfXSS(X)
    [Rb, pb, ~, Rlf, plf, ~, ~, ~, ~]  = Estimation.DLGEKF_RIE.State.extract(X);
    hinv = [Rlf'*Rb -Rlf'*(plf - pb); zeros(1, 3) 1];
end

function H = calcHSS(X)
    [~, ~, ~, Zl, pl, ~, ~, ~, ~]  = Estimation.DLGEKF_RIE.State.extract(X);
    ZlT_S_of_dl = (Zl')*Utils.skew(pl);
             
             H = [-Zl'    ZlT_S_of_dl   zeros(3)      Zl'  -ZlT_S_of_dl zeros(3) zeros(3); ...
                 zeros(3)        -Zl'   zeros(3) zeros(3)           Zl' zeros(3) zeros(3)];

end

function X_rand = getRandomGroupElement(estimate_bias)
Rb = getRandomRotation();
pb = getRandomVector(-10, 10);
vb = getRandomVector(-10, 10);
Rlf = getRandomRotation();
plf = getRandomVector(-10, 10);
Rrf = getRandomRotation();
prf = getRandomVector(-10, 10);
ba = [0 0 0]';
bg = [0 0 0]';
X_rand = Estimation.DLGEKF_RIE.State.construct(Rb, pb, vb, Rlf, plf, Rrf, prf, ba, bg, estimate_bias);
end

function R = getRandomRotation()
R = Utils.rpy2rot(Utils.getRandomDouble(-pi, pi), ...
                 Utils.getRandomDouble(-pi, pi), ...
                 Utils.getRandomDouble(-pi, pi));
end

function p = getRandomVector(lb, ub)
p = [Utils.getRandomDouble(lb, ub), ...
     Utils.getRandomDouble(lb, ub), ...
     Utils.getRandomDouble(lb, ub)]';
end