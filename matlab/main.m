%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% copyright 2021 Istituto Italiano di Tecnologia (IIT). This software may be modified and
% distributed under the terms of the GNU Lesser General Public License v2.1 or any later version.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% clear workspace
clc
clear 
close all

%% flags
enable_ocekf = false;
enable_invekf = false;
enable_swa = false;
enable_diligent = true;
enable_diligent_rie = true;

%% load dataset
experiment_name = 'com-sinusoid'; % walking | com-sinusoid
experiment.mat_file = ['./resources/' experiment_name '.mat'];

load(experiment.mat_file);

experiment.robot_name='iCubGenova04'; %% Name of the robot
experiment.model_path = './resources/';  %% Path to the robot model
experiment.model_file = 'model.urdf';

experiment.joints_list = ["neck_pitch", "neck_roll", "neck_yaw", ...
                          "torso_pitch", "torso_roll", "torso_yaw", ...
                          "l_shoulder_pitch", "l_shoulder_roll", "l_shoulder_yaw", ...
                          "l_elbow", ...
                          "r_shoulder_pitch", "r_shoulder_roll", "r_shoulder_yaw", ...
                          "r_elbow", ...
                          "l_hip_pitch", "l_hip_roll", "l_hip_yaw", ...
                          "l_knee", "l_ankle_pitch", "l_ankle_roll", ...
                          "r_hip_pitch", "r_hip_roll", "r_hip_yaw", ...
                          "r_knee", "r_ankle_pitch", "r_ankle_roll"];
experiment.root_link = 'root_link'; 
%% Timing iterators
estimatorJointsTime = estimatorJointsTime - estimatorJointsTime(1);
dt = mean(diff(estimatorJointsTime));

startIter = 1;
endIter = length(estimatorJointsTime);
nrIters = endIter - startIter + 1;

%% load iDynTree KinDynComputations
experiment.kindyn_debug = false;
KinDynModel = iDynTreeWrappers.loadReducedModel(experiment.joints_list, ...
                                                experiment.root_link, ...
                                                experiment.model_path, ...
                                                experiment.model_file, ...
                                                experiment.kindyn_debug);
kinDyn = KinDynModel.kinDynComp;

%% Load ModelComputations
lfootLink = 'l_foot';
rfootLink = 'r_foot';
baseLink = 'root_link';
baseLinkIMU = 'root_link_imu_acc';
lfVertexNames = {'l_sole'};
rfVertexNames = {'r_sole'};

LFVertexIDs = -ones(length(lfVertexNames), 1);
for idx = 1:length(lfVertexNames)
    LFVertexIDs(idx) = kinDyn.model().getFrameIndex(lfVertexNames{idx});
end

RFVertexIDs = -ones(length(lfVertexNames), 1);
for idx = 1:length(lfVertexNames)
    RFVertexIDs(idx) = kinDyn.model().getFrameIndex(rfVertexNames{idx});
end

LSoleIds = [LFVertexIDs(1)];
RSoleIds = [RFVertexIDs(1)];
SoleIds = [LSoleIds; RSoleIds];

modelComp = Model.ModelComputations(kinDyn, baseLink, baseLinkIMU, ...
                                    lfootLink, rfootLink, ...
                                    LFVertexIDs, RFVertexIDs);

%% Load Estimator options
options = Estimation.Proprioception.EstimatorOptions();
options.nr_joints_est = nr_joints_est;
options.enable_bias_estimation = true;

%% Setup initial states
rpy0 = linkBaseRot(startIter ,:);
R0 = Utils.rpy2rot(rpy0(1), rpy0(2), rpy0(3));
pos0 = linkBasePos(startIter ,:)';
pose0 = LieGroups.SE3.constructSE3(R0, pos0);

s = estimatorJointsPos(1, :)';
sDot = zeros(size(estimatorJointsPos(1, :)'));
modelComp.setRobotState(pose0, zeros(6, 1), s, sDot);

b_H_w0 = modelComp.kindyn.getWorldBaseTransform().inverse();
pW0 = b_H_w0.getPosition().toMatlab();
qW0 = b_H_w0.getRotation().asQuaternion().toMatlab();

w_H_IMU0 = modelComp.kindyn.getWorldTransform(modelComp.base_link_imu_idx);
pIMU0 = w_H_IMU0.getPosition().toMatlab();
qIMU0 = w_H_IMU0.getRotation().asQuaternion().toMatlab();

w_H_LF0 = modelComp.kindyn.getWorldTransform(LFVertexIDs(1));
pLF0 = w_H_LF0.getPosition().toMatlab();
qLF0 = w_H_LF0.getRotation().asQuaternion().toMatlab();

w_H_RF0 = modelComp.kindyn.getWorldTransform(RFVertexIDs(1));
pRF0 = w_H_RF0.getPosition().toMatlab();
qRF0 = w_H_RF0.getRotation().asQuaternion().toMatlab();

Nl = length(LFVertexIDs);
Nr = length(RFVertexIDs);

for idx = 1:Nl
    frameName = modelComp.kindyn.model.getFrameName(LFVertexIDs(idx));
    dl = modelComp.kindyn.getWorldTransform(frameName).getPosition().toMatlab();
end

for idx = 1:Nr
    frameName = modelComp.kindyn.model.getFrameName(RFVertexIDs(idx));
    dr = modelComp.kindyn.getWorldTransform(frameName).getPosition().toMatlab();
end

ba0 = zeros(3, 1);
bg0 = zeros(3, 1);

if (options.enable_bias_estimation)
    initial_state = zeros(30, 1);
else
    initial_state = zeros(24, 1);
end
initial_v = zeros(3, 1);
initial_state(1:4) = qIMU0; % imu orientation
initial_state(5:7) = pIMU0; % imu position
initial_state(8:10) = initial_v; % imu linear velocity
initial_state(11:14) = qLF0;  % left foot orientation
initial_state(15:17) = pLF0; % left foot position
initial_state(18:21) = qRF0;  % right foot orientation
initial_state(22:24) = pRF0; % right foot position


if (options.enable_bias_estimation)
    initial_state(25:27) = ba0; % acc bias
    initial_state(28:30) = bg0;  % gyro bias
end

[initial_invekf_state, initial_theta] = Estimation.InvEKF.State.construct(w_H_IMU0.getRotation().toMatlab(), initial_v, initial_state(5:7), ...
        initial_state(22:24), initial_state(15:17), bg0, ba0, []);
    
initial_dlgekf_state = Estimation.DLGEKF.State.construct(w_H_IMU0.getRotation().toMatlab(), initial_state(5:7), initial_v, ...
        w_H_LF0.getRotation().toMatlab(), initial_state(15:17), ...
        w_H_RF0.getRotation().toMatlab(), initial_state(22:24), ...
        ba0, bg0, options.enable_bias_estimation);

%% setup priors
priors = Config.configMatlabPriors();

%% setup sensor noise
sensors_dev = Config.configMatlabSensorDev(options.nr_joints_est);

%% setup schmitt trigger
[leftSchmittParams, rightSchmittParams] = Config.configSchmittParams;
leftFootContactStateMachine = Estimation.ContactHandler.ContactStateMachine(leftSchmittParams);
rightFootContactStateMachine = Estimation.ContactHandler.ContactStateMachine(rightSchmittParams);

%% Estimator collectors - collect outputs for plotting
outMap = containers.Map('KeyType','char','ValueType','any');

gt = 'gt';
outMap(gt) = EstCollector(gt);
outMap(gt).resizeBuffers(nrIters, Nl, Nr);

if enable_ocekf
    ocekfmat = 'ocekfmat';
    outMap(ocekfmat) = EstCollector(ocekfmat);
    outMap(ocekfmat).resizeBuffers(nrIters, Nl, Nr);
end

if enable_invekf
    invekfmat = 'invekfmat';
    outMap(invekfmat) = EstCollector(invekfmat);
    outMap(invekfmat).resizeBuffers(nrIters, Nl, Nr);
end

if enable_swa
    swamat = 'swamat';
    outMap(swamat) = EstCollector(swamat);
    outMap(swamat).resizeBuffers(nrIters, Nl, Nr);
end

if enable_diligent
    diligentmat = 'diligentmat';
    outMap(diligentmat) = EstCollector(diligentmat);
    outMap(diligentmat).resizeBuffers(nrIters, Nl, Nr);
end


if enable_diligent_rie
    diligentriemat = 'diligentriemat';
    outMap(diligentriemat) = EstCollector(diligentriemat);
    outMap(diligentriemat).resizeBuffers(nrIters, Nl, Nr);
end

%% Initialize filters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ocekf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if enable_ocekf
    ocekf = Estimation.RotellaEstimator.Filter();
    ocekf.setup(priors, sensors_dev, modelComp, options);
    ocekf.initialize(initial_state, ba0, bg0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% invekf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if enable_invekf
    invekfm = Estimation.InvEKF.Filter();
    invekfm.setup(priors, sensors_dev, modelComp, options);
    invekfm.initialize(initial_invekf_state, initial_theta);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% swa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if enable_swa
    attest_x0 = [initial_state(1:4); zeros(6, 1)];
    attest_type='qekf';
    fuse_position = true;
    primary_foot = 'right';
    flat_floor = false;
    switching_pattern = 'alternate';
    qekfParams = Config.configQEKFParams(dt);
    swa = Estimation.SimpleBipedEstimator.LO_IMU(attest_type);
    swa.setup(modelComp, primary_foot, flat_floor, switching_pattern, qekfParams, fuse_position);
    swa.initialize('base_link', b_H_w0.asHomogeneousTransform().toMatlab(), s, attest_x0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diligent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if enable_diligent
    diligentm = Estimation.DLGEKF.Filter();
    diligentm.setup(priors, sensors_dev, modelComp, options);
    diligentm.initialize(initial_dlgekf_state);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diligent_rie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if enable_diligent_rie
    diligentriem = Estimation.DLGEKF_RIE.Filter();
    diligentriem.setup(priors, sensors_dev, modelComp, options);
    diligentriem.initialize(initial_dlgekf_state);
end

%% init buffers
lfC = zeros(length(estimatorJointsTime), 1);
rfC = zeros(length(estimatorJointsTime), 1);

jPosDyn = iDynTree.JointPosDoubleArray(modelComp.kindyn.model);
XPrev = [];


%% Run filters (main loop)
for iter = startIter : endIter
    t = estimatorJointsTime(iter);
    s = estimatorJointsPos(iter, :)';
    sdot = estimatorJointsVel(iter, :)';
    fzl = lfForceZ(iter);
    fzr = rfForceZ(iter);
    acc = imuAcc(iter, :);
    gyro = imuOmega(iter, :);
       
    % update contact states using Schmitt Trigger
    leftFootContactStateMachine.contactMeasurementUpdate(fzl, t);
    rightFootContactStateMachine.contactMeasurementUpdate(fzr, t);
    contacts = [leftFootContactStateMachine.contactState(); ...
        rightFootContactStateMachine.contactState()];            
            
    lfC(iter) = contacts(1);
    rfC(iter) = contacts(2);    
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ocekf
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if enable_ocekf
        [Xocekf, ~, ocekfBase, ~] = ocekf.advance(t, gyro', acc', s, contacts);
        w_H_B_ocekf = Utils.idynPose(Utils.quat2rot(ocekfBase.q), ocekfBase.I_p);
        outMap(ocekfmat).updateBasePose(iter, w_H_B_ocekf.getRotation().asRPY().toMatlab(), ...
            ocekfBase.I_p, ocekfBase.I_pdot, ocekfBase.I_omega);
        [~,~,~, ...
            qLF_ocekf,p_LF_ocekf, ...
            qRF_ocekf,p_RF_ocekf, ...
            ba_ocekf, bg_ocekf] = Estimation.RotellaEstimator.State.extract(Xocekf, ba0, bg0);
        outMap(ocekfmat).updateFootRotation('left',iter, Utils.rot2rpy(Utils.quat2rot(qLF_ocekf)));
        outMap(ocekfmat).updateFootRotation('right',iter, Utils.rot2rpy(Utils.quat2rot(qRF_ocekf)));
        outMap(ocekfmat).updateFootPosition('left',iter, 1, p_LF_ocekf);
        outMap(ocekfmat).updateFootPosition('right',iter, 1, p_RF_ocekf);
        outMap(ocekfmat).updateBias(iter, ba_ocekf, bg_ocekf);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % invekfmat
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if enable_invekf
        [Xinvekf, thetainekf, ~, invekfmBase, ~] = invekfm.advance(t, gyro', acc', s, contacts);
        w_H_B_invekfm = Utils.idynPose(Utils.quat2rot(invekfmBase.q), invekfmBase.I_p);
        outMap(invekfmat).updateBasePose(iter, w_H_B_invekfm.getRotation().asRPY().toMatlab(), ...
            invekfmBase.I_p, invekfmBase.I_pdot, invekfmBase.I_omega);
        [~, ~, ~, ...
            p_RF_invekf, p_LF_invekf, ...
            ba_invekf, bg_invekf, ~] = Estimation.InvEKF.State.extract(Xinvekf, thetainekf);
        outMap(invekfmat).updateFootPosition('left',iter, 1, p_LF_invekf);
        outMap(invekfmat).updateFootPosition('right',iter, 1, p_RF_invekf);
        outMap(invekfmat).updateBias(iter, ba_invekf, bg_invekf);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % swa
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if enable_swa
        mag = [0;0;0]; % unused
        [~, swamBase, ~, ~] = swa.advance(t, acc', gyro', mag,...
            s, sDot,...
            contacts, zeros(6, 1), zeros(6, 1));
        
        w_H_B_swam = Utils.idynPose(Utils.quat2rot(swamBase.q), swamBase.I_p);
        outMap(swamat).updateBasePose(iter, w_H_B_swam.getRotation().asRPY().toMatlab(), ...
            swamBase.I_p, swamBase.I_pdot, swamBase.I_omega);
    end
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % diligentmat
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if enable_diligent
        [Xdil, Pk, diligentmBase, ~] = diligentm.advance(t, gyro', acc', s, contacts);
        w_H_B_diligentm = Utils.idynPose(Utils.quat2rot(diligentmBase.q), diligentmBase.I_p);
        outMap(diligentmat).updateBasePose(iter, w_H_B_diligentm.getRotation().asRPY().toMatlab(), ...
            diligentmBase.I_p, diligentmBase.I_pdot, diligentmBase.I_omega);
        [~, ~, ~, ...
            Rdil_LF, pdil_LF, ...
            Rdil_RF, pdil_RF, ...
            bias_acc_dil, bias_gyro_dil] = Estimation.DLGEKF.State.extract(Xdil);
        outMap(diligentmat).updateFootRotation('left',iter, Utils.rot2rpy(Rdil_LF));
        outMap(diligentmat).updateFootRotation('right',iter, Utils.rot2rpy(Rdil_RF));
        outMap(diligentmat).updateFootPosition('left',iter, 1, pdil_LF);
        outMap(diligentmat).updateFootPosition('right',iter, 1, pdil_RF);
        outMap(diligentmat).updateBias(iter, bias_acc_dil, bias_acc_dil);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % diligentriemat
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if enable_diligent_rie
        [Xdilrie, Pkrie, diligentriemBase, ~] = diligentriem.advance(t, gyro', acc', s, contacts);
        w_H_B_diligentriem = Utils.idynPose(Utils.quat2rot(diligentriemBase.q), diligentriemBase.I_p);
        outMap(diligentriemat).updateBasePose(iter, w_H_B_diligentriem.getRotation().asRPY().toMatlab(), ...
            diligentriemBase.I_p, diligentriemBase.I_pdot, diligentriemBase.I_omega);
        [~, ~, ~, ...
            Rdilrie_LF, pdilrie_LF, ...
            Rdilrie_RF, pdilrie_RF, ...
            bias_acc_dilrie, bias_gyro_dilrie] = Estimation.DLGEKF_RIE.State.extract(Xdilrie);
        outMap(diligentriemat).updateFootRotation('left',iter, Utils.rot2rpy(Rdilrie_LF));
        outMap(diligentriemat).updateFootRotation('right',iter, Utils.rot2rpy(Rdilrie_RF));
        outMap(diligentriemat).updateFootPosition('left',iter, 1, pdilrie_LF);
        outMap(diligentriemat).updateFootPosition('right',iter, 1, pdilrie_RF);
        outMap(diligentriemat).updateBias(iter, bias_acc_dilrie, bias_acc_dilrie);
    end
   
%     pause
    if (mod(iter, 50) == 0)
        disp(['========> Iter: ' num2str(iter) '/' num2str(endIter-startIter)])
    end
end

%% Plot base rotation and position
estBaseTime = estimatorJointsTime;

figure
color = {'r', 'g', 'b'};
color2 = {'g', 'b', 'r'};

title1 = {'Roll (deg)', 'Pitch (deg)', 'Yaw (deg)'};
title2 = {'x (m)', 'y (m)', 'z (m)'};
for idx = 1:2
    subplot(2, 1, idx)
    simp = plot(estBaseTime(startIter:endIter), rad2deg(linkBaseRot(startIter:endIter, idx)), 'LineWidth', 4, 'Color', 'k', 'LineStyle', '--');
    legend_line = [];
    legendtex = {};
    legend_line = [legend_line simp];
    legendtex = [legendtex 'Vicon'];
    hold on
    if enable_ocekf
        ocekfmatline = plot(estBaseTime(startIter:endIter), rad2deg(outMap(ocekfmat).baseRPY(startIter:endIter, idx)),  'Color',[0.9290,0.6940, 0.1250], 'LineWidth', 2, 'LineStyle', '--');
        legend_line = [legend_line ocekfmatline];
        legendtex = [legendtex 'OCEKF'];
        
    end
    if enable_invekf
        invekfmatline = plot(estBaseTime(startIter:endIter), rad2deg(outMap(invekfmat).baseRPY(startIter:endIter, idx)), 'r', 'LineWidth', 2, 'LineStyle', '--');
        legend_line = [legend_line invekfmatline];
        legendtex = [legendtex 'InvEKF'];
    end
    if enable_swa
        swamatline = plot(estBaseTime(startIter:endIter), rad2deg(outMap(swamat).baseRPY(startIter:endIter, idx)),  'LineWidth', 2, 'LineStyle', '--', 'Color', [0,0.4470, 0.7410]);
        legend_line = [legend_line swamatline];
        legendtex = [legendtex 'SWA'];
    end
    if enable_diligent
        diligentmatline = plot(estBaseTime(startIter:endIter), rad2deg(outMap(diligentmat).baseRPY(startIter:endIter, idx)),  'g', 'LineWidth', 2);    
        legend_line = [legend_line diligentmatline];
        legendtex = [legendtex 'DILIGENT-KIO'];
    end
    if enable_diligent_rie
        diligentriematline = plot(estBaseTime(startIter:endIter), rad2deg(outMap(diligentriemat).baseRPY(startIter:endIter, idx)),  'b', 'LineWidth', 2);    
        legend_line = [legend_line diligentriematline];
        legendtex = [legendtex 'DILIGENT-RIE-KIO'];
    end
    xlabel('Time(s)', 'FontSize', 18)
    ylabel(title1{idx}, 'FontSize', 24) 
    legend(legend_line, legendtex,'FontSize', 24)    
    xlim([0.03 estBaseTime(end)]) % remove uninitialized part of SWA from plots for proper visualization
    set(gca,'FontSize',20)
end
sgtitle('Base Orientation: Roll and Pitch', ...
    'FontSize', 24);
drawnow;

figure
for idx = 1:3
    subplot(3, 1, idx)
    simp = plot(estBaseTime(startIter:endIter), linkBasePos(startIter:endIter, idx), 'LineWidth', 4, 'Color', 'k', 'LineStyle', '--');
    legend_line = [];
    legendtex = {};
    legend_line = [legend_line simp];
    legendtex = [legendtex 'Vicon'];
    hold on
    if enable_ocekf
        ocekfmatline = plot(estBaseTime(startIter:endIter), outMap(ocekfmat).basePos(startIter:endIter, idx), 'Color',[0.9290,0.6940, 0.1250], 'LineWidth', 2, 'LineStyle', '--');
        legend_line = [legend_line ocekfmatline];
        legendtex = [legendtex 'OCEKF'];
    end
    if enable_invekf
        invekfmatline = plot(estBaseTime(startIter:endIter), outMap(invekfmat).basePos(startIter:endIter, idx), 'r', 'LineWidth', 2, 'LineStyle', '--');
        legend_line = [legend_line invekfmatline];
        legendtex = [legendtex 'InvEKF'];
    end
    if enable_swa
        swamatline = plot(estBaseTime(startIter:endIter), outMap(swamat).basePos(startIter:endIter, idx), 'LineWidth', 2, 'LineStyle', '--', 'Color', [0,0.4470, 0.7410]);
        legend_line = [legend_line swamatline];
        legendtex = [legendtex 'SWA'];
    end
    if enable_diligent
        diligentmatline = plot(estBaseTime(startIter:endIter), outMap(diligentmat).basePos(startIter:endIter, idx),   'g', 'LineWidth', 2);
        legend_line = [legend_line diligentmatline];
        legendtex = [legendtex 'DILIGENT-KIO'];
    end
    if enable_diligent_rie
        diligentriematline = plot(estBaseTime(startIter:endIter), outMap(diligentriemat).basePos(startIter:endIter, idx),   'b', 'LineWidth', 2);
        legend_line = [legend_line diligentriematline];
        legendtex = [legendtex 'DILIGENT-RIE-KIO'];
    end
    
    xlabel('Time(s)', 'FontSize', 18)
    ylabel(title2{idx}, 'FontSize', 24)    
    xlim([0.03 estBaseTime(end)]) % remove uninitialized part of SWA from plots for proper visualization
    legend(legend_line,legendtex,'FontSize', 24)    
    set(gca,'FontSize',20)
end
sgtitle('Base Position', ...
    'FontSize', 24);
drawnow;

%% Compute Vicon linear velocity by numerical differentiatio
linkBaseLinVel = [];
for idx = 1:3    
    linkBaseLinVel(:, idx) = [0; smooth(diff(linkBasePos(:, idx))/dt, 'sgolay',  1)];
end

%% Plot Base Linear Velocity
figure
for idx = 1:3
    subplot(3, 1, idx)
    simp = plot(estBaseTime(startIter:endIter), linkBaseLinVel(startIter:endIter, idx), 'LineWidth', 4, 'Color', 'k', 'LineStyle', '--');
    legend_line = [];
    legendtex = {};
    legend_line = [legend_line simp];
    legendtex = [legendtex 'Vicon'];
    hold on
    if enable_ocekf
        ocekfmatline = plot(estBaseTime(startIter:endIter), outMap(ocekfmat).baseLinVel(startIter:endIter, idx), 'Color',[0.9290,0.6940, 0.1250], 'LineWidth', 2, 'LineStyle', '--');
        legend_line = [legend_line ocekfmatline];
        legendtex = [legendtex 'OCEKF'];
    end
    if enable_invekf
        invekfmatline = plot(estBaseTime(startIter:endIter), outMap(invekfmat).baseLinVel(startIter:endIter, idx), 'r', 'LineWidth', 2, 'LineStyle', '--');
        legend_line = [legend_line invekfmatline];
        legendtex = [legendtex 'InvEKF'];
    end
    if enable_swa
        swamatline = plot(estBaseTime(startIter:endIter), outMap(swamat).baseLinVel(startIter:endIter, idx), 'LineWidth', 2, 'LineStyle', '--', 'Color', [0,0.4470, 0.7410]);
        legend_line = [legend_line swamatline];
        legendtex = [legendtex 'SWA'];
    end
    if enable_diligent
        diligentmatline = plot(estBaseTime(startIter:endIter), outMap(diligentmat).baseLinVel(startIter:endIter, idx),   'g', 'LineWidth', 2);
        legend_line = [legend_line diligentmatline];
        legendtex = [legendtex 'DILIGENT-KIO'];
    end
    if enable_diligent_rie
        diligentriematline = plot(estBaseTime(startIter:endIter), outMap(diligentriemat).baseLinVel(startIter:endIter, idx),   'b', 'LineWidth', 2);
        legend_line = [legend_line diligentriematline];
        legendtex = [legendtex 'DILIGENT-RIE-KIO'];
    end
    
    xlabel('Time(s)', 'FontSize', 18)
    ylabel(title2{idx}, 'FontSize', 24)    
    xlim([0.03 estBaseTime(end)]) % remove uninitialized part of SWA from plots for proper visualization
    legend(legend_line,legendtex,'FontSize', 24)    
    set(gca,'FontSize',20)
end
sgtitle('Base Linear Velocity', ...
    'FontSize', 24);
drawnow;

%% Compute errors

RPEgap = 100;
align_yaw = false;
% [rotError, posError, velError, ATErot, ATEpos, ATEvel, RPErot, RPEpos]
if enable_ocekf
[errors.ocekfL.rotError, errors.ocekfL.posError, errors.ocekfL.velError, ...
    errors.ocekfL.ATErot, errors.ocekfL.ATEpos, errors.ocekfL.ATEvel,  ...
    errors.ocekfL.RPErot, errors.ocekfL.RPEpos] = Estimation.Utils.getLeftInvariantErrorMetrics(linkBasePos, linkBaseRot, linkBaseLinVel, ...
                                               outMap(ocekfmat).basePos, outMap(ocekfmat).baseRPY, outMap(ocekfmat).baseLinVel, RPEgap, align_yaw);
end
if enable_invekf
[errors.invekfL.rotError, errors.invekfL.posError, errors.invekfL.velError, ...
    errors.invekfL.ATErot, errors.invekfL.ATEpos, errors.invekfL.ATEvel, ...
    errors.invekfL.RPErot, errors.invekfL.RPEpos] = Estimation.Utils.getLeftInvariantErrorMetrics(linkBasePos, linkBaseRot, linkBaseLinVel, ...
                                               outMap(invekfmat).basePos, outMap(invekfmat).baseRPY, outMap(invekfmat).baseLinVel, RPEgap, align_yaw);
end
if enable_swa
[errors.swaL.rotError, errors.swaL.posError, errors.swaL.velError, ... 
    errors.swaL.ATErot, errors.swaL.ATEpos, errors.swaL.ATEvel, ... 
    errors.swaL.RPErot, errors.swaL.RPEpos] = Estimation.Utils.getLeftInvariantErrorMetrics(linkBasePos, linkBaseRot, linkBaseLinVel, ...
                                               outMap(swamat).basePos, outMap(swamat).baseRPY, outMap(swamat).baseLinVel, RPEgap, align_yaw);                                           
end
if enable_diligent      
[errors.diligentL.rotError, errors.diligentL.posError, errors.diligentL.velError, ...
    errors.diligentL.ATErot, errors.diligentL.ATEpos, errors.diligentL.ATEvel, ... 
    errors.diligentL.RPErot, errors.diligentL.RPEpos] = Estimation.Utils.getLeftInvariantErrorMetrics(linkBasePos, linkBaseRot, linkBaseLinVel, ...
                                               outMap(diligentmat).basePos, outMap(diligentmat).baseRPY, outMap(diligentmat).baseLinVel, RPEgap, align_yaw);
end
if enable_diligent_rie
[errors.diligentrieR.rotError, errors.diligentrieR.posError, errors.diligentrieR.velError, ...
    errors.diligentrieR.ATErot, errors.diligentrieR.ATEpos, errors.diligentrieR.ATEvel, ... 
    errors.diligentrieR.RPErot, errors.diligentrieR.RPEpos] = Estimation.Utils.getRightInvariantErrorMetrics(linkBasePos, linkBaseRot, linkBaseLinVel, ...
                                               outMap(diligentriemat).basePos, outMap(diligentriemat).baseRPY, outMap(diligentriemat).baseLinVel, RPEgap, align_yaw);
end