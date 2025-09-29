% This Matlab script can be used to generate Fig. 5(a) in the paper:
% R. Liu, M. Li, M. Zafari, B. Ottersten, and A. L. Swindlehurst, 
% ``Multi-domain optimization framework for ISAC: From electromagnetic shaping to network cooperation,''
% IEEE Wireless Commun., under revision.
% Download this paper at: https://arxiv.org/abs/2506.16011v1
% More information can be found at https://rangliu0706.github.io/publications
% Last edited by Rang Liu (rangl2@uci.edu) in 2025-09-28

clear; clc; rng('shuffle');

% ================ System Parameters ================
% Network Configuration
params.J = 6;                    % Number of Access Points (APs)
params.K = 8;                    % Number of User Equipment (UEs)
params.M = 4;                    % Number of antennas per AP (dual-port: 2M RF ports)
params.L = 256;                  % Number of OFDM symbols

% Frequency and Bandwidth
params.fc = 28e9;                % Carrier frequency [Hz]
params.lambda = physconst('LightSpeed') / params.fc;  % Wavelength [m]
params.deltaf = 120e3;           % Subcarrier spacing [Hz]
params.Nsub = 32;                % Number of resource blocks
params.Ns = 4 * params.Nsub;     % Total number of subcarriers (128)
params.BW = params.Ns * params.deltaf;  % Total bandwidth: 15.36 MHz

% Quality of Service
params.Rc = 20;                  % Communication rate requirement [bps/Hz]

% Noise Parameters
params.NF_dB = 7;                % Noise figure for gNB class-A [dB]
P_N_dBm = -174 + 10*log10(params.BW) + params.NF_dB;  % Noise power [dBm]
params.sigma2_r = 10^((P_N_dBm - 30)/10);  % Radar noise variance [W]
params.sigma2_c = params.sigma2_r;         % Communication noise variance [W]

% Radar Target Parameters
params.sigma_b2 = 10;            % Mean Radar Cross Section (RCS) [m²]
params.nPat = 2;                 % RCS directivity exponent

% Polarization Parameters
chi_ant = 0.1;                   % Antenna cross-polarization discrimination
XPD_dB_main = 20;                % Main polarization XPD [dB]
XPD_cross_dB = 20;               % Cross polarization XPD [dB]

params.XPD_ant = (1/sqrt(1 + chi_ant)) * [1, sqrt(chi_ant); sqrt(chi_ant), 1];
params.XPD_env = [1, 1/sqrt(10^(XPD_cross_dB/10));
    1/sqrt(10^(XPD_cross_dB/10)), 1/10^(XPD_dB_main/10)];
params.XPD = sqrtm(params.XPD_env .* params.XPD_ant);
params.Rsqrt_rx = chol(params.XPD, 'lower');           % Receive correlation
params.Rsqrt = kron(eye(params.M), chol(params.XPD, 'lower'));  % Transmit correlation

% Path Loss Parameters
params.FSPL_1m = 32.4 + 20*log10(params.fc/1e9);  % Free space path loss at 1m [dB]
params.nPLE_LOS = 2.0;           % Path loss exponent for LoS links
params.nPLE_NLOS = 3.0;          % Path loss exponent for NLoS links
params.sigma_SF_AP = 2;          % Shadow fading std for AP-Target & AP-AP [dB]
params.sigma_SF_UE = 5;          % Shadow fading std for AP-UE [dB]

% Network Geometry
params.Dradius = 50;             % Network coverage radius [m]

% ================ Simulation Settings ================
N_sim = 100;                     % Number of Monte Carlo trials
Power_dBm = 30:4:50;             % Transmit power range [dBm]

SINR_prop = zeros(numel(Power_dBm), N_sim);   % Proposed method
SINR_radar = zeros(numel(Power_dBm), N_sim);  % Radar-only benchmark
SINR_EM = zeros(numel(Power_dBm), N_sim);     % Baseband and EM-domain optimization
SINR_BP = zeros(numel(Power_dBm), N_sim);     % Baseband-domain only
SINR_NC = zeros(numel(Power_dBm), N_sim);     % Baseband and Network-domain optimization

% Generate AP positions uniformly on a circle
theta_AP_pos = (0:params.J-1).' * 2*pi / params.J;
ap_positions = params.Dradius * [cos(theta_AP_pos), sin(theta_AP_pos)];

% ================ Monte Carlo Simulation ================
alg.maxIter = 50; alg.res_thr = 1e-2;
% cvx_solver Mosek
for sim_idx = 1:N_sim
    tic;

    % Generate random target position
    theta_t = 2*pi * rand();   r_t = params.Dradius * sqrt(rand());
    target_pos = [r_t * cos(theta_t), r_t * sin(theta_t)];

    % Generate random UE positions
    theta_ue = 2*pi * rand(params.K, 1);
    r_ue = params.Dradius * sqrt(rand(params.K, 1));
    ue_positions = [r_ue .* cos(theta_ue), r_ue .* sin(theta_ue)];

    % Calculate angles from APs to target
    theta_AP = zeros(params.J, 1);
    for j = 1:params.J
        delta = target_pos - ap_positions(j, :);
        theta_AP(j) = mod(atan2(delta(2), delta(1)), 2*pi);
    end

    % Generate channels
    Gijn = generate_ap_to_ap_channels(ap_positions, params);
    Hjnk = generate_ap_to_ue_channels(ap_positions, ue_positions, params);
    Phiij = generate_target_scattering(ap_positions, target_pos, params);

    for p_idx = 1:numel(Power_dBm)     
        fprintf('Simulation %d/%d, Power level %d/%d (%.1f dBm)\n', ...
                sim_idx, N_sim, p_idx, numel(Power_dBm), Power_dBm(p_idx));

        P_linear = 10^((Power_dBm(p_idx) - 30)/10);

        % Radar-only
        v_radar = opt_radar(theta_AP, Gijn, Phiij, P_linear, params.sigma2_r, params.XPD, params.Nsub, alg);
        SINR_radar(p_idx,sim_idx) = SINR_radar(p_idx,sim_idx) + v_radar;

        % Proposed cross-domain joint optimization
        v_prop = opt_prop(theta_AP, Hjnk, Gijn, Phiij, P_linear, ...
            params.sigma2_r, params.sigma2_c, params.Nsub, params.XPD, params.Rc, alg);
        SINR_prop(p_idx,sim_idx) = SINR_prop(p_idx,sim_idx) + v_prop;

        % EM and baseband domains only
        v_eml = opt_BP_EM(theta_AP, Hjnk, Gijn, Phiij, P_linear, ...
            params.sigma2_r, params.sigma2_c, params.Nsub, params.XPD, params.Rc, alg);
        SINR_EM(p_idx,sim_idx) = SINR_EM(p_idx,sim_idx) + v_eml;

        % Baseband domain only
        v_spl = opt_BP(theta_AP, Hjnk, Gijn, Phiij, P_linear, ...
            params.sigma2_r, params.sigma2_c, params.Nsub, params.XPD, params.Rc, alg);
        SINR_BP(p_idx,sim_idx) = SINR_BP(p_idx,sim_idx) + v_spl;

        % Network and baseband domains
        v_ncl = opt_BP_NC(theta_AP, Hjnk, Gijn, Phiij, P_linear, ...
            params.sigma2_r, params.sigma2_c, params.Nsub, params.XPD, params.Rc, alg);
        SINR_NC(p_idx,sim_idx) = SINR_NC(p_idx,sim_idx) + v_ncl;
    end
    toc
    save('data_power.mat');
end

processing_gain = params.L;
N_sim = sim_idx;

SINR_prop = processing_gain * sum(SINR_prop,2) / N_sim;
SINR_radar = processing_gain * sum(SINR_radar,2) / N_sim;
SINR_EM = processing_gain * sum(SINR_EM,2) / N_sim;
SINR_BP = processing_gain * sum(SINR_BP,2) / N_sim;
SINR_NC = processing_gain * sum(SINR_NC,2) / N_sim;

figure; hold on; box on;
plot(Power_dBm, 10*log10(SINR_radar),'-s','Color',[0.15 0.15 0.15],'LineWidth',1.5);
plot(Power_dBm, 10*log10(SINR_prop),'-o','Color',[0.8 0 0],'LineWidth',1.5);
plot(Power_dBm, 10*log10(SINR_NC),'-^','Color',[0 0.45 0.74],'LineWidth',1.5);
plot(Power_dBm, 10*log10(SINR_EM),'-d','Color',[0.4 0.6 0.1],'LineWidth',1.5);
plot(Power_dBm, 10*log10(SINR_BP),'-x','Color',[0.49 0.18 0.56],'LineWidth',1.5);
xlabel('Total transmit power (dBm)');
ylabel('Radar SINR (dB)');
grid on;
xticks(Power_dBm);
legend('Radar-only','Proposed','BP-NC','BP-EM','BP-only');
axis([30 50 -17 13])
