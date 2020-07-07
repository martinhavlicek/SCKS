% SCKS_hdm
%--------------------------------------------------------------------------
% Example of inverison of hemodynamic model with SCKF-SCKS algorithm
%
clear all; close all;
addpath(genpath('C:\Backup\spm8\'));
% it needs paths to SPM8 toolbox
%==========================================================================
f       = spm_figure('GetWin','Graphics');
%------------------------------------------------------------------
M(1).E.linear = 0;                          % linear model
M(1).E.s      = 1;                          % smoothness
M(1).E.dt     = 0.1;                        % integration step for simulation 
% level 1
%------------------------------------------------------------------
[pE pC] = spm_hdm_priors(1,3);              % parameters
pE(6) = 0.02;
pE(7) = 0.54;
M(1).n  = 4;
M(1).f  = 'spm_fx_hdm';
M(1).g  = 'spm_gx_hdm';
M(1).pE = pE;                               % prior expectation
M(1).V  = exp(6);                           % error precision
M(1).W  = exp(10);                          % error precision

% level 2
%------------------------------------------------------------------
M(2).l  = exp(0);                                % inputs
M(2).V  = exp(0);                                % with shrinkage priors

M       = spm_DEM_M_set(M);

% free parameters
%--------------------------------------------------------------------------
P       = M(1).pE;                                % true parameters
ip      = [1:length(P)];                          % free parameters
pE      = spm_vec(P);
np      = length(pE);
pE      = spm_unvec(pE,P);
M(1).pE = pE;
M(1).pC = pC;
M(1).ip = ip;

% generate data
%==========================================================================
% DEM-steps
N         = 40/M(1).E.dt;                      % length of data sequence ()
U         = (exp(-([1:15] - 8).^2/(1.^2)));      % this is the Gaussian cause
input     = zeros(1,N);
input([21 37 54 82 99 200 217]) = [1 1 1 1 1 1 1];
U         = conv(conv(U,ones(1,2)),input);
U         = U(1:N)/max(U(:));
SCKS      = spm_DEM_generate(M,U,{P},{2,8},{3});
spm_DEM_qU(SCKS.pU)

SCKS.Y = SCKS.Y(1:1/M(1).E.dt:end) + 0.05*randn(1,40); 

pE([1 2 3])  = P([1 2 3])+ 1/12*randn(3,1);  % randomize first 3 parameters
SCKS.M(1).pE = pE;
SCKS.M(1).E.dt   = 1;    
SCKS.M(1,1).E.nD = 10;

% option to specify constrains on paramters (example, test)
cb(1,:) = [0.5,0.8];    % low and high bound
cb(2,:) = [0.3,0.5];    % low and high bound
cb(3,:) = [0.7,1.3];    % low and high bound
cb(4,:) = [0.29,0.35];  % low and high bound
cb(5,:) = [0.32,.4];   % low and high bound
cb(6,:) = [0.01,0.04]; % low and high bound
cb(7,:) = [0.5,0.6];   % low and high bound

SCKS.M(1).ip = ip;  % indices of model parameters to be estimated
SCKS.M(1).cb = [];  % option to specify constrain on parameters values [min max]
SCKS.M(2).v  = 0;   % input initial condition
SCKS.M(2).V  = 50;   % inpud noise precison (fixed)
SCKS.M(1).xP = eye(4)*1e-2^2;   % state error covariance matrix
SCKS.M(1).uP = eye(1)*1e-2^2;   % input error covariance matrix
SCKS.M(1).wP = eye(np)*1e-10^2;  % parameter error covariance matrix
SCKS.M(1).pC = diag([1e-4; 1e-4; 1e-5; 1e-8; 1e-8; 1e-8; 1e-8]);  % covarinace matrix of paramters noise
SCKS.M(1).f  = 'spm_fx_hdm_all2';  % state equations rewriten for matrix operations
SCKS.M(1).g  = 'spm_gx_hdm_all2';  % observation equations rewriten for matrix operations
SCKS.M(1).Q  = {speye(M(1).l,M(1).l)}; % if Q is specified then algorithm performs
% estimation of measurement noise covariance
%SCKS.M(1).Q  = [];                     % if presion on measurement noise
                                        % is known then Q = [];

SCKS.M(1).Qf      = 'all';  % form of estimation of measurement noise covariance
% (after online VB estimatin); options:
% [auto,all,min,mean]
SCKS.M(1).E.nN    = 10;    % max number of iteration of SCKF-SCKS algorithm
SCKS.M(1).E.Itol  = 1e-5;  % convergence tolerance value for SCKF_SCKS algorithm
SCKS.M(1).E.RM    = [1e2 1e6];  % scaling parmater for Robbins-Monro approximation of
% parameter noise covariance [scaling
% parameter, max-limit]
SCKS.M(1).VB.N    = 5;      % max number of VB iteration during one SCKF-SCKS run
SCKS.M(1).VB.Itol = -inf;   % convergence tolerance value for VB algorithm
SCKS.M(1).VB.l    = 1 - exp(-4.5);    % scaling parameter for VB algorithm,
                                    % controls dynamics of covarince (smaller = faster)
SCKS0 = SCKS;

%--------------------------------------------------------------------------
% SCKF-SCKS Estimation:
%[SCKS] = spm_CKS(SCKS0);         % AUGMENTED STATE SCKS
% [SCKS] = spm_CKS2(SCKS0);       % NON-AUGMENTED STATE SCKS
[SCKS] = spm_CKS2new(SCKS0);      % new NON-AUGMENTED STATE SCKS


% Visualisation:
nD = SCKS.M(1,1).E.nD;
x  = SCKS.pU.x{1}(:,1:nD:end);
v  = SCKS.pU.v{2}';
r  = SCKS.pU.v{1}(:,1:(1/M(1).E.dt):end);

% Display State results:
f1 = spm_figure('Create','Graphics','SCKS estimates');
set(f1,'RendererMode','auto','Renderer','painter');
clf(f1);
for p = 1:2
    subplot(2,1,p),
    hax = gca;
    si  = spm_invNcdf(1 - 0.05);

    if p == 1,
        xxfig = SCKS.qU.x{2}(:,1:nD:end);
        Sfig  = SCKS.qU.S{2}(:,1:nD:end);
        tit   = 'Forward Estimate';
    else
        xxfig = SCKS.qU.x{1}(:,1:nD:end);
        Sfig  = SCKS.qU.S{1}(:,1:nD:end);
        tit   = 'Backward Estimate';
    end

    s = abs(Sfig);

    % conditional covariances
    %------------------------------------------------------------------
    j       = [1:size(xxfig,1)];
    ss      = si*s(j,:);
    s(j,:)  = [];
    [ill indss] = sort(mean(ss,2),'descend');
    pf = plot(1:size(xxfig,2),xxfig,'linewidth',1.5);
    set(hax,'xlim',[1,N*M(1).E.dt],'nextplot','add')
    for ic = 1:size(xxfig,1)
        col0 = get(pf(indss(ic)),'color');
        col  = (ones(1,3)-col0)*0.65 + col0;
        fill([(1:size(xxfig,2)) fliplr(1:size(xxfig,2))],[(xxfig(indss(ic),:) + ss(indss(ic),:)) fliplr((xxfig(indss(ic),:) - ss(indss(ic),:)))],...
            'r',...
            'FaceColor',col,...
            'EdgeColor',col);
        hold on
        COL{ic} = col0;
    end
    for ic = 1:size(xxfig,1)
        plot(xxfig(indss(ic),:),'color',COL{ic},'linewidth',0.75);
    end
    h1 = plot(x','Color',[0 0 0],'linewidth',1);

    title(tit);
    grid(hax,'on')
    axis(hax,'tight')
    set(hax,'box','on','Layer','top');
    set(hax,'tickdir','out')

end

% Display input  results:
f2 = spm_figure('Create','Graphics','SCKS Input estimate');
set(f2,'RendererMode','auto');
clf(f2);
for p = 1:2
    subplot(2,1,p),
    hax = gca;
    si    = spm_invNcdf(1 - 0.05);
    if p == 1,
        xxfig = SCKS.qU.v{2};
        Sfig  = SCKS.qU.C{2};
        tit   = 'Forward Estimate';
    else
        xxfig = SCKS.qU.v{1};
        Sfig  = SCKS.qU.C{1};
        tit   = 'Backward Estimate';
    end

    s = abs(Sfig);
    % conditional covariances
    %------------------------------------------------------------------
    j       = [1:size(xxfig,1)];
    ss      = si*s(j,:);
    s(j,:)  = [];

    pf = plot(1:size(xxfig,2),xxfig,'linewidth',1.5);
    set(hax,'xlim',[1,N],'nextplot','add')

    for ic = 1:size(xxfig,1)
        col0 = get(pf(ic),'color');
        col = (ones(1,3)-col0)*0.65 + col0;
        fill([(1:size(xxfig,2)) fliplr(1:size(xxfig,2))],[(xxfig(ic,:) + ss(ic,:)) fliplr((xxfig(ic,:) - ss(ic,:)))],...
            'r',...
            'FaceColor',col,...
            'EdgeColor',col);%,...
        %'FaceAlpha',0.25);
        hold on
        plot(1:size(xxfig,2),xxfig(ic,:),'color',col0,'linewidth',0.75);
    end

    h1 = plot(v','Color',[0 0 0],'linewidth',1);

    hold off
    title(tit); %drawnow;
    grid(hax,'on')
    axis(hax,'tight')
    set(hax,'box','on','Layer','top','tickdir','out');
    legend([h1(1) pf(1)],{'True','Input'});
end


% Display response results:
f3 = spm_figure('Create','Graphics','SCKS response estimate');
set(f3,'RendererMode','auto');
clf(f3);
subplot(211)
pf = plot(1:size(r,2),SCKS.qU.r{1}(1:nD:end),'r','linewidth',1.5);
set(gca,'xlim',[1,N],'nextplot','add')
pf2 = plot(1:size(r,2),SCKS.qU.z{1}(1:nD:end),'.-b','linewidth',1.5);
h1 = plot(r','Color',[0 0 0],'linewidth',1);
title(tit); %drawnow;
grid(gca,'on')
axis(gca,'tight')
set(gca,'box','on','tickdir','out');
legend([h1(1) pf(1) pf2(1)],{'True','Estimate','Residuals'});
