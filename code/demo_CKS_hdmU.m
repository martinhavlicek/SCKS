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

% level 1
%------------------------------------------------------------------
[pE pC] = spm_hdm_priors(1,3);              % parameters
pE(6) = 0.02;
pE(7) = 0.54;
M(1).n  = 4;
M(1).f  = 'spm_fx_hdm';
M(1).g  = 'spm_gx_hdm';
M(1).pE = pE;                               % prior expectation
M(1).V  = exp(0);                                % error precision
M(1).W  = exp(12);                          % error precision
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
N         = 60;                                % length of data sequence
U         = exp(-([1:11] - 6).^2/(2.^2));      % this is the Gaussian cause
input     = zeros(1,N);
input([5 10 20 34 43 50]) = [1 0.8 1 0.2 .9 0.4];
U         = conv(U,input);
U         = U(1:N);
V         = 2;   %!!!! precision on observation noise
SCKS      = spm_DEM_generate(M,U,{P},{V,6},{6});
spm_DEM_qU(SCKS.pU)


pE([1 2 3])  = P([1 2 3]) + 1/12*randn(3,1);  % randomize first 3 parameters



% MODEL INVERSION SPECIFICATION:
%================================

% Specification of model parameters and its covariances: 
%--------------------------------------------------------------------------
np = 8;                   % number of model parameters
SCKS.M(1).pE = [pE;0];    % initialization values of model parameters (8th one is additional for controling exogenous input U)
SCKS.M(1).ip = [1,2,3,8];  % indices of model parameters to be estimated (ignoring 4-6, because these can not be well estimated)
SCKS.M(1).wP = diag([1e-4; 1e-4; 1e-5; 1e-6; 1e-6; 1e-6; 1e-6; 1e-4;]).^2;  % parameter error covariance matrix
SCKS.M(1).pC = diag([1e-3; 1e-4; 1e-6; 1e-4; 1e-6; 1e-4; 1e-6; 1e-3;]);     % for parameter noise covarince
                                                                            % specified are all 8, but used are only those in "ip" vector

SCKS.M(1).E.RMw  = [1e2 1e8];  % scaling parmater for Robbins-Monro approximation of
                                 % parameter noise covariance [scaling
                                 % parameter, max-limit]                                                                            
SCKS.M(1).cb = [];               % option to specify constrain on parameters values [min max]
                                                                            


% Specification of model states and its covariances: 
%--------------------------------------------------------------------------
SCKS.M(1).n  = 5;     % number of hidden states
SCKS.M(1).x  = zeros(SCKS.M(1).n ,1);   % initial conditions
SCKS.M(1).xP = diag([1e-2^2, 1e-4^2,1e-4^2,1e-4^2,1e-4^2]);   % state error covariance matrix
W            = exp(spm_vec(sparse(1,1,(6 - 16),1,5) + 16));
SCKS.M(1).W  = diag(W)*1e-1;%    % noise covariance for hidden states (the first value is for the neuronal state)
SCKS.M(1).E.RMq  = [5e2 1e8];     % similarly as for parameter noise covariance (Robbins-Monro)


% Specification of exogenous input:
%--------------------------------------------------------------------------
SCKS.M(1).U  = [0 0 0 input(1:end-3)]; % or [] % if no known input


% Specifiaction of measurement noise estimation scheme:
%----------------------------------------------------------------
SCKS.M(1).V  = SCKS.M(1).V + exp(V); % this is for case that we assume the noise covariace to be known  
SCKS.M(1).Q  = {speye(M(1).l,M(1).l)}; % if Q is specified (in this way) then algorithm performs
                                       % estimation of measurement noise covariance
%SCKS.M(1).Q  = [];                    % if presion on measurement noise
                                       % is known then Q = [];

SCKS.M(1).VB.N    = 5;              % max number of VB iteration during one SCKF-SCKS run
SCKS.M(1).VB.Itol = -inf;           % convergence tolerance value for VB algorithm
SCKS.M(1).VB.l    = 1 - exp(-6);    % scaling parameter for VB algorithm,
                                    % controls dynamics of covarince
                                    % (smaller = faster)


                                    
% Specification of input level (NOT used here!, neuronal level is part of state equations)
%--------------------------------------------------------------------------
SCKS.M(2).v  = 0;   % input initial condition %
SCKS.M(2).V  = 1;   % noise varinace of the input  
SCKS.M(1).uP = eye(1)*1e-5^2;   % input error covariance matrix


% general setting
%--------------------------------------------------------------------------
SCKS.M(1,1).E.nN = 15; % number of iterations
SCKS.M(1,1).E.nD = 3;  % number of integration steps between regular sampling (e.g. step size TR/3)
SCKS.M(1).E.Itol = 1e-3;  % convergence tolerance value for SCKF_SCKS algorithm
SCKS.M(1).f  = 'spm_fx_hdm_all2U';  % state equations rewriten for matrix operations
SCKS.M(1).g  = 'spm_gx_hdm_all2U';  % observation equations rewriten for matrix operations



SCKS0 = SCKS;

%--------------------------------------------------------------------------
% SCKF-SCKS Estimation:

[SCKS]  = spm_CKS2newU(SCKS0);      % NON-AUGMENTED STATE SCKS

% Visualisation:
nD = SCKS.M(1,1).E.nD;
x  = SCKS.pU.x{1};
v  = interp1(SCKS.pU.v{2}',[1:1/nD:N]);
r  = SCKS.pU.v{1};

% Display State results:
f1 = spm_figure('Create','Graphics','SCKS estimates');
set(f1,'RendererMode','auto','Renderer','painter');
clf(f1);
for p = 1:2
    subplot(2,1,p),
    hax = gca;
    si    = spm_invNcdf(1 - 0.05);

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
    pf = plot(1:N,xxfig,'linewidth',1.5);
    set(hax,'xlim',[1,N],'nextplot','add')
    for ic = 1:size(xxfig,1)
        col0 = get(pf(indss(ic)),'color');
        col = (ones(1,3)-col0)*0.65 + col0;
        fill([(1:N) fliplr(1:N)],[(xxfig(indss(ic),:) + ss(indss(ic),:)) fliplr((xxfig(indss(ic),:) - ss(indss(ic),:)))],...
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
    legend([h1(1) pf(1)],{'True','x1'});
end


% Display response results:
f3 = spm_figure('Create','Graphics','SCKS response estimate');
set(f3,'RendererMode','auto');
clf(f3);
subplot(211)
pf = plot(1:N,SCKS.qU.r{1}(1:nD:end),'r','linewidth',1.5);
set(gca,'xlim',[1,N],'nextplot','add')
pf2 = plot(1:N,SCKS.qU.z{1}(1:nD:end),'.-b','linewidth',1.5);
h1 = plot(r','Color',[0 0 0],'linewidth',1);
title(tit); %drawnow;
grid(gca,'on')
axis(gca,'tight')
set(gca,'box','on','tickdir','out');
legend([h1(1) pf(1) pf2(1)],{'True','Estimate','Residuals'});
