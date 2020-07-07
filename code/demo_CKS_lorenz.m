% SCKS_lorenz
%--------------------------------------------------------------------------
% Example of inverison of lorenz system with SCKF-SCKS algorithm
%
clear all; close all;
addpath(genpath('C:\Backup\spm8\'));
% it needs paths to SPM8 toolbox
% 
%==========================================================================
f       = spm_figure('GetWin','Graphics');

% get model
%--------------------------------------------------------------------------
M       = spm_CKS_M('Lorenz');

% create data
%==========================================================================
% create innovations & add causes
%--------------------------------------------------------------------------
N       = 120;
U       = sparse(1,N);
%--------------------------------------------------------------------------
M(1).V  = exp(0);
M(1).W  = exp(6);

% level 2 precisions
%--------------------------------------------------------------------------
M(2).v  = 0;
M(2).V  = exp(16);

SCKS     = spm_DEM_generate(M,U);
SNR      = 5;
Y0          = SCKS.Y(:,:)-SCKS.pU.z{1}(:,:);
stdY        = mean(std(Y0,[],2));
stdZ        = mean(std(SCKS.pU.z{1}(:,:),[],2));
SCKS.pU.z{1} = SCKS.pU.z{1}*((stdY/SNR)/stdZ);
SCKS.Y       = Y0 + SCKS.pU.z{1};
SCKS.pU.r{1} = SCKS.Y;

disp(mean(std(SCKS.pU.z{1}(:,:),[],2)))

spm_DEM_qU(SCKS.pU)

% Initialization:
%--------------------------------------------------------------------------
SCKS.M(1).x = [rand(2,1)*10;25]; % partly random initialization of the states
%SCKS.M(1).x = [0.9;0.8;30];     % True state values: 

pE = [2;-12;46.92;];   % initialization of parameters (True values: pE =[18;-4;46.92;])
%pE =[18;-4;46.92;];
SCKS.M(1).pE = pE;
np = length(pE);
ip = 1:np;             

% option to specify constrains on paramters (example, not used further)
cb = [0, 30;
    -10, 10;
     40, 50;];   % [min max]

SCKS.M(1).cb = [];    % option to specify constrain on parameters values [min max]
SCKS.M(1).ip = ip;    % indices of model parameters to be estimated
SCKS.M(2).v  = [];    % input initial condition
SCKS.M(2).voff  = 1;  % to switch off input
SCKS.M(2).V  = [];    % inpud noise precison (fixed)
SCKS.M(1).xP = eye(3)*1e-1^2;     % state error covariance matrix
SCKS.M(1).uP = eye(1)*1e-1^2;    % input error covariance matrix

SCKS.M(1).wP = eye(np)*5e-1^2;    % parameter error covariance matrix
SCKS.M(1).pC = diag([1e-1; 1e-1; 1e-3;]);  % covarinace matrix of paramters noise
SCKS.M(1).Q  = {speye(M(1).l,M(1).l)};     % if Q is specified like this then algorithm performs
                                           % estimation of measurement
                                           % noise covariance  via
                                           % Variational Bayesian method
%SCKS.M(1).Q  = [];                        % if presion on measurement noise
                                           % is known then Q = [];
SCKS.M(1,1).E.nD  = 1;
SCKS.M(1).E.nN    = 30;         % max number of iteration of SCKF-SCKS algorithm
SCKS.M(1).E.Itol  = 1e-6;       % convergence tolerance value for SCKF_SCKS algorithm
SCKS.M(1).E.RM    = [1e3 1e6];  % scaling parmater for Robbins-Monro approximation of 
                                % parameter noise covariance [scaling parameter, max-limit]
SCKS.M(1).VB.N    = 5;          % max number of VB iteration during one SCKF-SCKS run
SCKS.M(1).VB.Itol = -inf;
SCKS.M(1).VB.l    = 1 - exp(-5);    % scaling parameter for VB algorithm, 
                                    % controls dynamics

SCKS0 = SCKS;

%--------------------------------------------------------------------------
% SCKF-SCKS Estimation:
%[SCKS]  = spm_CKSdec(SCKS0);     % AUGMENTED STATE SCKS
[SCKS] = spm_CKS2(SCKS0);     % NON-AUGMENTED STATE SCKS
%[SCKS] = spm_CKS2newU(SCKS0);   % new version NON-AUGMENTED STATE SCKS
%--------------------------------------------------------------------------
% Visualisation:
nD = SCKS.M(1,1).E.nD;
x  = SCKS.pU.x{1};
v  = interp1(SCKS.pU.v{2}',[1:1/nD:N]);
r  = SCKS.pU.v{1};

% Display State results:
f1 = figure;%spm_figure('Create','Graphics','CKF estimates');
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
    legend([h1(1) pf(1) pf(2) pf(3)],{'True','x1','x2','x3',});
end
% Display response results:
f3 = figure;%spm_figure('Create','Graphics','SCKS response estimate');
set(f3,'RendererMode','auto');
clf(f3);
subplot(211)


pf1 = plot(1:N,r' - SCKS.qU.r{1}(1:nD:end)','-b','linewidth',1);
set(gca,'xlim',[1,N],'nextplot','add')
h1 = plot(SCKS.Y,'Color',[0 0 0],'linewidth',1);
pf2 = plot(1:N,SCKS.Y-SCKS.pU.z{1}(1:nD:end),'k','linewidth',2);
pf3 = plot(1:N,SCKS.qU.r{1}(1:nD:end),'r','linewidth',1);

title(tit); %drawnow;
grid(gca,'on')
axis(gca,'tight')
set(gca,'box','on','tickdir','out');
%legend([h1(1) pf1(1) pf3(1)],{'True','Estimate','Residuals'});



% Display State results:
f3 = figure;%spm_figure('Create','Graphics','CKF estimates');
set(f3,'RendererMode','auto','Renderer','painter');
clf(f3);
for p = 1
    subplot(2,1,p),
    hax = gca;
    si    = spm_invNcdf(1 - 0.05);


    xxfig = SCKS.qU.ww{1}(:,1:nD:end);
    Sfig  = SCKS.qU.ss{1}(:,1:nD:end);
 

    s = abs(Sfig);
    N = size(Sfig,2);
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
    %h1 = plot(x','Color',[0 0 0],'linewidth',1);

    title(tit);
    grid(hax,'on')
    axis(hax,'tight')
    set(hax,'box','on','Layer','top');
    set(hax,'tickdir','out')
end


f4 = figure;%spm_figure('Create','Graphics','SCKS response estimate');
set(f4,'RendererMode','auto');
clf(f4);

pf1 = plot(1:50,SCKS.F,'k','linewidth',1);
grid(gca,'on')
axis(gca,'tight')
set(gca,'box','on','tickdir','out');

f5 = figure;%spm_figure('Create','Graphics','SCKS response estimate');
set(f5,'RendererMode','auto');
clf(f5);

pf1 = plot(1:50,diff([SCKS.F,SCKS.F(end)]),'k','linewidth',1);
grid(gca,'on')
axis(gca,'tight')
set(gca,'box','on','tickdir','out');