% SCKS_convolution
%--------------------------------------------------------------------------
% Example of inverison of linear convolution system with SCKF-SCKS algorithm
%

clear all; close all;
addpath(genpath('C:\Backup\spm8\'));
% it needs paths to SPM8 toolbox
% 
%==========================================================================
f       = spm_figure('GetWin','Graphics');
M       = spm_CKS_M('convolution model');
% SCKS needs to be able to evaluate model equations in matrix form, i.e. for all
% cubature point at once! (in order to avaid FOR loops and speed up the process)
M(1).f  = inline('[P(1,:).*x(1,:)+P(3,:).*x(2,:)+P(13,:).*v(:,:);P(2,:).*x(1,:)+P(4,:).*x(2,:)+P(14,:).*v(:,:);]','x','v','P');
M(1).g  = inline('[P(5,:).*x(1,:)+P(9,:).*x(2,:);P(6,:).*x(1,:)+P(10,:).*x(2,:);P(7,:).*x(1,:)+P(11,:).*x(2,:);P(8,:).*x(1,:)+P(12,:).*x(2,:);]','x','v','P');
M(1).V  = exp(8);                            % error precision
M(1).W  = exp(12);                           % error precision
M(2).v  = 0;
M(2).V  = exp(16);
% free parameters
%--------------------------------------------------------------------------
P       = spm_vec(M(1).pE);                                % true parameters
ip      = [2 4 5 9];   %6 9                                % free parameters
pE      = P;
np      = length(pE);
pE(ip)  = 0;
pC      = sparse(ip,ip,exp(4),np,np);
M(1).pE = pE;
M(1).pC = pC;

% generate data and invert
%==========================================================================
N         = 32;                                % length of data sequence
U         = exp(-([1:N] - 12).^2/(2.^2));      % this is the Gaussian cause
SCKS      = spm_DEM_generate(M,U,{P},{8,32},{32});
spm_DEM_qU(SCKS.pU);


% Initialization:
%--------------------------------------------------------------------------
SCKS.M(1).cb = [];   % option to specify constrain on parameters values [min max]
SCKS.M(1).ip = ip;   % indices of model parameters to be estimated
SCKS.M(2).v  = 0;    % input initial condition
SCKS.M(2).voff = 0;
SCKS.M(2).V  = 1;    % inpud noise precison (fixed)
SCKS.M(1).xP = eye(2)*1e-10^2;   % state error covariance matrix
SCKS.M(1).uP = eye(1)*1e-10^2;   % input error covariance matrix
SCKS.M(1).wP = sparse(ip,ip,10e-1^2,np,np); % parameter error covariance matrix
SCKS.M(1).pC = sparse(ip,ip,1e-12,np,np);   % covarinace matrix of paramters noise
SCKS.M(1).Q  = {speye(M(1).l,M(1).l)};      % if Q is specified then algorithm performs
                                            % estimation of measurement noise covariance 
%SCKS.M(1).Q  = [];                         % if presion on measurement noise
                                            % is known then Q = [];
SCKS.M(1).E.nN    = 50;         % max number of iteration of SCKF-SCKS algorithm
SCKS.M(1).E.Itol  = 1e-3;       % convergence tolerance value for SCKF_SCKS algorithm
SCKS.M(1).E.RM    = [1e2 1e6];  % scaling parmater for Robbins-Monro approximation of 
                                % parameter noise covariance [scaling
                                % parameter, max-limit]
SCKS.M(1).VB.N    = 5;      % max number of VB iteration during one SCKF-SCKS run
SCKS.M(1).VB.Itol = -inf;   % convergence tolerance value for VB algorithm
SCKS.M(1).VB.l    = 1 - exp(-4);   % scaling parameter for VB algorithm, 
                                   % controls dynamics

SCKS0 = SCKS;                                   
%--------------------------------------------------------------------------
% SCKF-SCKS Estimation:
%[SCKS] = spm_CKS(SCKS0);          % AUGMENTED STATE SCKS
%[SCKS] = spm_CKS2(SCKS0);        % NON-AUGMENTED STATE SCKS
%[SCKS] = spm_CKS2new(SCKS0);     % new version NON-AUGMENTED STATE SCKS
[SCKS] = spm_CKS_fl2Y(SCKS0);
%--------------------------------------------------------------------------
% Visualisation:
nD = SCKS.M(1,1).E.nD;
x  = SCKS.pU.x{1};
v  = interp1(SCKS.pU.v{2}',[1:1/nD:N]);
r  = SCKS.pU.v{1};

% Display State results:
f1 = spm_figure('Create','Graphics','CKF estimates');
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

% Display Input results:
f2 = spm_figure('Create','Graphics','CKF Input estimate');
set(f2,'RendererMode','auto');
clf(f2); 
for p = 1:2
subplot(2,1,p), 
hax = gca;
si    = spm_invNcdf(1 - 0.05);
if p == 1,
    xxfig = SCKS.qU.v{2};
    Sfig  = SCKS.qU.S{2};
    tit   = 'Forward Estimate';
else
    xxfig = SCKS.qU.v{1};
    Sfig  = SCKS.qU.S{1};
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
pf = plot(1:N,SCKS.qU.r{1}(:,1:nD:end),'r','linewidth',1.5);
set(gca,'xlim',[1,N],'nextplot','add')
pf2 = plot(1:N,SCKS.qU.z{1}(:,1:nD:end),'.-b','linewidth',1.5);
h1 = plot(r','Color',[0 0 0],'linewidth',1);
title(tit); %drawnow;
grid(gca,'on')
axis(gca,'tight')
set(gca,'box','on','tickdir','out');
legend([h1(1) pf(1) pf2(1)],{'True','Estimate','Residuals'});
