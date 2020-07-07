function [f] = spm_fx_hdm_allU(x,u,P,M)
% state equation for the hemodynamic model
% FORMAT [f] = spm_fx_hdm(x,u,P,M)
% x      - state vector
%   x(1) - vascular signal                                    s
%   x(2) - rCBF                                           log(f)
%   x(3) - venous volume                                  log(v)
%   x(4) - dHb                                            log(q)
% u      - input (neuronal activity)                      (u)
% P      - free parameter vector
%   P(1) - signal decay                                   d(ds/dt)/ds)
%   P(2) - autoregulation                                 d(ds/dt)/df)
%   P(3) - transit time                                   (t0)
%   P(4) - exponent for Fout(v)                           (alpha)
%   P(5) - resting oxygen extraction                      (E0)
%   P(6) - ratio of intra- to extra-vascular components   (epsilon)
%          of the gradient echo signal   
%
%   P(6 + 1:m)   - input efficacies                       d(ds/dt)/du)
%
% y      - dx/dt
%__________________________________________________________________________
%
% Ref Buxton RB, Wong EC & Frank LR. Dynamics of blood flow and oxygenation
% changes during brain activation: The Balloon model. MRM 39:855-864 (1998)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_fx_hdm.m 2495 2008-11-27 12:18:33Z karl $

% exponentiation of hemodynamic state variables
%--------------------------------------------------------------------------

% neuronal equation:

 f(1,:) = -P(7,:).*x(1,:) + 0.5*exp(P(8,:)).*repmat(u,1,size(x,2));

% hemodynamic states
x(3:end,:) = exp(x(3:end,:)); 

% Fout = f(v) - outflow
%--------------------------------------------------------------------------
fv(:,:)       = x(4,:).^(1./P(4,:));

% e = f(f) - oxygen extraction
%--------------------------------------------------------------------------
ff(:,:)       = (1 - (1 - P(5,:)).^(1./x(3,:)))./P(5,:);

% implement differential state equations
%--------------------------------------------------------------------------
f(2,:)     = x(1,:) - P(1,:).*x(2,:) - P(2,:).*(x(3,:) - 1);
f(3,:)     = x(2,:)./x(3,:);
f(4,:)     = (x(3,:) - fv)./(P(3,:).*x(4,:));
f(5,:)     = (ff.*x(3,:) - fv.*x(5,:)./x(4,:))./(P(3,:).*x(5,:));
%f        = f(:);
% f(2:end,:)=exp(f(2:end,:));
% f(3:end,:) = f(3:end,:)+1;
% adjust motion for DEM (that uses time-bins as units of time)
%--------------------------------------------------------------------------
try, global dt, f  = f*dt; end

return