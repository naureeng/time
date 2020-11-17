function y = bimodal_hazard_model(param,x)
%% BIMODAL_HAZARD_MODEL     Bimodal subjective hazard rate model
%       RT = BIMODAL_HAZARD_MODEL(PARAM, ITI) calculates model results (RT,
%       reaction time) at given points (ITI, foreperiod) based on the
%       following subjective hazard rate model:
%               RT(t) = we + wb * Ab(t)
%       where 
%               we    = constant parameter
%               wb    = weight of bimodal anticipation function
%               Ab    = biomodal anticipation function
%                     = a mixture of two Gaussians (means = 0.3 and 2s; 
%                       SD = 0.15s) and a uniform random variable
%                       constrained between 0.1 and 3 seconds
%                       (mixing proportions: 0.35, 0.35, 0.3)
%       Subjective hazard rate is computed from the hazard rate by
%       convolving it with a Gaussian pdf with variability linearly
%       increasing in time (phi = slope parameter of this increase)
%
%  REFERENCE:
%  Janssen P, Shadlen MN (2005) A representation of the hazard rate of
%  elapsed time in macaque area LIP. Nat Neurosci 8(2):234-41.
%       

%% handle both vector and struct input
if ~isstruct(param)
    st.we  = param(1);
    st.wb  = param(2);
    st.phi = param(3);
    param = st;
end

%% call the model
y = param.we + param.wb * subjective_hazard_rate(param.phi, x);

% -------------------------------------------------------------------------

function v = subjective_hazard_rate( phi, ITI )

%% paramaters of the foreperiod distribution

%  for the uniform distribution:
ITIMin = 0.1; 
ITIMax = 3;

%  for the Gaussians:
mng1 = 0.3;
mng2 = 2;
sdg  = 0.15;
%  mixing probabilities:
pmx1 = 0.35; 
pmx2 = 0.35;
pmx3 = 1 - pmx1 - pmx2 ;

%  time
dt = 0.01; % resolution
t  = ITIMin: dt : ITIMax; % time

% failure density function
ft = pmx3 / (ITIMax - ITIMin) + ...
    pmx1 * normpdf( t, mng1, sdg ) + ...
    pmx2 * normpdf( t, mng2, sdg ); 

% sum normalization (necessary for estimating cdf)
ft = ft / sum(ft);

% plot
figure;
plot(t, ft, 'Linewidth', 2);
xlabel(' Time of GO signal [sec] ' );
ylabel(' f(t) ');
title( ' Probability Density Function ' );
set( gca, 'LineWidth', 1);
set( gca, 'fontname', 'Te X Gyre Heros'); % due to Linux compatability issue with Helvetica font


% hazard rate
Ft = cumsum(ft);
Rt = 1 - Ft;
ht = ( Rt(1:end-1) - Rt(2:end) ) ./ (dt * Rt(1:end-1));

% plot
figure;
plot(t(1:end-1), ht, 'Linewidth', 2);
xlabel(' Time of GO signal [sec] ' );
ylabel(' h(t) ');
title( ' Hazard Rate ' );
set( gca, 'LineWidth', 1);
set( gca, 'fontname', 'Te X Gyre Heros'); % due to Linux compatability issue with Helvetica font

% subjective hazard rate (Fwt) = 
% convolution of the failure density function with a time-dependent
% Gaussian
intfg = zeros(size(t));
tau = t;
cntr = 1; % counter

% we need a longer time range to better estimate Fwt
t2 = [fliplr(t(1)-dt: -dt: dt) t t(end)+dt: dt: 10*t(end)]; 

for tm = t2
    fg = ft .* exp(-(tau-tm).^2/(2*phi^2*tm^2));
    intfg(cntr) = sum( fg * dt );
    cntr = cntr + 1;
end

fwt = 1./ (phi * t2 * sqrt(2*pi) ) .* intfg;
fwt = fwt / sum(fwt);

% subjective hazard rate
Fwt = cumsum(fwt);
Rwt = 1 - Fwt;
hwt = (Rwt(1:end-1) - Rwt(2:end)) ./ (dt * Rwt(1:end-1));

% plot
figure;
plot( t2(1:end-1), hwt, 'Linewidth', 2);
xlabel(' Time of GO signal [sec] ' );
ylabel(' subjective h(t) ');
title( ' Subjective Hazard Rate (Anticipation Function) ' );
set( gca, 'LineWidth', 1);
set( gca, 'fontname', 'Te X Gyre Heros'); % due to Linux compatability issue with Helvetica font

% values at requested times
v = interp1( t2(1:end-1), hwt, ITI );

end

end
















