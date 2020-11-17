%% Michael Shadlen bimodal failure distribution

NumTrials = 1e4;
% RANDOM    generic MATLAB function for many distributions
ITIs1 = random('Rayleigh', 1/6, 1, NumTrials) + 0.1;
ITIs2 = random('Rayleigh', 1/sqrt(30), 1, NumTrials) + 1.75;
rr = round( rand(1,NumTrials) );
ITIs = rr .* ITIs1 + (1 - rr) .* ITIs2;

figure;
hist(ITIs, 50);

xlabel('ITI Distribution [Rayleigh]');
ylabel('Counts');
title('Bimodal Failure Distribution [Rayleigh]');
set( gca, 'LineWidth', 2);
set( gca, 'fontname', 'Te X Gyre Heros'); % due to Linux compatability issue with Helvetica font

%% Hangya Balazs failure distribution [bimodal, Normal]

% define parameters 
NumTrials = 1e4;
ITIMin = 0.1;
ITIMax = 3;
% parameters for the Gaussians
mng1 = 0.3;
mng2 = 2;
sdg = 0.15;
% mixing probabilities
pmx1 = 0.35;
pmx2 = 0.35;
pmx3 = 1 - pmx1 - pmx2;

% use parameters to create distributions
ITIs1 = random('Normal', mng1, sdg, 1, NumTrials);
while any(ITIs1>ITIMax) || any(ITIs1<ITIMin)
    inx = ITIs1 > ITIMax | ITIs1 < ITIMin;
    ITIs1(inx) = random('Normal', mng1, sdg, 1, sum(inx) );
end

ITIs2 = random('Normal', mng2, sdg, 1, NumTrials);
while any(ITIs2>ITIMax) || any(ITIs2<ITIMin)
    inx = ITIs2 > ITIMax | ITIs2 < ITIMin;
    ITIs2(inx) = random('Normal', mng2, sdg, 1, sum(inx) );
end

ITIs3 = random('Uniform', ITIMin, ITIMax, 1, NumTrials);
prr = rand(1, NumTrials);
rr  = zeros(3, NumTrials);
rr(1, prr<pmx1) = 1;
rr(2, prr>=pmx1&prr<(pmx1+pmx2)) = 1;
rr(3, prr>=(pmx1+pmx2)) = 1;
ITIs = rr(1,:) .* ITIs1 + rr(2,:) .* ITIs2 + rr(3,:) .* ITIs3;
ITIMin = min(ITIs);
ITIMax = max(ITIs);

figure;
hist(ITIs, 50);

xlabel('ITI Distribution [Normal]');
ylabel('Counts');
title('Bimodal Failure Distribution [Normal]');
set( gca, 'LineWidth', 2);
set( gca, 'fontname', 'Te X Gyre Heros'); % due to Linux compatability issue with Helvetica font

%% Hangya Balazs failure distribution [unimodal, Normal]

% define parameters
NumTrials = 1e4;
ITIMin = 0.1;
ITIMax = 3;
% parameters for the Gaussians
mng = 1.4;
sdg = 0.25;
% mixing probabilities
pmx1 = 0.65;
pmx2 = 1 - pmx1; 

% make distribution
ITIs1 = random('Normal', mng, sdg, 1, NumTrials);
while any(ITIs1>ITIMax) || any(ITIs1<ITIMin)
    inx = ITIs1 > ITIMax | ITIs1 < ITIMin;
    ITIs1(inx) = random('Normal', mng, sdg, 1, sum(inx));
end

ITIs2 = random('Uniform', ITIMin, ITIMax, 1, NumTrials);
prr = rand(1, NumTrials);
rr = zeros(2, NumTrials);
rr(1,prr<pmx1) = 1;
rr(2,prr>=pmx1) = 1;

ITIs = rr(1,:) .* ITIs1 + rr(2,:) .* ITIs2;
ITIMin = min(ITIs);
ITIMax = max(ITIs);

figure;
hist(ITIs, 50);

xlabel('ITI Distribution [Normal]');
ylabel('Counts');
title('Unimodal Failure Distribution [Normal]');
set( gca, 'LineWidth', 2);
set( gca, 'fontname', 'Te X Gyre Heros'); % due to Linux compatability issue with Helvetica font

%% Hangya Balazs failure distribution [bimodal, Exponential]

% define parameters
ITIMean = 1.4;
ITIMin = 0.1;
ITIMax = 5;
NumTrials = 1e4;

% make distribution
ITIMean2 = ITIMean - ITIMin;
ITIMax2 = ITIMax - ITIMin;

temp = exprnd(ITIMean2,1,NumTrials);
while any(temp>ITIMax2)
    inx = temp > ITIMax2;
    temp(inx) = exprnd(ITIMean2,1, sum(inx));
end

temp = temp + ITIMin;
ITIs = temp;

%% failure density function

dt = 0.05;
times = ITIMin-dt : dt : ITIMax+dt;
cnts = (times(1:end-1) + times(2:end)) / 2;
ft = histc(ITIs, times);
ft = ft(1:end-1);
ft = ft / sum(ft);
figure;
bar(cnts, ft);

xlabel('Counts');
ylabel('Failure Rate');
title('Cumulative Failure Density');
set( gca, 'LineWidth', 1);
set( gca, 'fontname', 'Te X Gyre Heros'); % due to Linux compatability issue with Helvetica font

%% hazard rate

Ft = cumsum(ft);
Rt = 1 - Ft;
ht = (Rt(1:end-1) - Rt(2:end)) ./ (dt * Rt(1:end-1));
figure;
plot(cnts(1:end-1), ht, 'LineWidth', 2);

xlabel('Counts');
ylabel('Hazard Rate');
title('Cumulative Hazard');
set( gca, 'LineWidth', 1);
set( gca, 'fontname', 'Te X Gyre Heros'); % due to Linux compatability issue with Helvetica font

%% subjective hazard rate

% convolution of the failure density function with a time-dependent
% Gaussian
phi = 0.26; % Shadlen paper [2005]
intfg = zeros(size(cnts));
tau = cnts;
cntr = 1; % counter
cnts2 = [cnts cnts(end)+dt : dt : cnts(end)+5e3*dt];
for t = cnts2
    fg = ft .* exp(-(tau-t).^2 / (2*phi^2*t^2));
    intfg(cntr) = sum(fg*dt);
    cntr = cntr + 1;
end

fwt = 1 ./ (phi * cnts2 * sqrt(2*pi)) .* intfg;
fwt = fwt / sum(fwt);

Fwt = cumsum(fwt);
Rwt = 1 - Fwt; 
hwt = (Rwt(1:end-1) - Rwt(2:end)) ./ (dt * Rwt(1:end-1));
figure; 
plot(cnts2(1:end-1), hwt, 'LineWidth', 2);
xlim([0 5]);

xlabel('Counts');
ylabel('Subjective Hazard Rate');
title('Subjective Hazard');
set( gca, 'LineWidth', 1);
set( gca, 'fontname', 'Te X Gyre Heros'); % due to Linux compatability issue with Helvetica font

































