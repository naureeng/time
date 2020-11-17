%% failure distribution 

ITIMean = 1.5;
ITIMin  = 0.1;
ITIMax  = 5;
NumTrials = 1e4;

ITIMean2 = ITIMean - ITIMin;
ITIMax2  = ITIMax  - ITIMin;

temp = exprnd( ITIMean2, NumTrials, 1);
while any(temp>ITIMax2)
    inx = temp > ITIMax2;
    temp(inx) = exprnd(ITIMean2, sum(inx), 1)';
end

temp = temp + ITIMin;
ITIs = temp;

%% failure density function [pdf]

dt = 0.1;
times = ITIMin-dt : dt : ITIMax+dt;
cnts  = (times(1:end-1) + times(2:end)) / 2;
ft = histcounts( ITIs, times );
ft = ft(1:end-1);
ft = ft/ sum(ft); % sum normalization
figure; 
bar(cnts(1:end-1), ft);

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
plot(cnts(1:end-2), ht, 'LineWidth', 2);

xlabel('Counts');
ylabel('Hazard Rate');
title('Cumulative Hazard');
set( gca, 'LineWidth', 1);
set( gca, 'fontname', 'Te X Gyre Heros'); % due to Linux compatability issue with Helvetica font

%% subjective hazard rate
% convolution of the failure density function with a time-dependent
% Gaussian

phi = 0.26; % value that works best in Shadlen, 2005 paper
intfg = zeros(size(cnts));
tau = cnts;
cntr = 1; % counter
cnts2 = [cnts(1)-95*dt: dt : cnts(1)-dt cnts cnts(end)+dt: dt : cnts(end)+195*dt];

for t = cnts2
    fg = ft' .* exp(-(tau-t).^2/(2*phi^2*t^2));
    fg = ft( : , 1 : end-1 );
    intfg(cntr) = sum(fg*dt);
    cntr = cntr + 1;
end

fwt = 1 ./ (phi * cnts2 * sqrt(2*pi)) .* intfg;
fwt = fwt/ sum(fwt);

% compute subjective hazard rate
Fwt = cumsum(fwt);
Rwt = 1 - Fwt;
hwt = (Rwt(1:end-1) - Rwt(2:end)) ./ (dt * Rwt(1:end-1));
figure;
plot(cnts2(1:end-1), hwt, 'LineWidth', 2);

xlabel('Counts');
ylabel('Subjective Hazard Rate');
title('Subjective Hazard');
set( gca, 'LineWidth', 1);
set( gca, 'fontname', 'Te X Gyre Heros'); % due to Linux compatability issue with Helvetica font













