clear
figure(1)
clf
colormap jet

fParams = fopen('params.txt', 'r');
params = fscanf(fParams, '%lf');
%disp(params);
fclose(fParams);

sizeX = params(1);
sizeY = params(2);
nX = params(3);
nY = params(4);
P0 = params(5);
coh = params(6);
loadValue = params(7)

dX = sizeX / (nX - 1);
dY = sizeY / (nY - 1);

x = -0.5 * sizeX : dX : 0.5 * sizeX;
y = -0.5 * sizeY : dY : 0.5 * sizeY;
[x, y] = ndgrid(x, y);
sinsin = y ./ (sqrt(x.*x + y.*y));
coscos = x ./ (sqrt(x.*x + y.*y));

% DATA READING

fPres = fopen('pressure.txt', 'r');
P = fscanf(fPres, '%lf');
fclose(fPres);
P = reshape(P, nX, nY);

fTauXX = fopen('tau_xx.txt', 'r');
tauXX = fscanf(fTauXX, '%lf');
fclose(fTauXX);
tauXX = reshape(tauXX, nX, nY);

fTauYY = fopen('tau_yy.txt', 'r');
tauYY = fscanf(fTauYY, '%lf');
fclose(fTauYY);
tauYY = reshape(tauYY, nX, nY);

fPlast = fopen('plast.txt', 'r');
iPlast = fscanf(fPlast, '%lf');
fclose(fPlast);
iPlast = reshape(iPlast, nX, nY);

% ANALYTIC SOLUTION
Sanrr = zeros(nX, nY);
Sanrr(iPlast > 0) = -P0 + sign(loadValue) * 2.0 * coh * log(sqrt(sqrt(x(iPlast > 0) .* x(iPlast > 0) + y(iPlast > 0) .* y(iPlast > 0) ))) / sqrt(2);
Sanrr(sqrt(x.*x + y.*y) < 1.0) = 0.0;

Sanff = zeros(nX, nY);
Sanff(iPlast > 0) = -P0 + sign(loadValue) * 2.0 * coh * (log(sqrt(sqrt(x(iPlast > 0) .* x(iPlast > 0) + y(iPlast > 0) .* y(iPlast > 0)))) + 1.0) / sqrt(2);
Sanff(sqrt(x.*x + y.*y) < 1.0) = 0.0;

% NUMERIC SOLUTION
Snurr = zeros(nX, nY);
Snurr(iPlast > 0) = (tauXX(iPlast > 0) - P(iPlast > 0)); % .* coscos(iPlast > 0) .* coscos(iPlast > 0) + ...
                    %2.0 * tauxyAv(iPlast > 0) .* sinsin(iPlast > 0) .* coscos(iPlast > 0) + ...
                    %(tauYY(iPlast > 0) - P(iPlast > 0)) .* sinsin(iPlast > 0) .* sinsin(iPlast > 0);
Snurr(sqrt(x.*x + y.*y) < 1.0) = 0.0;

Snuff = zeros(nX, nY);
Snuff(iPlast > 0) = tauYY(iPlast > 0) - P(iPlast > 0);
%Snuff(iPlast > 0) = (tauxx(iPlast > 0) - P(iPlast > 0)) .* sinsin(iPlast > 0) .* sinsin(iPlast > 0) - ...
                    %2.0 * tauxyAv(iPlast > 0) .* sinsin(iPlast > 0) .* coscos(iPlast > 0) + ...
                    %(tauyy(iPlast > 0) - P(iPlast > 0)) .* coscos(iPlast > 0) .* coscos(iPlast > 0);
Snuff(sqrt(x.*x + y.*y) < 1.0) = 0.0;

% FIGURES
subplot(2, 1, 1)
plot(x, 0.5 * (Sanrr(:, nY/2) + Sanrr(:, nY/2 - 1)), 'g', x, 0.5 * (Snurr(:, nY/2) + Snurr(:, nY/2 - 1)), 'r')
title("\sigma_{rr}")

subplot(2, 1, 2)
plot(x, 0.5 * (Sanff(:, nY/2) + Sanff(:, nY/2 - 1)), 'g', x, 0.5 * (Snuff(:, nY/2) + Snuff(:, nY/2 - 1)), 'r')
%plot(x, 0.5 * (Sanrr(:, nY/2) - Sanff(:, nY/2)), 'g', x, 0.5 * (Snurr(:, nY/2) - Snuff(:, nY/2)), 'r')
title("\sigma_{rr} - \sigma_{\phi \phi}")

drawnow