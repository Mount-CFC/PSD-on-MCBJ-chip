function [fitresult, gof] = createFit(X_mincor, Y_mincor, N)
%CREATEFIT(X_MINCOR,Y_MINCOR,N)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : X_mincor
%      Y Input : Y_mincor
%      Z Output: N
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 02-Jul-2019 20:03:40 自动生成


%% Fit: 'untitled fit 1'.
[xData, yData, zData] = prepareSurfaceData( X_mincor, Y_mincor, N );

% Set up fittype and options.
ft = fittype( 'exp(a*x.^2 + x.*y.*b + x.*c + y.*d + e*y.^2 + f)', 'independent', {'x', 'y'}, 'dependent', 'z' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.647617630172684 0.679016754093202 0.635786710514084 0.945174113109401 0.208934922426023 0.709281702710545];

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, [xData, yData], zData );
% legend( h, 'untitled fit 1', 'N vs. X_mincor, Y_mincor', 'Location', 'NorthEast' );
% Label axes
% xlabel X_mincor
% ylabel Y_mincor
% zlabel N
% grid on


