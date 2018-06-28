function [RMSE R2] = goodnessFit(y,f,varargin)

% y is the measure
% f is the fit model

RMSE = sqrt(mean((y-f).^2));

% SStot = 0;
% SSres = 0;
% for i = 1:length(y)
%     SStot = SStot+(y(i)-mean(y))^2;
%     SSres = SSres+(y(i)-f(i))^2;
% end
% R2 = 1-SSres/SStot;
R2 = corr2(y,f)^2;