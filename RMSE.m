function [rmse] = RMSE(y1,y2)
% [rmse] = RMSE(y1,y2) computes the root-mean square error (rmse) between
% two time series y1 and y2
%
% Input:
% y1: [1xN] or [Nx1] double
% y1: [1xN] or [Nx1] double
%
% Output
% rmse: [1x1] double: rmse value betwen y1 and y2
% 
% Author: E. Cheynet - UiB - 03-08-2020

rmse = sqrt(nanmean((y1(:)-y2(:)).^2));

end