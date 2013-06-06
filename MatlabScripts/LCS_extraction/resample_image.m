function [resampled_Y_grid,resampled_X_grid,resampled_Image] = resample_image( Image, X, Y, n )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%Compute Interpolant
[X_grid,Y_grid]=ndgrid(X,Y);
F=griddedInterpolant(X_grid,Y_grid,Image); % FIXME why cubic is not accepted
%% Resample axis
resampled_X=linspace(min(X),max(X),n);
resampled_Y=linspace(min(Y),max(Y),n);
[resampled_X_grid,resampled_Y_grid]=ndgrid(resampled_X,resampled_Y);

resampled_Image=F(resampled_X_grid,resampled_Y_grid);
end

