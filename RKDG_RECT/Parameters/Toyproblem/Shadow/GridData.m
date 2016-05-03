function [ out ] = GridData()
%GRIDDATA Summary of this function goes here
%   Detailed explanation goes here

out.xspan = [0,12];
out.yspan = [0,6];
out.Nx = 2*round(sqrt(29931/2/2));
out.Ny = round(sqrt(29931/2/2));
out.CFL = 0.1;
out.Time = [0,40];
out.N = 500000;
end

