function [ out ] = GridData()
%GRIDDATA Summary of this function goes here
%   Detailed explanation goes here

out.xspan = [-10,10];
out.yspan = [-10,10];
out.Nx = round(sqrt(67548/2));
out.Ny = out.Nx;
out.CFL = 0.1;
out.Time = [0,6];
out.N = 500000;
end

