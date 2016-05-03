function [ out ] = GridData()
%GRIDDATA Summary of this function goes here
%   Detailed explanation goes here

out.xspan = [0,7];
out.yspan = [0,7];
out.Nx = 112;
out.Ny = 112;
out.CFL = 0.1;
out.Time = [0,2.8];
out.N = 500000;
end

