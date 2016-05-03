function [ out ] = GridData()
%GRIDDATA Summary of this function goes here
%   Detailed explanation goes here

out.xspan = [0,1];
out.yspan = [0,1];
out.Nx = 20;
out.Ny = 20;
out.CFL = 0.1;
out.Time = [0,1.2];
out.N = 500000;
end

