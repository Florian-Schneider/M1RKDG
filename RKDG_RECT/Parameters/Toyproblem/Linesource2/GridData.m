function [ out ] = GridData()
%GRIDDATA Summary of this function goes here
%   Detailed explanation goes here

out.xspan = [-0.5,0.5];
out.yspan = [-0.5,0.5];
out.Nx = 200;
out.Ny = 200;
out.CFL = 0.1;
out.Time = [0,0.45];
out.N = 500000;
end

