function [ out ] = GridData()
%GRIDDATA Summary of this function goes here
%   Detailed explanation goes here

out.xspan = [-5,5];
out.yspan = [-5,5];
out.Nx = 100;
out.Ny = 100;
out.CFL = 0.1;
out.Time = [0,0.45];
out.N = 500000;
end

