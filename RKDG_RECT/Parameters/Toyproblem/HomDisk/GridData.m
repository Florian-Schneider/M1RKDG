function [ out ] = GridData()
%GRIDDATA Summary of this function goes here
%   Detailed explanation goes here

out.xspan = [-5,5];
out.yspan = [-5,5];
out.Nx = round(sqrt(47310/2));
out.Ny = out.Nx;
out.CFL = 0.1;
out.Time = [0,3.75];
out.N = 500000;
end

