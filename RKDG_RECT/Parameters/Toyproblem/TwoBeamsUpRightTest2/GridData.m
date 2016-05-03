function [ out ] = GridData()
%GRIDDATA Summary of this function goes here
%   Detailed explanation goes here

x = 1;

out.xspan = [0,7]*x;
out.yspan = [0,7]*x;
out.Nx = 40;
out.Ny = 40;
out.CFL = 0.1;
out.Time = [0,7]*x;
out.N = 500000;
end

