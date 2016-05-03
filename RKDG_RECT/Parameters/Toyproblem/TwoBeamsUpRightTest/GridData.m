function [ out ] = GridData()
%GRIDDATA Summary of this function goes here
%   Detailed explanation goes here

x = 1/7;

out.xspan = 3*[-3.5,3.5]*x;
out.yspan = 3*[-3.5,3.5]*x;
out.Nx = 20;
out.Ny = 20;
out.CFL = 0.1;
out.Time = [0,7]*x;
out.N = 500000;
end

