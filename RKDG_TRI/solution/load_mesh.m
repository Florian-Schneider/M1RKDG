function [p,t] = load_mesh(filename)
fid = fopen(filename);
n_p = fscanf(fid,'%e',1);         % number of points
p = fscanf(fid,'%e %e',[2 n_p])'; % coordinates of points
n_t = fscanf(fid,'%e',1);         % number of triangles
t = fscanf(fid,'%e %e',[3 n_t])'; % edge indices of triangles
fclose(fid);
