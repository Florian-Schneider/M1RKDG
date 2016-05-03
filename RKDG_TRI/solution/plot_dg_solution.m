function vax = plot_dg_solution(filename,parameters)
%PLOT_DG_SOLUTION    Visualize true DG solution.
%   The true solution (given in polynomial order up to 2) is evaluated
%   at 6 points in each triangluar elements, and it is plotted in a
%   refined triangulation (each element split into four sub-triangles).
%   If an output argument is given, then the solution is not
%   plotted, but just its range computed: vax = [min value,max value].

% (C) 2012/03/14 by Benjamin Seibold

%-----------------------------------------------------------------------
% Parameters
%-----------------------------------------------------------------------
try parameters.moment_nr; catch, parameters.moment_nr = 0; end
try parameters.plot_style; catch, parameters.plot_style = 2; end
try parameters.flag_plotmesh; catch, parameters.flag_plotmesh = 0; end
try parameters.flag_average; catch, parameters.flag_average = 1; end
%try parameters.rescale = @rescale_log; end
%-----------------------------------------------------------------------

% Load solution
fid = fopen(filename);
filename_mesh = fscanf(fid,'%s',1); % filename of mesh file
time = fscanf(fid,'%e',1); % current time
sol = fscanf(fid,'%d %d %e',[3 inf])'; % coefficients of solution
fclose(fid);

% Load mesh

% Load mesh
filename_mesh = strcat('../meshes/',filename_mesh);
[p,t] = load_mesh(filename_mesh);


% Obtain solution for one specific moment
sol = sol(sol(:,2)==parameters.moment_nr,:);

% Check for NaN entries
nanind = isnan(sol(:,3));
if any(nanind)
    nanelements = unique(sol(nanind,1));
    fprintf('Warning: NaN values in element')
    if length(nanelements)==1
        fprintf(': %d.\n',nanelements)
    else
        fprintf('s: %d',nanelements(1))
        fprintf(',%d',nanelements(2:end))
        fprintf('.\n')
    end
end

n_fct = sum(sol(:,1)==sol(1,1)); % number of basis functions

% Compute solution on vertices
switch n_fct
    case 1 % piecewise constant functions
        v{1} = sol(1:1:end,3); % first vertex
        v{2} = v{1}; % second vertex
        v{3} = v{1}; % third vertex
        v{4} = v{1}; % vertex between first and second
        v{5} = v{1}; % vertex between second and third
        v{6} = v{1}; % vertex between third and first
    case 3 % piecewise linear functions
        v{1} = sol(1:3:end,3); % first vertex
        v{2} = v{1}+sol(2:3:end,3); % second vertex
        v{3} = v{1}+sol(3:3:end,3); % third vertex
        v{4} = v{1}+sol(2:3:end,3)*.5; % vertex between first and second
        v{5} = v{1}+sol(2:3:end,3)*.5+...
            sol(3:3:end,3)*.5; % vertex between second and third
        v{6} = v{1}+sol(3:3:end,3)*.5; % vertex between third and first
    case 6 % piecewise quadratic functions
        v{1} = sol(1:6:end,3); % first vertex
        v{2} = v{1}+sol(2:6:end,3)+sol(4:6:end,3); % second vertex
        v{3} = v{1}+sol(3:6:end,3)+sol(6:6:end,3); % third vertex
        v{4} = v{1}+sol(2:6:end,3)*.5+...
            sol(4:6:end,3)*.25; % vertex between first and second
        v{5} = v{1}+sol(2:6:end,3)*.5+sol(3:6:end,3)*.5+...
            sol(4:6:end,3)*.25+sol(5:6:end,3)*.25+...
            sol(6:6:end,3)*.25; % vertex between second and third
        v{6} = v{1}+sol(3:6:end,3)*.5+...
            sol(6:6:end,3)*.25; % vertex between third and first
end

% Modify solution
try % Only rescale if rescale function is defined
    v{1} = feval(parameters.rescale,v{1});
    v{2} = feval(parameters.rescale,v{2});
    v{3} = feval(parameters.rescale,v{3});
    v{4} = feval(parameters.rescale,v{4});
    v{5} = feval(parameters.rescale,v{5});
    v{6} = feval(parameters.rescale,v{6});
end

% If output present, then merely report range of values
if nargout>0
    values = [v{1};v{2};v{3};v{4};v{5};v{6}];
    vax = [min(values),max(values)];
    return
end

if parameters.flag_average
    np = size(p,1); % number of vertices
    % Make solution continuous by averaging at vertices
    for j = 1:np % loop over all vertices (v{1},v{2},v{3})
        i1 = find(t(:,1)==j);
        i2 = find(t(:,2)==j);
        i3 = find(t(:,3)==j);
        val = mean([v{1}(i1);v{2}(i2);v{3}(i3)]); % average values
        v{1}(i1) = val; v{2}(i2) = val; v{3}(i3) = val;
    end
    n = size(t,1); % number of elements
    edges = [sort(t(:,[1 2]),2)*[2*np;1],(1:n)',ones(n,1)*4;...
        sort(t(:,[2 3]),2)*[2*np;1],(1:n)',ones(n,1)*5;...
        sort(t(:,[3 1]),2)*[2*np;1],(1:n)',ones(n,1)*6];
    [~,si] = sort(edges(:,1));
    edges = edges(si,:); % sort edges
    j = 1;
    while j<size(edges,1)
        if edges(j,1)==edges(j+1,1) % edge is shared
            val = (v{edges(j,3)}(edges(j,2))+...
                v{edges(j+1,3)}(edges(j+1,2)))/2;
            v{edges(j,3)}(edges(j,2)) = val;
            v{edges(j+1,3)}(edges(j+1,2)) = val;
            j = j+2;
        else % boundary edge
            j = j+1;
        end
    end
end

% Create new triangulation (each element split into four)
n = size(t,1); % number of elements (can be reduced if desired)
p1 = p(t(1:n,:)',:); % vertices of elements
a = reshape(1:3*n,3,[]); ind = reshape(a([2 3 1],:),1,[]);
p2 = (p1(1:3*n,:)+p1(ind,:))/2; % verticles at edge centers
% new list of vertices
p_new = [p1(1:3:end,:);p2(1:3:end,:);p2(3:3:end,:);...
    p1(2:3:end,:);p2(2:3:end,:);p2(1:3:end,:);...
    p1(3:3:end,:);p2(3:3:end,:);p2(2:3:end,:);...
    p2(1:3:end,:);p2(2:3:end,:);p2(3:3:end,:)];
p_new = p_new(reshape(reshape(1:3*4*n,n,[])',1,[]),:);
% new list of triangles
t_new = reshape(1:3*4*n,3,[])';
% new function values
c_new = [v{1}(1:n);v{4}(1:n);v{6}(1:n);...
    v{2}(1:n);v{5}(1:n);v{4}(1:n);...
    v{3}(1:n);v{6}(1:n);v{5}(1:n);...
    v{4}(1:n);v{5}(1:n);v{6}(1:n)];
c_new = c_new(reshape(reshape(1:3*4*n,n,[])',1,[]),:);

% Create coarse meshing (for mesh lines on 3d plot)
t_coarse = reshape(1:3*n,3,[])';
p_coarse = p1(1:3*n,:);
c_coarse = [v{1}(1:n);v{2}(1:n);v{3}(1:n)];
c_coarse = c_coarse(reshape(reshape(1:3*n,n,[])',1,[]),:);


% Visualize solution on mesh
%clf
switch parameters.plot_style
    case 1 % 3d surface plot
        trisurf(t_new,p_new(:,1),p_new(:,2),c_new,...
            'FaceColor','interp','EdgeColor','none')
%            'FaceColor','interp','EdgeColor',[1 1 1]*.8,'linewidth',.25)
        if parameters.flag_plotmesh
            hold on
            trisurf(t_coarse,p_coarse(:,1),p_coarse(:,2),c_coarse,...
                'FaceColor','none','EdgeColor','k')
            hold off
        end
        view(3)
    case 2 % 2d color plot
        patch('Faces',t_new,'Vertices',p_new,...
            'FaceVertexCData',c_new,...
            'FaceColor','interp','EdgeColor','none')
        if parameters.flag_plotmesh
            hold on
            %iptsetpref('ImshowBorder','tight')
            patch('Faces',t,'Vertices',p,...
                'FaceColor','none','EdgeColor','k')
            hold off
        end
        
        colorbar
        axis equal
        set(gcf,'renderer','zbuffer')
end
%colormap jet(250)
colormap parula
xlabel('x'), ylabel('y')
title(['All Limiters, t = ',num2str(time)])
%caxis([0,186]);
%caxis([-7,0]);
% Make 2d axis boundaries tight

set(gca,'XLim',[min(p(:,1)) max(p(:,1))],...
   'YLim',[min(p(:,2)) max(p(:,2))])
min(c_new)
max(c_new)
function w = rescale_log(v)
% Logarithmic rescaling
%w = max(min(log10(max(v,1e-12)),0),-3);
 w = log10(max(v,1e-7));

