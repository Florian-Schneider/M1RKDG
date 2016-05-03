function domain = generateDomain(x_limIn,y_limIn,n_xIn,n_yIn)

domain.x_lim = x_limIn;
domain.y_lim = y_limIn;
domain.n_x = n_xIn;
domain.n_y = n_yIn;
domain.n_int = domain.n_x*domain.n_y;
domain.dx = diff(x_limIn)/n_xIn;
domain.dy = diff(y_limIn)/n_yIn;
%[x,y] = meshgrid(linspace(x_limIn(1),x_limIn(2),n_xIn),linspace(y_limIn(1),y_limIn(2),n_yIn));
[x,y] = meshgrid(x_limIn(1) + domain.dx*(1/2+(0:(domain.n_x-1))),y_limIn(1) + domain.dy*(1/2+(0:(domain.n_y-1))));
%keyboard
domain.cell_centers = [x(:),y(:)];
%domain.dx = min(diff(linspace(x_limIn(1),x_limIn(2),n_xIn)));
%domain.dy = min(diff(linspace(y_limIn(1),y_limIn(2),n_yIn)));
domain.cell_neighbor_index = zeros(length(domain.cell_centers),4); %[x_l,x_r,y_l,y_r]
%keyboard
% Inner points
domain.cell_neighbor_index(:,1) = (1:domain.n_int)' - n_yIn; %left neighbor in x direction
domain.cell_neighbor_index(:,2) = (1:domain.n_int)' + n_yIn; %right neighbor in x direction
domain.cell_neighbor_index(:,3) = (1:domain.n_int)' - 1; %left neighbor in y direction = downwards
domain.cell_neighbor_index(:,4) = (1:domain.n_int)' + 1; %right neighbor in y direction = upwards

% Create boundary:
facet1 = [x(abs(x-x_limIn(1)-domain.dx/2)<domain.dx/100)-domain.dx,y(abs(x-x_limIn(1)-domain.dx/2)<domain.dx/100)];
facet2 = [x(abs(x-x_limIn(2)+domain.dx/2)<domain.dx/100)+domain.dx,y(abs(x-x_limIn(2)+domain.dx/2)<domain.dx/100)];
facet3 = [x(abs(y-y_limIn(1)-domain.dy/2)<domain.dy/100),y(abs(y-y_limIn(1)-domain.dy/2)<domain.dy/100)-domain.dy];
facet4 = [x(abs(y-y_limIn(2)+domain.dy/2)<domain.dy/100),y(abs(y-y_limIn(2)+domain.dy/2)<domain.dy/100)+domain.dy];

% Link facets
ntmp = domain.n_int;
domain.cell_neighbor_index((abs(x-x_limIn(1)-domain.dx/2)<domain.dx/100),1) = ntmp+(1:length(facet1))';
ntmp = ntmp + length(facet1);
domain.cell_neighbor_index((abs(x-x_limIn(2)+domain.dx/2)<domain.dx/100),2) = ntmp+(1:length(facet2))';
ntmp = ntmp + length(facet2);
domain.cell_neighbor_index((abs(y-y_limIn(1)-domain.dy/2)<domain.dy/100),3) = ntmp+(1:length(facet3))';
ntmp = ntmp + length(facet3);
domain.cell_neighbor_index((abs(y-y_limIn(2)+domain.dy/2)<domain.dy/100),4) = ntmp+(1:length(facet4))';
%ntmp = ntmp + length(facet4);   
domain.cell_centers = [domain.cell_centers;facet1;facet2;facet3;facet4];

% Quadrature points

Points = domain.cell_centers;
dx = domain.dx;
dy = domain.dy;
z = [-1,-sqrt(1/5),sqrt(1/5),1];
w = [1/6,5/6,5/6,1/6];
[Zx,Zy] = meshgrid(z);
W = w'*w;
W = W(:);
Zx = Zx(:);
Zy = Zy(:);
QuadraturePointsx = Zx;
QuadraturePointsy = Zy;
Quadrature.ReferencePointsx = QuadraturePointsx;
Quadrature.ReferencePointsy = QuadraturePointsy;
Quadrature.Weights = W/4*dx*dy;

Quadrature.LeftEdgex = 0*z-1;
Quadrature.LeftEdgey = z;
Quadrature.RightEdgex = 0*z+1;
Quadrature.RightEdgey = z;
Quadrature.LowerEdgex = z;
Quadrature.LowerEdgey = 0*z-1;
Quadrature.UpperEdgex = z;
Quadrature.UpperEdgey = 0*z+1;

Quadrature.PhysicalPointsx = bsxfun(@plus,dx/2*Quadrature.ReferencePointsx',x(:));
Quadrature.PhysicalPointsy = bsxfun(@plus,dy/2*Quadrature.ReferencePointsy',y(:));

domain.Quadrature = Quadrature;

% plotDomain();
    function plotDomain()
        clf
        plot(domain.cell_centers(:,1),domain.cell_centers(:,2),'x')
        hold on;
%         plot3(facet1(:,1),facet1(:,2),facet1(:,3),'or')
%         plot3(facet2(:,1),facet2(:,2),facet2(:,3),'ok')
%         plot3(facet3(:,1),facet3(:,2),facet3(:,3),'om')
%         plot3(facet4(:,1),facet4(:,2),facet4(:,3),'oc')
%         plot3(facet5(:,1),facet5(:,2),facet5(:,3),'og')
%         plot3(facet6(:,1),facet6(:,2),facet6(:,3),'oy')
        
        for i=1:domain.n_int
            h(1) = plot(domain.cell_centers(i,1),domain.cell_centers(i,2),'or');
            Style = {'or','ok','om','oc','og','ob'};
            for j=1:4
                h(j+1) = plot(domain.cell_centers(domain.cell_neighbor_index(i,j),1),domain.cell_centers(domain.cell_neighbor_index(i,j),2),Style{j});
            end
            drawnow;
            pause(2);
            delete(h)
        end
    end



end

