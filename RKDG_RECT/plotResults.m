function plotResults(Results,varargin)
% plotResults:
%   Input parameters:
%       ReferenceSolution: Fokker-Planck solution as provided by the code
%                          of Andreas Roth (roth@mathematik.uni-kl.de)
%       Scaling: cell array {string,[cmin cmax]}, string is either
%                   none: rescale in every time step, linear scale
%                   nlog: as none, but in log scale
%                   lin: scale as provided by [cmin,cmax], linear scale
%                   log: as lin, but in log scale
%       Time: [tmin,tmax], time interval in which the solution should be
%               plotted
%       FirstMoment: bool, shows norm(F)/E with black and white spots for
%                     loss of realizability
%       Movie: string, filename of the movie to be generated. The default
%               [] does not create a movie


close all
n = 3;

for i=1:2:length(varargin)
    if strcmp(varargin{i},'ReferenceSolution')
        FPSol = varargin{i+1};
        Results(end+1).tframe=FPSol.tframe;
        Results(end).E=FPSol.E;
        Results(end).points=FPSol.points;
        Results(end).Name = 'Fokker-Planck';
    end
    if strcmp(varargin{i},'Scaling')
        tmp = varargin{i+1};
        if ~iscell(tmp) || ~ischar(tmp{1}) || ~isnumeric(tmp{2}) || ~numel(tmp{2})==2
            error('Scaling must be a cell array {how to scale,[cmin cmax]}')
        end
        
        if strcmp(tmp{1},'lin') || strcmp(tmp{1},'log') || strcmp(tmp{1},'none') || strcmp(tmp{1},'nlog')
            Scaling = tmp{1};
        else
            error('This type of Scaling is not implemented or unknown!');
        end
        
        if numel(tmp{2})==2
            cmin = min(tmp{2});
            cmax = max(tmp{2});
        else
            cmin = [];
            cmax = [];
        end
        
    end
    if strcmp(varargin{i},'Time')
        Time = varargin{i+1};
    end
    
    if strcmp(varargin{i},'Cuts')
        Cuts = varargin{i+1};
    end
    
    if strcmp(varargin{i},'Movie')
        Movie = true;
        Filename = varargin{i+1};
    end
    
    if strcmp(varargin{i},'FirstMoment')
        FirstMoment = varargin{i+1};
    end
end

close all

if ~exist('Scaling','var')
    Scaling = 'lin';
end

if ~exist('cmin','var')
    cmin = [];
    cmax = [];
end

if ~exist('FirstMoment','var')
    FirstMoment = false;
end

if ~exist('Cuts','var')
    Cuts = {'XMean',[];'YMean',[];'Diagonal',[]};
end

PlotMode = '2D';

if ~exist('Time','var')
    Time = [Results(1).t(1),Results(1).t(end)];
end

if ~exist('Movie','var')
    Movie = false;
    Filename = [];
end

figure(24);
clf

scale_min = inf;
scale_max = -inf;

maplength = 500;


map = parula(maplength);

dx = Results(1).parameter.Domain.dx;
dy = Results(1).parameter.Domain.dy;

parameter = Results(1).parameter;
n_x = parameter.Domain.n_x;
n_y = parameter.Domain.n_y;

%p1 = @(x,y) 1;
p2 = @(x,y) x;
p3 = @(x,y) y;
p4 = @(x,y) p2(x,y).*p3(x,y);
p5 = @(x,y) p2(x,y).^2-1/3;
p6 = @(x,y) p3(x,y).^2-1/3;

if strcmp(PlotMode,'3D')
    z = 0.99*[NaN,-1,-sqrt(1/5),sqrt(1/5),1,NaN];
else
    z = [-1,-sqrt(1/5),sqrt(1/5),1];
end
[Zx,Zy] = meshgrid(z);
Zx = Zx(:);
Zy = Zy(:);


p2Val = p2(Zx,Zy);

p3Val = p3(Zx,Zy);
p4Val = p4(Zx,Zy);
p5Val = p5(Zx,Zy);
p6Val = p6(Zx,Zy);

PVal = [ones(1,length(p6Val));p2Val';p3Val';p4Val';p5Val';p6Val'];

if strcmp(PlotMode,'2D') %Shift the x/y-Data a bit
    z = 0.99*[-1,-sqrt(1/5),sqrt(1/5),1];
    [Zx,Zy] = meshgrid(z);
    Zx = Zx(:);
    Zy = Zy(:);
end

Px = bsxfun(@plus,parameter.Domain.cell_centers(1:parameter.Domain.n_int,1),dx/2*Zx');
Py = bsxfun(@plus,parameter.Domain.cell_centers(1:parameter.Domain.n_int,2),dy/2*Zy');

[Pxy,I] = sortrows([Px(:),Py(:)],[2; 1]);
Px = reshape(Pxy(:,1)',length(z)*n_x,length(z)*n_y);
Py = reshape(Pxy(:,2)',length(z)*n_x,length(z)*n_y);

LResults = length(Results);

for i=1:LResults
    
    if ~isfield(Results(i),'CallName')
        Results(i).CallName = Results(i).Name;
    end
    
    if isfield(Results(i),'Fx') && FirstMoment
        Results(end+1) = Results(i); %#ok<*AGROW>
        Results(end).FirstMoment = true;
        Results(i).FirstMoment = false;
    else
        Results(i).FirstMoment = false;
    end
    
    
    if isempty(cmin)
        for k=1:length(Results(i).E)
            scale_min = min(scale_min,min(min(min(Results(i).E{k}(:,1)))));
            scale_max = max(scale_max,max(max(max(Results(i).E{k}(:,1)))));
        end
        if strcmp(Scaling,'log')
            scale_min = max(scale_min,0);
        end
    else
        scale_min = cmin;
        scale_max = cmax;
    end
end

N = {};
for i=1:length(Results)
    N{i} = Results(i).CallName;
    
end
t_ind = 1;
t_length = length(Results(1).tframe);
for i=2:length(Results)
    if length(Results(i).tframe)<t_length
        t_ind = i;
        t_length = length(Results(i).tframe);
    end
end

[line1x,line1y,s1,s1lab,t1lab] = giveCuts(Cuts{1,1},Cuts{1,2});
[line2x,line2y,s2,s2lab,t2lab] = giveCuts(Cuts{2,1},Cuts{2,2});
[line3x,line3y,s3,s3lab,t3lab] = giveCuts(Cuts{3,1},Cuts{3,2});

tref = Results(t_ind).tframe;
Time = [max(Time(1),tref(1)),min(Time(2),tref(end))];
warning('off','all');
Specs = {'-','-.r',':k','--g','.m','-k'};
for i=1:length(Results)
    setaxis(2,n,i);
    
    if isfield(Results(i),'FirstMoment') && Results(i).FirstMoment
        
        rho = Results(i).E{1};
        rho = rho*PVal;
        rho = rho(I);
        rho = reshape(rho',length(z)*n_x,length(z)*n_y);
        
        Fx = Results(i).Fx{1};
        Fx = Fx*PVal;
        Fx = Fx(I);
        Fx = reshape(Fx',length(z)*n_x,length(z)*n_y);
        
        Fy = Results(i).Fy{1};
        Fy = Fy*PVal;
        Fy = Fy(I);
        Fy = reshape(Fy',length(z)*n_x,length(z)*n_y);
        
        rhosmall = rho<=1e-8;
        rho = sqrt(Fx.^2+Fy.^2)./rho;
        
        
        if strcmp(PlotMode,'3D')
            
            h(i) = surf(Px,Py,rho);
            
            shading interp
            
            axis equal tight
        else
            
            C1 = rho>1;
            C2 = rho<0;
            rho(rhosmall) = 0;
            rho(rho>1) = NaN;
            rho(rho<0) = NaN;
            
            ImageR = reshape(interp1(linspace(0,1,maplength),map(:,1),rho(:)),size(rho,1),size(rho,2));
            ImageG = reshape(interp1(linspace(0,1,maplength),map(:,2),rho(:)),size(rho,1),size(rho,2));
            ImageB = reshape(interp1(linspace(0,1,maplength),map(:,3),rho(:)),size(rho,1),size(rho,2));
            
            
            
            
            
            ImageR(C1) = 0;
            ImageR(C2) = 1;
            ImageG(C1) = 0;
            ImageG(C2) = 1;
            ImageB(C1) = 0;
            ImageB(C2) = 1;
            Image = cat(3,ImageR,ImageG,ImageB);
            
            
            
            h(i) = imagesc(Px(:,1),Py(1,:),Image,[0,1]); 
            
            axis equal tight
            
        end
        
        
    else
        
        rho = Results(i).E{1};
        rho = rho*PVal;
        rho = rho(I);
        rho = reshape(rho,length(z)*n_x,length(z)*n_y);
        
        if strcmp(PlotMode,'3D')
            
            if strcmp(Scaling,'none')
                h(i) = surf(Px,Py,rho);
                set(h(i),'CDataMapping','direct');
            else
                h(i) = surf(Px,Py,rho,(rho+scale_min)/(scale_min+scale_max));
                set(h(i),'CDataMapping','direct');
            end
            
            shading interp
            
            axis xy equal tight
        else
            
            
            if strcmp(Scaling,'log') || strcmp(Scaling,'nlog')
                rho(rho<1e-15) = 1e-15;
                
                h(i) = imagesc(Px(:,1),Py(1,:),log10(rho),[log10(scale_min),log10(scale_max)]); 
            else
                h(i) = imagesc(Px(:,1),Py(1,:),rho,[scale_min,scale_max]); 
            end
            
            
            
            axis xy equal tight
            
        end
        
    end
    
    drawnow
    
    title(Results(i).CallName)
    xlabel('x');
    ylabel('y');
    
    if ~(isfield(Results(i),'FirstMoment') && Results(i).FirstMoment)
        
        if strcmp(Scaling,'lin') || strcmp(Scaling,'none')
            setaxis(2,n,4);
            h1(i) = plot(s1,interp2(reshape(Results(i).points{1}(:),size(Results(i).points{1}')),reshape(Results(i).points{2}(:),size(Results(i).points{2}')),reshape(Results(i).E{1}(:,1),size(Results(i).points{1}')),line1x,line1y),Specs{i});
            xlabel(s1lab);
            hold on;
            legend(N);
            title(t1lab);
            setaxis(2,n,5);
            h2(i) = plot(s2,interp2(reshape(Results(i).points{1}(:),size(Results(i).points{1}')),reshape(Results(i).points{2}(:),size(Results(i).points{2}')),reshape(Results(i).E{1}(:,1),size(Results(i).points{1}')),line2x,line2y),Specs{i});
            xlabel(s2lab);
            hold on;
            legend(N);
            title(t2lab);
            setaxis(2,n,6);
            h3(i) = plot(s3,interp2(reshape(Results(i).points{1}(:),size(Results(i).points{1}')),reshape(Results(i).points{2}(:),size(Results(i).points{2}')),reshape(Results(i).E{1}(:,1),size(Results(i).points{1}')),line3x,line3y),Specs{i});
            xlabel(s3lab);
            hold on;
            legend(N);
            title(t3lab);
        else
            setaxis(2,n,4);
            h1(i) = semilogy(s1,interp2(reshape(Results(i).points{1}(:),size(Results(i).points{1}')),reshape(Results(i).points{2}(:),size(Results(i).points{2}')),reshape(Results(i).E{1}(:,1),size(Results(i).points{1}')),line1x,line1y),Specs{i});
            xlabel(s1lab);
            hold on;
            legend(N);
            title(t1lab);
            setaxis(2,n,5);
            h2(i) = semilogy(s2,interp2(reshape(Results(i).points{1}(:),size(Results(i).points{1}')),reshape(Results(i).points{2}(:),size(Results(i).points{2}')),reshape(Results(i).E{1}(:,1),size(Results(i).points{1}')),line2x,line2y),Specs{i});
            xlabel(s2lab);
            hold on;
            legend(N);
            title(t2lab);
            setaxis(2,n,6);
            h3(i) = semilogy(s3,interp2(reshape(Results(i).points{1}(:),size(Results(i).points{1}')),reshape(Results(i).points{2}(:),size(Results(i).points{2}')),reshape(Results(i).E{1}(:,1),size(Results(i).points{1}')),line3x,line3y),Specs{i});
            xlabel(s3lab);
            hold on;
            legend(N);
            title(t3lab);
        end
    end
end
warning('on','all');

[~,tminind] = min(abs(tref-Time(1)));
[~,tmaxind] = min(abs(tref-Time(2)));



screen_size = get(0, 'ScreenSize');
set(figure(24), 'Position', [0 0 min(screen_size(3:4)) min(screen_size(3:4)) ] );

if Movie
    fig1=figure(24);
    winsize = get(fig1,'Position');
    winsize(1:2) = [0 0];
    numframes=tmaxind-tminind;
    A=moviein(numframes,fig1,winsize);
    set(fig1,'NextPlot','replacechildren')
end


for j=tminind:tmaxind
    
    if strcmp(Scaling,'none')
        scale_min = inf;
        scale_max = -inf;
        [~,t_ind] = min(abs(tref(j)-Results(i).tframe));
        for i=1:LResults
            for k=t_ind
                scale_min = min(scale_min,min(min(min(Results(i).E{k}(:,1)))));
                scale_max = max(scale_max,max(max(max(Results(i).E{k}(:,1)))));
            end
        end
        if scale_min == scale_max
            scale_max = scale_max + 1e-9;
        end
        
        
    end
    
    for i=1:length(Results)
        setaxis(2,n,i);
        [~,t_ind] = min(abs(tref(j)-Results(i).tframe));
        
        if isfield(Results(i),'FirstMoment') && Results(i).FirstMoment
            
            rho = Results(i).E{t_ind};
            rho = rho*PVal;
            rho = rho(I);
            rho = reshape(rho,length(z)*n_x,length(z)*n_y);
            
            Fx = Results(i).Fx{t_ind};
            Fx = Fx*PVal;
            Fx = Fx(I);
            Fx = reshape(Fx',length(z)*n_x,length(z)*n_y);
            
            Fy = Results(i).Fy{t_ind};
            Fy = Fy*PVal;
            Fy = Fy(I);
            Fy = reshape(Fy',length(z)*n_x,length(z)*n_y);
            rhosmall = rho<=1e-8;
            rho = sqrt(Fx.^2+Fy.^2)./rho;
            
            if strcmp(PlotMode,'3D')
                
                h(i) = surf(Px,Py,rho);
                
                shading interp
                
                axis equal tight
            else
                
                C1 = rho>1;
                C2 = rho<0;
                rho(rhosmall) = 0;
                rho(C1) = NaN;
                rho(C2) = NaN;
                
                ImageR = reshape(interp1(linspace(0,1,maplength),map(:,1),rho(:)),size(rho,1),size(rho,2))';
                ImageG = reshape(interp1(linspace(0,1,maplength),map(:,2),rho(:)),size(rho,1),size(rho,2))';
                ImageB = reshape(interp1(linspace(0,1,maplength),map(:,3),rho(:)),size(rho,1),size(rho,2))';
                
                ImageR(C1') = 0;
                ImageR(C2') = 1;
                ImageG(C1') = 0;
                ImageG(C2') = 1;
                ImageB(C1') = 0;
                ImageB(C2') = 1;
                Image = cat(3,ImageR,ImageG,ImageB);
                
                set(h(i),'CData',Image);
                
                title([Results(i).CallName ' FM - t = ' num2str(Results(i).tframe(t_ind))])
                
                axis xy equal tight
                
            end
        else
            
            
            
            rho = Results(i).E{t_ind};
            rho = rho*PVal;
            rho = rho(I);
            rho = reshape(rho',length(z)*n_x,length(z)*n_y);
            
            if strcmp(PlotMode,'3D')
                
                set(h(i),'ZData',rho);
                if strcmp(Scaling,'lin')
                    set(h(i),'CData',(rho-scale_min)/(-scale_min+scale_max)*length(colormap));
                elseif strcmp(Scaling,'log')
                    set(h(i),'CData',(log10(rho)-log10(scale_min))/(-log10(scale_min)+log10(scale_max))*length(colormap));
                else
                    set(h(i),'CData',(rho-scale_min)/(-scale_min+scale_max)*length(colormap));
                    
                end
            else
                if strcmp(Scaling,'log') || strcmp(Scaling,'nlog')
                    rho(rho<1e-15) = 1e-15;
                    
                    h(i) = imagesc(Px(:,1),Py(1,:),log10(rho'),[log10(max(scale_min,1e-15)),log10(scale_max)]);
                else
                    h(i) = imagesc(Px(:,1),Py(1,:),rho',[scale_min,scale_max]); 
                end
                
                axis xy equal tight
                shading interp
            end
            
            title([Results(i).CallName ' - t = ' num2str(Results(i).tframe(t_ind))])
            
            setaxis(2,n,4)
            set(h1(i),'yData',interp2(reshape(Results(i).points{1}(:),size(Results(i).points{1}')),reshape(Results(i).points{2}(:),size(Results(i).points{2}')),reshape(Results(i).E{t_ind}(:,1),size(Results(i).points{1}')),line1x,line1y));
            
            setaxis(2,n,5)
            set(h2(i),'yData',interp2(reshape(Results(i).points{1}(:),size(Results(i).points{1}')),reshape(Results(i).points{2}(:),size(Results(i).points{2}')),reshape(Results(i).E{t_ind}(:,1),size(Results(i).points{1}')),line2x,line2y));
            
            setaxis(2,n,6)
            set(h3(i),'yData',interp2(reshape(Results(i).points{1}(:),size(Results(i).points{1}')),reshape(Results(i).points{2}(:),size(Results(i).points{2}')),reshape(Results(i).E{t_ind}(:,1),size(Results(i).points{1}')),line3x,line3y));
            
        end
        drawnow
        
        
        
        
        
    end
    
    
    if Movie
        A(:,j-tminind+1)=getframe(fig1,winsize);
    end
    
end

if Movie
    myVideo = VideoWriter(Filename);
    myVideo.FrameRate = 5;
    open(myVideo);
    writeVideo(myVideo, A);
    close(myVideo)
end

    function [x,y,s,slab,tlab] = giveCuts(name,options)
        
        
        xspan = Results(1).parameter.Domain.x_lim;
        yspan = Results(1).parameter.Domain.y_lim;
        
        switch name
            case 'YMean'
                yMitte = sum(Results(1).parameter.Domain.y_lim)/2;
                x = linspace(xspan(1),xspan(2));
                y = yMitte+0*linspace(xspan(1),xspan(2));
                s = x;
                slab = 'x';
                tlab = ['y = ' num2str(yMitte)];
                
            case 'XMean'
                xMitte = sum(Results(1).parameter.Domain.x_lim)/2;
                x = xMitte+0*linspace(yspan(1),yspan(2));
                y = linspace(yspan(1),yspan(2));
                s = y;
                slab = 'y';
                tlab = ['x = ' num2str(xMitte)];
                
            case 'Diagonal'
                x = linspace(xspan(1),xspan(2));
                y = linspace(yspan(1),yspan(2));
                s = x;
                slab = 'x=y';
                tlab = 'Diagonal';
                
            case 'Circle'
                if isempty(options)
                    a = [mean(xspan),mean(yspan)];
                    r = 0.8*min(diff(xspan),diff(yspan));
                else
                    a = options.a;
                    r = options.r;
                end
                
                
                
                x = a(1)+r*cos(linspace(0,2*pi));
                y = a(2)+r*sin(linspace(0,2*pi));
                s = linspace(0,2*pi);
                slab = '\alpha';
                tlab = ['Circle with r = ' num2str(r) ' around [' num2str(a(1)) ',' num2str(a(2)) ']'];
                
                
            case 'Ray'
                if isempty(options)
                    a = [0,0];
                    b = [0,1];
                else
                    a = options.a;
                    b = options.b;
                end
                
                sx1 = -(a(1) - xspan(2))/b(1);
                sy1 = -(a(2) - yspan(2))/b(2);
                sx2 = -(a(1) - xspan(1))/b(1);
                sy2 = -(a(2) - yspan(1))/b(2);
                
                S = [sx1,sy1,sx2,sy2];
                S(isinf(S)) = min(S(~isinf(S)));
                
                L = linspace(min(S),max(S));
                x = a(1)+b(1)*L;
                y = a(2)+b(2)*L;
                s = L;
                slab = 's';
                tlab = ['Ray [' num2str(a(1)) ',' num2str(a(2)) ']+s*[' num2str(b(1)) ',' num2str(b(2)) ']'];
                
                
            otherwise
                x = [];
                y = [];
                s = [];
                slab = [];
                tlab = [];
                
        end
    end

    function setaxis(n1,n2,i)
        subplot(n1,n2,i)
    end
end



