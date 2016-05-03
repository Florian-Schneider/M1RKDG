function plot_dg_movie
%PLOT_DG_MOVIE    Create sequence of frames from DG computation.
%   Calls PLOT_DG_SOLUTION with all files of one name, and creates
%   a sequence of plots for them.

% (C) 2011-06-07 by Benjamin Seibold
%     2011-11-20 video encoding added by Philipp Monreal

%-----------------------------------------------------------------------
% Parameters
%-----------------------------------------------------------------------
% filename_base = 's10unif0.01';
% filename_base = 's18cunif0.50';
% filename_base = 's18chkr0.10';
%filename_base = 's10mesh32768_2';
%filename_base = 's20chkr_unif3.msh';
% filename_base = 's9mesh_earfemale_0p18000';
% filename_base = 's17reflect7072';
%filename_base = 's27cunif0.03';
%filename_base= 's6unif0_00625';
filename_base='s27cunif0.03';
parameters.moment_nr = 0;
parameters.plot_style = 2;
parameters.flag_plotmesh = 0;
parameters.flag_average = 0; 
%parameters.rescale = @rescale_log;
parameters.flag_create_movie = 1;
%-----------------------------------------------------------------------
figure;

if parameters.flag_create_movie
    fid = fopen('filelist.txt','w');
end

% Find all data files of specified name
files = dir([filename_base,'_*']);
frame_nr = ones(1,length(files))*inf;
frame_suffix = cell(1,length(files));
for j = 1:length(files)
    if not(files(j).isdir)
        ind = find(files(j).name=='_',1,'last');
        frame_suffix{j} = files(j).name(ind+1:end);
        frame_nr(j) = str2double(frame_suffix{j});
    end
end
ind = not(isinf(frame_nr));
frame_nr = frame_nr(ind); frame_suffix = frame_suffix(ind);
[frame_nr,ind] = sort(frame_nr); frame_suffix = frame_suffix(ind);
n = length(frame_nr);

% Loop over all frames to find maximal data range
Vax = zeros(n,2);
for k = 1:n
    filename = [filename_base,'_',frame_suffix{k}];
    Vax(k,:) = plot_dg_solution(filename,parameters);
end
vax = [min(Vax(:,1)),max(Vax(:,2))];

% Loop over all frames
for k = 1:n
    filename = [filename_base,'_',frame_suffix{k}];
    fprintf('Frame %d  from  %s\n',k,filename)
    plot_dg_solution(filename,parameters)
    %set(gca,'ZLim',vax+[-1 1]*diff(vax)*.02)
    %caxis(vax)
    pause(.2)

    if parameters.flag_create_movie
        % Store .png
        pic_name = strcat('solution_',int2str(k),'.png');
        print('-dpng','-r300',pic_name,'-noui')
        
        % Save figure
        fig_name = strcat('solution_',int2str(k),'.jpg');
        saveas(gcf,fig_name)
        
        % Create file-list used to encode the video
        fprintf(fid, '%s\n', pic_name);
    end
end

if parameters.flag_create_movie
    fclose(fid);
    % Render movie using MEncoder
    % Change fps here if needed
    !mencoder mf://@filelist.txt -mf fps=5:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o solution.avi
end

%========================================================================

function w = rescale_log(v)
% Logarithmic rescaling
w = max(min(log10(max(v,1e-12)),0),-7);
% w = log10(max(v,1e-12));
