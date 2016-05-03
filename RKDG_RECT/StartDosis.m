function [Results]=StartDosis(name,varargin)
% STARTDOSIS: 
% Input parameters
% name [string]:
%   2BeamsUpRight
%   Checkerboard
%   Flash
%   HomDisk
%   Linesource
%   Shadow
%
% Methods [array of int]:
%   1 = First Order
%   Third order methods:
%   2 = Unlimited
%   3 = Scalar minmod, M = 0
%   4 = Scalar minmod, M = 50
%   5 = Realizability limiter
%   6 = Characteristic minmod, M = 0;
%   7 = Characteristic minmod, M = 50;
%   8 = 5+6
%   9 = 5+7
%   10 = 3+5
%   11 = 4+5
%
% Compression [bool]:
%   true/false discards/saves unnecessary quantities after calculation
%
% See also StartDosisScript for a wrapper which allows to change Nx and Ny
%
% Example: 
% Results = StartDosisScript('HomDisk',30,30,'Methods',[1 2 5 9])
% plotResults(Results,'FirstMoment',true,'Scaling',{'log',[1e-10,1e0]})

Equations = 'Toyproblem';

for i=1:2:length(varargin)
    if strcmp(varargin{i},'Methods')
        if isequal(union(varargin{i+1},1:11),1:11)
            Methods = varargin{i+1};
        else
            error('This type of Methods is not implemented or unknown!');
        end
    end
    if strcmp(varargin{i},'Compression')
        if strcmp(varargin{i+1},'true') || strcmp(varargin{i+1},'false')
            Compression = varargin{i+1};
        else
            error('This type of Compression is not implemented or unknown!');
        end
    end
    if strcmp(varargin{i},'GridData')
        if isstruct(varargin{i+1})
            GridData = varargin{i+1};
            if ~isfield(GridData,'xspan') || ~isfield(GridData,'yspan') ||~isfield(GridData,'Nx') || ~isfield(GridData,'Ny') || ~isfield(GridData,'CFL') || ~isfield(GridData,'Time')
                error('GridData has the wrong format!')
            end
        else
            error('This type of GridData is not implemented or unknown!');
        end
    end
end

if ~exist('Methods','var')
    Methods = 1;
end
if ~exist('Compression','var')
    Compression = 'true';
end
if ~exist('GridData','var')
    run(['Parameters/' Equations '/' name '/GridData'])
    GridData = ans; %#ok<NOANS>
end

cd(['Parameters/' Equations])
parameter = generateParameters(name,GridData);
cd(['Methods/' Equations])
methods = generateMethods(parameter,Methods);

for i=1:length(methods)
    clear functions
    
    cd(['Methods/' Equations]); 
    parameter.Method = methods(i);
    parameter.Name = name;
    
    fprintf('----------Starting Example %d - %s ----------\n',i,parameter.Method.CallName)
    
    [fields,vals] = run_calculation(parameter);
    for j=1:length(fields)
        Results(i).(fields{j}) = vals{j};
    end
    Results(i).parameter = parameter;
    Results(i).Name = methods(i).Name;
    Results(i).CallName = methods(i).CallName;
    Results(i).n_x = Results(i).parameter.Domain.n_x;
    Results(i).n_y = Results(i).parameter.Domain.n_y;
    Results(i).n_int = Results(i).parameter.Domain.n_int;
    Results(i).points = {reshape(Results(i).parameter.Domain.cell_centers(1:Results(i).n_int,1) ,Results(i).n_x,Results(i).n_y),reshape(Results(i).parameter.Domain.cell_centers(1:Results(i).n_int,2) ,Results(i).n_x,Results(i).n_y )};
    Results(i).E = Results(i).Y;
    Results(i).Fx = Results(i).Y;
    Results(i).Fy = Results(i).Y;
    for k=1:length(Results(i).E)
        Results(i).E{k} = Results(i).E{k}(:,:,1);
        Results(i).Fx{k} = Results(i).Fx{k}(:,:,2);
        Results(i).Fy{k} = Results(i).Fy{k}(:,:,3);       
    end
end




if strcmp(Compression,'true')
    R = rmfield(Results,{'Y','n_x','n_y','n_int'});
    Results = R;
end

if nargout==0
    tmp = cd;
    if ~(exist('Results','dir')==7)
       mkdir('Results')
    end
    cd('Results')
    warning off all
    mkdir(Equations);
    cd(Equations)
    mkdir(name)
    cd(name)
    warning on all
    
    
    str = [datestr(now,'yyyy-mm-dd_HH-MM-SS') ' - ' name ' - ' num2str(GridData.Nx) 'x' num2str(GridData.Ny)];
    for i=1:length(methods)
        str = [str '-' methods(i).CallName];
    end
    robustSave(str,'Results');
    cd(tmp);
end

end

function robustSave(varargin)
% try to save as given, and if get sizeTooBigForMATFile warning, try again in -v7.3 mode

    if any(strcmp(varargin, '-v7.3'));
        error('robustSave will try to save without v7.3, otherwise with. Please omit -v7.3 flag');
    end

    % clear last warning
    lastwarn('','')
    
    % try to save 
    cmdStr = sprintf('save(%s%s)', repmat(['''%s'','], 1, numel(varargin) - 1), '''%s''');
    cmd = sprintf(cmdStr, varargin{:});
    disp(cmd)
    evalin('caller', cmd);
    
    % if we get a sizeTooBigForMATFile warning, try -v7.3
    [msg, mi] = lastwarn();
    if strcmp(mi, 'MATLAB:save:sizeTooBigForMATFile')
        warning('ROBUSTSAVE:save:sizeTooBigForMATFile', ...
            'Big File Warning caught during save. trying -v7.3');
        delete(varargin{1});
        
        cmdStr = sprintf('save(%s%s)', repmat(['''%s'','], 1, numel(varargin)), '''%s''');
        cmd = sprintf(cmdStr, varargin{:}, '-v7.3');
        disp(cmd)
        evalin('caller', cmd);
        
    end
       
    % if we got a warning and it wasn't sizeTooBigForMATFile, throw an error
    if ~strcmp(mi, '') && ~strcmp(mi, 'MATLAB:save:sizeTooBigForMATFile')
        error(msg);
    end
    
end