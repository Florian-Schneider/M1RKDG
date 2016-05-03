function Results = StartDosisScript(name,Nx,Ny,varargin)
%STARTDOSISSCRIPT 
%   Same as StartDosis but allows to specify Nx and Ny
%   See also StartDosis

Equations = 'Toyproblem';

run(['Parameters/' Equations '/' name '/GridData'])
GridData = ans; %#ok<NOANS>
GridData.Nx = Nx;
GridData.Ny = Ny;
varargin{end+1} = 'GridData';
varargin{end+1} = GridData;
if nargout==1
    Results = StartDosis(name,varargin{:});
else
    StartDosis(name,varargin{:});
    Results = ans; %#ok<NOANS>
end
end

