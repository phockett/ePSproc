%% ePSproc - demo calculations
%  ePSproc: Post-processing code for ePolyScat calculations
%  https://github.com/phockett/ePSproc
%  Released under a GNU General Public License (v3)
%
%  ePSproc contains code to:
%   - Read raw matrix elements from ePS output files with "dumpIdy" segments
%   - Calculate MF PADs from the matrix elements
%   - Plot MF PADs
%   - Plot partial waves
%   - Plot X-sects
%
%  This demo script illustrates the use of the various functions, and reproduces the sample NO2 results included with the ePSproc suite.
%
%  13/04/16    Tidied up a little for first Github release
%  25/09/15    Fixed issue with z-axis defn.
%  18/08/15    ePSproc v1, consolidated from various old codes. 
%              Tested & verified against NO2 calculations from Toffoli et. al., JCP 126, 054307, 2007
%
%  ePSproc code:
%  Paul Hockett
%  paul.hockett@nrc.ca
%  femtolab.ca
%  github.com/phockett
%
%  For details about ePolyScat (ePS), a tool for computation of electron-molecule scattering, see:
%   - ePS website & manual, maintained by R.R. Lucchese
%       http://www.chem.tamu.edu/rgroup/lucchese/ePolyScat.E3.manual/manual.html)
%   - F. A. Gianturco, R. R. Lucchese, and N. Sanna, J. Chem. Phys. 100, 6464 (1994).  
%       http://dx.doi.org/10.1063/1.467237
%   - A. P. P. Natalense and R. R. Lucchese, J. Chem. Phys. 111, 5344 (1999). 
%       http://dx.doi.org/10.1063/1.479794
%
%% *** SETTINGS
%  Set up basic environment

% Name & path to ePS output file. In this version set full file name here, and working directory below.
fileName='no2_demo_ePS.out'

% Set paths for Linux or Win boxes (optional!)
if isunix
    dirSlash='/';
else
    dirSlash='\';
end

filePath=pwd;                       % Root to working directory, here set as current dir.
fileBase=[filePath dirSlash fileName];   % Full path to ePS results file, here set as current working direcory
    
path(path,[filePath]);   % Add path to ePSproc scrips to Matlab path list, here assumed to be in root dir
    
%% *** Read data
%  Variables:
%       rlAll contains matrix elements (from DumpIdy segments)
%       params contains various calculation parameters
%       getCro contains cross-section (from GetCro segments), if present

[rlAll, params, getCro]=ePSproc_read(fileBase);

params.fileBase=fileBase;
params.fileName=fileName;

%% Plot GetCro results for each symm & total

col=2;  % Select column from getCro output (see params.GetCroHeader)

figure('color',[1 1 1],'name','GetCro outputs');

for n=1:length(getCro)
    plot(getCro(n).GetCro(:,1)-params.IP,getCro(n).GetCro(:,col));
    hold on;
end

title({['NO_2 ePS resutls, files ' strrep(fileName,'_','\_')]; 'X-sects from ePS(GetCro) results'});
xlabel('eKE/eV');
ylabel('X-sect/Mb');

legend([params.symmList 'Sum']);

%% *** Calculate MFPADs - single polarization geometry, all energies and symmetries
%  Calculate for specified Euler angles (polarization geometry) & energies

% Set resolution for calculated I(theta,phi) surfaces
res=100;

% ip components to use from ePS output (1=length gauge, 2=velocity gauge)
ipComponents=1;

% it components to use from ePS output (for degenerate cases), set an array here for as many components as required, e.g. it=1, it=[1 2] etc.
it=1;

% Set light polarization and axis rotations LF -> MF
p=0;                % p=0 for linearly pol. light, +/-1 for L/R circ. pol.
eAngs=[0 0 0];      % Eugler angles for rotation of LF->MF, set as [0 0 0] for z-pol, [0 pi/2 0] for x-pol, [pi/2 pi/2 0] for y-pol

% Run calculation - outputs are D, full set of MFPADs (summed over symmetries); Xsect, calculated X-sects; calcsAll, structure with results for all symmetries.
[Xsect, calcsAll, pWaves]=ePSproc_MFPAD(rlAll,p,eAngs,it,ipComponents,res);
    

%% Plotting - MFPAD panel plots

% Set plot ranges
symmInd=1;     % Select symmetry (by index into calcsAll rows). Final symmetry state is set as sum over all symmetries
eRange=1;      % Select energies (by index into calcsAll cols)

% Additional options (optional)
sPlotSet=[1 2];             % Set [rows cols] for subplot panels. The final panel will be replaced with a diagram of the geometry
titlePrefix='NO2 testing';  % Set a title prefix for the figure

ePSproc_MFPAD_plot(calcsAll,eRange,symmInd,params,sPlotSet,titlePrefix);
