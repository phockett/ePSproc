%% function [rlAll, params, getCro]=ePSproc_read(fileBase)
%
% Function to read-in & parse ePolyScat results files
%
% INPUT fileBase   File to read in (string). Include full path to file if it is not in the current working directory.
%                  This function looks for certain flags in the ePS file:
%                  - DumpIdy segments for the raw matrix elements.
%                  - GetCro segments for X-sections
%                  If DumpIdy segments are not present in the ePS file, further processing with ePSproc based on the raw matrix elements will fail.
%
% OUTPUT structures     rlAll   with one page per DumpIdy segment in the ePS results
%                       params  which contains various global properties & indexes
%                       getCro  contains all cross-sections
%
% 13/04/16         ePSproc version for release, see notes below
% 29/03/16         Added error checking and correction for case of blank DumpIdy segments, common issue when clipping of ScatEng line occurs with ePS input
%                  Added some diagnostics to params output to track this.
%               
% 17/08/15         2015 version, streamlined to allow for reading of multiple DumpIdy segments from a single file.  Use find or grep to parse file, then scan each DumpIdy section as previously.
%                               Added additional output structure "params", which contains various global properties & indexes
%
% 08/11/10         2d, Following below, but read in header line which now provides necessary normalization values.
% 16/06/10         2c, small change for newer version of ePolyScat (E2 dated June 2009) which has an extra header line in output files.
% 11/03/10         2b, changed textscan lineskips which seemed inconsistent on windows machines due to win/unix EOL differences
% 14/01/10 Version 2, still ugly, but added option to read only certain files by setting N to a vector of file numbers
% 17/11/09 Version 1, ugly but it works!
%
% Known BUGS:
%       On Linux (Ubuntu 14 LTS) seems to be giving occasional issues with drop-out during file IO, not sure why.  Rerunning function seems to fix.
%
% *** NOTES
%
%  ePSproc: Post-processing code for ePolyScat calculations
%  https://github.com/phockett/ePSproc
%  Released under a GNU General Public License (v3)
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

function [rlAll, params, getCro]=ePSproc_read(fileBase)

%% Display message
disp('*** Reading ePS output file');
disp(fileBase);


%% *** Parse input file for DumpIdy segments
haystack = fileBase;
needle   = 'DumpIdy - dump';

if ispc % Windows
    [~,lines] = system(['find /n "' needle '" ' haystack]);
    % dumpSeg=textscan(lines,'%*[^[] [%n] %*[^[]')  % Grab line #s from output of above string search.  Output from find starts with full file name string, skip it, and line numbers defined by opening [.  Everything is returned as a single string.
    dumpSeg=textscan(lines,'%*[^[] [%n %*[^coefs]');
elseif isunix % Mac, Linux
    [~,lines] = system(['grep -n "' needle '" ' haystack]);     % This returns a line for each token found, with the line number first.
    % Scan result to find correct part of file
    dumpSeg=textscan(lines,'%n %*[^\n]');  % Grab line #s
else 
    error('Unknown operating system!');
end

dumpSeg=cell2mat(dumpSeg);

N=size(dumpSeg,1);
disp(['Found ' num2str(N) ' sets of matrix elements']);

Lmax=0; % Use to log max L

fid=fopen(fileBase);
% textscan(fid,'%s','HeaderLines',dumpSeg(1)),'Delimiter','\n'); %Skip some lines to get to first segment

% Loop over segments & read in matrix elements
for record=1:N
    
    eKE=textscan(fid,'%*s %*s %*s %*s %*c %n','HeaderLines',dumpSeg(record)+1); % Read in energy value. Set headerlines == line #, after segment read reset file position to beginning... maybe a little slow, but works OK.
    textscan(fid,'%s',1,'Delimiter','\n'); %Skip lines
    MbNorm=textscan(fid,'%*s %*s %*s %*s %*s %*s %*s %*c %n %n',1); % Read scale factor to sqrt(Mb)
    symm=textscan(fid,'%*s %*s %s',1);      %Read continuum symmetry
    textscan(fid,'%s',5,'Delimiter','\n'); %Skip some lines
    rlIn=textscan(fid,'%n %n %n %n %n %n %n'); %Read raw matrix elements
    
    frewind(fid);   % Set file pointer back to start for next read, might be faster to reindex dumpSeg (as per mol. geom code below)... but not trivial with multi-structured input
        
    % Compile into Matlab structure for use
    rlAll(record).eKE=eKE{1};
    rlAll(record).symm=symm{1}{1}(2:end);
    rlAll(record).MbNorm=[MbNorm{1} MbNorm{2}];
       
    index=(rlIn{1,4}==1);   % Assign form 1 (length) and also reformat into mag/phase.
    %rlAll(record).rawIdy1=[rlIn{1,2}(index) rlIn{1,1}(index) rlIn{1,6}(index)+rlIn{1,7}(index)*1i];
    rlAll(record).rawIdyHead={'Raw matrix elements, legnth (1) or velocity (2) form';'m,l,mu,it,rawIdy'};
    rlAll(record).rawIdy1=[rlIn{1,1}(index) rlIn{1,2}(index) rlIn{1,3}(index) rlIn{1,5}(index) rlIn{1,6}(index)+rlIn{1,7}(index)*1i];
    rlAll(record).rlnlHead={'Raw matrix elements, legnth (1) or velocity (2) form, reformatted into (real) magnitude and phase';'l,m,rlm,nlm,mu,it'};
    rlAll(record).rlnl1=[rlAll(record).rawIdy1(:,2) rlAll(record).rawIdy1(:,1) abs(rlAll(record).rawIdy1(:,5)) angle(rlAll(record).rawIdy1(:,5)) rlAll(record).rawIdy1(:,3) rlAll(record).rawIdy1(:,4)];
    
    index=(rlIn{1,4}==2);  % Assign form 2 (velocity) and also reformat into mag/phase.
    %rlAll(record).rawIdy2=[rlIn{1,2}(index) rlIn{1,1}(index) rlIn{1,6}(index)+rlIn{1,7}(index)*1i rlIn{1,5}(index)];
    rlAll(record).rawIdy2=[rlIn{1,1}(index) rlIn{1,2}(index) rlIn{1,3}(index) rlIn{1,5}(index) rlIn{1,6}(index)+rlIn{1,7}(index)*1i];
    rlAll(record).rlnl2=[rlAll(record).rawIdy1(:,2) rlAll(record).rawIdy1(:,1) abs(rlAll(record).rawIdy1(:,5)) angle(rlAll(record).rawIdy1(:,5)) rlAll(record).rawIdy1(:,3) rlAll(record).rawIdy1(:,4)];
    
    symmAll{record}=rlAll(record).symm;
    eAll(record)=rlAll(record).eKE;
    
    if max(rlAll(record).rlnl1(:,1))>Lmax   % Track global Lmax for use later
        Lmax=max(rlAll(record).rlnl1(:,1));
    end
    
    %*** Error checking. In some cases miss records due to, e.g. clipping in input file (Fortran fixed length restrictions etc.)
    if size(rlAll(record).rlnl1,1)<2
        errFlag(record)=true;
    else
        errFlag(record)=false;
    end
    
end

fclose(fid);

%% *** Remove null (erroneous) records
missingE=eAll(errFlag);
rlAll(errFlag)=[];
eAll(errFlag)=[];
symmAll(errFlag)=[];
N=N-sum(errFlag);
disp(['Read ' num2str(N) ' sets of matrix elements (' num2str(sum(errFlag)) ' blank records)']);

%% *** Scan and sort results
rInd=1;
symmInd=1;
while rInd<=N
    index=find(strcmp(rlAll(rInd).symm,symmAll));  % Take reference symm and check which records it corresponds to
    
    eInd(:,symmInd)=index.';  % Log index into energies
    symmList{symmInd}=rlAll(rInd).symm; % Log unique symmetries
    eKE=eAll(index);        % Set vector of energies
    
    rInd=rInd+length(index);
    symmInd=symmInd+1;
end

% Echo to screen
% size(symmList);
disp(['Found ' num2str(symmInd-1) ' symmetries']);
disp(symmList(:).');
disp(['Found ' num2str(length(eKE)) ' energies']);

% Log global properties & indexes to output structure
params.symmList=symmList;
params.symmList{end+1}='All';
params.eKE=eKE;
params.symmAll=symmAll;
params.eAll=eAll;
params.nRecords=N;         % Set nRecords & nEnergies for back-compatibility with old code
params.nEnergies=length(eKE);
params.nSymms=symmInd-1;
params.gLmax=Lmax;              % Global Lmax
params.blankRec=sum(errFlag);
params.missingE=missingE;

% Reshape rlAll to dimensions (E x symmetry)
rlAll=reshape(rlAll,length(eKE),symmInd-1).';     

%% *** Create re-indexed matrix elements, summed over some variables, for easy plotting later

%*** Generate index for sorted partial waves, set for all possible L,M (no mu in this case, will be summed over below)
LMallInd=[0 0];
for l=1:Lmax
    for m=-l:l
        LMallInd(end+1,:)=[l m];
    end
end

%**** Extract partial waves
ip=1;
pWaveAll=zeros(size(LMallInd,1),params.nEnergies,params.nSymms);
pWaveAllMb=pWaveAll;

for indSymm=1:params.nSymms  % Loop over symmetries
    for n=1:params.nEnergies                
        % Loop over terms and extract for each energy & symmetry, summed over all other indicies (mu, it)
        clear pWaveTemp
        for m=1:size(LMallInd,1)
            indSort=(rlAll(indSymm,n).rawIdy1(:,2)==LMallInd(m,1)) & (rlAll(indSymm,n).rawIdy1(:,1)==LMallInd(m,2));
            pWaveTemp(m)=sum(rlAll(indSymm,n).rawIdy1(indSort,5));
        end

        rlAll(indSymm,n).pWaveAll=[LMallInd pWaveTemp.'];    % Sorted pWaves for (symm, E)

        pWaveAll(:,n,indSymm)=pWaveAll(:,n,indSymm)+pWaveTemp.';

        pWaveAllMb(:,n,indSymm)=pWaveAllMb(:,n,indSymm)+pWaveTemp.'.*rlAll(indSymm,n).MbNorm(ip);

        % xSectAll(indSymm,n)=calcsAll(indSymm,n).Xsect;
        
    end
end

% Sum over symmetries
pWaveAll(:,:,end+1)=sum(pWaveAll,3);
pWaveAllMb(:,:,end+1)=sum(pWaveAllMb,3);

% Log to structure
params.pWaveAll=pWaveAll;
params.pWaveAllMb=pWaveAllMb;
params.LMallInd=LMallInd;

%% *** Read molecular structure

needle   = 'Atoms found';

if ispc % Windows
    [~,lines] = system(['find /n "' needle '" ' fileBase]);
    atoms=textscan(lines,'%*[^[] %*c %n %*s %*s %n %*[^\n]');  % Grab line # and # of atoms from output of above string search.  Output from find starts with full file name string, skip it, and line number defined by opening [
elseif isunix % Mac, Linux
    [~,lines] = system(['grep -n "' needle '" ' fileBase]);     % This returns a line for each token found, with the line number first.
    % Scan result to find correct part of file
    atoms=textscan(lines,'%n %*s %*s %n %*[^\n]');  % Grab line # and # of atoms from output of above string search
else 
    error('Unknown operating system!');
end

fid=fopen(fileBase);
% coords=textscan(fid,'%*c %*c %d %*c %*c %d %*c %*c %d %d %d',atoms{2},'HeaderLines',atoms{1});
coords=textscan(fid,'%*c %*c %d %*c %*c %*c %d %*c %*c %f %f %f',atoms{2},'HeaderLines',atoms{1});  % Read atoms{2} lines, which start at line number atoms{1}+1
fclose(fid);

params.coords=coords;

disp(['Found ' num2str(length(coords{1})) ' atoms']);

%% *** Read some other useful parameters - any other input records set 
% needle='+ End of input'
needle='+ Data Record';
if ispc % Windows
    [~,lines] = system(['find /n "' needle '" ' fileBase]);
    dumpSeg=textscan(lines,'%*[^[] [%n %*[^Record]');         % Grab line # from output of above string search.  Output from find starts with full file name string, skip it, and line number defined by opening [
elseif isunix % Mac, Linux
    [~,lines] = system(['grep -n "' needle '" ' fileBase]);     % This returns a line for each token found, with the line number first.
    % Scan result to find correct part of file
    % atoms=textscan(lines,'%n %*s %*s %n %*[^\n]')  % Grab line # and # of atoms from output of above string search
    dumpSeg=textscan(lines,'%n %*[^\n]');  % Grab line #s
else 
    error('Unknown operating system!');
end

dumpSeg=cell2mat(dumpSeg);

N=size(dumpSeg,1);
disp(['Found ' num2str(N) ' data records']);

if ~isempty(dumpSeg)
    fid=fopen(fileBase);

    % If records are contiguous can use this version, otherwise set loop below
    % temp=textscan(fid,'%*c %*s %*s %s %*c %n %*[^\n]','HeaderLines',dumpSeg{1});
    
    % Loop over found data records & assign record name and value pairs to output
    for seg=1:N
      % temp=textscan(fid,'%*c %*s %*s %s %*c %n %*[^\n]',1,'HeaderLines',dumpSeg(seg));    % OK for numerical values, but not otherwise
      temp=textscan(fid,'%*c %*s %*s %s %*c %[^\n]',1,'HeaderLines',dumpSeg(seg)-1);
      params.dataRecords(seg,:)={temp{1,1}{1} temp{1,2}{1}};
      
      frewind(fid);
      
      if strcmp(temp{1,1}{1},'IPot')    % Assign IP if found
          params.IP=str2num(temp{1,2}{1});
      end
      
    end
    
    fclose(fid);
end

%% *** Read any cross-sections (GetCro)
needle='COMPOSITE CROSS SECTIONS';      % This grabs position of table of GetCro outputs

if ispc % Windows
    [~,lines] = system(['find /n "' needle '" ' fileBase]);
    dumpSeg=textscan(lines,'%*[^[] [%n %*[^ENERGIES]');         % Grab line # from output of above string search.  Output from find starts with full file name string, skip it, and line number defined by opening [
elseif isunix % Mac, Linux
    [~,lines] = system(['grep -n "' needle '" ' fileBase]);     % This returns a line for each token found, with the line number first.
    % Scan result to find correct part of file
    % atoms=textscan(lines,'%n %*s %*s %n %*[^\n]')  % Grab line # and # of atoms from output of above string search
    dumpSeg=textscan(lines,'%n %*[^\n]');  % Grab line #s
else 
    error('Unknown operating system!');
end

dumpSeg=cell2mat(dumpSeg);

N=size(dumpSeg,1);
disp(['Found ' num2str(N) ' sets of cross sections']);

if ~isempty(dumpSeg)
    fid=fopen(fileBase);

    for seg=1:N
        % xSectAllePS(seg).xSectL=textscan(fid,'%f %f','HeaderLines',dumpSeg(seg)+9);   % Read X-sect from start of GetCro segment
        % temp=textscan(fid,'%*s %f %f %f %f %f %f %f','HeaderLines',dumpSeg(seg)+1,'CollectOutput',true);     % Read table of GetCro output

        GetCroHeader=textscan(fid,'%[^\n]',1,'HeaderLines',dumpSeg(seg));       % Grab headerline for reference.
        temp=textscan(fid,'%*s %f %f %f %f %f %f %f','CollectOutput',true);     % Read table of GetCro output

        % xSectAllePS(seg).xSectL=cell2mat(temp);
        getCro(seg).GetCro=cell2mat(temp);

        frewind(fid);

    end

    fclose(fid);
    
else
    getCro=[];
    
end

params.GetCroHeader=GetCroHeader;
    
