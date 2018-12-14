%% function [Xsect, calcsAll, pWaves]=ePSproc_MFPAD(rlAll, varargin)
% 
%   Function to calculate MF PADs from ePolyScat matrix elements
%
%   INPUTS (only rlAll is required):
%       rlAll   Stucture containing matrix elements, dimensions (Ns rows, Ne cols), where Ns is the number of symmetries and Ne the number of energies
%               This struture is usually generated from ePSproc_read.m, which in turn requires a valid ePS results file.
%       p       Light helicity, scalar.  Optional, default is linear, p=0
%       eAngs   Eugler angles defining rotation of light into MF, vector [phi theta chi].  Optional, default =[0 0 0]
%       it      Components to use in degenerate cases, vector. Optional, default =1
%       ipComponents    Set ip components to use (1=length, 2=velocity gauge). Optional, default =1
%       res     Resolution for angular grid, scalar.  Optional, default is 50.
%
%   OUTPUTS:
%       D       MFPADs, 3D array (theta,phi,E), where [theta,phi]=meshgrid(linspace(0,pi,res),linspace(0,2*pi,res));
%       Xsect   Partial Xsects for each symmetry & E, from sum|partial waves|^2
%       calcsAll    Full structure of calculation results, dimensions (Ns rows, Ne cols) + an extra row for sum over symmetries
%
%   13/04/16        ePSproc version for release, see notes below
%   26/01/16        Small modification to 'it' to allow variable size.
%   18/08/15    v1, consolidated from various old codes during NO2 calculations & verification
%                   Converted to 2D indexing (previously linear)
%                   NEEDS A TIDY UP!
%
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

function [Xsect, calcsAll, pWaves]=ePSproc_MFPAD(rlAll, varargin)

%% Parse input arguments, set defaults if not passed

if nargin>4
    res=varargin{5};
else
    res=50;
end

if nargin>3
    ipComponents=varargin{4};
else
    ipComponents=1;
end

if nargin>2
    it=varargin{3};
else
    it=1; %[1 0];
end

if nargin>1   
    eAngs=varargin{2};
else
    eAngs=[0 0 0];
end

if nargin>0
    p=varargin{1};
else
    p=0;
end
    

%% Set up frame
pRot=eAngs(1);         % Set Euler angles (phi, theta, chi).  For circ pol rotation of LF z-axis (propagation axis) to MF z-axis (mol axis) requires tRot=pi/2
tRot=eAngs(2);
cRot=eAngs(3);

%*** Set rotation matrices (anon. fun. form)

% Construct for the moment as anonymous functions
Rx= @(t) [1 0 0; 0 cos(t) -sin(t); 0 sin(t) cos(t)];
Ry= @(t) [cos(t) 0 sin(t); 0 1 0; -sin(t) 0 cos(t)];
Rz= @(t) [cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1];

R= @(p,t,c) Rz(p)*Ry(t)*Rz(c);

Rlf=R(pRot,tRot,cRot);  % Define frame rotations for use later.

%*** Init variables
nRecords=length(rlAll(:));
nEnergies=size(rlAll,2);
nSymms=size(rlAll,1);

Lmax=0; % Use this for checking Lmax later

[theta,phi]=meshgrid(linspace(0,pi,res),linspace(0,2*pi,res));

for indE=1:nEnergies 
    calcsAll(nSymms+1,indE).D=zeros(res,res);   % Set master outputs for symmetry-summed MFPADs to zero.
end

%% *** Calculate d=sum(d_l,m,mu*Y_lm*D) for each record

% for n=1:nRecords  % Loop over energies - NOW REPLACED BY 2D indexing below
%     
% %    Lmax=0;     % Use this to check Lmax when looping over mat. elements, should also be able to check vs. ePS input file
%     
%     % Select column (energy) to add result to for para and perp components of d, and correct for case where n==nEnergies and rem returns 0
%     col=rem(n,nEnergies);
%     if col==0
%        col=nEnergies;
%     end
%     
%     % for ip=ipComponents  % loop over components
%     ip=ipComponents;    % Case for single ip component
%     
%         % Switch on ip
%         if ip==1
%            rawIdy=rlAll(n).rawIdy1;
%         elseif ip==2
%            rawIdy=rlAll(n).rawIdy2;
%         end

for indSymm=1:nSymms
    for indE=1:nEnergies
        
        ip=ipComponents;    % Case for single ip component
        % Switch on ip
        if ip==1
           rawIdy=rlAll(indSymm,indE).rawIdy1;
        elseif ip==2
           rawIdy=rlAll(indSymm,indE).rawIdy2;
        end

        % index=((rawIdy(:,4)==it(1))|(rawIdy(:,4)==it(2)));  % Select component, explicit size set for it
        index=sum(repmat(rawIdy(:,4),1,length(it))==repmat(it,length(rawIdy(:,4)),1),2);  % Select component, varible length it

        Cind=1;
        clear C;

        for m=find(index)'  % Loop over selected terms & sum

            l=rawIdy(m,2);
            Mf=rawIdy(m,1);        
            % for Mlf=-l:l    % Removed loop and second D in order to calculate MF (mu) terms for given polarization state (currently doesn't include summation over different polarization geometries p)
                C(Cind,:)=[l Mf rawIdy(m,5).*ePSproc_wignerD(1,rawIdy(m,3),-p,pRot,tRot,cRot) rawIdy(m,3)];
                Cind=Cind+1;
            % end

        end

        indThres=abs(C(:,3))>1e-10;  % Threshold values to discard small terms
        Cthres=C(indThres,:);
        % Cthres=C;                   % No thresholding

        % Loop over C terms & sum terms with same indicies (in above code there are/can be duplicates since there is no looping over l, only a selection of l) - NEATER way to do this?  Somehow include above?
    %     CsortInd=1;
    %     clear Csort;
    %     for l=0:max(C(:,1))
    %         for m=-l:l
    %             indSort=(C(:,1)==l) & (C(:,2)==m);
    %             Csort(CsortInd,:)=[l m sum(C(indSort,3))];
    %             CsortInd=CsortInd+1;
    %         end
    %     end
        % C2sort=[Csort(:,1:2) abs(Csort(:,3)) angle(Csort(:,3))];    % mag, phase form     
        Csort=Cthres;    % NO SORTING REQUIRED for MF version!

        if exist('Csort','var')&&(~isempty(Csort));
%            Ylm=ePSproc_Ylm_calc(Csort,theta,phi);  % Calculate Ylm partial wave  
                
%             % Alternate version with loop to check conj(Ylm) - tic-toc test shows no real difference to the above in any case.
                Ylm=zeros(res,res); 
                for cInd=1:size(Csort,1);
                    Ylm=Ylm+Csort(cInd,3).*conj(ePSproc_Ylm_calc([Csort(cInd,1:2) 1],theta,phi));
                    % Ylm=Ylm+Csort(cInd,3).*(ePSproc_Ylm_calc([Csort(cInd,1:2) 1],theta,phi));
                    % Ylm=Ylm+Csort(cInd,3).*conj(ePSproc_Ylm_calc([Csort(cInd,1:2) 1],theta,phi).*conj(Y1mu(:,:,Csort(cInd,4)+2)));
                end


            %         if rawIdy(m,1)<0
            %             Ylm=-Ylm;
            %         end

            %*** Scale factor to sqrt(Mb)
            sf=sqrt(2)*rlAll(indSymm,indE).MbNorm(ip); % +1i*rlAll(n).MbNorm(2)];  % Factor of 2 seems to match RLLs values, but why?  Could be Ylm normalization, or messed up degeneracy factor - PROBABLY THIS, since I'm only using it=1 component... maybe ip component?
            % sf=1;

            % Dlmmu=conj(Ylm)*sf;  % Calculate d(l,m,mu) for selected (l,m,mu)
            % Dlmmu=conj(Ylm);
            Dlmmu=Ylm;
            %        Dlmmu=conj(Ylm).*ePSproc_wignerD(1,rawIdy(m,3),p,pRot,tRot,cRot).*ePSproc_wignerD(1,rawIdy(m,1),p,pRot,tRot,cRot);  % Calculate d(l,m,mu) for selected (l,m,mu)
            
            calcsAll(indSymm,indE).C=C;    % Log partial waves, no thresholding
            calcsAll(indSymm,indE).Cthres=Cthres;    % Log partial waves, after thresholding

            % size(Csort)
        else
            Dlmmu=zeros(res,res);
            
            % calcsAll(n).C=[];
            calcsAll(indSymm,indE).C=C;    % Log partial waves, no thresholding
            
            sf=0;
        end
        
        % Main calc. result is summed over symmetries, individual components also logged below to calcsAll structure
        % D(:,:,indE)=D(:,:,indE)+Dlmmu*sf;
        calcsAll(nSymms+1,indE).D=calcsAll(nSymms+1,indE).D+Dlmmu*sf;
    
        % Log parameters to output structure - NOTE CURRENTLY ASSUMES ONLY ONE IP COMPONENT
        % calcsAll(n).Lmax=Lmax;
        Xsect(indSymm,indE)=sum(abs(Csort(:,3)).^2)*rlAll(indSymm,indE).MbNorm(ip);  % Total X-sect as sum |partial waves|^2
        calcsAll(indSymm,indE).eKE=rlAll(indSymm,indE).eKE;
        calcsAll(indSymm,indE).symm=rlAll(indSymm,indE).symm;
        calcsAll(indSymm,indE).euler=[pRot tRot cRot];
        calcsAll(indSymm,indE).Xsect=Xsect(indSymm,indE);
        calcsAll(indSymm,indE).D=Dlmmu*sf;
        calcsAll(indSymm,indE).XsectD=sum(sum(abs(Dlmmu*sf).^2))./(res.^2);      % Total X-sect as sum over MFPAD, renormalized by number of points.
        calcsAll(indSymm,indE).Rlf=Rlf;
        calcsAll(indSymm,indE).p=p;
        
        % Log global Lmax for later
        if max(calcsAll(indSymm,indE).C(:,1))>Lmax
            Lmax=max(calcsAll(indSymm,indE).C(:,1));
        end
    end
    
end

%% Resort output & sum over symmetries

% calcsAll=reshape(calcsAll,nEnergies,nSymms).';
% Xsect=reshape(Xsect,nSymms,nEnergies);

%*** Generate index for sorted partial waves, set for all possible L,M (no mu in this case, will be summed over below)
LMallInd=[0 0];
for l=1:Lmax
    for m=-l:l
        LMallInd(end+1,:)=[l m];
    end
end

%**** Extract partial waves
pWaveAll=zeros(size(LMallInd,1),nEnergies,nSymms);
pWaveAllMb=pWaveAll;

for indSymm=1:nSymms  % Loop over symmetries
    for n=1:nEnergies                
        % Loop over terms and extract for each energy & symmetry, summed over all other indicies (mu, it)
        clear pWaveTemp
        for m=1:size(LMallInd,1)
            indSort=(calcsAll(indSymm,n).C(:,1)==LMallInd(m,1)) & (calcsAll(indSymm,n).C(:,2)==LMallInd(m,2));
            pWaveTemp(m)=sum(calcsAll(indSymm,n).C(indSort,3));
        end

        calcsAll(indSymm,n).Cind=[LMallInd pWaveTemp.'];    % Sorted C for (symm, E)

        pWaveAll(:,n,indSymm)=pWaveAll(:,n,indSymm)+pWaveTemp.';

        pWaveAllMb(:,n,indSymm)=pWaveAllMb(:,n,indSymm)+pWaveTemp.'.*rlAll(indSymm,n).MbNorm(ip);

        xSectAll(indSymm,n)=calcsAll(indSymm,n).Xsect;
        xSectAllD(indSymm,n)=calcsAll(indSymm,n).XsectD;
        
    end
end

% Sum over symmetries
pWaveAll(:,:,end+1)=sum(pWaveAll,3);
pWaveAllMb(:,:,end+1)=sum(pWaveAllMb,3);
xSectAll(nSymms+1,:)=sum(xSectAll,1);
xSectAllD(nSymms+1,:)=sum(xSectAllD,1);

%*** Assign to output structure
pWaves.pWaveAll=pWaveAll;
pWaves.pWaveAllMb=pWaveAllMb;
pWaves.xSectAll=xSectAll;
pWaves.xSectAllD=xSectAllD;
pWaves.euler=[pRot tRot cRot];
pWaves.p=p;
pWave.LMallInd=LMallInd;

%*** Set final row in output to sum over symmetries
for n=1:nEnergies
    
    % calcsAll(nSymms+1,n).Dlmmu=D(:,:,n);
    
    calcsAll(nSymms+1,n).eKE=calcsAll(1,n).eKE;
    calcsAll(nSymms+1,n).symm='Sum';
    calcsAll(nSymms+1,n).euler=[pRot tRot cRot];
    calcsAll(nSymms+1,n).Xsect=sum(Xsect(:,n));  % Total X-sect as sum |partial waves|^2
    calcsAll(nSymms+1,n).XsectD=sum(abs(calcsAll(nSymms+1,n).D).^2)./(res.^2);      % Total X-sect as sum over MFPAD, renormalized by number of points.
    calcsAll(nSymms+1,n).Rlf=Rlf;
    calcsAll(nSymms+1,n).p=p;
    
end
    
