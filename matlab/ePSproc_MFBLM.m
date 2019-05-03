%% function calcsAll=ePSproc_MFBLM_2019(rlAll, p, eAngs, comment)
% 
%   Function to calculate MF BLMs from ePolyScat matrix elements
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
%       calcsAll    Full structure of calculation results, dimensions (Ns rows, Ne cols) + an extra row for sum over symmetries
%
%   25/04/19        Development version for ePSproc.
%                   Based on old code family "pad_simple_MF.m", and recent updates/testing for ePS matrix elements "pad_simple_MF_ePS_2019.m", for direct BLM calculation with external gamma calcs.
%                   Implement here using loop over ePS matrix elements, as per "ePSproc_MFPAD.m", with gamma calcualtion included.
%
%   STATUS:         Tested for NO2 matrix elements vs. numerics (ePSproc_MFPAD.m). Working, but quite slow - could do with some more optimisation & tidy up.
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


function calcsAll=ePSproc_MFBLM_2019(rlAll, varargin)

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

%% Init variables
nRecords=length(rlAll(:));
nEnergies=size(rlAll,2);
nSymms=size(rlAll,1);

Lmax=0; % Use this for checking Lmax later

for indE=1:nEnergies 
    calcsAll(nSymms+1,indE).D=zeros(res,res);   % Set master outputs for symmetry-summed MFPADs to zero.
end

%% Set up frame
pRot=eAngs(1);         % Set Euler angles (phi, theta, chi).  For circ pol rotation of LF z-axis (propagation axis) to MF z-axis (mol axis) requires tRot=pi/2
tRot=eAngs(2);
cRot=eAngs(3);

[theta,phi]=meshgrid(linspace(0,pi,res),linspace(0,2*pi,res));


%% *** Calculate BLM for each record

for indSymm=1:nSymms
    for indE=1:nEnergies
        
        ip=ipComponents;    % Case for single ip component
        % Switch on ip
        if ip==1
           rawIdy=rlAll(indSymm,indE).rawIdy1;
        elseif ip==2
           rawIdy=rlAll(indSymm,indE).rawIdy2;
        end
        
        %*** Scale factor to sqrt(Mb)
        % sf=sqrt(2)*rlAll(indSymm,indE).MbNorm(ip); % +1i*rlAll(n).MbNorm(2)];  % Factor of 2 seems to match RLLs values, but why?  Could be Ylm normalization, or messed up degeneracy factor - PROBABLY THIS, since I'm only using it=1 component... maybe ip component?    
        sf=2*rlAll(indSymm,indE).MbNorm(ip);
        
        % index=sum(repmat(rawIdy(:,4),1,length(it))==repmat(it,length(rawIdy(:,4)),1),2)  % Select component, varible length it
        index=sum((repmat(rawIdy(:,4),1,length(it))==repmat(it,length(rawIdy(:,4)),1)).*(abs(rawIdy(:,5))>1e-3),2);  % Select component, varible length it + thresholding

        Cind=1;
        clear C;

        for m=find(index)'  % Loop over selected terms & sum
            % Set (l,m)
            l=rawIdy(m,2);
            Mf=rawIdy(m,1); 
            mu=rawIdy(m,3);
            
            for mp=find(index)'  % Loop over selected terms & sum
                % Set (l',m')
                lp=rawIdy(mp,2);
                Mfp=rawIdy(mp,1);  
                mup=rawIdy(mp,3);
                          
                for L=0:(l+lp)  % Loop over allowed B(L,M) terms
                    % for M=-L:L
                    M=(-Mf+Mfp); % *** Might be a phase issue here... TBC
                                 % UPDATE: seems to be correct, with -M in the 3j term below and in the output.  Switching to M=-(-Mf+Mfp) here produces cylindrically symmetric output.
                    
                    % Calculate associated gamma term, (L,M) values
                    gammaLM=ePSproc_3j(l,lp,L,0,0,0)*ePSproc_3j(l,lp,L,-Mf,Mfp,-M);
                    
                    % Gamma terms for polarization - sum over P
                    gammaP=0;
                    matEle = rawIdy(m,5)*conj(rawIdy(mp,5));  % Product of matrix elements for this set of QNs - TO DO CHECK INDEXING HERE
                    % matEle = conj(rawIdy(m,5))*(rawIdy(mp,5));
                    for P=0:2
                        Rp=mup-mu;
                        if abs(Rp)<= P  % Check for allowed terms, otherwise wignerD() will throw an error
                            
                            % Sum over R,R' projections
                            % WORKS with full summation, but slow (lots of expensive zeros).
                            % NOW working with specific R set only.
                            R=0; % For single pol state, i.e. p-p'=0
                            % for R=-P:P
                               % for Rp=-P:P
                               % Rp=mup-mu;
                                    gammaP = gammaP + (2*P+1)*(-1)^(Rp-R)*ePSproc_3j(1,1,P,mu,-mup,Rp)*ePSproc_3j(1,1,P,-p,p,R)*ePSproc_wignerD(P,-Rp,-R,pRot,tRot,cRot)*matEle;   % Note polarization terms, here mu is MF, and p is LF.
                               % end
                            % end
                            % gammaP = gammaP + (2*P+1)*ePSproc_3j(1,1,P,mu,-mup,mu-mup)*ePSproc_3j(1,1,P,-p,p,0)*ePSproc_wignerD(P,-mu+mup,0,pRot,tRot,cRot)*matEle;   % Note polarization terms, here mu is MF, and p is LF.
                            % gammaP = gammaP + (2*P+1)*ePSproc_3j(1,1,P,mu,-mup,mu-mup)*ePSproc_3j(1,1,P,-p,p,0)*ePSproc_wignerD(P,0,-mu+mup,pRot,tRot,cRot)*matEle;   % Test order switching - no change.
                        end
                    end
                    
                    %*** Final set of terms
                    % phase = (-1)^M*(-1)^Mf*(-1)^(mup+p)*(-1)^(mu-mup);
                    % phase = (-1)^Mf*(-1)^(mup+p)*(-1)^(mu-mup);
                    phase = (-1)^M*(-1)^Mf*(-1)^(mup+p);
                    degen = sqrt(((2*l+1)*(2*lp+1)*(2*L+1))/(4*pi));
                    
                    gamma = gammaLM*gammaP*phase*degen;
                    
                    LMterm = gamma*(sf^2);
        
                    % Log term in C list
                    % [l Mf lp Mfp L M LMterm gammaLM gammaP phase degen] % Echo for debug
                    % C(Cind,:)=[l Mf lp Mfp L -M LMterm gammaLM gammaP phase degen];   % With M>-M
                    C(Cind,:)=[l Mf lp Mfp L M LMterm gammaLM gammaP phase degen];      % for +M
                    Cind=Cind+1;
                    % end
                end
            end
        end

        indThres=abs(C(:,7))>1e-10;  % Threshold values to discard small terms
        % Cthres=C(indThres,:);
        Cthres = C;
        
        % Sum over (L,M) terms
        clear BLM;
        BLMind = 1;
        for L=0:max(Cthres(:,5))
            for M=-L:L
                maskLM = (Cthres(:,5)==L).*(Cthres(:,6)==M);
                BLM(BLMind,:) = [L M sum(maskLM.*Cthres(:,7))];
                BLMind = BLMind+1;
            end
        end
        
        % Log parameters to output structure - NOTE CURRENTLY ASSUMES ONLY ONE IP COMPONENT
        % calcsAll(n).Lmax=Lmax;
        % Xsect(indSymm,indE)=sum(abs(Csort(:,3)).^2)*rlAll(indSymm,indE).MbNorm(ip);  % Total X-sect as sum |partial waves|^2
        calcsAll(indSymm,indE).eKE=rlAll(indSymm,indE).eKE;
        calcsAll(indSymm,indE).symm=rlAll(indSymm,indE).symm;
        calcsAll(indSymm,indE).euler=[pRot tRot cRot];
        % calcsAll(indSymm,indE).Xsect=BLM(1,3);              % X-sect from B00 term
        % calcsAll(indSymm,indE).D=Dlmmu*sf;
        % calcsAll(indSymm,indE).XsectD=sum(sum(abs(Dlmmu*sf).^2))./(res.^2);      % Total X-sect as sum over MFPAD, renormalized by number of points.
        % calcsAll(indSymm,indE).Rlf=Rlf;
        calcsAll(indSymm,indE).p=p;
        calcsAll(indSymm,indE).BLM=BLM;
        calcsAll(indSymm,indE).Cthres=Cthres;
    end
end

%% Sum over symmetries

%*** Set final row in output to sum over symmetries: MAY NOT BE correct, since some interferences may be missing in BLM summation.
% for n=1:nEnergies
%     % calcsAll(nSymms+1,n).Dlmmu=D(:,:,n);    
%     calcsAll(nSymms+1,n).eKE=calcsAll(1,n).eKE;
%     calcsAll(nSymms+1,n).symm='Sum';
%     calcsAll(nSymms+1,n).euler=[pRot tRot cRot];
%     calcsAll(nSymms+1,n).p=p;
%      
%     % Sum BLMs - currently set assuming consistent size of symmetries, may not be the case (CHECK OLD CODES HERE for better routines)
%     BLMsum = calcsAll(1,n).BLM;
%     BLMsum(:,3)=0;
%     for m=1:nSymms
%         BLMsum(:,3) = BLMsum(:,3) + calcsAll(nSymms,n).BLM(:,3);
%     end
%     
%     calcsAll(nSymms+1,n).BLM = BLMsum;
%     
% end