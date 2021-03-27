% function D=DJMM(J,M1,M2,phi,theta,chi)
%
% Calculate Wigner D function for (quantized) frame rotations.  
%                  D=e^(-i*phi*M1)*dJMM*e(-i*chi*M2)
% See Angular Momentum by Zare, eqn 3.54
%
% INPUTs: J, M1, M2 are scalars
%         (phi, theta, chi) are single valued, or arrays of the same size.
%
% 13/04/16  ePSproc version for release, see notes below
% 24/06/10  Version 1, dJMM function adapted from old C code.
% Oct. 2005 C version, adapted from original Zare group Fortran code, 3jlib.for
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
%       https://epolyscat.droppages.com)
%   - F. A. Gianturco, R. R. Lucchese, and N. Sanna, J. Chem. Phys. 100, 6464 (1994).  
%       http://dx.doi.org/10.1063/1.467237
%   - A. P. P. Natalense and R. R. Lucchese, J. Chem. Phys. 111, 5344 (1999). 
%       http://dx.doi.org/10.1063/1.479794
%

function D=ePSproc_wignerD(J,M1,M2,phi,theta,chi)

d=dJMM(J,M1,M2,theta);

D=exp(-1i*phi*M1).*d.*exp(-1i*chi*M2);




%*********************************
%**
%**  dJMM(J,M1,M2,theta)
%**
%**  Calculate (reduced) rotation matrix element d(J)M1M2(theta)
%**  Theta is in radians
%**
%**  See Zare Angular Momentum eqn. 3.57
%**
%**  Adapted into C from Zare group Fortran code, 3jlib.for, P. Hockett Oct. 2005
%**  Adapted into Matlab from C code, P. Hockett 24/06/10
%**
%*********************************

function d=dJMM(J,M1,M2,theta)

sumV=zeros(size(theta));
%condition=1;

prefactor=(factorial(J+M2)*factorial(J-M2)*factorial(J+M1)*factorial(J-M1))^0.5;	% Calculate prefactor

if ((M1-M2)<0) 
    startV=abs(M1-M2);	%Set summation limits
else
    startV=0;
end

if ((J-M1)>(J+M2)) 
    endV=J+M2;
else
    endV=J-M1;
end

%loops=0;

for V=startV:endV
    	
	t1=J-M1-V;		%Calculate factorial terms
  	t2=J+M2-V;
  	t3=V+M1-M2;

%     if ((t1<0)||(t2<0)||(t3<0)) 
%         condition=0;  %Stop when one of the terms goes negative, defines upper limit for V
%     end
% 
%     if(condition)

      % Calculate term in V (factorials)
      termV=((-1)^V)/(factorial(t1)*factorial(t2)*factorial(t3)*factorial(V));

      % Cosine term
      termC=cos(theta/2).^(2*J+M2-M1-2*V);
            
	  %Sine term
      termS=(-sin(theta/2)).^(M1-M2+2*V);
            
      sumV=sumV+termV.*termC.*termS;
      
%    end
    
   	%condition=1;
    %loops=loops+1;
end

  d=prefactor*sumV;

