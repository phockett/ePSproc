% function I=ePSproc_Ylm_calc(a,theta,phi)
% Function to calculate spherhical harmonic funcions YLM.
% Based on older code from "velocity_dist_3Df.m" and notes "photoelectron wavefunction plotting 16 02 09.m", also "Ylm_sum.m" and "pad_sph.m".
% 
% INPUT: a=[L M A; ...] to define harmonic(s) to calculate (and sum), where A is the (complex) magnitude
%        (theta,phi) are points, vectors or arrays defining co-ord system for calculation.  Theta and phi must be same size and dimensionality.
%
% OUTPUT: I 2D array of calculation result (complex)
%         
% 13/04/16  ePSproc version for release, see notes below
% 24/06/10  adapted from Ylm_sum to take generalized input from ePolyScat results, and all forms of angular input
% 01/09/09  Version 1 - note that plotting functions are included, but only implemented by uncommenting.
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


function I=ePSproc_Ylm_calc(a,theta,phi)

I=zeros(size(theta));  %Init I

for n=1:size(a,1)
    L=a(n,1);
    M=abs(a(n,2));
    Ms=a(n,2);  %Ms preserves sign of M.
    t=legendre(L,cos(theta)); 
    if L~=0
       t2=squeeze(t(M+1,:,:)); %M+1 for array offset
    else
       t2=t;        %Don't squeeze for L=0
    end
    a1=((2*L+1)/(4*pi));
    a2=factorial(L-M)/factorial(L+M);
    C=sqrt(a1*a2)*a(n,3);
    if Ms<0                 %Phase for M<0, from Yl,-|m|=(-1)^|m|*Ylm. Yl|m| only calculated above due to defn of Legendre function, not Condon-Shortley phase (-1)^|m| already incorporated in Legendre function.
        C=C*(-1)^M;
    end
    
    I=I+C*t2.*exp(1i*(Ms*phi)); % Calculate Ylm and add to sum
    t=[];
    t2=[];
end