function [s,eFlag]=ePSproc_3j(j1,j2,j3,m1,m2,m3)
%3j code
%adapted from zare, angular momentum, fortran code appendix a/theory p49
%complete rewrite using matlab code only
%
%Version 3 4/11/08 - as version 2b, but renamed to prevent confusion with older copies of v2b in matlab path
%Version 2b 1/12/06 - added error checking to set results to zero for out of bounds combinations
%Version 2 15/06/06 - all working nicely
%
%Equation for 3j (Zare p49 & p325) is split into 3 parts:
%	- polarity
%	- delta(j1,j2,j3)
%	- w3j(j1,j2,j3,m1,m2,m3) which includes summation term
%
% *** NOTES
%  This version distributed as part of ePSproc suite.
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

%error code adapted from Wigner3j.m script
if ( j1 - m1 ~= floor ( j1 - m1 ) )
    eFlag=1;
elseif ( j2 - m2 ~= floor ( j2 - m2 ) )
    eFlag=1;
elseif ( j3 - m3 ~= floor ( j3 - m3 ) )
    eFlag=1;
elseif j3 > j1 + j2 | j3 < abs(j1 - j2)
    eFlag=1;
elseif abs(m1) > j1
    eFlag=1;
elseif abs(m2) > j2
    eFlag=1;
elseif abs(m3) > j3
    eFlag=1;
elseif (m1+m2+m3)~=0    %triangle condition, could also be set as m1+m2-m3~=0 depending on form of inputs
    eFlag=1;
else
    eFlag=0;
end

if(eFlag)   %calculate (or not) based on eFlag value
    s=0;
else
    polarity=(-1)^(j1-j2-m3);

    delta=(factorial(j1+j2-j3)*factorial(j1-j2+j3)*factorial(-j1+j2+j3)/factorial(j1+j2+j3+1))^(1/2);
    w3ja=(factorial(j1+m1)*factorial(j1-m1)*factorial(j2+m2)*factorial(j2-m2)*factorial(j3+m3)*factorial(j3-m3))^(1/2);

    s=0;

    %find min of summation, zero if both other terms are negative

    sumMin = max( 0, max( -1*(j3-j1-m2), -1*(j3-j2+m1) ) );

    %find max of summation
    sumMax = min( j1+j2-j3, min( j2+m2, j1-m1) );


    for i=sumMin:sumMax
        w3jb=factorial(i)*factorial(j1+j2-j3-i)*factorial(j1-m1-i)*factorial(j2+m2-i)*factorial(j3-j2+m1+i)*factorial(j3-j1-m2+i);
        w3jb=((-1)^i)/w3jb;
        s=s+w3jb;
    end

    s=polarity*delta*w3ja*s;
end



