%% WignerD and 3j definition tests
%   09/08/19
%
%   For numerically testing my usual (Zare) convenctions (hence ePSproc definitions) against Moble's python code.
%
%   Defns....
%
%     See:
%         - https://en.wikipedia.org/wiki/Wigner_D-matrix#Definition_of_the_Wigner_D-matrix
%         - https://moble.github.io/spherical_functions/WignerDMatrices.html
%         - https://moble.github.io/spherical_functions/#euler-angles
%         - https://moble.github.io/spherical_functions/SWSHs.html
%         - https://github.com/moble/spherical_functions/blob/master/Notes/conventions.ipynb
%             (Comparison with Mathematica and Sympy definitions)
% 
%     Looks like:
%         - WignerD defn. matches wikipedia (with matches Zare, eqns. 3.54-3.55)
%         - Euler angle defns. also consistent, (alpha, beta, gamma) == (phi, theta, chi)
%         - Calcs should match Matlab ePSproc_wignerD.m
%
%% Path to ePSproc scripts

ePSprocPath='D:\code\ePSproc\distro_120416\ePSproc-master\ePSproc-master'

path(path,ePSprocPath);   % Add path to ePSproc scrips to Matlab path list


%% WignerD tests & benchmarking
% Set QNs and angles

wDbench = zeros(1,4);

for Lmax = 6
    QNs=[0 0 0];
    for l=1:Lmax
        for m=-l:l
            for mp= -l:l
                QNs(end+1,:)=[l m mp];
            end
        end
    end

    % Spherical coord system
    res=50;
    [theta,phi]=meshgrid(linspace(0,pi,res),linspace(0,2*pi,res));

    % Set a range of Eugler angles for testing - define rotation of LF->MF, set as [0 0 0] for z-pol, [0 pi/2 0] for x-pol, [pi/2 pi/2 0] for y-pol
    Nangs = 1000
    pRot = linspace(0,180,Nangs);
    tRot = linspace(0,90,Nangs);
    cRot = linspace(0,180,Nangs);
    eAngs = [pRot' tRot' cRot'].*pi./180;


    % Test wignerD

    % Set of angles, all QNs
    tStart = tic;
    wD_QNs = zeros(Nangs,5);

    r = Nangs;
    N = 1;
    for n = 1:size(QNs,1)
        wD_QNs(N:(N+r-1),1:4) = [repmat(QNs(n,:),r,1) ePSproc_wignerD(QNs(n,1), QNs(n,2), QNs(n,3), eAngs(:,1), eAngs(:,2), eAngs(:,3))];
        tEnd = toc;
        wD_QNs(N:(N+r-1),5) = repmat(tEnd,r,1);
        N = N+r;
    end

    tTotal = wD_QNs(end,5)-wD_QNs(1,5)

    wDbench(end+1,:) = [Nangs Lmax size(wD_QNs,1) tTotal]

    % Running on Stimpy... doesn't even max a single core.
    % ~ 100ms over multiple runs for Lmax = 6, Nangs = 5 > 2280 rows in wD_QNs.
    % ~ 2000ms over multiple runs for Lmax = 6, Nangs = 1000 > 455454 rows in wD_QNs.
    % ~ 29100ms over multiple runs for Lmax = 10, Nangs = 1000 > 1772770 rows in wD_QNs.
end

%% Plot benchmarks
figure; 
subplot(1,2,1);
plot(wDbench(:,3),wDbench(:,4));
xlabel('Output rows');
ylabel('t/s');

subplot(1,2,2);
plot(wDbench(:,2),wDbench(:,4));
xlabel('Lmax');
ylabel('t/s');

title('WignerD ePSproc function benchmarks, 09/08/19');

%% Plot values
pRange = 1:10000;
figure; 
plot([real(wD_QNs(pRange,4)) imag(wD_QNs(pRange,4))]);

legend('Re','Im');
title('Wigner D values vs. (l,m,\Omega)');


%% 3j tests
%   Note ePSproc version is not vecorised, single set of QNs only

w3jbench = zeros(1,3);

for Lmax = 6 % 0:10
    QNs_3j = zeros(1,6);
    for l=0:Lmax
        for lp=0:Lmax
            for m = -l:l
                for mp = -lp:lp
                    for L = 0:(l+lp)
                        M = -(m+mp);
                        QNs_3j(end+1,:)=[l lp L m mp M];
                    end
                end
            end
        end
    end

    QNs_3j(1,:)=[]  % Remove extraneous 1st row
    
    % Set of angles, all QNs
    tStart = tic;
    wD_QNs = zeros(1,5);

    % N = 1;
    for n = 1:size(QNs_3j,1)
        QNs_3j(n,7) = ePSproc_3j(QNs_3j(n,1), QNs_3j(n,2), QNs_3j(n,3), QNs_3j(n,4), QNs_3j(n,5), QNs_3j(n,6));
        QNs_3j(n,8) = toc;
        % N = N+1;
    end

    tTotal = QNs_3j(end,8)- QNs_3j(1,8);
    
    w3jbench(end+1,:) = [Lmax size(QNs_3j,1) tTotal]

    % Running on Stimpy
    % ~3ms Lmax = 1 > 28 rows in output
    % ~2500ms Lmax = 6 > 21793 rows
    % ~28000ms Lmax = 10 > 212401 rows
end

%% Plot
figure; 

subplot(1,2,1);
plot(w3jbench(:,2),w3jbench(:,3));
xlabel('Output rows');
ylabel('t/s');
title('Wigner 3j ePSproc function benchmarks, 09/08/19');

subplot(1,2,2);
plot(w3jbench(:,1),w3jbench(:,3));
xlabel('Lmax');
ylabel('t/s');
title('Wigner 3j ePSproc function benchmarks, 09/08/19');


