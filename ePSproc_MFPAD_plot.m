%%  function ePSproc_MFPAD_plot(calcsAll,eRange,symmInd,params,plotLayout,titlePrefix,saveFlag,axisType)
%   Function to plot MFPADs from ePS results stucture
%  
%  INPUTS:
%   calcsAll        Structure of calculation results
%   eRange          Vector of energy indicies to plot
%   symmInd         Scalar defining symmetry index to plot
%   params          (optional) Stucture of calculation parameters, used to add molecular structure to plots
%   plotLayout      (optional) Vector [row col] defining subplot layout (default is 2x5)
%   titlePrefix     (optional) Text string defining plot title prefix - params.fileName, symmetry etc. added automatically.
%   saveFlag        (optional) Set to 's' or 'y' to save, default is figure not saved.
%   axisType        (optional) Set axis style
%
%  13/04/16     ePSproc version for release, see notes below
%  23/10/15     Added additional axis styles for plotting
%  09/10/15 v1  Code ported from "ePS_MFPADs_scripts_2015.m", tidied up a little bit
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
%

function ePSproc_MFPAD_plot(calcsAll,eRange,symmInd,params,plotLayout,titlePrefix,saveFlag,aType)

%% Parse input arguments, set defaults if not passed

% Set axis type, grid is default
if nargin<8
    aType='grid';
end

% Set saveflag
if nargin<7
    saveFlag=0;
end

% Set title text
if nargin<6
    titlePrefix='MFPADs';
end


%% Set-up plot

% Set number of panels
if exist('plotLayout','var')&&(~isempty(plotLayout))
    pRow=plotLayout(1);
    pCol=plotLayout(2);
else
    pRow=2;
    pCol=5;
end
pFig=pRow*pCol-1;     % Panels per figure, subtract 1 to use for molecule plot

% Set grid
res=size(calcsAll(symmInd,1).D,1);
[theta,phi]=meshgrid(linspace(0,pi,res),linspace(0,2*pi,res));

figTot=ceil(length(eRange)./pFig);      % Set total number of figures, == eRange./(plots per figure)

%% Loop over results and plot

for figInd=1:figTot
    panel=1;
    
    % Set eRange for figure, based on no. of panels. 
    if figInd==figTot
        ePlot=eRange(((figInd-1)*pFig+1):end);      % For last figure, set to end of range
    else
        ePlot=eRange(((figInd-1)*pFig+1):(figInd*pFig));
    end
    
    % Create figure
    figure('color',[1 1 1],'name',['MFPADs, eke=' num2str(calcsAll(symmInd,ePlot(1)).eKE,'%4.2f') '-' num2str(calcsAll(symmInd,ePlot(end)).eKE,'%4.2f') ' eV']);

    for n=ePlot     % Loop over results and plot
        
        d_plot=calcsAll(symmInd,n).D.*conj(calcsAll(symmInd,n).D);

        %[x,y,z]=sph2cart(phi, theta+pi/2,abs(d_plot.*conj(d_plot)));  % Convert to (x,y,z) form, note +pi/2 offset to allow for different axis definitions in Matlab's sph2cart function!
        %pMap=angle(d_plot.*conj(d_plot));

        [x,y,z]=sph2cart(phi, theta+pi/2,abs(d_plot));  % Convert to (x,y,z) form, note +pi/2 offset to allow for different axis definitions in Matlab's sph2cart function!
        % pMap=angle(d_plot);   % Set colour map to phase
        % pMap=imag(d_plot);    % Set colour map to imag
        pMap=abs(d_plot);       % Set colour map to abs

       % if (length(ePlot)>1)
            subplot(pRow,pCol,panel);  % Select panel, note number of panels in figure is also set here
       % end

        % surf(x,y,z,pMap);   % Plot d and use phase for colour map
        % surf(x,y,z,z);   % Plot d and use Z for colour map
        surf(x,y,z,pMap,'FaceAlpha',0.8,'FaceColor','interp','EdgeAlpha',0);   % Plot d and use abs(D*D) for colour map
        % colormap('jet');
        % This version should have alpha mapping, but currently doesn't seem to work?
        % surf(x,y,z,abs(d_plot),'AlphaDataMapping','direct','AlphaData',abs(d_plot),'FaceColor','interp','EdgeAlpha',0.13);

        axis equal;


        % if n>size(rlAll,2);
        %     title('Sum over symmetries');
        % elseif length(ePlot)<6
        if pFig<5
            title({[titlePrefix ', MF code, files ' strrep(params.fileName,'_','\_')]; ['\mu_0=' num2str(calcsAll(symmInd,n).p) ', D^1_{\mu,-\mu_0}(' num2str(calcsAll(symmInd,n).euler(1)*180/pi) ',' num2str(calcsAll(symmInd,n).euler(2)*180/pi) ',' num2str(calcsAll(symmInd,n).euler(3)*180/pi) '), \Gamma_e=' calcsAll(symmInd,n).symm ', eKE = ' num2str(calcsAll(symmInd,n).eKE,'%4.2f') ' eV']});  % Title with eKE
        else
            % title({['\Gamma_e=' rlAll(n).symm ', eKE = ' num2str(rlAll(n).eKE,'%4.2f') ' eV']});  % Title with symm & eKE only
            title({[num2str(calcsAll(symmInd,n).eKE,'%4.2f') ' eV']});  % Title with eKE only, other info added globally.
        end
        %title({['H' num2str(round((rlAll(n).eKE+IP)/1.55),'%3.0f')]});  % Title with harmonic number
        %colorbar;

        %*** Add coord frames
        % set limits & add axes
        limFactor=1.2;      % Set to give a little extra space
        if sum(sum(d_plot))>0
            lims=[-max(abs(d_plot(:))) max(abs(d_plot(:)))];
            lims=lims*limFactor;
        else
            lims=[-1 1];
        end
        
        zlim(lims);
        xlim(lims);
        ylim(lims);
        
        a=[lims(1) 0 lims(2)];  % Set points for plotting axes
        X=[a; zeros(2,3)];
        Y=circshift(X,1);
        Z=circshift(X,2);
        
        hold on;
        plot3(X(1,:),X(2,:),X(3,:),'k');
        plot3(Y(1,:),Y(2,:),Y(3,:),'k');
        plot3(Z(1,:),Z(2,:),Z(3,:),'k');
        
        %*** Add axis labels if set
        
        % Set bounds and labels, round from limits above
        limsRound=floor((lims.*(1./limFactor))./2).*2;    % Round (down) to nearest multiple of 2
        
        a=[limsRound(1) limsRound(2)];  % Set points for plotting axes
        X=[a; zeros(2,2)];
        Y=circshift(X,1);
        Z=circshift(X,2);
        
        % Set axis type - MIGHT HAVE SOME BETTER OLD CODE FOR THIS STUFF?
        switch aType
            case 'grid'     % Default, just leave grid as is
                
            case 'aLabels'  % Turn off main grid, add labels to axes
                axis off
                text(X(:),Y(:),Z(:),num2str(limsRound(2)));
                
            case 'circles'
                axis off
                               
                lineColors=gray(8);
                               
                for cInd=1:2
                    circX=(1./cInd).*limsRound(2).*cos(linspace(0,2*pi,res));
                    circY=(1./cInd).*limsRound(2).*sin(linspace(0,2*pi,res));
                    
                    ccInd=2*cInd+2;
                    plot3(circX,circY,zeros(1,res),'color',lineColors(ccInd,:));
                    plot3(circX,zeros(1,res),circY,'color',lineColors(ccInd,:));
                    plot3(zeros(1,res),circX,circY,'color',lineColors(ccInd,:));
                    
                    % Add text, z-axis only
                    text(0,0,(1./cInd).*limsRound(2),num2str((1./cInd).*limsRound(2)));
                    
                end
                
            case 'cbar'
                axis off
                
                ticks=[0 limsRound(2)./2 limsRound(2)];
                
                % colorbar('east','Ticks',ticks);
                % colorbar('Position',[0.9 0.1 0.1 0.5]);   % Set for FIGURE, not subplot
                
%                 p=get(gca,'position'); % save position
%                 colorbar;
%                 set(gca,'position',p); % restore position
                
                % Add colorbar, but get & reset figure position first to avoid resizing
                % Code adapted from http://www.mathworks.com/matlabcentral/newsreader/view_thread/20750
                % Get axis settings
                p = get(gcf,'pos');
                ax = gca;    % findobj(gcf,'type','axes');
                p0 = get(ax,'pos');
                % Add colourbar
                cb = colorbar;
                pcb = get(cb,'pos');
                % Reset axis
                set(ax,'pos',p0);
                % set(cb,'pos',[p0(1)+p0(3)+0.02 pcb(2:4)]) % move
                % set([ax cb],'units','pixels')
                % set(gcf,'pos',[p(1) p(2) p(3)*0.5 p(4)])
                % Rescale colourbar
                cbScale=0.3;
                set(cb,'pos',[pcb(1)+2*pcb(3) pcb(2)+0.5*(pcb(4)-pcb(4).*cbScale) pcb(3).*0.5 pcb(4).*cbScale])
                
        end
            
        
        %*** Add molecule & laser pol (local function) - set here to include with every plot
        % molPlotOverlay(params,lims,calcsAll(symmInd,n).Rlf,Z);


        view([50 25]);
        caxis([0 lims(2)]);     % Reset caxis to prevent molecule plotting above overriding main MFPAD surface plot too much.

        panel=panel+1;
        
%       % In case plots exceed # of panels, stop plotting (otherwise throws an error and skips plot labels below)        
%         if panel>(pRow*pCol)
%             break;
%         end
        
    end
    
    % Add info label for large panel cases (not included above)
    if length(eRange)>=6
        % annotation('textbox',[0.1 0.45 0.1 0.1],'String',{['NO_2 ePS testing, MF code, files ' strrep(params.fileName,'_','\_')]; ['\mu_0=' num2str(calcsAll(plotSymm,n).p) ', D^1_{\mu,-\mu_0}(' num2str(calcsAll(plotSymm,n).euler(1)*180/pi) ',' num2str(calcsAll(plotSymm,n).euler(2)*180/pi) ',' num2str(calcsAll(plotSymm,n).euler(3)*180/pi) '), \Gamma_e=' calcsAll(plotSymm,n).symm ', eKE = ' num2str(calcsAll(plotSymm,n).eKE,'%4.2f') ' eV']; date},'LineStyle','none');  % Title with eKE
        % Use textarrow type (normalized fig units), or text() (current axis), instead of textbox to allow for rotated text
        annotation('textarrow',[0.05 0.05], [0.9 0.9],'String',{titlePrefix; ['File: ' strrep(params.fileName,'_','\_')]; ['\mu_0=' num2str(calcsAll(symmInd,n).p) ', D^1_{\mu,-\mu_0}(' num2str(calcsAll(symmInd,n).euler(1)*180/pi) ',' num2str(calcsAll(symmInd,n).euler(2)*180/pi) ',' num2str(calcsAll(symmInd,n).euler(3)*180/pi) '), pol: ' calcsAll(symmInd,n).polLabel ', \Gamma_e=' calcsAll(symmInd,n).symm]; date},'LineStyle','none','HeadStyle','none','TextRotation',90);  % Title with eKE
        % text(0.1,0.1,{['NO_2 ePS testing, MF code, files ' strrep(params.fileName,'_','\_')]; ['\mu_0=' num2str(calcsAll(plotSymm,n).p) ', D^1_{\mu,-\mu_0}(' num2str(calcsAll(plotSymm,n).euler(1)*180/pi) ',' num2str(calcsAll(plotSymm,n).euler(2)*180/pi) ',' num2str(calcsAll(plotSymm,n).euler(3)*180/pi) '), \Gamma_e=' calcsAll(plotSymm,n).symm ', eKE = ' num2str(calcsAll(plotSymm,n).eKE,'%4.2f') ' eV']; date},'LineStyle','none','Rotation',90);
    end
    
    %*** Add molecule (local function) - set here to include as last plot only
    subplot(pRow,pCol,panel);
    molPlot(params,calcsAll(symmInd,n).Rlf,Z);
    view([50 25]);
    axis equal;
    axis off;
    
    if (saveFlag=='s')||(saveFlag=='y')
        figFileName=[params.fileName(1:end-8) '_MFPADs_' calcsAll(symmInd,n).symm '_' calcsAll(symmInd,n).polLabel '-pol_' num2str(calcsAll(symmInd,ePlot(1)).eKE,'%4.2f') '-' num2str(calcsAll(symmInd,ePlot(end)).eKE,'%4.2f') 'eV_' date '.fig']
        savefig(figFileName);
    end
end

end

%% Additional functions
%

%% Add molecule & light to plot, version for overlay on MFPADs plot
function molPlotOverlay(params, lims, Rlf, Z)
% Rotated axes
%     Xr=Rlf*X;
%     Yr=Rlf*Y;
%     Zr=Rlf*Z;

% plot3(Xr(1,:),Xr(2,:),Xr(3,:),'--');
% plot3(Yr(1,:),Yr(2,:),Yr(3,:),'--');
% plot3(Zr(1,:),Zr(2,:),Zr(3,:));

%*** Add molecule
sf=0.2*lims(2);
offset=0.5*lims(2); %-0.5*lims(2);
XYZ=sf.*cell2mat(params.coords(3:5));  % Set Cart. coords, scale to current axes
% XYZ(:,2:3)=XYZ(:,2:3)+offset;       % Add offset
% XYZ=XYZ+offset;
% XYZ(:,1)=XYZ(:,1)-offset;   % Plot on (YZ) plane
XYZ(:,3)=XYZ(:,3)-offset;       % Plot on XY plane
S=params.coords{2}; %*10;             % Set atomic number, use for marker size & colour

% plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'o');
scatter3(XYZ(:,1),XYZ(:,2),XYZ(:,3),S,S./lims(2),'filled');

% Add molecular plane - use atomic positions directly (OK for small molecules)
% patch(XYZ(1,:),XYZ(2,:),XYZ(3,:),0.5*[1 1 1])
% patch(XYZ(:,1),XYZ(:,2),XYZ(:,3),0.8*[0 0 1],'FaceAlpha',0.2);

% Add molecular plane - define planes
% Find maxima and use to define planes
Xm=[max(XYZ(:,1)) min(XYZ(:,1))];
Ym=[max(XYZ(:,2)) min(XYZ(:,2))];
Zm=[max(XYZ(:,3)) min(XYZ(:,3))];

% Plot XY plane
patch([Xm fliplr(Xm)], [Ym(1) Ym(1) Ym(2) Ym(2)], repmat(Zm(1),1,4),0.8*[0 0 1],'FaceAlpha',0.2);
% Plot XZ
patch([Xm fliplr(Xm)], repmat(Ym(1),1,4), [Zm(1) Zm(1) Zm(2) Zm(2)],0.8*[0 0 1],'FaceAlpha',0.2);
% Plot YZ
patch(repmat(Xm(1),1,4), [Ym fliplr(Ym)], [Zm(1) Zm(1) Zm(2) Zm(2)],0.8*[0 0 1],'FaceAlpha',0.2);
% axis equal;

% Add light polarization vector
Zr=Rlf*(sf*Z./(max(Z(:))));     % Rescale z-vector & rotate
Zr=Zr-offset;
% plot3(Zr(1,:)-(2.8*offset),Zr(2,:),Zr(3,:),'LineWidth',3,'Color',[1 0 0]);
plot3(Zr(1,:)+offset,Zr(2,:)+0.5*offset,Zr(3,:),'LineWidth',3,'Color',[1 0 0]);

end

%% Add molecule & light to plot, version for separate plot containing molecule only
function molPlot(params, Rlf, Z)
      
%*** Add molecule

% Set limits
lims=[-1.2 1.2].*max(max(cell2mat(params.coords(3:5))));
if abs(max(lims(:)))<=0.1   % Reset in case of null limits (for atoms)
    lims=[-1 1];
end

sf=1; %0.2*lims(2);
offset=0; %.5*lims(2); %-0.5*lims(2);
XYZ=sf.*cell2mat(params.coords(3:5));  % Set Cart. coords, scale to current axes
% XYZ(:,2:3)=XYZ(:,2:3)+offset;       % Add offset
% XYZ=XYZ+offset;
% XYZ(:,1)=XYZ(:,1)-offset;   % Plot on (YZ) plane
XYZ(:,3)=XYZ(:,3)-offset;       % Plot on XY plane
S=params.coords{2}; %*10;             % Set atomic number, use for marker size & colour

% plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'o');
scatter3(XYZ(:,1),XYZ(:,2),XYZ(:,3),S*10,S,'filled');
hold on;

% Add molecular plane - use atomic positions directly (OK for small molecules)
% patch(XYZ(1,:),XYZ(2,:),XYZ(3,:),0.5*[1 1 1])
% patch(XYZ(:,1),XYZ(:,2),XYZ(:,3),0.8*[0 0 1],'FaceAlpha',0.2);

% Add molecular plane - define planes
% Find maxima and use to define planes
Xm=[max(XYZ(:,1)) min(XYZ(:,1))];
Ym=[max(XYZ(:,2)) min(XYZ(:,2))];
Zm=[max(XYZ(:,3)) min(XYZ(:,3))];

% Plot XY plane
patch([Xm fliplr(Xm)], [Ym(1) Ym(1) Ym(2) Ym(2)], repmat(Zm(1),1,4),0.8*[0 0 1],'FaceAlpha',0.2);
% Plot XZ
patch([Xm fliplr(Xm)], repmat(Ym(1),1,4), [Zm(1) Zm(1) Zm(2) Zm(2)],0.8*[0 0 1],'FaceAlpha',0.2);
% Plot YZ
patch(repmat(Xm(1),1,4), [Ym fliplr(Ym)], [Zm(1) Zm(1) Zm(2) Zm(2)],0.8*[0 0 1],'FaceAlpha',0.2);
% axis equal;

%*** Add axes
% Set cart axes
a=[lims(1) 0 lims(2)];  % Set points for plotting axes
X=[a; zeros(2,3)];
Y=circshift(X,1);
Z=circshift(X,2);

% Plot cart axes, replace any existing plot content
plot3(X(1,:),X(2,:),X(3,:),'k');
plot3(Y(1,:),Y(2,:),Y(3,:),'k');
plot3(Z(1,:),Z(2,:),Z(3,:),'k');

text(lims(1),lims(1).*0.2,0,'x');
text(lims(1)*0.2,lims(1),0,'y');
text(0,lims(1).*0.2,lims(2),'z');


%*** Add light polarization vector
Zr=Rlf*(sf*Z./(max(Z(:))));     % Rescale z-vector & rotate
Zr=Zr-offset;
% plot3(Zr(1,:)-(2.8*offset),Zr(2,:),Zr(3,:),'LineWidth',3,'Color',[1 0 0]);
plot3(Zr(1,:)+offset,Zr(2,:)+0.5*offset,Zr(3,:),'LineWidth',3,'Color',[1 0 0]);




end
