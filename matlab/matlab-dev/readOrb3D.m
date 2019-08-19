% [Orb, X,Y,Z]=readOrb3D(file)
% Function to read in grid data from ePolyScat "ViewOrb" & sort to f(x,y,z)
%
% Input:    file    path for reading
% Outputs:  Orb     structure of fields from file
%           [X,Y] cartesian grid
%
% 17/04/17  3D  For 3D ePS outputs. [NOTE TO SELF - check back on original readMO.m for alternatives?]
% 14/02/11  v1  Based on "readMO.m" script for reading MacMolPlot outputs
%
% Paul Hockett
% paul.hockett.nrc.ca

function [Orb,X,Y,Z]=readOrb3D(file)

% Read in data file
fid=fopen(file);

Orb.file=file;

for page=1:3    % Repeat read in and processing for 3 components of data from ePS - Re, Im and Abs
    title=textscan(fid,'%s',1,'delimiter','\n');                 % Title row
    Orb(page).title=title{1};
    
    grid1=(textscan(fid,'%n',1))   % Dimensions
    grid2=textscan(fid,'%n',1)   % Data parameters
    grid3=textscan(fid,'%n',1)   % Data parameters
    grid4=textscan(fid,'%n',1)   % Data parameters
    
    % grid1=grid1{1};
    % Orb(page).grid=[abs(grid1{1}) grid2{1} grid3{1} grid4{1}];
    Orb(page).grid=[grid2{1} grid3{1} grid4{1}];
    
    gridX=textscan(fid,'%f',grid2{1});   % x-axis
    gridY=textscan(fid,'%f',grid3{1});   % y-axis
    gridZ=textscan(fid,'%f',grid4{1});   % y-axis
    Orb(page).gridX=gridX{1};
    Orb(page).gridY=gridY{1};
    Orb(page).gridZ=gridZ{1};
    
    dataIn=textscan(fid,'%f',grid2{1}*grid3{1}*grid4{1});  % Read data
    % temp=textscan(fid,'%f',1);                                        % Skip to next line ready for next segment of file
    %fclose(fid);

    % Re-sort data into 3D matrix form
    C=1;
    %Orb.data=zeros(Orb.gridSize{3},Orb.gridSize{2},Orb.gridSize{1});
    %Orb(page).data=zeros(grid2{1},grid3{1},grid4{1});
    Orb(page).data=zeros(grid2{1},grid3{1},grid4{1});
    % Orb(page).data=zeros(grid3{1},grid2{1},grid4{1});

    % data=[dataIn{:}].'; % Transpose for correct ordering of linear indexing by column
   


    % Variation A
    for z=1:grid4{1}
        for y=1:grid3{1}   
         Orb(page).data(1:grid2{1},y,z)=dataIn{1}(C:(C+grid2{1}-1));  % Read in x-data line by line using linear indexing.  Is there a neater way to do this?
         C=C+grid2{1}; 
        end
    end
    
%     % Variation B
%     for y=1:grid3{1}
%         for z=1:grid4{1}   
%          Orb(page).data(1:grid2{1},y,z)=dataIn{1}(C:(C+grid2{1}-1));  % Read in x-data line by line using linear indexing.  Is there a neater way to do this?
%          C=C+grid2{1}; 
%         end
%     end

%     % Variation C
%     for z=1:grid4{1}
%         for y=1:grid3{1}   
%          Orb(page).data(y,1:grid2{1},z)=dataIn{1}(C:(C+grid2{1}-1));  % Read in x-data line by line using linear indexing.  Is there a neater way to do this?
%          C=C+grid2{1}; 
%         end
%     end

    textscan(fid,'%s',1,'delimiter','\n');  % Shift pointer to next line, seems to be necessary here for some reason!
    
end

fclose(fid);

% Reformat structure into complex form
Orb(1).data=Orb(1).data+1i*Orb(2).data;
Orb(1).dataAbs=Orb(3).data;

Orb(3)=[];  % Delete unnecessary components
Orb(2)=[];

% Set X,Y grid data values
% [X,Y,Z]=meshgrid(Orb.gridX,Orb.gridY,Orb.gridZ);  % Set values with meshgrid, (x,y,z) convention to match sort code above.
% [Y,X,Z]=meshgrid(Orb.gridY,Orb.gridX,Orb.gridZ);  % Set values with meshgrid, (x,y,z) convention to match sort code above.
  
% Set grid, convert to Cart if necessary, assuming that grid won't be larger than 10 Angs
if max(Orb.gridY>50)
   [THETA,RHO,PHI]=meshgrid(Orb.gridY,Orb.gridX,Orb.gridZ);  % Set values with meshgrid, (x,y,z) convention to match sort code above.
   [X,Y,Z] = sph2cart(PHI.*pi./180,(THETA.*pi./180)-pi/2,RHO);  % Convert to Cart, sph2cart uses convention (az,el,R).
    
else
   % [Y,X,Z]=meshgrid(Orb.gridY,Orb.gridX,Orb.gridZ);  % Set values with meshgrid, (x,y,z) convention to match sort code above.
      
   % (X,Y,Z) gridding causes an issue in some cases, with error "Input grid is not a valid MESHGRID.", so regenerate coords
   %  AH - problem actually seemed to be with switching X and Y coords... although no idea why.  Now fixed below, although might need switching again later!
   [X,Y,Z]=meshgrid(Orb.gridX,Orb.gridY,Orb.gridZ);  % Set values with meshgrid.
end
    