clear all
clc

%% Header
% This program creates a blockMeshDict file to create an axisymmetric wedge
% mesh for reproducing the isentropic expansion V&V case in Stuart and
% Lumpkin (2012).
%
% Jonathan S. Pitt, Aegis Aerospace, Inc., for LENS
% Created: 21 Sep 2021
% Updated: 10 Nov 2021
%

%% Inputs
% The domain geometry
% Note, the domain is in the x-y plane with z being the wedge direction
rstar = 1.05;   % Radius of the inflow surface,  [m]
L = 20.0;         % Radius of the exterior of the domain, [m]
thetaW = 5.0;   % Wedge angle thickness, [deg]

% The mesh parameters - there are only 3
dr = 0.5;        % Radial grid spacing [m]
dt = 0.05;       % Inflow grid spacing [m]
rexp = 10.0;      % Expansion factor for radial mesh - use the eye test

% compute the parameters for the block definition from the mesh inputs
Nr = floor(L/dr)+1;             % Number of radial points
Nt = floor(2*pi*rstar/4/dt)+1;  % Number of angular points


%% Compute the vertex points

% Chamber Points
x(1) = rstar;   y(1) = 0;   z(1) = 0;
x(2) = 0; y(2) = rstar*cosd(thetaW/2);   z(2) = -rstar*sind(thetaW/2); 
x(3) = 0; y(3) = rstar*cosd(thetaW/2);   z(3) = rstar*sind(thetaW/2); 
x(4) = 0; y(4) = L*cosd(thetaW/2);   z(4) = -L*sind(thetaW/2); 
x(5) = 0; y(5) = L*cosd(thetaW/2);   z(5) = L*sind(thetaW/2); 
x(6) = L;   y(6) = 0;   z(6) = 0;



%% Write the blockMeshDict file
fid = fopen('blockMeshDict','w');

% write the header
fprintf(fid,'/*--------------------------------*- C++ -*----------------------------------*\\');
fprintf(fid,'\n');
fprintf(fid,'| =========                 |                                                 |\n');
fprintf(fid,'| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n');
fprintf(fid,'|  \\    /   O peration     | Version:  v2106                                 |\n');
fprintf(fid,'|   \\  /    A nd           | Website:  www.openfoam.com                      |\n');
fprintf(fid,'|    \\/     M anipulation  |                                                 |\n');
fprintf(fid,'\\*---------------------------------------------------------------------------*/');
fprintf(fid,'\n');
fprintf(fid,'FoamFile\n');
fprintf(fid,'{\n');
fprintf(fid,'\t version     2.0;\n');
fprintf(fid,'\t format      ascii;\n');
fprintf(fid,'\t class       dictionary;\n');
fprintf(fid,'\t object      blockMeshDict;\n');
fprintf(fid,'}\n');
fprintf(fid,'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n');
fprintf(fid,'scale 1;\n\n');
fprintf(fid,'convertToMeters 1.0;\n\n');
fprintf(fid,'vertices\n(\n');

% write the points
for n=1:length(x)    
    fprintf(fid,'\t(%e\t%e\t%e)\n',[x(n) y(n) z(n)]);
end
fprintf(fid,');\n\n');

% write the block - there is only one block
fprintf(fid,'blocks\n(\n');

% chamber block
fprintf(fid,['\thex (0 1 2 0 5 3 4 5) (' int2str(Nt) ' 1 ' int2str(Nr) ') simpleGrading (1 1 ' sprintf('%.4f', rexp) ')\n']);
fprintf(fid,');\n\n');

%fprintf(fid,['\thex (0 2 1 0 3 5 4 3) (' int2str(Nr) ' 1 ' int2str(LcNz) ') simpleGrading (' sprintf('%.4f', expR) ' 1 1)\n']);

% % write the edges
fprintf(fid,'edges\n(\n');
% 

% vacuum arc points
arc1x = L*cosd(45);
arc1y = L*sind(45)*cosd(thetaW/2);
arc1z = L*sind(45)*sind(thetaW/2);

fprintf(fid,['\tarc 4 5 (' sprintf('%.4f', arc1x) ' ' sprintf('%.4f', arc1y) ' ' sprintf('%.4f', arc1z) ' )\n']);
fprintf(fid,['\tarc 3 5 (' sprintf('%.4f', arc1x) ' ' sprintf('%.4f', arc1y) ' ' sprintf('%.4f', -arc1z) ' )\n']);

% inflow arc points
arc2x = rstar*cosd(45);
arc2y = rstar*sind(45)*cosd(thetaW/2);
arc2z = rstar*sind(45)*sind(thetaW/2);

fprintf(fid,['\tarc 2 0 (' sprintf('%.4f', arc2x) ' ' sprintf('%.4f', arc2y) ' ' sprintf('%.4f', arc2z) ' )\n']);
fprintf(fid,['\tarc 1 0 (' sprintf('%.4f', arc2x) ' ' sprintf('%.4f', arc2y) ' ' sprintf('%.4f', -arc2z) ' )\n']);

fprintf(fid,');\n\n');
% 
% % write the boundaries
fprintf(fid,'boundary\n(\n');
% 
% inflow
fprintf(fid,'\tinflow\n\t{\n\t\ttype patch;\n\t\tfaces\n\t\t(\n\t\t\t');
fprintf(fid,'(0 1 2 0)');
fprintf(fid,'\n\t\t);\n\t}\n');
% 
% vacuum
fprintf(fid,'\tvacuum\n\t{\n\t\ttype patch;\n\t\tfaces\n\t\t(\n\t\t\t');
fprintf(fid,'(5 3 4 5)');
fprintf(fid,'\n\t\t);\n\t}\n');
% 
% x min symmetry plane
fprintf(fid,'\tsym\n\t{\n\t\ttype symmetryPlane;\n\t\tfaces\n\t\t(\n\t\t\t');
fprintf(fid,'(1 3 4 2)');
fprintf(fid,'\n\t\t);\n\t}\n');
%
% wedge positive z
fprintf(fid,'\tfront\n\t{\n\t\ttype wedge;\n\t\tfaces\n\t\t(\n\t\t\t');
fprintf(fid,'(0 2 4 5)');
fprintf(fid,'\n\t\t);\n\t}\n');
% 
% wedge negative z
fprintf(fid,'\tback\n\t{\n\t\ttype wedge;\n\t\tfaces\n\t\t(\n\t\t\t');
fprintf(fid,'(0 1 3 5)');
fprintf(fid,'\n\t\t);\n\t}\n');
% 
% wedge axis pole
fprintf(fid,'\taxis\n\t{\n\t\ttype empty;\n\t\tfaces\n\t\t(\n\t\t\t');
fprintf(fid,'(0 5 5 0)');
fprintf(fid,'\n\t\t);\n\t}\n');
% 
fprintf(fid,');\n\n');
% 
% % write the merged patch pairs
fprintf(fid,'mergePatchPairs\n(\n');
% 
% 
fprintf(fid,');\n\n');
% 
fprintf(fid,'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //');
fclose(fid);


