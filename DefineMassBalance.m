


function [UserVar,as,ab,dasdh,dabdh]=DefineMassBalance(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)
% Defines mass balance along upper and lower ice surfaces.
%
%   [UserVar,as,ab]=DefineMassBalance(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)
%
%   [UserVar,as,ab,dasdh,dabdh]=DefineMassBalance(UserVar,CtrlVar,MUA,CtrlVar.time,s,b,h,S,B,rho,rhow,GF);
%
%   as        mass balance along upper surface 
%   ab        mass balance along lower ice surface
%   dasdh     upper surface mass balance gradient with respect to ice thickness
%   dabdh     lower surface mass balance gradient with respect to ice thickness
%  
% dasdh and dabdh only need to be specified if the mass-balance feedback option is
% being used. 
%
% In Úa the mass balance, as returned by this m-file, is multiplied internally by the local ice density. 
%
% The units of as and ab are water equivalent per time, i.e. usually
% as and ab will have the same units as velocity (something like m/yr or m/day).
%
% Examples:  
%
% To set upper surface mass balance to zero, and melt along the lower ice
% surface to 10 over all ice shelves:
%
%   as=zeros(MUA.Nnodes,1);
%   ab=-(1-GF.node)*10 
%
%
% To set upper surface mass balance as a function of local surface elevation and
% prescribe mass-balance feedback for the upper surface:
%
%   as=0.1*h+b;
%   dasdh=zeros(MUA.Nnodes,1)+0.1;
%   ab=s*0;
%   dabdh=zeros(MUA.Nnodes,1);
%
% 
%%



as=zeros(MUA.Nnodes,1);
%ab=-(1-GF.node)*10;
dasdh=0;
dabdh=0;

%% Mass balance measurements
%
%   xm (km)            ym (km)             asm  (m/yr)
%  -98               -100                -1
%    0                 0                 +1
%   -67                0                  0
%    89                50                 3
%%

% x=MUA.coordinates(:,1) ;  
% y=MUA.coordinates(:,2) ; 


% Fas=scatteredInterpolant(xm,ym,asm) ;
% as=Fas(x,y) ; 




% time=CtrlVar.time; 
%MUA=load('AntarcticMesh_msh2d_msh2d.mat');

 x=MUA.coordinates(:,1) ;  
 y=MUA.coordinates(:,2) ; 
% 
% ftime=string(time);
% fpre=string("Interpolant_racmo_") ;
% adj=string(".mat")
% file=append(fpre,ftime)
% file=append(file,adj)
% 
 S=load('RACMO_climate_1995_2014_interpolant.mat');
 Fsmb=S.Fsmb;

 
as=Fsmb(x,y);
 
% figure; PlotMeshScalarVariable([],MUA,as);
%%%%%%%%%%%%%%%%%
%Define AB

%MUA=data.MUA;
%AGlen=data.AGlen;
%GF=data.GF;
%CtrlVar=data.CtrlVar;
GLgeo=GLgeometry(MUA.connectivity,MUA.coordinates,GF,CtrlVar);
%TRI=[]; DT=[]; 
xGL=[] ; yGL=[] ;
x=MUA.coordinates(:,1);  y=MUA.coordinates(:,2);

  CtrlVar.GLsubdivide=1
  [GF,GLgeo,GLnodes,GLele]=IceSheetIceShelves(CtrlVar,MUA,GF);
   
 % figure
 % PlotMuaMesh(CtrlVar,MUA,GF.ElementsUpstreamOfGroundingLines,'color','k')
 % hold on
 % PlotMuaMesh(CtrlVar,MUA,GF.ElementsDownstreamOfGroundingLines,'color','b')
 % PlotMuaMesh(CtrlVar,MUA,GF.ElementsCrossingGroundingLines,'color','r')

 %ab=-(1-GF.node)*10;
 
 %%apply constant melt rate to the ice shelves
 
 UnifrmBslMltRt=''; %'-apply-;
 if contains(UnifrmBslMltRt,'-apply-')
 %ab=-GF.NodesDownstreamOfGroundingLines*10;
 end
 
 
 %apply meltrate from PICO
 
 
 BslPICO='-apply-';
 if contains(BslPICO,'-apply-')
load MeshBoundaryCoordinatesForAntarcticaBasedOnBedmachine.mat; %MeshBoundaryCoordinates.mat
load BasinsInterpolant.mat;

%% Run PICO:

% the PICO_opts structure can be defined by the user to change from default
% values (which are set in PICO_DefaultParameters)

PICO_opts = struct;
PICO_opts.algorithm = 'watershed';%'polygon','oneshelf';
PICO_opts.C1 = 1e6; 
PICO_opts.gamTstar = 2e-5;
PICO_opts.nmax = 5;
PICO_opts.minArea = 2e9;
PICO_opts.minNumShelf = 20;
PICO_opts.SmallShelfMelt = 0;
PICO_opts.PICOres = 10000; % resolution in km (for watershed algorithm only)
PICO_opts.BasinsInterpolant = Fbasins;
PICO_opts.FloatingCriteria = 'GLthreshold'; %'GLthreshold' or 'StrictDownstream'
PICO_opts.persistentBC = 0;
PICO_opts.InfoLevel = 0; %100; % 0,1,10,100

% these two vectors are the salinity and temperature of each basin in the BasinsInterpolant file - provided by Ronja
PICO_opts.Sbasins = [34.6505;34.5273;34.3222;34.3259;34.3297;34.5315;34.4819;34.5666;34.5766;34.6677;34.7822;34.6254;34.4107;34.5526;34.6902;34.6668;34.5339;34.5849;34.6644];
PICO_opts.Tbasins = [-1.75725;-1.65931;-1.58212;-1.54757;-1.51301;-1.72267;-1.69117;-0.67609;-1.61561;-1.30789;-1.83764;-1.5798;-0.368398;0.45505;1.04046;1.17196;0.229878;-1.23091;-1.79334];

PICO_opts.MeshBoundaryCoordinates = Boundary; %MeshBoundaryCoordinates;

tic

%[Mk,ShelfID,T0,S0,Tkm,Skm,q,PBOX,Ak] = PICO_driver(1,CtrlVarInRestartFile,MUA,GF,F.h,median(F.rho),F.rhow,PICO_opts);
[Mk,ShelfID,T0,S0,Tkm,Skm,q,PBOX,Ak] = PICO_driver(1,CtrlVar,MUA,GF,h,median(rho),rhow,PICO_opts);

ab=Mk;
toc
 end

end