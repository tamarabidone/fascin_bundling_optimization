% FIND OPTIMUM ALPHA VALUES FOR EACH FILAMENT
% All distances in nanometers
% Model Parameters (MP) --------------------------------------------------------------
clear
clc

MP.TwistTheta  =  -166.7;  % degrees
MP.MonLength   =    2.70;  % nm
MP.FilSep      =   12.15; %12.15; % nm  STD = 0.06  Filament axis separation distance
MP.R.ASC  = 1.58;         % nm  STD = 0.06  Actin subunit centroid radius
MP.R.S1C  = 2.63;         % nm  STD = 0.06  Subdomain 1 centroid radius
MP.R.D24  = 3.62;
MP.R.R95a = 3.64;
MP.R.S350 = 3.66;
MP.R.R95b = 3.64;

MP.nMonomers = 17;  % monomers
MP.BondThreshold = 0.4; % 0.4

% Binding probability functions ---------------------------------------------------------------------------------------
SDM = 4; % Standard Deviation Multiplier
MP.BPF.ASC_im1_j_bp   = @(D) exp( -0.5*((D-14.97)/(2)).^2 ); % Actin Subunit Centroid (ASC) i-1 to j binding probability function 
MP.BPF.ASC_ip3_jp4_bp = @(D) exp( -0.5*((D-15.36)/(2)).^2 ); % Actin Subunit Centroid (ASC) i+3 to j+4 

% Gaussian distribution bond probability functions ----------------------------------------------------
MP.BPF.R95a_R95b_bp   = @(D) exp( -0.5*((D- 5.80)/(SDM*0.069)).^2 ); % R95a to R95b
MP.BPF.R95a_S350_bp   = @(D) exp( -0.5*((D- 8.08)/(SDM*0.133)).^2 ); % R95a to S350
MP.BPF.D24_S350_bp    = @(D) exp( -0.5*((D- 6.07)/(SDM*0.170)).^2 ); % D24 to S350
MP.BPF.D24_R95b_bp    = @(D) exp( -0.5*((D- 5.34)/(SDM*0.078)).^2 ); % D24 to R95b
MP.BPF.FascinCent_bp  = @(D) exp( -0.5*((D- 11.287)/(SDM*3.0)).^2 ); % Distance of fascin centroid to plane defined by two fascin 

 
% FILAMENT ARRANGEMENTS ===============================================================
%===== Create adjacent bundle array of points ------------------------
%         [~,X,Y] = AttachHexagons01(4,MP.FilSep);
%         Y = Y(X<87);
%         X = X(X<87);
%         Z = MP.MonLength*(0:MP.nMonomers-1)';
%         N = length(X);
%         MP.FirstFilamentIndex = 13; % Select which filament to start with. 
%         PlotPoints01(X,Y); drawnow;

%===== Create hexagonal array of points ------------------------------
        nPts = 7; % Nummber of filaments in hexagon arrangment
        [X,Y] = CreateHexagonArray(nPts,MP.FilSep);
        Z = MP.MonLength*(0:MP.nMonomers-1)';
        N = length(X); 
        MP.FirstFilamentIndex = 1; % Select which filament to start with.
        PlotPoints01(X,Y); drawnow;

%===== Create line of points -----------------------------------------
%         nPts = 2;
%         [X,Y] = CreateLineOfPoints(nPts,MP.FilSep);
%         N = length(X); 
%         Z = MP.MonLength*(0:MP.nMonomers-1)';
%         MP.FirstFilamentIndex = 1; % Select which filament to start with.
%         MP.FirstPairIndex = 1;
%         PlotPoints01(X,Y); drawnow;
%
%=========================================================================================

%-------------------------------------------------------------------------------------------------------------
MP.nFilaments  = N;
MP.FilCoords.X = X;
MP.FilCoords.Y = Y;
MP.FilCoords.Z = Z;
%-------------------------------------------------------------------------------------------------------------  
MP.FilamentsStart = repmat( mod((0:MP.TwistTheta:(MP.nMonomers-1)*MP.TwistTheta)',360)', [1, MP.nFilaments]);

% Alpha steps for the two filaments
MP.Alpha1 = 0:1:359;
MP.Alpha2 = 0:1:359;
% Alpha steps for each additional filament
MP.Alpha   = 0:0.1:359.9;

% Create all point pair combinations
c = nchoosek(1:MP.nFilaments,2); 
% Now filter out the combinations that do are not a distance FilSep apart
Dist = sqrt( ( X(c(:,2))-X(c(:,1)) ).^2 + ( Y(c(:,2))-Y(c(:,1)) ).^2 );
ind = find( Dist < (MP.FilSep + 0.001) );
C = c(ind,:);
nPairs = size(C,1);

MP.FilPairCom = C;
MP.nPairs = nPairs;


%=============================================================================================
% RUN ALPHA OPTIMIZATION MODEL FOR FIRST PAIR ================================================
%=============================================================================================
% Run with filament arrangement: "Create line of points" nPts = 2 (lines 52-58 above)
%     CreatePairPlot = true;
%     HeatMapPlot = true;
%     PMO = CalculateOptimalAphasForFilamentPair( MP, CreatePairPlot, HeatMapPlot);

%=============================================================================================
% RUN ADD-ONE-FILAMENT-AT-A-TIME MODEL =======================================================
%=============================================================================================
% Run with filament arrangement: "Adjacent bundle" (lines 35-41 above) or "hexagonal array" (lines 44-49 above)

for SA = 1:10 
    MP.StartAlpha = randi([0,30],1);
    CreateFilamentPlot = true;
    MO = OptimizingOneFilamentAtATime02( MP, CreateFilamentPlot); % For rotating FirstFilamentIndex a fixed alpha
end






