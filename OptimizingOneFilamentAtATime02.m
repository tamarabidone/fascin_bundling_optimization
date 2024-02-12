function MO = OptimizingOneFilamentAtATime02(MP,CreateFilamentPlot)


% Model parameters ----------------------------------------------------------------------------------------------------------
MO = [];

TwistTheta      = MP.TwistTheta;  % degrees
MonLength       = MP.MonLength;   % nm
FilSep          = MP.FilSep;      % nm  STD = 0.06, Filament axis separation distance
R.ASC  = MP.R.ASC;    % nm  STD = 0.06  Actin subunit centroid radius
R.S1C  = MP.R.S1C;    % nm  STD = 0.06  Subdomain 1 centroid radius
R.D24  = MP.R.D24;
R.R95a = MP.R.R95a;
R.S350 = MP.R.S350;
R.R95b = MP.R.R95b;

nMonomers = MP.nMonomers;      % monomers
BondThreshold = MP.BondThreshold;
% Grab Pair points -----------------------------------------------------------------------------------

N = MP.nFilaments;
nPts = N;
AllPairs = MP.FilPairCom;
nPairs = MP.nPairs;
nMonomers = MP.nMonomers;
X = MP.FilCoords.X;
Y = MP.FilCoords.Y;
Z = MP.FilCoords.Z;

%-----------------------------------------------------------------------------------------------------------------------------  
                      
                % Create initial filaments orientations ----------------------------------------------------------------------  
                    FilamentsStart = repmat( mod((0:MP.TwistTheta:(MP.nMonomers-1)*MP.TwistTheta)',360), 1, N);
                % Initialized filament alphas to zero
                    Alphas = MP.Alpha'; 
                    AlphaIdx = NaN(nPairs,1);
                    FinalAlphas = zeros(size(X));
                    nAlpha = length(Alphas);
                % Preallocate space to record binding probability for each pair and all     
                    PairBondStrength  = NaN(nAlpha,nPairs);
                    BondStrengthUp    = zeros(nPairs,nMonomers,nAlpha);
                    BondStrengthDn    = zeros(nPairs,nMonomers,nAlpha);
                    MonomersTaken     = false(nMonomers,N);
                    FilamentOptimized = false(size(X));
                    FilamentAlphas    = zeros(size(X));
                    FilamentOrder     = NaN(N,1);
                    Iteration = 0;
                    
                    
        
                    while any(~FilamentOptimized) % While there are still filaments that are not optimized, keep looping through filament pairs in proximity to FilIdx
                            Iteration = Iteration + 1;
                            if Iteration == 1 
                                StartFilament = MP.FirstFilamentIndex;       % <<<<<<<<<<  set start filament
                                FinalAlphas(StartFilament,1) = MP.StartAlpha; %P MO.Peakvalues(1); %0;     % <<<<<<<<<<  set start filament's initial orientation (from CalculateOptimalalphasForFilamentPair)
                                FilamentOptimized(StartFilament,1) = true;
                                FilamentOrder(1,1) = StartFilament;
                            end
        
                            % Get distance between optimzed and non-optimzed filaments (This is the start of calculating the closest/adjacent non-optimzed filaments to the optimized filaments
                            OptimX = X;     OptimX(~FilamentOptimized,1) = NaN;
                            OptimY = Y;     OptimY(~FilamentOptimized,1) = NaN;
                            NonOpX = X;     NonOpX( FilamentOptimized,1) = NaN;
                            NonOpY = Y;     NonOpY( FilamentOptimized,1) = NaN;
                            
                            % Calculate distances between optimized and non-optimzed filaments, and set distances for non-optimzed filaments that are more 
                            % than a distance FilSep away to NaN (since we want the filaments that are next to the current group of optimzed filaments).
                            Distances  = sqrt( (OptimX-NonOpX').^2 + (OptimY-NonOpY').^2  ); 
                            Distances( Distances > (FilSep + 0.1)) = NaN;  
                            % Find the indices of all the filaments that are adjacent to the current optimzed bundle
                            [~,AdjacentFilaments] = find( any(~isnan(Distances),1) );
                            
        
                   %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                            
                            nAdjPairs = 1;
                            nLoops = 0;

                                while nAdjPairs < 2 
                                    % Randomly grab one of the adjacent filaments so we can optimze its alpha relative to any optimzed filaments in proximity to it.
                                    FilIdx = AdjacentFilaments(1, randi([1,length(AdjacentFilaments)]) );
                                    FilamentOrder(Iteration+1,1) = FilIdx;
                                    % Find the optimzed filaments in immediate proximity to FilIdx
                                    [AdjFil,~] = find( ~isnan(Distances(:,FilIdx)) ); 
                                    % Create filament pairs that need to be optimized based on filaments adjacent (AdjFil) to the current filament selected for optimizing (FilIdx)
                                    CurrentPairs = [repmat(FilIdx,length(AdjFil),1), AdjFil];
                                    nAdjPairs = size(CurrentPairs,1);
                                    if length(find(~FilamentOptimized))==1
                                        nAdjPairs = 2;
                                    end
                                    % Find the indices in the AllPairs varible that match CurrentPairs
                                    [~,idxAP,~] = intersect( sort(AllPairs,2), sort(CurrentPairs,2), 'rows');
                                    if nLoops > 100
                                        nAdjPairs = 2;
                                    end
                                    nLoops = nLoops + 1;
                                end


                   %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  
        
                            for pair = idxAP'
                                        % UP BOND calculations ============================================================================================================
                                            F1 = AllPairs(pair,1);
                                            F2 = AllPairs(pair,2);
    
                                            if F1 == FilIdx      
                                                A1 = Alphas;
                                                A2 = repmat(FinalAlphas(F2),length(Alphas),1);
                                            else 
                                                A1 = repmat(FinalAlphas(F1),length(Alphas),1);
                                                A2 = Alphas;
                                            end
    
                                            Filaments = cat(3, FilamentsStart(:,F1)+A1', FilamentsStart(:,F2)+A2' ); % Offset filments and angle Alpha 
                                
                                            
                                            [D_up, PTS_up] = CalculateXYZForCentroidsAndResidues(X,Y,Z,Filaments,F1,F2,R);
%                                             BPa = MP.BPF.ASC_im1_j_bp(D_up.ASC_im1_j);      BPa(isnan(BPa)) = 1;
%                                             BPb = MP.BPF.ASC_ip3_jp4_bp(D_up.ASC_ip3_jp4);  BPb(isnan(BPb)) = 1;
                                            BP1 = MP.BPF.R95a_R95b_bp(D_up.R95a_R95b);         BP1(isnan(BP1)) = 0;  %BP1( BP1 < BondThreshold) = 0;
                                            BP2 = MP.BPF.R95a_S350_bp(D_up.R95a_S350);         BP2(isnan(BP2)) = 0;  %BP2( BP2 < BondThreshold) = 0;
                                            BP3 = MP.BPF.D24_S350_bp(D_up.D24_S350);           BP3(isnan(BP3)) = 0;  %BP3( BP3 < BondThreshold) = 0;
                                            BP4 = MP.BPF.D24_R95b_bp(D_up.D24_R95b);           BP4(isnan(BP4)) = 0;  %BP4( BP4 < BondThreshold) = 0;
                                            BP5 = MP.BPF.FascinCent_bp(D_up.Plane2Centroid);   BP5(isnan(BP5)) = 0;  %BP5( BP5 < BondThreshold) = 0;
                                            BondProbabilityUp = BP1 .* BP2 .* BP3 .* BP4 .* BP5;             % Calculate bond probabilities
                                            % Calculate binding probability -----------------------------------------------------------------    
                                            BondProbabilityUp(isnan(BondProbabilityUp)) = 0;
                                            BondProbabilityUp(BondProbabilityUp < BondThreshold) = 0;   % Threshold bond probabilities
                                            % Remove probabilities of monomers that are already used in another pair binding ---------------
                                            BondProbabilityUp( MonomersTaken(1:end-3,F1), : ) = 0;
                                            BondProbabilityUp( MonomersTaken(2:end-2,F2), : ) = 0;
                                        
                                            
                                        % DOWN BOND calculation ===========================================================================================================
                                        
                                            F1 = AllPairs(pair,2); % Switch which filament is F1 in order to compute DOWN bonds of the same pair
                                            F2 = AllPairs(pair,1);
    
                                            if F1 == FilIdx      
                                                A1 = Alphas;
                                                A2 = repmat(FinalAlphas(F2),length(Alphas),1);
                                            else 
                                                A1 = repmat(FinalAlphas(F1),length(Alphas),1);
                                                A2 = Alphas;
                                            end
    
                                            Filaments = cat(3,FilamentsStart(:,F1)+A1', FilamentsStart(:,F2)+A2'); % Offset filments and angle Alpha 
                                        
                                            % Calculate distances and probabilites ----------------------------------------------------------
                                            [D_dn,PTS_dn] = CalculateXYZForCentroidsAndResidues(X,Y,Z,Filaments,F1,F2,R);
%                                             BPa = MP.BPF.ASC_im1_j_bp(D_dn.ASC_im1_j);      BPa(isnan(BPa)) = 1;
%                                             BPb = MP.BPF.ASC_ip3_jp4_bp(D_dn.ASC_ip3_jp4);  BPb(isnan(BPb)) = 1;
                                            BP1 = MP.BPF.R95a_R95b_bp(D_dn.R95a_R95b);         BP1(isnan(BP1)) = 0;   %BP1( BP1 < BondThreshold) = 0;
                                            BP2 = MP.BPF.R95a_S350_bp(D_dn.R95a_S350);         BP2(isnan(BP2)) = 0;   %BP2( BP2 < BondThreshold) = 0;
                                            BP3 = MP.BPF.D24_S350_bp(D_dn.D24_S350);           BP3(isnan(BP3)) = 0;   %BP3( BP3 < BondThreshold) = 0;
                                            BP4 = MP.BPF.D24_R95b_bp(D_dn.D24_R95b);           BP4(isnan(BP4)) = 0;   %BP4( BP4 < BondThreshold) = 0;
                                            BP5 = MP.BPF.FascinCent_bp(D_dn.Plane2Centroid);   BP5(isnan(BP5)) = 0;   %BP5( BP5 < BondThreshold) = 0;
                                            % Calculate binding probability -----------------------------------------------------------------    
                                            BondProbabilityDn = BP1 .* BP2 .* BP3 .* BP4 .* BP5;             % Calculate bond probabilities
    
                                            % Calculate binding probability -----------------------------------------------------------------  
                                            BondProbabilityDn(isnan(BondProbabilityDn)) = 0;
                                            BondProbabilityDn(BondProbabilityDn < BondThreshold) = 0;   % Threshold bond probabilities
                        
                                            % Remove probabilities of monomers that are already used in another pair binding ----------------
                                            BondProbabilityDn( MonomersTaken(1:end-3,F1), : ) = 0;
                                            BondProbabilityDn( MonomersTaken(2:end-2,F2), : ) = 0;
                        
                            
                                        
                                        % Now make sure Up and Down bonds don't overlap---------------------------------------------------------------
                                        [~,ColUp] = find(BondProbabilityUp >= BondThreshold);   % Find all monomers with bonds above threshold
                                        [~,ColDn] = find(BondProbabilityDn >= BondThreshold);   % Find all monomers with bonds above threshold
                                        Col = unique([ColUp;ColDn]);
                                        % Run through all the up and down bond combinations and find ones that overalap. If the overlap, take the whichever one is stronger.
                                        for col = Col'
                                                    IdxUp2 = find(BondProbabilityUp(:,col) >= BondThreshold);
                                                    IdxDn2 = find(BondProbabilityDn(:,col) >= BondThreshold);      
                                                    for m1 = 1:numel(IdxUp2)
                                                        for m2 = 1:numel(IdxDn2) 
                                                            I = intersect( [IdxUp2(m1):IdxUp2(m1)+3], [IdxDn2(m2):IdxDn2(m2)+3] );
                                                            if ~isempty(I)
                                                                if BondProbabilityUp(IdxUp2(m1),col) <= BondProbabilityDn(IdxDn2(m2),col)
                                                                    BondProbabilityUp(IdxUp2(m1),col) = 0;
                                                                else
                                                                    BondProbabilityDn(IdxDn2(m2),col) = 0;
                                                                end
                                                            end
                                                        end
                                                    end
                                        end

                            
                                        % Now record the monomers that are in another pair combination -----------------------------------------------
                                        % Start with UP bonds
                                        F1 = AllPairs(pair,1);
                                        F2 = AllPairs(pair,2);
                                        [~,ColUp] = find(BondProbabilityUp >= BondThreshold);      % Find all alpha columns having non-zero values
                                        for col = ColUp'  
                                            IdxUp3 = find(BondProbabilityUp(:,col) >= BondThreshold);
                                            for ii = 1:numel(IdxUp3) % For each possible UpBond above th
                                                MonomersTaken(IdxUp3(ii)  ,F1) = true;
                                                MonomersTaken(IdxUp3(ii)+1,F2) = true;
                                            end
                                        end
                            
                                        % Now record DOWN bond monomers
                                        F1 = AllPairs(pair,2);
                                        F2 = AllPairs(pair,1);
                                        [~,ColDn] = find(BondProbabilityDn >= BondThreshold);      % Find all alpha columns having non-zero values
                                        for col = ColDn'  
                                            IdxDn3 = find(BondProbabilityDn(:,col) >= BondThreshold);
                                            for ii = 1:numel(IdxDn3) % For each possible UpBond above th
                                                MonomersTaken(IdxDn3(ii)  ,F1) = true;
                                                MonomersTaken(IdxDn3(ii)+1,F2) = true;
                                            end
                                        end
                                        % ------------------------------------------------------------------------------------------------------------
                                    
                                         
                                        BondStrengthUp(pair,1:end-3,:) = BondProbabilityUp;    
                                        BondStrengthDn(pair,1:end-3,:) = BondProbabilityDn;  
                                        
                                        Sum = sum( BondProbabilityUp + BondProbabilityDn, 1); 
                                        Sum(Sum<BondThreshold)  = 0;
                                        SumBP = Sum; 
                                        SumBP(isinf(SumBP)) = 10;
                                        SumBP(isnan(SumBP)) = 0;

                                        PairBondStrength(:,pair) = SumBP';
                            end
        
                            % Record optimal alpha for FilIdx based on the =========================
                            BondSums = sum(PairBondStrength(:,idxAP),2); % Combine the probability sum of all adjacent filaments
                            [~,aIdx] = max(BondSums);
                            AlphaIdx(idxAP,1) = aIdx;
                            FinalAlphas(FilIdx,1) = Alphas(aIdx);
                            FilamentOptimized(FilIdx,1) = true;
                    end
        
        
                    OptimizedFilaments = FilamentsStart + FinalAlphas';
                    OrderAndAlpha = [FilamentOrder, FinalAlphas(FilamentOrder)];
                    

                    if CreateFilamentPlot
                
                        figure(6); clf
                        set(gcf,'Color','w')
                        TL = tiledlayout(3,1,"TileSpacing","tight","Padding","loose");
                        nexttile(TL,1,[2,1])
                            
                        ShowAlphas = true;
                        ShowColorbar = false;
                        ShowFourPoints = false;
                        
                        %SumProbArray = zeros(7*nMonomers,7*nMonomers); % Uncomment for contact map

                        PlotOutput3
                        % Grab XYZ of Actin Subunit centroids and Fascin bonds
                        MO.FBCxyz = FBCxyz; 
                        MO.ASCxyz = ASCxyz;

                        if exist('SumProbArray')
                            MO.ContactMap = SumProbArray;
                        end

                        view(0,0)
                        %ylim([-2,2])

                        nexttile(TL,3)
                        ShowAlphas = true;
                        ShowColorbar = true;
                        ShowFourPoints = false;

                        PlotOutput3
                        view(0,90)
                        drawnow
                    end


                    MO.nBonds = nBonds;
                    MO.ProbabilitySum = BondSums(aIdx);

        %==================================================================

