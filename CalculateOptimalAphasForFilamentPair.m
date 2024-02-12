function MO = CalculateOptimalAphasForFilamentPair(MP, PairPlot, HeatMapPlot)


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

% Grab Pair points ------------------------------------------------------------------------------------
x = MP.FilCoords.X(MP.FilPairCom(MP.FirstPairIndex,:));
y = MP.FilCoords.Y(MP.FilPairCom(MP.FirstPairIndex,:)); 
z = MonLength*(0:nMonomers-1)';
N = length(x); 
%-------------------------------------------------------------------------------------------------------------  
FilamentsStart = repmat( mod((0:TwistTheta:(nMonomers-1)*TwistTheta)',360), 1, N);

% For the first two filaments, very gamma and alpha of F2 two find its optimal binding location
Alpha1   = MP.Alpha1;
Alpha2   = MP.Alpha2;
SumBP    = zeros(length(Alpha1),length(Alpha2));
% % Create all point pair combinations
%     c = nchoosek(1:N,2); 
% % Now filter out the combinations that do are not a distance FilSep apart
%     Dist = sqrt( ( x(c(:,2))-x(c(:,1)) ).^2 + ( y(c(:,2))-y(c(:,1)) ).^2 );
%     ind = find( Dist < (FilSep + 0.001) );
%     C = c(ind,:);
%     nPairs = size(C,1);

p = 1; %MP.FirstPairIndex; % select which filament pair to optimize. Max value is size(C,1)
C = MP.FilPairCom;
BondStrengthUp = zeros( length(Alpha1), length(Alpha2), nMonomers);
BondStrengthDn = zeros( length(Alpha1), length(Alpha2), nMonomers);

BondThreshold = MP.BondThreshold;


% MODEL START ==================================================================================================================================================

            [~, Alpha1grid, Alpha2grid] = ndgrid( 1:nMonomers, Alpha1, Alpha2);

            % UP BOND calculation ---------------------------
                F1 = C(p,1);
                F2 = C(p,2);

                Filament1 =  mod( repmat(FilamentsStart(:,F1),[1,size(Alpha1grid,[2,3])]) + Alpha1grid,  360);
                Filament2 =  mod( repmat(FilamentsStart(:,F2),[1,size(Alpha2grid,[2,3])]) + Alpha2grid,  360);

                % Calculate distances and binding probabilities
                [D_up,PTS_up] = CalculateXYZForCentroidsAndResidues( MP.FilCoords.X, MP.FilCoords.Y, MP.FilCoords.Z ,cat(4,Filament1,Filament2), F1, F2, R);
%                 BPa = MP.BPF.ASC_im1_j_bp(  D_up.ASC_im1_j);      BPa(isnan(BPa)) = 1;
%                 BPb = MP.BPF.ASC_ip3_jp4_bp(D_up.ASC_ip3_jp4);    BPb(isnan(BPb)) = 1;
                BP1 = MP.BPF.R95a_R95b_bp(  D_up.R95a_R95b);       BP1(isnan(BP1)) = 0;   %BP1( BP1 < BondThreshold) = 0;
                BP2 = MP.BPF.R95a_S350_bp(  D_up.R95a_S350);       BP2(isnan(BP2)) = 0;   %BP2( BP2 < BondThreshold) = 0;
                BP3 = MP.BPF.D24_S350_bp(   D_up.D24_S350);        BP3(isnan(BP3)) = 0;   %BP3( BP3 < BondThreshold) = 0;
                BP4 = MP.BPF.D24_R95b_bp(   D_up.D24_R95b);        BP4(isnan(BP4)) = 0;   %BP4( BP4 < BondThreshold) = 0;
                BP5 = MP.BPF.FascinCent_bp(D_up.Plane2Centroid);   BP5(isnan(BP5)) = 0;  %BP5( BP5 < BondThreshold) = 0;
                % Calculate binding probability ----------------------------------------------------------------- 
                ProdWeightsUp = BP1.*BP2.*BP3.*BP4.*BP5;
                %DiffWeightsUp = repmat(  std(ProdWeightsUp,[],1),[],1),  [size(BP1,1),1,1] );
                BondProbabilityUp = ProdWeightsUp; %./DiffWeightsUp; % Calculate bond probabilities
                BondProbabilityUp(isnan(BondProbabilityUp)) = 0;
                BondProbabilityUp(BondProbabilityUp < BondThreshold) = 0;   % Threshold bond probabilities


            % DOWN BOND calculation ------------------------- 
                F1 = C(p,2); % just switch indices on C 
                F2 = C(p,1);

                Filament1 =  mod( repmat(FilamentsStart(:,F1),[1,size(Alpha2grid,[2,3])]) + Alpha2grid,  360);
                Filament2 =  mod( repmat(FilamentsStart(:,F2),[1,size(Alpha1grid,[2,3])]) + Alpha1grid,  360);

                % Calculate distances and binding probabilities
                [D_dn,PTS_dn] = CalculateXYZForCentroidsAndResidues( MP.FilCoords.X, MP.FilCoords.Y, MP.FilCoords.Z ,cat(4,Filament1,Filament2), F1, F2, R);
%                 BPa = MP.BPF.ASC_im1_j_bp(D_dn.ASC_im1_j);      BPa(isnan(BPa)) = 1;
%                 BPb = MP.BPF.ASC_ip3_jp4_bp(D_dn.ASC_ip3_jp4);  BPb(isnan(BPb)) = 1;
                BP1 = MP.BPF.R95a_R95b_bp(D_dn.R95a_R95b);         BP1(isnan(BP1)) = 0;   %BP1( BP1 < BondThreshold) = 0;
                BP2 = MP.BPF.R95a_S350_bp(D_dn.R95a_S350);         BP2(isnan(BP2)) = 0;   %BP2( BP2 < BondThreshold) = 0;
                BP3 = MP.BPF.D24_S350_bp(D_dn.D24_S350);           BP3(isnan(BP3)) = 0;   %BP3( BP3 < BondThreshold) = 0;
                BP4 = MP.BPF.D24_R95b_bp(D_dn.D24_R95b);           BP4(isnan(BP4)) = 0;   %BP4( BP4 < BondThreshold) = 0;
                BP5 = MP.BPF.FascinCent_bp(D_dn.Plane2Centroid);   BP5(isnan(BP5)) = 0;   %BP5( BP5 < BondThreshold) = 0;
                % Calculate binding probability -----------------------------------------------------------------  
                ProdWeightsDn = BP1.*BP2.*BP3.*BP4.*BP5;
                %DiffWeightsDn = repmat( max(ProdWeightsDn,[],1) - min(ProdWeightsDn,[],1), [size(BP1,1),1,1]);
                %DiffWeightsDn = repmat( std(ProdWeightsDn,[],1), [size(BP1,1),1,1]);
                BondProbabilityDn = ProdWeightsDn; %./DiffWeightsDn; % Calculate bond probabilities
                BondProbabilityDn(isnan(BondProbabilityDn)) = 0;
                BondProbabilityDn(BondProbabilityDn < BondThreshold) = 0;   % Threshold bond probabilities


                % Make sure Up and Down bonds don't overlap ------------------------------------------------------------------------------------------
                for d2 = 1:size(BondProbabilityUp,2)
                    for d3 = 1:size(BondProbabilityUp,3)
                                IdxUp = find(BondProbabilityUp(:,d2,d3) >= BondThreshold);   % Find all monomers with bonds above threshold
                                IdxDn = find(BondProbabilityDn(:,d2,d3) >= BondThreshold);   % Find all monomers with bonds above threshold
                                % Run through all the up and down bond combinations and find ones that overalap. If the overlap, take the whichever one is stronger.
                                if ~isempty(IdxUp) && ~isempty(IdxDn)
                                            for m1 = 1:numel(IdxUp)
                                                for m2 = 1:numel(IdxDn) 
                                                    I = intersect( IdxUp(m1):IdxUp(m1)+3, IdxDn(m2):IdxDn(m2)+3 );
                                                    if ~isempty(I)
                                                        if BondProbabilityUp(IdxUp(m1),d2,d3) <= BondProbabilityDn(IdxDn(m2),d2,d3)
                                                            BondProbabilityUp(IdxUp(m1),d2,d3) = 0;
                                                        else
                                                            BondProbabilityDn(IdxDn(m2),d2,d3) = 0;
                                                        end
                                                    end
                                                end
                                            end
                                end
                    end
                end

%                 DiffWeightsUp = repmat(  max(BondProbabilityUp,[],1) - min(BondProbabilityUp,[],1),  [size(BP1,1),1,1] );
%                 DiffWeightsDn = repmat(  max(BondProbabilityDn,[],1) - min(BondProbabilityDn,[],1),  [size(BP1,1),1,1] );
%                 SumWeightsUp = repmat( sum(BondProbabilityUp,1), [size(BP1,1),1,1]);
%                 SumWeightsDn = repmat( sum(BondProbabilityDn,1), [size(BP1,1),1,1]);

                DiffWeightsUp = squeeze( std(BondProbabilityUp,[],1) );
                DiffWeightsDn = squeeze( std(BondProbabilityDn,[],1) );
                SumWeightsUp = squeeze( sum(BondProbabilityUp,1) );
                SumWeightsDn = squeeze( sum(BondProbabilityDn,1) );

% MODEL END ==================================================================================================================================================
               
    % Record final Up and Down bonds and probability sum -----------------------------------------
    for m = 1:nMonomers-3 
        BondStrengthUp(:,:,m) = BondProbabilityUp(m,:,:);    
        BondStrengthDn(:,:,m) = BondProbabilityDn(m,:,:);  
    end

    Sum = squeeze( sum( BondProbabilityUp + BondProbabilityDn, 1) ); 
    Sum(Sum<BondThreshold)  = 0;
    stdBP = squeeze( std( BondProbabilityUp + BondProbabilityDn,[],1) );
    SumBP = Sum; %./stdBP;
    SumBP(isinf(SumBP)) = 10;
    SumBP(isnan(SumBP)) = 0;
    
%     MaxVal = max([BondStrengthUp(:);BondStrengthDn(:)]);
%     BondStrengthUp = BondStrengthUp/MaxVal;
%     BondStrengthDn = BondStrengthDn/MaxVal;
%     SumBP = sum( BondStrengthUp + BondStrengthDn, 3 ); 
%     SumBP = SumBP./max(SumBP(:));
      
% WeightsUp = BondProbabilityUp./DiffWeightsUp;      
%     WeightsUp(isnan(WeightsUp) ) = 0;
% WeightsDn = BondProbabilityDn./DiffWeightsDn;      
%     WeightsDn(isnan(WeightsDn) ) = 0;
% SumBP = WeightsUp + WeightsDn;


% Record final Up and Down bonds and probability sum -----------------------------------------
MO.BondStrengthDn = BondStrengthDn;  
MO.BondStrengthUp = BondStrengthUp;    
MO.SumBP = SumBP; 



% Find peaks in ALpah heat map
BW = imregionalmax(SumBP);
%BW( SumBP < 0.99*max(SumBP(:)) ) = false;
%[A1,A2] = find(BW);
Aindex = find(BW);
PV = SumBP(find(BW));
[A1,A2] = ind2sub(size(BW),Aindex);
B = sortrows([PV,A1,A2],1,'descend'); % Sort peaks and Alpha indices by peak height

% grab at most the first ten highest peaks
if size(B,1) > 20
    A1 = B(1:20,2);
    A2 = B(1:20,3);
else
    A1 = B(:,2);
    A2 = B(:,3);
end



 figure(5); clf; 
TL = tiledlayout(1,3,'TileSpacing','compact','Padding','compact'); %subplot(1,2,1)
    set(gcf,'Color','w')
    nexttile(TL,1,[1,1])
    imagesc(Alpha2,Alpha1,SumBP); 
    set(gca,'FontSize',18)
    xlabel('\it\alpha_{2}','FontSize',32); ylabel('\alpha_{1}','FontSize',32); 
%     xticks(-60:30:60)
%     yticks(-60:30:60)
    xticks(0:30:360)
    yticks(0:30:360)
    colormap(turbo); 
    cm = colorbar; 
    cm.Label.String = {'Sum monomer binding probabilities'};
    cm.Label.FontSize = 20;
    
    grid on
    set(gca,'GridAlpha',0.75,'GridColor',[0.6,0.6,0.6],'YDir','normal')
    title(['nMonomers = ', num2str(nMonomers)],'FontWeight','normal')
    axis equal tight
    hold on
    plot(Alpha2(A2),Alpha1(A1),'.w','MarkerSize',20)
    hold off

% Randomly grab one of the peaks as the alpha
Aind = randi(length(A1),1);
A1 = A1(Aind);
A2 = A2(Aind);
PeakValues = [Alpha1(A1)',Alpha2(A2)'];
MO.Peakvalues = PeakValues;
MO.A1 = A1;
MO.A2 = A2;


    nexttile(TL,2,[1,2])
    stem((1:nMonomers)', squeeze(BondStrengthUp(A1(1),A2(1),:)),'filled','LineWidth',3,'MarkerEdgeColor','k','MarkerFaceColor','r','Marker','o','MarkerSize',8)
    hold on
    stem((1:nMonomers)', squeeze(BondStrengthDn(A1(1),A2(1),:)),'filled','LineWidth',3,'MarkerEdgeColor','k','MarkerFaceColor','g','Marker','o','MarkerSize',8)
    hold off
    set(gca,'FontSize',24,'LineWidth',2)%,'PlotBoxAspectRatio',[2,1,1])  
    xlabel('Monomer i_{1}','FontSize',30); ylabel('Binding probability','FontSize',30); 
    grid on
    xlim([0,nMonomers-2])
    title(['For \alpha_{1} = ',sprintf('%4.1f',Alpha1(A1(1))),'^{o}  and  \alpha_{2} = ',...
           sprintf('%4.1f',Alpha2(A2(1))),'^{o}       \delta = ',sprintf('%4.1f',mode(abs( Alpha1(A1(1))- Alpha2(A2(1)) ),360))],'FontWeight','normal' )

    Filaments = mod( [FilamentsStart(:,F1)+Alpha1(A1(1)), FilamentsStart(:,F2)+Alpha2(A2(1))], 360);
    
    if PairPlot
        PlotOutput
    end



end