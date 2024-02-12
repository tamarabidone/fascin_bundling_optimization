% Run   /Users/keithcarney/Dropbox/ENM/Actin Bundle Assembly/ActinBundleModel_07.m   first
% Model parameters --------------------------------------------------------------
TwistTheta      = MP.TwistTheta; %-166.7;  % degrees
MonLength       = MP.MonLength;     % nm
FilSep          = MP.FilSep;   % nm  STD = 0.06, Filament axis separation distance
FILAMENTradius  = 3.5;     % nm  Filament external face structure radius
ASCradius       = MP.R.ASC;    % nm  STD = 0.06  Actin subunit centroid radius
%S1Cradius       = 2.63;    % nm  STD = 0.06  Subdomain 1 centroid radius
% figure(6); clf
% set(gcf,'Color','w')
z1 = 0:0.01:((nMonomers-1)*MonLength);
%[X1GammaVec, Y1GammaVec] = pol2cart(Gamma21(G)*(pi/180),FilSep);
%Filaments = repmat( mod((0:TwistTheta:(nMonomers-1)*TwistTheta)',360), 1, 2);
ASCxyz = []; % Actin Subunit Centroid XYZ locations

for f = 1:N
    
    xp = X(f);
    yp = Y(f);

    [xm1,ym1,zm1] = pol2cart( (pi/180)*(z1.*(TwistTheta/MonLength) + OptimizedFilaments(1,f)), ASCradius, z1);
    xFil = xp + xm1;
    yFil = yp + ym1;
    plot3(xFil,yFil,zm1,'-','LineWidth',2,'Color',[0.5,0.5,0.5]); 
    hold on

    [xm2,ym2,zm2] = pol2cart( (pi/180)*(OptimizedFilaments(1:nMonomers,f)), ASCradius, MonLength*(0:(nMonomers-1)) );
    plot3( [xp.*ones(nMonomers,1),(xm2+xp).*ones(nMonomers,1)]',...
           [yp.*ones(nMonomers,1),(ym2+yp).*ones(nMonomers,1)]',...
           [zm2',zm2']', '-b','LineWidth',2)
   
    plot3( (xm2+xp).*ones(nMonomers,1), (ym2+yp).*ones(nMonomers,1), zm2', '.k', 'MarkerSize', 20 )
    plot3( [xp,xp], [yp,yp], [z1(1),z1(end)], ':k' ,'LineWidth',2)

    ASCxyz = [ASCxyz; [(xm2+xp).*ones(nMonomers,1), (ym2+yp).*ones(nMonomers,1), zm2']];

%     text((xm2+xp).*ones(nMonomers,1),...
%               (ym2+yp).*ones(nMonomers,1),...
%               zm2'+0.5,compose('%1.0f',OptimizedFilaments(1,f)),'Fontsize',12,'HorizontalAlignment','center','Color','b','BackgroundColor','none'); 
end
hold off
axis equal
%set(gca,'CameraViewAngle',8)
% Now add N+1 connections
nPairs = size(AllPairs,1);
nColors = 10;
%Colors = turbo(nColors);
[Colors,~] = plasma(nColors);
%[Colors,~] = colorGradient([1,0,0],[0,1,0.2],nColors);
Cedges = linspace(0,1,nColors+1);
cIdx = 0;
%BondThreshold = 0.8;
%OptimizedFilaments = FilamentsStart + Alphas(idx,:);
nBonds = 0;
TotalBondStrength = 0;
nBondsPerFilament = zeros(N,1);
nBondsPerMonomer  = zeros(N,nMonomers);
ProbabilitesPerMonomer = zeros(N,nMonomers);
RecordBondProbabilities = [];
FBCmono = cell(1,nPairs); % fascin connected monomers for each fascin
FBCxyz  = []; % 
FascinTightness = 1;

for pair = 1:nPairs
    % Create graphics for UP bonds --------------------------------------------------------------------------------------------------
    F1 = AllPairs(pair,1);
    F2 = AllPairs(pair,2);
    FilamentsNew = mod( [OptimizedFilaments(:,F1), OptimizedFilaments(:,F2)], 360); % Offset filments and angle Alpha 
   [F1_ASCxi, F1_ASCyi, F2_ASCxjp1, F2_ASCyjp1, F1_ASCxip2, F1_ASCyip2, F2_ASCxjp3, F2_ASCyjp3] = CalculateXYForMonomerPositionsUp(X,Y,FilamentsNew,F1,F2,ASCradius);

    for m = 1:nMonomers-3
            if BondStrengthUp(pair,m,AlphaIdx(pair)) >= BondThreshold
                P = BondStrengthUp(pair,m,AlphaIdx(pair));
                X1 = [F1_ASCxi(m), F1_ASCxip2(m), F2_ASCxjp3(m), F2_ASCxjp1(m), F1_ASCxi(m)];
                Y1 = [F1_ASCyi(m), F1_ASCyip2(m), F2_ASCyjp3(m), F2_ASCyjp1(m), F1_ASCyi(m)];
                Z1 = [zm2(m),  zm2(m+2), zm2(m+3), zm2(m+1), zm2(m)];

                FBCmono{1,pair} = [FBCmono{1,pair}; [m,m+1]];

                % Move corners in for XYZ file.....
                FascinCenter = [mean(X1(1:4)), mean(Y1(1:4)), mean(Z1(1:4))];

                X1new = mean([X1(1:4); repmat(FascinCenter(1),FascinTightness,4)]);
                Y1new = mean([Y1(1:4); repmat(FascinCenter(2),FascinTightness,4)]);
%                 X1 = [X1new,X1new(1)];
%                 Y1 = [Y1new,Y1new(1)];

                % Record the four fascin corners
                FBCxyz  = [FBCxyz;  [X1new', Y1new', Z1(1:4)', repmat(P,4,1) ] ];

                if exist('SumProbArray')
                    F1b = FilamentOrientationConversion(F1);
                    F2b = FilamentOrientationConversion(F2);
                    index_m_F1   = m + (F1b-1)*nMonomers;
                    index_mp1_F2 = m+1+(F2b-1)*nMonomers;
                    SumProbArray(index_m_F1, index_mp1_F2) = SumProbArray(index_m_F1,index_mp1_F2) + P;
                    SumProbArray(index_mp1_F2, index_m_F1) = SumProbArray(index_mp1_F2,index_m_F1) + P;
                end

                cIdx =  discretize(BondStrengthUp(pair,m,AlphaIdx(pair)),Cedges) ;
                hold on
                    patch(X1,Y1,Z1,Colors(cIdx,:),'FaceAlpha',0.8)
                if ShowFourPoints

                    % Overlap R95, D24, S350, R95b  points
                    PTS = CalculateXYZforResidue(X,Y,Z,FilamentsNew,F1,F2,R,'Up');
                    Xres = [PTS.R95a_x(m); PTS.D24_x(m+2); PTS.S350_x(m+3); PTS.R95b_x(m+1); PTS.R95a_x(m)];
                    Yres = [PTS.R95a_y(m); PTS.D24_y(m+2); PTS.S350_y(m+3); PTS.R95b_y(m+1); PTS.R95a_y(m)];
                    Zres = [PTS.R95a_z(m); PTS.D24_z(m+2); PTS.S350_z(m+3); PTS.R95b_z(m+1); PTS.R95a_z(m)];
                    patch(Xres,Yres,Zres,0.5*ones(1,3))
                
                    plot3(PTS.R95a_x(m),   PTS.R95a_y(m),   PTS.R95a_z(m),  '.r','MarkerSize',30)
                    plot3(PTS.D24_x(m+2),  PTS.D24_y(m+2),  PTS.D24_z(m+2), '.b','MarkerSize',30)
                    plot3(PTS.S350_x(m+3), PTS.S350_y(m+3), PTS.S350_z(m+3),'.g','MarkerSize',30)
                    plot3(PTS.R95b_x(m+1), PTS.R95b_y(m+1), PTS.R95b_z(m+1),'.k','MarkerSize',30)
                end

                hold off

                 nBonds = nBonds + 1;
%                 nBondsPerFilament(F1,1) = nBondsPerFilament(F1,1) + 1;
%                 nBondsPerFilament(F2,1) = nBondsPerFilament(F2,1) + 1;
%                 nBondsPerMonomer(F1,m:m+2)  = nBondsPerMonomer(F1,m:m+2) + 1;
%                 nBondsPerMonomer(F2,m+1:m+3)  = nBondsPerMonomer(F2,m+1:m+3) + 1;
%                 ProbabilitesPerMonomer(F1,m:m+2) = ProbabilitesPerMonomer(F1,m:m+2) + BondStrengthUp(pair,m,AlphaIdx(pair));
%                 ProbabilitesPerMonomer(F2,m+1:m+3) = ProbabilitesPerMonomer(F2,m+1:m+3) + BondStrengthUp(pair,m,AlphaIdx(pair));
%                 RecordBondProbabilities = [RecordBondProbabilities; BondStrengthUp(pair,m,AlphaIdx(pair))];
            end
    end

    % Create graphics for DOWN bonds --------------------------------------------------------------------------------------------------
    F1 = AllPairs(pair,2);
    F2 = AllPairs(pair,1);
   [F1_ASCxi, F1_ASCyi, F2_ASCxjp1, F2_ASCyjp1, F1_ASCxip2, F1_ASCyip2, F2_ASCxjp3, F2_ASCyjp3] = CalculateXYForMonomerPositionsDn(X,Y,FilamentsNew,F1,F2,ASCradius);

    for m = 1:nMonomers-3            
            if BondStrengthDn(pair,m,AlphaIdx(pair)) >= BondThreshold
                P = BondStrengthDn(pair,m,AlphaIdx(pair));

                X1 = [F1_ASCxi(m), F1_ASCxip2(m), F2_ASCxjp3(m), F2_ASCxjp1(m), F1_ASCxi(m)];
                Y1 = [F1_ASCyi(m), F1_ASCyip2(m), F2_ASCyjp3(m), F2_ASCyjp1(m), F1_ASCyi(m)];
                Z1 = [zm2(m),  zm2(m+2), zm2(m+3), zm2(m+1), zm2(m)];

                FBCmono{1,pair} = [FBCmono{1,pair}; [m+1,m]];

                % Move corners in for XYZ file.....
                FascinCenter = [mean(X1(1:4)), mean(Y1(1:4)), mean(Z1(1:4))];

                X1new = mean([X1(1:4); repmat(FascinCenter(1),FascinTightness,4)]);
                Y1new = mean([Y1(1:4); repmat(FascinCenter(2),FascinTightness,4)]);
%                 X1 = [X1new,X1new(1)];
%                 Y1 = [Y1new,Y1new(1)];

                % Record the four fascin corners
                FBCxyz  = [FBCxyz;  [X1new', Y1new', Z1(1:4)', repmat(P,4,1) ] ];

                if exist('SumProbArray')
                    F1b = FilamentOrientationConversion(F1);
                    F2b = FilamentOrientationConversion(F2);
                    index_m_F1   = m + (F1b-1)*nMonomers;
                    index_mp1_F2 = m+1+(F2b-1)*nMonomers;
                    SumProbArray(index_m_F1, index_mp1_F2) = SumProbArray(index_m_F1,index_mp1_F2) + P;
                    SumProbArray(index_mp1_F2, index_m_F1) = SumProbArray(index_mp1_F2,index_m_F1) + P;
                end

                cIdx = discretize(BondStrengthDn(pair,m,AlphaIdx(pair)),Cedges) ;
                hold on
                patch(X1,Y1,Z1,Colors(cIdx,:),'FaceAlpha',0.8)

                % Overlap R95, D24, S350, R95b  points
                if ShowFourPoints
                    PTS = CalculateXYZforResidue(X,Y,Z,FilamentsNew,F1,F2,R,'Dn');
                    Xres = [PTS.R95a_x(m); PTS.D24_x(m+2); PTS.S350_x(m+3); PTS.R95b_x(m+1); PTS.R95a_x(m)];
                    Yres = [PTS.R95a_y(m); PTS.D24_y(m+2); PTS.S350_y(m+3); PTS.R95b_y(m+1); PTS.R95a_y(m)];
                    Zres = [PTS.R95a_z(m); PTS.D24_z(m+2); PTS.S350_z(m+3); PTS.R95b_z(m+1); PTS.R95a_z(m)];
                    patch(Xres,Yres,Zres,0.5*ones(1,3))

                    plot3(PTS.R95a_x(m),   PTS.R95a_y(m),   PTS.R95a_z(m),  '.r','MarkerSize',30)
                    plot3(PTS.D24_x(m+2),  PTS.D24_y(m+2),  PTS.D24_z(m+2), '.b','MarkerSize',30)
                    plot3(PTS.S350_x(m+3), PTS.S350_y(m+3), PTS.S350_z(m+3),'.g','MarkerSize',30)
                    plot3(PTS.R95b_x(m+1), PTS.R95b_y(m+1), PTS.R95b_z(m+1),'.k','MarkerSize',30)

                end
                hold off

                 nBonds = nBonds + 1;
%                 nBondsPerFilament(F1,1) = nBondsPerFilament(F1,1) + 1;
%                 nBondsPerFilament(F2,1) = nBondsPerFilament(F2,1) + 1;
%                 nBondsPerMonomer(F1,m:m+2)  = nBondsPerMonomer(F1,m:m+2) + 1;
%                 nBondsPerMonomer(F2,m+1:m+3)  = nBondsPerMonomer(F2,m+1:m+3) + 1;
%                 ProbabilitesPerMonomer(F1,m:m+2) = ProbabilitesPerMonomer(F1,m:m+2)     + BondStrengthDn(pair,m,AlphaIdx(pair));
%                 ProbabilitesPerMonomer(F2,m+1:m+3) = ProbabilitesPerMonomer(F2,m+1:m+3) + BondStrengthDn(pair,m,AlphaIdx(pair));
%                 RecordBondProbabilities = [RecordBondProbabilities; BondStrengthDn(pair,m,AlphaIdx(pair))];
            end
    end
end
%set(gca,'CameraViewAngle',10)
%view(0,0)
if ShowAlphas
    T = text(X,Y,z1(end)*ones(N,1)+2,compose('%1.1f',FinalAlphas'),'Fontsize',18,'HorizontalAlignment','center','Color','w','BackgroundColor','k'); 
    
    for idx = 1:size(AllPairs,1)
        ADx = mean(X(AllPairs(idx,:)));
        ADy = mean(Y(AllPairs(idx,:)));
        ADidx1 = AllPairs(idx,1);
        ADidx2 = AllPairs(idx,2);
        if X(ADidx1) > X(ADidx2)
            AlphaDiff =  FinalAlphas(AllPairs(idx,1)) - FinalAlphas(AllPairs(idx,2)); % abs(diff(FinalAlphas(AllPairs(idx,:))))
        else
            AlphaDiff =  FinalAlphas(AllPairs(idx,2)) - FinalAlphas(AllPairs(idx,1));
        end
        T2 = text(ADx,ADy,z1(end)+2,compose('%1.1f',AlphaDiff),'Fontsize',18,'HorizontalAlignment','center','Color','k','BackgroundColor',[1,1,1],'EdgeColor','k'); 
    end
    SF = MP.FirstFilamentIndex;
    delete(T(SF))
    T3 = text(X(SF),Y(SF),z1(end)+2,compose('%1.0f',FinalAlphas(SF)),'Fontsize',18,'HorizontalAlignment','center','Color','w','BackgroundColor',[0,0,1],'FontWeight','bold');
    
end
hold off

% xlabel('X (nm)','Fontsize',25)
% ylabel('Y (nm)','Fontsize',25)
% zlabel('Z (nm)','Fontsize',25)
axis off

if ShowColorbar
    CB = colorbar;
    colormap(CB,Colors)
    CB.FontSize = 20;
    CB.Limits = [0,1];
    CB.Ticks = 0:0.1:1;
    CB.TickLabels = compose('%0.1f',CB.Ticks);
    CB.TickDirection = 'both';
    CB.Label.String = ['Fascin bond probability'];
    CB.Label.FontSize = 30;
    CB.LineWidth = 2;
end
