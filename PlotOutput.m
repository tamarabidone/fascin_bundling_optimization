% Run   /Users/keithcarney/Dropbox/ENM/Actin Bundle Assembly/ActinBundleModel_07.m   first
% Model parameters --------------------------------------------------------------
TwistTheta      = MP.TwistTheta; %-166.7;  % degrees
MonLength       = MP.MonLength;     % nm
FilSep          = MP.FilSep;   % nm  STD = 0.06, Filament axis separation distance
FILAMENTradius  = 3.5;     % nm  Filament external face structure radius
ASCradius       = MP.R.ASC;    % nm  STD = 0.06  Actin subunit centroid radius
%S1Cradius       = 2.63;    % nm  STD = 0.06  Subdomain 1 centroid radius
figure(4); clf
set(gcf,'Color','w')
% Plot Filaments
%nMonomers = 100;
z1 = 0:0.01:((nMonomers-1)*MonLength);
%[X1GammaVec, Y1GammaVec] = pol2cart(Gamma21(G)*(pi/180),FilSep);
%Filaments = repmat( mod((0:TwistTheta:(nMonomers-1)*TwistTheta)',360), 1, 2);
ASCxyz = []; % Actin Subunit Centroid XYZ locations



for n = 1:2 %N
    F = C(p,n);
    xp = x(F);
    yp = y(F);

    [xm1,ym1,zm1] = pol2cart( (pi/180)*(z1.*(TwistTheta/MonLength) + Filaments(1,n)), ASCradius, z1);
    xFil = xp + xm1;
    yFil = yp + ym1;
    plot3(xFil,yFil,zm1,'-b','LineWidth',2); 
    hold on

    [xm2,ym2,zm2] = pol2cart( (pi/180)*(Filaments(1:nMonomers,n)), ASCradius, MonLength*(0:(nMonomers-1)) );
    plot3( [xp.*ones(nMonomers,1),(xm2+xp).*ones(nMonomers,1)]',...
           [yp.*ones(nMonomers,1),(ym2+yp).*ones(nMonomers,1)]',...
           [zm2',zm2']', '-k','LineWidth',2)

    plot3( (xm2+xp).*ones(nMonomers,1), (ym2+yp).*ones(nMonomers,1), zm2,'.k', 'MarkerSize', 20)
    plot3(  [xp,xp], [yp,yp], [z1(1),z1(end)],'-k' )

    ASCxyz = [ASCxyz; [(xm2+xp).*ones(nMonomers,1), (ym2+yp).*ones(nMonomers,1), zm2']]; % Record monomer XYZ
end
hold off
axis equal


% Now add N+1 connections
nPairs = size(C,1);
nColors = 10;
[Colors,~] = plasma(nColors);
Cedges = linspace(0,1,nColors+1);
cIdx = 0;
RecordBondProbabilities = [];
FBCxyz  = []; % 

    for m = 1:nMonomers-1
            if BondStrengthUp(A1(1),A2(1),m) >= BondThreshold
                P = BondStrengthUp(A1(1),A2(1),m);
                cIdx = cIdx + 1;
                if cIdx > 7; cIdx = 1; end
                F1 = C(p,1);
                F2 = C(p,2);
                [F1_ASCxi, F1_ASCyi, F2_ASCxjp1, F2_ASCyjp1, F1_ASCxip2, F1_ASCyip2, F2_ASCxjp3, F2_ASCyjp3] = CalculateXYForMonomerPositionsUp(x,y,Filaments,F1,F2,ASCradius);
                X1 = [F1_ASCxi(m), F1_ASCxip2(m), F2_ASCxjp3(m), F2_ASCxjp1(m), F1_ASCxi(m)];
                Y1 = [F1_ASCyi(m), F1_ASCyip2(m), F2_ASCyjp3(m), F2_ASCyjp1(m), F1_ASCyi(m)];
                Z1 = [zm2(m),  zm2(m+2), zm2(m+3), zm2(m+1), zm2(m)];
                
                cIdx =  discretize(BondStrengthUp(A1(1),A2(1),m) ,Cedges) ;
                RecordBondProbabilities = [RecordBondProbabilities; BondStrengthUp(A1(1),A2(1),m)];
                % Record the four fascin corners
                FBCxyz  = [FBCxyz;  [X1(1:4)', Y1(1:4)', Z1(1:4)', repmat(P,4,1) ] ];

                hold on
                patch(X1,Y1,Z1,Colors(cIdx,:))
                hold off
            end
 
            if BondStrengthDn(A1(1),A2(1),m) >= BondThreshold
                P = BondStrengthDn(A1(1),A2(1),m);
                cIdx = cIdx + 1;
                if cIdx > 7; cIdx = 1; end
                F1 = C(p,2);
                F2 = C(p,1);
                [F1_ASCxi, F1_ASCyi, F2_ASCxjp1, F2_ASCyjp1, F1_ASCxip2, F1_ASCyip2, F2_ASCxjp3, F2_ASCyjp3] = CalculateXYForMonomerPositionsDn(x,y,Filaments,F1,F2,ASCradius);
                X1 = [F1_ASCxi(m), F1_ASCxip2(m), F2_ASCxjp3(m), F2_ASCxjp1(m), F1_ASCxi(m)];
                Y1 = [F1_ASCyi(m), F1_ASCyip2(m), F2_ASCyjp3(m), F2_ASCyjp1(m), F1_ASCyi(m)];
                Z1 = [zm2(m),  zm2(m+2), zm2(m+3), zm2(m+1), zm2(m)];

                cIdx =  discretize(BondStrengthDn(A1(1),A2(1),m) ,Cedges) ;
                RecordBondProbabilities = [RecordBondProbabilities; BondStrengthDn(A1(1),A2(1),m)];
                % Record the four fascin corners
                FBCxyz  = [FBCxyz;  [X1(1:4)', Y1(1:4)', Z1(1:4)', repmat(P,4,1) ] ];

                hold on
                patch(X1,Y1,Z1,Colors(cIdx,:))
                hold off
            end

    end

%disp(RecordBondProbabilities)
set(gca,'LineWidth',2,'FontSize',18)    
xlabel('X (nm)','Fontsize',25)
ylabel('Y (nm)','Fontsize',25)
zlabel('Z (nm)','Fontsize',25)
hold off
CB = colorbar;
colormap(CB,Colors)
%CB.Position(1) = 0.80;
CB.FontSize = 20;
CB.Limits = [0,1];
%CB.TickLabels = compose('%0.1f',CB.Ticks);
CB.TickDirection = 'both';
%CB.Label.String = ['Fascin bond probability  (bond threshold = ',sprintf('%0.2f',BondThreshold),')'];
CB.Label.String = ['Fascin bond probability'];
CB.Label.FontSize = 30;
CB.LineWidth = 2;


