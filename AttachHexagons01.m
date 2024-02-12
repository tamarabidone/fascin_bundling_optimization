function [nPts,X,Y] = AttachHexagons01(nHex,FilSep)

    
    nPts = 7; 
    [x,y] = CreateHexagonArray(nPts,FilSep);

    X = x;
    Y = y;

    x = x([1,2,3,4,6,7]);
    y = y([1,2,3,4,6,7]);

    Dx = 2*FilSep;
    Dy = 0;

    for n = 1:nHex-1
        X = [X; (x+n*Dx)];
        Y = [Y; (y+n*Dy)];
    end
        
    idx1 = find(Y > 5);
    SP1 = sortrows([X(idx1),Y(idx1)],1,'ascend');
    idx2 = find(Y <5 & Y > -5);
    SP2 = sortrows([X(idx2),Y(idx2)],1,'ascend');
    idx3 = find(Y < -5);
    SP3 = sortrows([X(idx3),Y(idx3)],1,'ascend');
    X = [SP1(:,1); SP2(:,1); SP3(:,1)];
    Y = [SP1(:,2); SP2(:,2); SP3(:,2)];
    X = X + abs(min(X));

    nPts = length(X);
%     figure(1)
%     plot(X,Y,'.','MarkerSize',20)
%     axis equal
end
