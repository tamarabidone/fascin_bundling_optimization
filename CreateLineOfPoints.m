function [X,Y] = CreateLineOfPoints(nPts,FilSep)

    if nPts > 1
        X = ( 0:FilSep:(nPts-1)*FilSep )';
        Y = zeros(nPts,1);
    else
        X = 0;
        Y = 0;
    end


end