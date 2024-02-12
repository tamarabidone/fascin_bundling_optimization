function [F1_xi,F1_yi, F2_xjp1,F2_yjp1, F1_xip2,F1_yip2, F2_xjp3,F2_yjp3] = CalculateXYForMonomerPositionsUp(x,y,Filaments,F1,F2,radius)

                % Calculate position of F1 i
                [xi, yi]     = pol2cart( (pi/180)*Filaments(1:end-3,1), radius);
                 F1_xi   = x(F1) + xi;
                 F1_yi   = y(F1) + yi;

                % Calculate position of F2 j+1
                [xjp1, yjp1] = pol2cart( (pi/180)*Filaments(2:end-2,2), radius);
                 F2_xjp1 = x(F2) + xjp1;
                 F2_yjp1 = y(F2) + yjp1;

                % Calculate position of F1 i+2
                [xip2, yip2] = pol2cart( (pi/180)*Filaments(3:end-1,1), radius);
                 F1_xip2 = x(F1) + xip2;
                 F1_yip2 = y(F1) + yip2;
                 
                % Calculate position of F2 j+3
                [xjp3, yjp3] = pol2cart( (pi/180)*Filaments(4:end ,2),  radius);
                 F2_xjp3 = x(F2) + xjp3;
                 F2_yjp3 = y(F2) + yjp3;