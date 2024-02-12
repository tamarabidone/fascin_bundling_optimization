function [D,PTS] = CalculateXYZForCentroidsAndResidues(x,y,z,Filaments,F1,F2,R)
    
    nDIMS = ndims(Filaments);
        
                
                % ACTIN SUBUNIT CENTROID (ASC), calculate positions of F1_i 
                    switch nDIMS
                        case 4; [xi, yi] = pol2cart( (pi/180)*Filaments(1:end-3,:,:,1), R.ASC);
                        case 3; [xi, yi] = pol2cart( (pi/180)*Filaments(1:end-3,:,1), R.ASC);
                        case 2; [xi, yi] = pol2cart( (pi/180)*Filaments(1:end-3,1), R.ASC);
                    end
                    
                    PTS.ASC1_xi   = x(F1) + xi;
                    PTS.ASC1_yi   = y(F1) + yi;
                    PTS.ASC1_zi   = z(1:end-3,1);
                    
                    XYZsize = size(xi);
                % ACTIN SUBUNIT CENTROID, calculate positions of F2_j+1
                    switch nDIMS
                        case 4; [xjp1, yjp1] = pol2cart( (pi/180)*Filaments(2:end-2,:,:,2), R.ASC);
                        case 3; [xjp1, yjp1] = pol2cart( (pi/180)*Filaments(2:end-2,:,2), R.ASC);
                        case 2; [xjp1, yjp1] = pol2cart( (pi/180)*Filaments(2:end-2,2), R.ASC);    
                    end

                    PTS.ASC2_xjp1 = x(F2) + xjp1;
                    PTS.ASC2_yjp1 = y(F2) + yjp1;
                    PTS.ASC2_zjp1 = z(2:end-2,1);   

    % Calculate distance between F1_i and F2_j+1
    D.ASC_i_jp1 = sqrt( (PTS.ASC2_xjp1 - PTS.ASC1_xi).^2 +...
                        (PTS.ASC2_yjp1 - PTS.ASC1_yi).^2 +...
                        (PTS.ASC2_zjp1 - PTS.ASC1_zi).^2 );  

                % ACTIN SUBUNIT CENTROID, calculate positions of F1_i+2
                    switch nDIMS
                        case 4; [xip2, yip2] = pol2cart( (pi/180)*Filaments(3:end-1,:,:,1), R.ASC);
                        case 3; [xip2, yip2] = pol2cart( (pi/180)*Filaments(3:end-1,:,1), R.ASC);
                        case 2; [xip2, yip2] = pol2cart( (pi/180)*Filaments(3:end-1,1), R.ASC);    
                    end

                    PTS.ASC1_xip2 = x(F1) + xip2;
                    PTS.ASC1_yip2 = y(F1) + yip2;
                    PTS.ASC1_zip2 = z(3:end-1,1);

                % ACTIN SUBUNIT CENTROID, calculate positions of F2_j+3
                    switch nDIMS
                        case 4; [xjp3, yjp3] = pol2cart( (pi/180)*Filaments(4:end,:,:,2),  R.ASC);
                        case 3; [xjp3, yjp3] = pol2cart( (pi/180)*Filaments(4:end,:,2),  R.ASC);
                        case 2; [xjp3, yjp3] = pol2cart( (pi/180)*Filaments(4:end,2),  R.ASC);    
                    end

                    PTS.ASC2_xjp3 = x(F2) + xjp3;
                    PTS.ASC2_yjp3 = y(F2) + yjp3;
                    PTS.ASC2_zjp3 = z(4:end,1);

    % Calculate distance between F1_i+2 and F2_j+3 
    D.ASC_ip2_jp3 = sqrt( (PTS.ASC2_xjp3 - PTS.ASC1_xip2).^2 +...
                          (PTS.ASC2_yjp3 - PTS.ASC1_yip2).^2 +...
                          (PTS.ASC2_zjp3 - PTS.ASC1_zip2).^2 );  


                % R95a relative to ASC F1_i 
                    switch nDIMS
                        case 4; [xR95a, yR95a] = pol2cart( (pi/180)*(Filaments(1:end-3,:,:,1) + 22.6), R.R95a); 
                        case 3; [xR95a, yR95a] = pol2cart( (pi/180)*(Filaments(1:end-3,:,1)   + 22.6), R.R95a); 
                        case 2; [xR95a, yR95a] = pol2cart( (pi/180)*(Filaments(1:end-3,1)     + 22.6), R.R95a);     
                    end

                    PTS.R95a_x   = x(F1) + xR95a;
                    PTS.R95a_y   = y(F1) + yR95a;
                    PTS.R95a_z   = z(1:end-3,1) + 0.85;

                % D24 relative to ASC F1_i+2 
                    switch nDIMS
                        case 4; [xD24, yD24]  = pol2cart( (pi/180)*(Filaments(3:end-1,:,:,1) + 8.8), R.D24);
                        case 3; [xD24, yD24]  = pol2cart( (pi/180)*(Filaments(3:end-1,:,1)   + 8.8), R.D24);
                        case 2; [xD24, yD24]  = pol2cart( (pi/180)*(Filaments(3:end-1,1)     + 8.8), R.D24);    
                    end

                    PTS.D24_x   = x(F1) + xD24;
                    PTS.D24_y   = y(F1) + yD24;
                    PTS.D24_z   = z(3:end-1,1) - 0.61;

                % R95b relative to ASC F2_j+1
                    switch nDIMS
                        case 4; [xR95b, yR95b] = pol2cart( (pi/180)*(Filaments(2:end-2,:,:,2) + 22.6), R.R95b);
                        case 3; [xR95b, yR95b] = pol2cart( (pi/180)*(Filaments(2:end-2,:,2)   + 22.6), R.R95b);
                        case 2; [xR95b, yR95b] = pol2cart( (pi/180)*(Filaments(2:end-2,2)     + 22.6), R.R95b);    
                    end

                    PTS.R95b_x = x(F2) + xR95b;
                    PTS.R95b_y = y(F2) + yR95b;
                    PTS.R95b_z = z(2:end-2,1) + 0.85;   

                % S350 relative to ASC F2_j+3
                    switch nDIMS
                        case 4; [xS350, yS350] = pol2cart( (pi/180)*(Filaments(4:end,:,:,2) + 16.4),  R.S350);
                        case 3; [xS350, yS350] = pol2cart( (pi/180)*(Filaments(4:end,:,2)   + 16.4),  R.S350);
                        case 2; [xS350, yS350] = pol2cart( (pi/180)*(Filaments(4:end,2)     + 16.4),  R.S350);    
                    end

                    PTS.S350_x = x(F2) + xS350;
                    PTS.S350_y = y(F2) + yS350;
                    PTS.S350_z = z(4:end,1) - 1.74;
                
                % ACTIN SUBUNIT CENTROID (ASC), calculate positions of F1_i-1 
                     switch nDIMS
                         case 4; [xim1, yim1] = pol2cart( (pi/180)*Filaments(1:end-4,:,:,1), R.ASC);
                                    PTS.ASC1_xim1   = x(F1) + cat(1,NaN([1,XYZsize(2:end)]), xim1);
                                    PTS.ASC1_yim1   = y(F1) + cat(1,NaN([1,XYZsize(2:end)]), yim1);
                         case 3; [xim1, yim1] = pol2cart( (pi/180)*Filaments(1:end-4,:,1),   R.ASC);
                                    PTS.ASC1_xim1   = x(F1) + cat(1,NaN([1,XYZsize(2:end)]), xim1);
                                    PTS.ASC1_yim1   = y(F1) + cat(1,NaN([1,XYZsize(2:end)]), yim1);
                         case 2; [xim1, yim1] = pol2cart( (pi/180)*Filaments(1:end-4,1),     R.ASC);  
                                    PTS.ASC1_xim1   = x(F1) + cat(1,NaN, xim1);
                                    PTS.ASC1_yim1   = y(F1) + cat(1,NaN, yim1);
                     end
                     
                     PTS.ASC1_zim1   = [NaN; z(1:end-4,1)];

                % ACTIN SUBUNIT CENTROID, calculate positions of F2_j
                     switch nDIMS
                         case 4; [xj, yj] = pol2cart( (pi/180)*Filaments(1:end-3,:,:,2), R.ASC);
                         case 3; [xj, yj] = pol2cart( (pi/180)*Filaments(1:end-3,:,2),   R.ASC);
                         case 2; [xj, yj] = pol2cart( (pi/180)*Filaments(1:end-3,2),     R.ASC);    
                     end

                     PTS.ASC2_xj = x(F2) + xj;
                     PTS.ASC2_yj = y(F2) + yj;
                     PTS.ASC2_zj = z(1:end-3,1);   
     
                % ACTIN SUBUNIT CENTROID, calculate positions of F1_i+
                     switch nDIMS
                         case 4; [xip3, yip3] = pol2cart( (pi/180)*Filaments(4:end,:,:,1), R.ASC);
                         case 3; [xip3, yip3] = pol2cart( (pi/180)*Filaments(4:end,:,1),   R.ASC);
                         case 2; [xip3, yip3] = pol2cart( (pi/180)*Filaments(4:end,1),     R.ASC);    
                     end

                     PTS.ASC1_xip3 = x(F1) + xip3;
                     PTS.ASC1_yip3 = y(F1) + yip3;
                     PTS.ASC1_zip3 = z(4:end,1);

                % ACTIN SUBUNIT CENTROID, calculate positions of F2_j+4
                    switch nDIMS
                        case 4; [xjp4, yjp4] = pol2cart( (pi/180)*Filaments(5:end,:,:,2),  R.ASC);
                                    PTS.ASC2_xjp4 = x(F2) + cat(1,xjp4,NaN([1,XYZsize(2:end)]));
                                    PTS.ASC2_yjp4 = y(F2) + cat(1,yjp4,NaN([1,XYZsize(2:end)]));
                        case 3; [xjp4, yjp4] = pol2cart( (pi/180)*Filaments(5:end,:,2),    R.ASC);
                                    PTS.ASC2_xjp4 = x(F2) + cat(1,xjp4,NaN([1,XYZsize(2:end)]));
                                    PTS.ASC2_yjp4 = y(F2) + cat(1,yjp4,NaN([1,XYZsize(2:end)]));
                        case 2; [xjp4, yjp4] = pol2cart( (pi/180)*Filaments(5:end,2),      R.ASC);  
                                    PTS.ASC2_xjp4 = x(F2) + cat(1,xjp4,NaN);
                                    PTS.ASC2_yjp4 = y(F2) + cat(1,yjp4,NaN);
                    end

                    PTS.ASC2_xjp4 = x(F2) + cat(1,xjp4,NaN([1,XYZsize(2:end)]));
                    PTS.ASC2_yjp4 = y(F2) + cat(1,yjp4,NaN([1,XYZsize(2:end)]));
                    PTS.ASC2_zjp4 = [z(5:end,1); NaN];


    % Calculate distance between R95a and R95b
    D.R95a_R95b = sqrt( (PTS.R95b_x - PTS.R95a_x).^2 +...
                        (PTS.R95b_y - PTS.R95a_y).^2 +...
                        (PTS.R95b_z - PTS.R95a_z).^2 );  

    % Calculate distance between R95a and S350
    D.R95a_S350 = sqrt( (PTS.S350_x - PTS.R95a_x).^2 +...
                        (PTS.S350_y - PTS.R95a_y).^2 +...
                        (PTS.S350_z - PTS.R95a_z).^2 );  

    % Calculate distance between D24 and S350
    D.D24_S350  = sqrt( (PTS.D24_x - PTS.S350_x).^2 +...
                        (PTS.D24_y - PTS.S350_y).^2 +...
                        (PTS.D24_z - PTS.S350_z).^2 ); 

    % Calculate distance between D24 and R95b
    D.D24_R95b  = sqrt( (PTS.D24_x - PTS.R95b_x).^2 +...
                        (PTS.D24_y - PTS.R95b_y).^2 +...
                        (PTS.D24_z - PTS.R95b_z).^2 );  

    % Calculate distance between F1_i-1 and F2_j
    D.ASC_im1_j  = sqrt( (PTS.ASC1_xim1 - PTS.ASC2_xj).^2 +...
                         (PTS.ASC1_xim1 - PTS.ASC2_yj).^2 +...
                         (PTS.ASC1_xim1 - PTS.ASC2_zj).^2 );  

    % Calculate distance between F1_i+3 and F2_j+4
    D.ASC_ip3_jp4 = sqrt( (PTS.ASC1_xip3 - PTS.ASC2_xjp4).^2 +...
                          (PTS.ASC1_xip3 - PTS.ASC2_yjp4).^2 +...
                          (PTS.ASC1_xip3 - PTS.ASC2_zjp4).^2 );   


        
    % CALCULATE DISTANCE BETWEEN THE FOUR FASCIN POINT'S CENTROID AND THE PLANE DEFINED BY THE TWO FILAMENT AXES
    % Calculate normal vector, n
        F1_mdpt = [x(F1), y(F2), mean(z)];
        R1 = [x(F2)-F1_mdpt(1), y(F2)-F1_mdpt(2),  z(end)-F1_mdpt(3)];
        R2 = [x(F2)-F1_mdpt(1), y(F2)-F1_mdpt(2),  z(1) - F1_mdpt(3)];
        n = cross(R1,R2);
        n = n./vecnorm(n);
    % Get plane equation coefficients.
        A = n(1);
        B = n(2);
        C = n(3);
        d = -( A*F1_mdpt(1) + B*F1_mdpt(2) + C*F1_mdpt(3) );
    % Get center of 4 fascin points 
        F_center_x = (PTS.R95a_x + PTS.D24_x + PTS.R95b_x + PTS.S350_x)/4;
        F_center_y = (PTS.R95a_y + PTS.D24_y + PTS.R95b_y + PTS.S350_y)/4;
        F_center_z = (PTS.R95a_z + PTS.D24_z + PTS.R95b_z + PTS.S350_z)/4;
        
        switch nDIMS
            case 4
                F_center_z = repmat(F_center_z, [1,size(F_center_x,[2,3])] );
            case 3
                F_center_z = repmat(F_center_z, [1,size(F_center_x,[2])] );
        end

        D.Plane2Centroid = abs(A.*F_center_x + B.*F_center_y + C.*F_center_z + d) ./ sqrt(A^2+B^2+C^2);


