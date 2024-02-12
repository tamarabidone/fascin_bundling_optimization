function [x,y] = CreateHexagonArray(nPts,R)
    
    %R = 12.16;
    % Create hexagram first layer (7 points);
    theta1   = (0:60:300)';
    radius1  = R*ones(length(theta1),1);
    [x1,y1] = pol2cart(theta1*(pi/180),radius1); % Points 2-7

    theta2   = (0:30:330)';
    radius2  = ones(length(theta2),1);
    radius2(1:2:end) = 2*R*radius2(1:2:end);
    radius2(2:2:end) = (sqrt(3)/2)*radius2(1:2:end); 
    [x2,y2] = pol2cart(theta2*(pi/180),radius2);
    % add center point
    x = [0;x1;x2];
    y = [0;y1;y2]; 

    % Build 3rd layer of points (cound't figure out a form so I just built it manually)
    xc = R/2;
    yc = R*sqrt(3)/2;
    x = [x; x(8)+2*xc; x(8)+xc; x(9)+xc; x(10)+xc; x(10)-xc; x(11)-xc; x(12)-xc; x(13)-xc; x(14)-xc; x(14)-2*xc;...
            x(15)-2*xc; x(15)-xc; x(16)-xc; x(16)+xc; x(17)+xc; x(18)+xc; x(19)+xc; x(8)+xc ];
    y = [y; y(8); y(8)+yc; y(9)+yc; y(10)+yc; y(10)+yc; y(11)+yc; y(12)+yc; y(13)+yc; y(14)+yc; y(14);...
            y(15) ; y(15)-yc; y(16)-yc; y(16)-yc; y(17)-yc; y(18)-yc; y(19)-yc; y(8)-yc ];
    
    % Add 4th layer manually
    x = [x; x(20)+2*xc; x(20)+xc;   x(21)+xc;   x(22)+xc;   x(23)+xc; x(24)+xc; x(25)+xc; x(26)+xc; x(26)-xc;...
            x(26)-2*xc; x(27)-2*xc; x(28)-2*xc; x(29)-2*xc; x(30)-2*xc; x(31)-2*xc; x(32)-2*xc; x(32)-xc;...
            x(33)-xc  ; x(33)+xc;   x(34)+xc;   x(35)+xc;   x(35)+2*xc; x(36)+2*xc; x(37)+2*xc];
    y = [y; y(20);      y(20)+yc;   y(21)+yc;   y(22)+yc;   y(23)+yc; y(24)+yc; y(25)+yc; y(26)+yc; y(26)+yc;...
            y(26);      y(27);      y(28);      y(29);      y(30);      y(31);      y(32);      y(32)-yc;...
            y(33)-yc;   y(33)-yc;   y(34)-yc;   y(35)-yc;   y(35);      y(36);      y(37)];

    x = x(1:nPts);
    y = y(1:nPts);
    
%     Plot and label locations
%     figure(1); clf
%     set(gcf,'Color','w')
%     plot(x,y,'-r','LineWidth',2);    hold on
%     plot(x,y,'.b','MarkerSize',120); hold off
%     N = length(x);
%     axis equal
%     axis(55*[-1,1,-1,1])
%     text(x,y,compose('%0.0f',(1:N)'),'Fontsize',20,'HorizontalAlignment','center','Color','w'); hold off
%     set(gca,'FontSize',30,'LineWidth',2); grid on; xlabel('X (nm)'); ylabel('Y (nm)')
%     yticks(-60:10:60)
%     xticks(yticks)
%     xtickangle(45)


