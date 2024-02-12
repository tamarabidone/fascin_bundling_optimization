function PlotPoints01(x,y)

    N = length(x); 
    % Plot and label locations
    figure(1); clf
    set(gcf,'Color','w')

    L1 = 1;
    L2 = 2:7;
    L3 = 8:19;
    L4 = 20:37;
    L5 = 38:61;
    
    ColorMap = lines(5);
    ColorMap(3,:) = [];
    ColorMap = [[0,0,0];ColorMap];

    for n = 1:N
        if any(ismember(n,L1))
            Color = ColorMap(1,:);
        elseif any(ismember(n,L2))
            Color = ColorMap(2,:);
        elseif any(ismember(n,L3))
            Color = ColorMap(3,:);
        elseif any(ismember(n,L4))
            Color = ColorMap(4,:);
        else
            Color = ColorMap(5,:);
        end
        plot(x(n),y(n),'.','MarkerSize',120,'Color','b'); hold on
    end

    pad = 5;
    axis equal
    axis([min(x)-pad,max(x)+pad,min(y)-pad,max(y)+pad])

    text(x,y,compose('%0.0f',(1:N)'),'Fontsize',20,'HorizontalAlignment','center','Color','w'); hold off
 
    set(gca,'FontSize',20,'LineWidth',2); grid on; 
    xlabel('X (nm)'); 
    ylabel('Y (nm)')
    %xticks(yticks)
    xtickangle(45)

end