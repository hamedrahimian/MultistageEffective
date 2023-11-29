function [ gcf ] = PlotPathEff( pathEff, NodeLocation )

    numStage=size(NodeLocation,1);
    numPath=size(pathEff,2);
    
    E=NodeLocation;
    k=0;
    for i=1:numPath
        if pathEff(i)==0
           k=k+1;
           I(:,k)=NodeLocation(:,i);
           d(k)=i; 
        end
    end
    if (k>0)
        E(:,d)=[];
    end

    figure('Visible','off');
    gcf=plot(E, 'r', 'LineWidth', 2);
    hold on;
    if (k>0)
        gcf=plot(I, '--b', 'LineWidth', 2);
    end
    ax=gca;
    ax.XTick=linspace(1,numStage,numStage);
    ax.XTickLabel=linspace(1,numStage,numStage);
    ax.YTick=[];
    axis tight;
    ax.DataAspectRatio=[1 numPath/(2*(numStage-1)) 1];
    box off;
    xlabel('Stage');


end

