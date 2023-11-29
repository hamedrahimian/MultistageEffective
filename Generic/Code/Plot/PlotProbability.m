function [ gcf ] = PlotProbability( worst_Prob, gamma, axis2 )

    %figure('Visible','off');
    subplot(axis2);
    gcf=barh(worst_Prob);
    ax=gca;
    ax.XLim=[0, 0.8];
    ax.YTick=[];
    ax.TickDir='in';
    ax.YLim=[0.5, size(worst_Prob,2)+0.5];
    ax.DataAspectRatio=[1 10 1];
    box on;
    name=strcat('$p^{*}_{\xi_{[T]}}: \quad \gamma =', num2str(gamma), '$');
    title(name, 'Interpreter', 'Latex');
    
end

