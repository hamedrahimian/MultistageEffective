function [ gcf ] = PlotCondEff( numScen, numScenNode, numScenNode_sum, condEff, NodeLocation, axis1 )

    numStage=size(NodeLocation,1);
    numPath=size(NodeLocation,2);
    
    subplot(axis1);
    
    for s=numStage-1:-1:1
        colorspec = cell(1,numPath);
        for i=1:numScenNode(s)
           k=(i-1)*numScen^(numStage-s);
           if (s>1) 
               w=i+numScenNode_sum(s-1);    
           else
               w=1;
           end

           for j=1:numScen
               kk=(j-1)*numScen^(numStage-s-1);
                if (condEff(w,j)==1)
                    for t=1:numScen^(numStage-s-1)
                        colorspec{k+kk+t}=[1,0,0];
                    end
                elseif (condEff(w,j)==0)
                    for t=1:numScen^(numStage-s-1)
                        colorspec{k+kk+t}=[0,1,0];
                    end
                elseif (condEff(w,j)==2)
                    for t=1:numScen^(numStage-s-1)
                        colorspec{k+kk+t}=[0,0,0];
                    end
                end
           end

        end
        
       x=NodeLocation(s:s+1,:);
       y=[repmat(s,[1,numScenNode(numStage)]);repmat(s+1,[1,numScenNode(numStage)])];
       for i=1:numScenNode(numStage)
           gcf=plot(y(:,i),x(:,i),'Color',colorspec{i}, 'LineWidth', 2);
           hold on;
       end
    end
    ax=gca;
    ax.XTick=linspace(1,numStage,numStage);
    ax.XTickLabel=linspace(1,numStage,numStage);
    ax.YTick=[];
    axis tight;
    ax.DataAspectRatio=[1 numPath/(2*(numStage-1)) 1];
    box off;
    xlabel('Stage');
    %legend('Ineffective', 'Effective', 'Location', 'northwestoutside');
   
end

