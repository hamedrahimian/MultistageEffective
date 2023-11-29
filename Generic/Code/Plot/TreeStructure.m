function [ numScenNode,numScenNode_sum,descendant,ancestor ] = TreeStructure( numStage, numScen )
    
    numScenNode=zeros(1,numStage);
    for s=1:numStage
        numScenNode(s)=numScen^(s-1);
    end

    numScenNode_sum=zeros(1,numStage);
    numScenNode_sum(1)=1;
    for s=2:numStage
        numScenNode_sum(s)=numScenNode_sum(s-1)+numScenNode(s);
    end
    
    descendant=zeros(numScenNode_sum(numStage - 1),numScen);
    for w = 1:numScenNode_sum(numStage - 1)
        for j = 1:numScen
            descendant(w,j) = 1 + numScen + numScen*(w - 2) + j;
        end
    end
    
    ancestor=zeros(numScenNode_sum(numStage - 1),numScenNode_sum(numStage));
    for t = 1:numStage-1
		for w = 1:numScenNode_sum(numStage - t)
            k=numScenNode_sum(t+1) + numScenNode(t+1)*(w - 2);
			for j = k+1:k+numScenNode(t+1)
				ancestor(w,j) = t; %w is ancestor j with "t" step // s=0 means w is not ancestor j
            end
        end
    end
    
    
end

