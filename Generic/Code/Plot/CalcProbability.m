function [ prob ] = CalcProbability( numStage, numPath, numScenNode_sum, ancestor, P )
    
    prob=ones(1,numPath);
    k=numScenNode_sum(numStage - 1);
    for w2=k+1:numScenNode_sum(numStage)
        w_node=w2;
        for s=numStage-1:-1:1
            prob(w2-k)=prob(w2-k)*P(w_node);
            w_ancestor=FindAncestor(w_node, s+1, ancestor, numScenNode_sum, 1);
            w_node=w_ancestor;
        end
        
    end
    
end

