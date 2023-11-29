function [w]=FindAncestor(scenario, stage, ancestor, numScenNode_sum, s)
    w2=0;
    if (stage>0 && stage >= s)
		IsAncestor = false;
		if (stage == 2)
			w2 = 1;
        else
            w2 = numScenNode_sum(stage - s -1 )+1;
        end
		
		while (w2<=numScenNode_sum(stage - s) && ~IsAncestor)
			if (ancestor(w2,scenario) == s)
				IsAncestor = true;
				w=w2;
            end
			w2=w2+1;
        end
    else
        w=1;
    end

end

