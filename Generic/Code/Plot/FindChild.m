function [ w ] = FindChild( scenario,ancestor, descendant, stage )

    j = 1;
	IsChild = false;
	w_ancestor = FindAncestor(scenario, stage, ancestor,numScenNode_sum,1);
	while (j<=numScen && ~IsChild)
		if (descendant(w_ancestor,j) == scenario)
			IsChild = true;
			w=j;
        end
		j=j+1;
    end
    
end

