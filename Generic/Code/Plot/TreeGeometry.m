function [ NodeLocation ] = TreeGeometry( numStage, numScen, numScenNode )
    X=[];
    X(numStage,:)=linspace(1,numScenNode(numStage),numScenNode(numStage));
    a(1:numScenNode(numStage))=X(numStage,:);
    for s=numStage-1:-1:1
       for i=1:numScenNode(s)
           x=a((i-1)*numScen+1:i*numScen);
           b(i)=mean(x);
           k=(i-1)*numScen^(numStage-s);
           A(k+1:k+numScen^(numStage-s))=repmat(b(i),[1,numScen^(numStage-s)]); 
       end
       a=b;
       b=[];
       X(s,:)=A;
       A=[];
    end
    NodeLocation=X;
end

