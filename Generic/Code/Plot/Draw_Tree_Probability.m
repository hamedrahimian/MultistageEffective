%%
numStage=3;
numScen=5;

%% gamma
gamma=[0,0.25,0.5,0.75,1];

%% Worst Probability in the order of nodes
% SGPF3Y3; gamma=0
Prob=[1, 0.216065883783934, 0.2152000957848, 0.175817481824183, 0.177053121822947, 0.215863416784137, ...
    0.21606588294152, 0.215200096096615, 0.175817483522757, 0.17705312051948,...
0.215863416919628, 0.216065884097003, 0.215200094520404, 0.175817481977331, 0.177053122690057, 0.215863416715204, 0.216065885871349, ...
0.215200095972254, 0.175817482131839, 0.17705312148652,...
0.215863414538038, 0.216065882193255, 0.215200096838733, 0.175817481490103, 0.17705312194382, 0.215863417534089, 0.216065883919553, ...
0.21520009571608, 0.175817479994769, 0.177053122438065,...
0.215863417931534];

% SGPF3Y3; gamma=0.5
Prob(2,:)=[1, 0.466065883783934, 0.2152000957848, 0.102870603647129, 0, 0.215863416784137, ...
    0.21606588294152, 0.215200096096615, 0.102870604042237, 0,...
0.465863416919628, 0.466065884097003, 0.215200094520404, 0.102870604667388, 0, 0.215863416715204, ...
0.216065885871349, 0.465200095972254, 0.175817482131839, 0,...
0.142916536024557, 0.216065882193255, 0.465200096838733, 0.175817481490103, 0, 0.142916539477909, ...
0.216065883919553, 0.46520009571608, 0.175817479994769, 0, 0.142916540369599];

% SGPF3Y3; gamma=1.0
Prob(3,:)=[1, 0.716065883783934, 0.0680706994319293, 0, 0, 0.215863416784137, 0.0689364869837573, ...
    0.215200096096615, 0, 0, 0.715863416919628, 0.216065884097003, 0.0680706991877922, 0, ...
    0, 0.715863416715204, 0.0689364894897085, 0.715200095972254, 0, 0, ...
0.215863414538038, 0.068936485627178, 0.715200096838733, 0, 0, 0.215863417534089, 0.0689364863523866, 0.715200095716079, 0, 0,...
0.215863417931534];

% SGPF3Y3; gamma=1.5
Prob(4,:)=[1, 0.784136583215863, 0, 0, 0, 0.215863416784137, 0, 0.0341365830803719, 0, 0,...
0.965863416919628, 0.784136583284796, 0, 0, 0, 0.215863416715204, 0.0347999040277462, 0.965200095972254, 0, 0,...
0, 0, 0.784136582465911, 0, 0, 0.215863417534089, 0.0347999042839205, 0.965200095716079, 0, 0, 0];

% SGPF3Y3; gamma=2
Prob(5,:)=[1, 0.784136583215863, 0, 0, 0, 0.215863416784137, 0, 0, 0, 0,...
1, 0.784136583284796, 0, 0, 0, 0.215863416715204, 0, 1, 0, 0,...
0, 0, 0.784136582465911, 0, 0, 0.215863417534089, 0, 1, 0, 0,0];

%% Effectiveness of Paths in the order of paths
pathEff= [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]; %gamma=0.0
pathEff(2,:)= [1,1,1,0,1,1,1,1,0,1,1,1,1,0,1,0,0,0,0,0,1,1,1,0,1]; %gamma=0.5
pathEff(3,:)= [1,1,0,0,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1]; %gamma=1.0
pathEff(4,:)= [0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0]; %gamma=1.5
pathEff(5,:)= [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0]; %gamma=2.0

%% Conditional Effectiveness: each row is a parent node
condEff=[1,1,1,1,1;1,1,1,1,1;1,1,1,1,1;1,1,1,1,1;1,1,1,1,1;1,1,1,1,1]; %gamma=0
condEff(:,:,2)= [1,1,1,0,1;1,1,1,0,1;1,1,1,0,1;1,1,1,0,1;1,1,1,0,1;1,1,1,0,1]; %gamma=0.5
condEff(:,:,3)= [1,1,0,0,1;1,1,0,0,1;1,1,0,0,1;1,1,0,0,1;1,1,0,0,1;1,1,0,0,1]; %gamma=1.0
condEff(:,:,4)= [1,0,0,0,1;1,0,0,0,1;1,0,0,0,1;1,1,0,0,0;0,1,0,0,1;1,1,0,0,0]; %gamma=1.5
condEff(:,:,5)= [1,0,0,0,1;0,0,0,0,1;1,0,0,0,1;0,1,0,0,0;0,1,0,0,1;0,1,0,0,0]; %gamma=2.0

%% Tree Structure
[ numScenNode,numScenNode_sum,descendant,ancestor ] = TreeStructure( numStage, numScen );

%% Node Locations
[ NodeLocation ] = TreeGeometry( numStage, numScen, numScenNode );

% for i=7:16
%     NodeLocation(i,:)=[1:1:32];
% end
% figure;
% plot(NodeLocation, 'b', 'LineWidth',2);
% ax=gca;
% ax.XTick=linspace(1,numStage+10,numStage+10);
% ax.XTickLabel=linspace(1,numStage+10,numStage+10);
% ax.YTick=[];
% axis tight;
% ax.DataAspectRatio=[1 size(NodeLocation,2)/(2*(numStage+10-1)) 1];
% box off;
% xlabel('Stage');
%% 
save_name=strcat('SGPF3Y3_');
for i=1:length(gamma)
    P=Prob(i,:);
    cEff=condEff(:,:,i);
    pEff=pathEff(i,:);
    
    worst_probability=CalcProbability( numStage, numScenNode(numStage), numScenNode_sum, ancestor, P );        
    
    figure('Visible','off');
    gca_Tree=subplot(1,2,1);
    gca_Prob=subplot(1,2,2);
    
    gcf_CondTree= PlotCondEff( numScen, numScenNode, numScenNode_sum, cEff, NodeLocation, gca_Tree );
%     save_name_f1=strcat('T', save_name, num2str(i), '_', num2str(gamma(i)), '_','tree.png');
%     save_name_f2=strcat('T', save_name, num2str(i), '_', num2str(gamma(i)), '_','tree.eps');
%     saveas(gcf_CondTree, save_name_f1);
%     saveas(gcf_CondTree, save_name_f2);
    
    gcf_Prob = PlotProbability( worst_probability, gamma(i), gca_Prob );
    
    save_name_f1=strcat(save_name, num2str(i), '_', num2str(gamma(i)), '_','tree_prob.png');
    save_name_f2=strcat(save_name, num2str(i), '_', num2str(gamma(i)), '_','tree_prob.eps');
    saveas(gcf_Prob, save_name_f1);
    saveas(gcf_Prob, save_name_f2);
    
   
    
    %gcf_PathTree=PlotPathEff( pEff, NodeLocation );
    %save_name_f1=strcat('03', save_name, num2str(gamma(i)), '_','tree.png');
    %save_name_f2=strcat('03', save_name, num2str(gamma(i)), '_','tree.eps');
    %saveas(gcf_PathTree, save_name_f1);
    %saveas(gcf_PathTree, save_name_f2);
    
end


