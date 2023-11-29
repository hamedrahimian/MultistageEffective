
figure ('Visible', 'off');
%X=[0; 0.25;0.45; 0.5; 0.55; 0.6; 0.65; 0.7;0.75;1.0];
X=[0; 0.25; 0.5; 0.75; 1.0];
% E= [Ineffective, Effective]
E= [0	3125	0;
    2101	1024	0;
    2829	214	82;
    3089	36	0; 
    3116 9 0];
gcf = bar(X, E,'FaceColor','flat');
labels = arrayfun(@(value) num2str(value),E(:,2),'UniformOutput',false);
h=text(X,E(:,2),labels,'HorizontalAlignment','left','VerticalAlignment','top');
set(h,'Rotation',90);
map=[0,1,0;1,0,0;0,0,0];
colormap(map);
ax=gca;
legend('Ineffective', 'Effective', 'Unknown', 'Location', 'North');
ax.XLim=[-0.1, 1.1];
ax.YLim=[0, 3700];
ax.XTick=X;
ax.XTickLabel=X;
ax.TickDir='out';
%box off;
xlabel('$\gamma$', 'Interpreter', 'Latex');

save_name=strcat('SGPF3Y6_Eff');
save_name_f1=strcat(save_name, '.png');
save_name_f2=strcat(save_name, '.eps');
saveas(ax, save_name_f1);
saveas(ax, save_name_f2);