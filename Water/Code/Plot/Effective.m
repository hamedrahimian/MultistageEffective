
figure ('Visible', 'off');
X=[0; 0.25; 0.5; 0.75; 1];
% E= [Ineffective, Effective]
E= [0, 125000;70128,54872;109374,15626;122803,2197;124999,1];
gcf = bar(X, E,'FaceColor','flat');
labels = arrayfun(@(value) num2str(value),E(:,2),'UniformOutput',false);
h=text(X,E(:,2),labels,'HorizontalAlignment','left','VerticalAlignment','top');
set(h,'Rotation',90);
map=[0,1,0;1,0,0];
colormap(map);
ax=gca;
legend('Ineffective', 'Effective', 'Location', 'North');
ax.XLim=[-0.1, 1.1];
ax.YLim=[0, 15*10^4];
ax.XTick=X;
ax.XTickLabel=X;
ax.TickDir='out';
%box off;
xlabel('$\gamma$', 'Interpreter', 'Latex');

save_name=strcat('WATER_Eff');
save_name_f1=strcat(save_name, '.png');
save_name_f2=strcat(save_name, '.eps');
saveas(ax, save_name_f1);
saveas(ax, save_name_f2);