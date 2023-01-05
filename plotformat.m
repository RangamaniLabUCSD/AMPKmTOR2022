set(gca,'fontname','cmss10')

set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','bold')
set(0,'defaultAxesFontSize', 2)
set(findall(gca, 'Type', 'Line'),'LineWidth',2)
set(gcf, 'Position',  [100, 100, 900, 600])
set(gca,'linewidth',2)

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',14,'FontWeight','bold')