clear
close all
x1=100;
y1=100;
width=800;
height=350;
ut = 0.0499;
dnu=1.006e-3;
jc=[47 54 71 105 130]
f=figure('OuterPosition',[x1,y1,width,height]);

subplot(1,2,1)
hold on
lt={'-','-.','--',':','-'}
pdfozp=[];
for i=1:5
    jcond=jc(i);
    mn=sprintf("../../data/pdfoz_j_%03d.mat",jcond)
    m=matfile(mn);
bc=m.edges(1,1:end-1)+ (m.edges(1,4)-m.edges(1,3))*0.5;
dbc=bc(2)-bc(1);
pdfoz=m.pdfoz;
 % pdfoz(pdfoz>1)=NaN;
pdfozn(i)=sum(pdfoz(bc<-0.005*m.edges(1,end)))*dbc;
pdfozp(i)=sum(pdfoz(bc>0.005*m.edges(1,end)))*dbc;
plot(bc./(ut/dnu),pdfoz,'LineStyle',lt(i),'LineWidth',2)


end
hold off
set(gca,'Yscale','log')
xline(0)
legend('y^+=39','y^+=53','y^+=93','y^+=197','y^+=297','Location','northwest')
legend boxoff
ylabel('pdf(\omega_z^+)')
xlabel('\omega_z^+')
ycp=[39,53,93,197,297]
ylim([0.003 10])
xlim([-0.4 0.3])
xticks([-0.4:0.1:0.3])
subplot(1,2,2)
plot(ycp,100*pdfozp,'-ok','LineWidth',2)
ylabel('Percent retrograde \omega_z')
xlabel('y^+')
ylim([0 40])

saveas(f,'pdfoz.fig')