clear
close all

ut = 0.0499;
dnu=1.006e-3;
jc=[47 54 71 105 130]

figure
subplot(1,2,1)
hold on
lt={'-','-.','--',':','-'}
pdfozp=[];
for i=1:5
    jcond=jc(i);
    mn=sprintf("../data/pdfoz_j_%03d.mat",jcond)
    m=matfile(mn);
bc=m.edges(1,1:end-1)+ (m.edges(1,4)-m.edges(1,3))*0.5;
dbc=bc(2)-bc(1);
pdfoz=m.pdfoz;
pdfozn(i)=sum(pdfoz(bc<-0.025))*dbc
pdfozp(i)=sum(pdfoz(bc>0.006))*dbc
plot(bc,m.pdfoz,'LineStyle',lt(i),'LineWidth',2)


end
hold off
set(gca,'Yscale','log')
xline(0)
legend('y_c^+=39','y_c^+=53','y_c^+=93','y_c^+=197','y_c^+=297','Location','northeast')
legend boxoff
ylabel('pdf(\omega_z)')
xlabel('\omega_z')
ycp=[39,53,93,197,297]
ylim([0.003 1])

subplot(1,2,2)
plot(ycp,100*pdfozp,'-o','LineWidth',2)
ylabel('% retrograde \omega_z')
xlabel('y_c^+')
