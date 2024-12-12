dnu=1e-3;
ut=0.0499;
yp=(yv+1)./dnu;

% %%
% eloc=-sqrt(edges(1:end-1).*edges(2:end));
% [peak,ind]=max(pconv(1:100,:));
% 
% %%
% %f=fit(eloc',peak','poly1');
% yscale=peak;
% lscale=abs(((yscale).^(-1)));
% for i =1:19
% yl(:,i)=yp./yp(ind(i));
% end
% pcn=pconv.*lscale;
% subplot(1,2,1)
% plot(yl(2:200,6:3:18),pcn(2:200,6:3:18))
% set(gca,'Xscale','log')
% %xlim([1 1000])
% 
%  subplot(1,2,2)
%  plot(yp(2:200),pconv(2:200,6:3:18))
% set(gca,'Xscale','log')
% xlim([1 1000])
% %%
