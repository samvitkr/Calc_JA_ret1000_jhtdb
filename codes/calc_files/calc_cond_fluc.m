close all
clear
baseDir =  '/data/geyink1/skumar67/Calc_JA_ret1000_jhtdb/data';
load(fullfile(baseDir,"bsplinedata.mat"),'yv')
y=yv(1:256)+1;
yp=y.*1000;
jcset=[47 54 71 105 130];
mp = matfile(fullfile(baseDir,'mean_profiles.mat'));
jj=3;
%for jj=4:4
    jcond=jcset(jj);
     fvgn=fullfile(baseDir,sprintf("conditionaln_jcond_%03d.mat",jcond));
     fvgp=fullfile(baseDir,sprintf("conditionalp_jcond_%03d.mat",jcond));

     mvgp = matfile(fvgp,'Writable',true);
     mvgn = matfile(fvgn,'Writable',true);
     dudyp=mvgp.dudy;
     dudyn=mvgn.dudy;
     up=mvgp.u;
     un=mvgn.u;


     % plot(yp,mp.U)
     % set(gca,'Xscale','log')
     % xlim([10 1000])

%      dudy=mvg.dudy;
% 
%      figure
%      hold on
%      plot(mp.dUdy,'-k')
%      plot(dudymc,'-.b')
%      xline(jcond)
%      hold off