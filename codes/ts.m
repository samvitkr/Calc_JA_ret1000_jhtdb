for time=2:4
ft=sprintf("Transfer_%03d.mat",time);
mt=matfile(ft,'Writable',true) ;
mt.convective_x=single(mt.convective_x);
mt.viscous_x=single(mt.viscous_x);
mt.total=single(mt.total);
end
