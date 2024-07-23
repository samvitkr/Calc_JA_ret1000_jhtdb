matObj = matfile('myFile.mat');
%m = matfile('myFile2.mat','Writable',true);
save('myFile.mat','-v7.3')
m.y(81:100,81:100,3) = magic(20);
