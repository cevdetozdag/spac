%  This file is to draw the synthesis 
%  seismogram of flat-layer models

clear;
load ur.dat;
load uf.dat;
load uz.dat;
uf=uf;
ur=ur;
uz=uz;
time=uz(:,1);
tmax=max(time);
rmax=max(ur(:,2));
rmin=min(ur(:,2));
fmax=max(uf(:,2));
fmin=min(uf(:,2));
zmax=max(uz(:,2));
zmin=min(uz(:,2));

figure;
theta=0;
u1=ur(:,2)*cos(theta)-uf(:,2)*sin(theta);
u2=ur(:,2)*sin(theta)+uf(:,2)*cos(theta);
u3=uz(:,2);

subplot(3,1,1);
plot(time,ur(:,2),'b');
ylabel('Ur');
grid;
title('dt=0.0064  m=11  r=1.0')

subplot(3,1,2);
plot(time,uf(:,2),'b');
ylabel('Uf');
%axis([0 tmax 1.5*fmin 3.5*fmax]);
grid;

subplot(3,1,3);
plot(time,uz(:,2),'b');
ylabel('Uz');
%axis([0 tmax 1.5*zmin 3.5*zmax]);
grid;



