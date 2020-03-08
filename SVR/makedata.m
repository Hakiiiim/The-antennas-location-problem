% création de données en 2D dans le carré [0, 10]^2

close all
clear all
xmin=0;
ymin=0;
xmax=10;
ymax=10;


nb1=1; % nombre de points d'apprentissage du premier groupe (vous pouvez modifier)
nb2=1; % nombre de points d'apprentissage du second groupe
  figure( 1) 
  axis([xmin xmax ymin ymax]);
  grid
  hold on

for i=1:nb1
   X1(:,i)=ginput(1)
   plot(X1(1,i),X1(2,i),'x','linewidth',2,'markersize',12)
   axis([xmin xmax ymin ymax])
   pause(0.1)
   hold on
end
for i=1:nb2
   X2(:,i)=ginput(1)
   plot(X2(1,i),X2(2,i),'o','linewidth',2,'markersize',12)
   axis([xmin xmax ymin ymax])
   pause(0.1)
   hold on
end
axis equal
X=[X1 X2];
lab=[-ones(nb1,1);ones(nb2,1)];

save 'data00' X lab
