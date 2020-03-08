
function Dessin_Modele
global ikernel
global maxy
ikernel=2;
%tracé de la fonction obtenue par régression. On dessine la Loire au dessus

load Modele X Y alpha b

xmin=min(X(1,:));
xmax=max(X(1,:));
ymin=min(X(2,:));
ymax=max(X(2,:));
% delx=xmax-xmin
% dely=ymax-ymin
% xmax=xmax+delx/10
% ymax=ymax+dely/10;
% xmin=xmin-delx/10;
% ymin=ymin-dely/10;
ikx=0;
ik=0;

for ix=xmin:0.5:xmax
     ikx=ikx+1;
     iky=0;
     for iy=ymin:0.3:ymax
         ik=ik+1;
         iky=iky+1;
         XS(1,ik)=ix;
         XS(2,ik)=iy;
         YS(ik)= prodwr(X,alpha,XS(:,ik))+b;
         TXS1(ikx,iky)=ix;
         TXS2(ikx,iky)=iy;
         TYS(ikx,iky)=YS(ik);
     end
 end

grid
surf(TXS1,TXS2,TYS) 
shading interp

%Optimisation locale
p0 = [0,0];
p1 = [-2,2];
p2 = [-1,0.3];

fun = @(p)pred(p); %regarder prochain paragraphe pour la fonction pred
[P,fval] = fminsearch(fun,p0);
[P1,fval1] = fminsearch(fun,p1);
[P2,fval2] = fminsearch(fun,p2);

disp(fval)
disp(fval1)
disp(fval2)

hold on 
load Loiremini
np=length(xs);
xs(np+1)=xs(1);
ys(np+1)=ys(1);
maxy=max(YS);
zs=maxy*ones(1,np+1);
axis([ xmin xmax ymin ymax])
plot3(xs,ys,zs,'linewidth',3)
axis equal
% view (2);
hold on

plot3([p0(1),P(1)],[p0(2),P(2)],[pred(p0),fval],'r--','linewidth',3)

hold on 

plot3(P(1),P(2),fval,'r*','linewidth',3)

hold on

plot3(p0(1),p0(2),pred(p0),'k*','linewidth',3)
hold on

plot3([p1(1),P1(1)],[p1(2),P1(2)],[pred(p1),fval1],'r--','linewidth',3)

hold on 

plot3(P1(1),P1(2),fval1,'r*','linewidth',3)

hold on

plot3(p1(1),p1(2),pred(p1),'k*','linewidth',3)
hold on

plot3([p2(1),P2(1)],[p2(2),P2(2)],[pred(p2),fval2],'r--','linewidth',3)

hold on 

plot3(P2(1),P2(2),fval2,'r*','linewidth',3)

hold on

plot3(p2(1),p2(2),pred(p2),'k*','linewidth',3)
