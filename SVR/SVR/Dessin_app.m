function Dessin_apprentissage
close all

% fonction écrite pour les SVR
% dessin des points d'apprentissage et du modèle "appris"
% Suivant la dimension de l'espace de dessin 

%données X et Y Les points Xi Yi pour la regression
% Xi est une colonne de Rn ( f cherchée de Rn dans R)
% Les Xi sont stockés dans X en ligne
global ikernel
ikernel=2;
load Modele X Y alpha b

ndim=size(X,1);
npoints=size(X,2);

if ndim==1
    plot(X(1,:),X(2,:),Y,'o')
else
plot3(X(1,:),X(2,:),Y,'o')
axis equal
hold on
end

% dessin des points 

%tracé de la fonction obtenue par régression
if ndim==1
    
else
xmin=min(X(1,:));
xmax=max(X(1,:));
ymin=min(X(2,:));
ymax=max(X(2,:));
delx=xmax-xmin;
dely=ymax-ymin;
xmax=xmax+delx/10;
ymax=ymax+dely/10;
xmin=xmin-delx/10;
ymin=ymin-dely/10;
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
         TXS1(ikx,iky)=XS(1,ik);
         TXS2(ikx,iky)=XS(2,ik);
         TYS(ikx,iky)=YS(ik);
     end
 end

grid
 surf(TXS1,TXS2,TYS) 
grid

hold on 
load Loiremini
np=length(xs);
xs(np+1)=xs(1);
ys(np+1)=ys(1);
maxy=max(YS);
zs=maxy*ones(1,np+1);
axis([ xmin xmax ymin ymax])
plot3(xs,ys,zs,'linewidth',3)
end