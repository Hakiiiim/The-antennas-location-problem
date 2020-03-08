% dessin Loire

load Loiremini
np=length(xs);
xs(np+1)=xs(1);
ys(np+1)=ys(1);
maxy=30;
%maxy=max(YS);
zs=maxy*ones(1,np+1);
%axis([ xmin xmax ymin ymax])
plot3(xs,ys,zs,'linewidth',3)
axis equal
view (2);
hold on