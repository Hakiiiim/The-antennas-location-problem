% Séparateurs à Vaste Marge pour la Classification
% Cas linéaire, non linéaire, marges "souples"
clear all
close all
global ikernel  % type de noyau choisi.  1: produit scalaire dans Rp (linéaire)
                                         % 2: noyau gaussien
prediction = false;
Analytic = false;
load data4 X lab     % chargement des données 
ikernel=2;             % choix du type de noyau
marge_souple = true;  % choix marge souple ou pas
C_souple= 10;          % coef souplesse de la marge souple

% bornes pour le dessin 2D
xmin=min(X(1,:));
ymin=min(X(2,:));
xmax=max(X(1,:));
ymax=max(X(2,:));

na=length(X); % nombre de points (d'apprentissage)

% dessin des points dans R^2
subplot(1,2,1)
for i=1:na
    if lab(i)==1
        plot(X(1,i),X(2,i),'o','linewidth',2)
        hold on
    else
        plot(X(1,i),X(2,i),'x','linewidth',2,'markersize',12)
        hold on
    end
end
axis([xmin xmax ymin ymax])
grid
axis equal
hold on

% Coeur du Programme - Méthode d'UZAWA
% assemblage matrice de la forme quadratique du problème dual

A = zeros(na);
%Ne fait intervenir que le produit scalaire (kernel) des points
for i=1:na
    for j=1:na
        A(i,j) = lab(i)*lab(j)* kernel(X(:,i),X(:,j));
    end
end
        

% gradient à pas constant pour problème dual
% On vous laisse quelques valeurs pour les paramettres reglables... à vous de voir
alpha0=0.5*ones(na,1); % point de départ (0, ou autre): réglable
%1000
pasgrad=5e-3; % pas du gradient : parametre réglable
%9e-2
u=ones(na,1);         % vecteur de 1
%crit_arret=1;         % initialisation critere d'arret
npas=0;               % comptage du nombre de pas
npas_max=100000;      % garde-fou convergence : on arrete si le nombre d'iterations est trop grand
epsi=1e-5; % seuil convergence
% 1e-1

% boucle de gradient projeté
% A vous de jouer!

%Initialisation des différentes grandeurs d'intérêt
alpha = alpha0 + pasgrad*(-A*alpha0+u);
alpha = alpha - (dot(lab/norm(lab),alpha))*(lab/norm(lab));
alpha = max(0,alpha);

crit_arret = norm(alpha-alpha0);

vit_conv = zeros(1,npas_max);

k=1;
while crit_arret >= epsi
    vit_conv(k) = crit_arret; %traçage de l'évolution du critère d'arrêt
    alpha0=alpha; %actualisation de la valeur de alpha0
    %On bouge suivant l'opposé du gradient
    alpha = alpha0 + pasgrad*(-A*alpha0+u);
    %on se projete sur l'hyperplan recherché
    alpha = alpha - (dot(lab/norm(lab),alpha))*(lab/norm(lab));
    alpha = max(0,alpha); %on se projete sur le demi espace positif
    %Introduction de la marge souple si demandé
    if marge_souple == true 
        alpha = min(alpha,C_souple);
    end
    %Actualisation du critère d'arrêt
    crit_arret = norm(alpha-alpha0);
    k=k+1; 
    if k==npas_max
        break %Vérification du nombre d'itérations maximal
    end
end
       

% recherche des points supports
epsi=1e-5;
Points_support = find(alpha>epsi)' ; % voir help find

% on ne se sert plus directement de w, mais de son produit scalaire 
% avec un vecteur, defini par un noyeau : fonction prodw

% Calcul de b (on le calcule pour chaque pt support, puis on moyenne)
% chaque pt support devrait fournit le même b (mais aux erreurs
% numériques près, c'est pouquoi on moyenne)

moy_b = zeros(1,length(Points_support)); %initialisation
for i=1:length(Points_support)
    indice = Points_support(i);
    l = lab(indice);
    moy_b(i) = (1/l) - prodw(X,lab,alpha,X(:,indice));
end
b=mean(moy_b); %on moyenne les b des points supports

% Fin du coeur du programme

%calcul et tracé des isovaleurs
xp=xmin:0.2:xmax;   % création d'une grille pour les besoins de contour
yp=ymin:0.2:ymax;
npx=length(xp);
npy=length(yp);
for i=1:npx
    for j=1:npy
        ps=prodw(X,lab,alpha,[xp(i),yp(j)]'); % calcul de <w,x> + b sur une grille
        V(i,j)=ps + b;   % on n'a pas besoin explicitement de w, mais de son 
                         % produit scalaire avec tout vecteur qu'on
                         % encapsule dans prodw (utilisation noyau si cas non
                         % linéaire
    end
end
contour(xp,yp,V',[-1 0 1],'linewidth',2,'color','r')
axis([xmin xmax ymin ymax])
title('Support Vector Machine (two points to predict)')
grid

if Analytic == true
    %Solution analytique dans le cas de deux points
    xk1=X(1,1);xk2=X(1,2);yk1=X(2,1);yk2=X(2,2);
    %Point milieu
    x01=(xk1+xk2)/2;y01=(yk1+yk2)/2;
    %équation de la droite médiatrice
    a=(yk2-yk1)/(xk2-xk1);
    b=yk1-a*xk1;
    b=y01+x01/a;
    
    x=xmin:0.1:xmax;
    y=-1/a*x+b;

    line([xk1 xk2],[yk1 yk2],'LineStyle','--')
    hold on
    plot(x,y,'k','linewidth',1)
    legend('classe 1','classe 2','SVM','segment','Solution analytique');
end

%Deux points à prédire
if prediction == true
    for i=1:2
        X_test=ginput(i);
        ps_test=prodw(X,lab,alpha,X_test');
        V_test=ps_test + b;
        if V_test >= 0
            plot(X_test(1),X_test(2),'o','linewidth',2,'markersize',12,'Color','k')
        else
            plot(X_test(1),X_test(2),'x','linewidth',2,'markersize',12,'Color','k')
        end
    end
end

subplot(1,2,2)
plot(1:k,log(vit_conv(1:k)),'linewidth',2)
xlabel(' nombre itérations')
ylabel(' log (ecart)')
title('comportement convergence')
grid

