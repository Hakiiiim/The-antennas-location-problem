% Séparateurs à Vaste Marge pour régression

clear all
close all
global ikernel  % type de noyau choisi. 1: produit scalaire dans Rp
                                       % 2: noyau gaussien
                                       
% Les parametres des SVR 
 epsilon=0.5;           % demi largeur de la bande
 ikernel=2;             % choix du type de noyau
 C_ecart= 50;           % coef pour l'écart permis 
 
 
%données X et Y Les points Xi Yi pour la regression
% Xi est  une colonne de Rn ( fonction cherchée de Rn dans R)
% Les Xi sont stockés dans X 
% Les Y sont stockés en ligne dans Y 

load antennes_2d_train  C  S  % Dans le fichier antenne, C contient X' et 
                              % S contient Y'
X=C';
Y=S';
npoints=size(X,2);    % nombre de points d'apprentissage

% dessin des points en 3D ( à modifier si essais avex X dans R   

plot3(X(1,:),X(2,:),Y, 'o')
axis equal
hold on

%
% -----------------------------------------------------------
% Partie que les élèves doivent programmer

% Méthode UZAWA gradient à pas constant pour problème dual
na=2*npoints; % nombre de multiplicateurs alpha 
u1=epsilon*ones(na,1);  % vecteurs colonne
u2=[Y'; -Y'];
u3=[ones(npoints,1); -ones(npoints,1)];

% assemblage matrice de la forme quadratique du problème dual, A

B = zeros(npoints);

for i=1:npoints
    for j=1:npoints
        B(i,j) = kernel(X(:,i),X(:,j));
    end
end

A = [B -B ; -B B];


% A faire....

% quelques initialisations pour la boucle de gradient
alpha0=2.5*abs(rand(na,1)); % point de départ (0, ou autre): réglable
pasgrad=1e-2;         % pas du gradient : parametre réglable
crit_arret=1;         % initialisation critere d'arret
npas=0;               % comptage du nombre de pas
npas_max=100000;      % garde-fou convergence
epsi=1e-4;            % seuil convergence

vit_conv = [];
while and(crit_arret>epsi,npas<npas_max)% boucle de gradient
    npas=npas+1;
    % gradient pas constant 
    % projection d' ABORD sur l'hyperplan
    alpha = alpha0 - pasgrad*(A*alpha0-u1+u2);
    % projection sur  alpha > 0
    alpha = alpha - (dot(u3/norm(u3),alpha))*(u3/norm(u3));
    % projection sur  alpha<C
    alpha = max(0,alpha);
    alpha = min(alpha,C_ecart);

    % stockage dans vit_conv de la distance entre deux alphas successifs par
    % exemple (donne une idée de la vitesse à laquelle la suite des alpha se
    % tasse (critère de Cauchy)
    crit_arret = norm(alpha-alpha0);
    
    vit_conv(npas) = crit_arret;
    
    alpha0 = alpha;
end

fprintf(' nb iterations: %d critère: %e ',npas,crit_arret)

Points_support=find(alpha>epsi)';

% on ne se sert pas directement de w, mais de son produit scalaire 
% avec un vecteur, defini par un noyeau : fonction prodw

% Calcul de b (on le calcule pour chaque pt support, puis on moyenne)
% on devrait vérifier que chaque pt support fournit le même b (aux erreurs
% numériques près, c'est pouquoi on prend la moyenne)

% calcul de la constante b
% a faire 

moy_b = [];

k=1;
for i=1:length(Points_support)
    indice = Points_support(i);
    alph = alpha(indice);
    
    if indice > npoints
        indice = indice-npoints;
    end
    
    if alph ~= C
        if indice <= npoints
            moy_b(k) = -epsilon + Y(indice) - prodwr(X,alpha,X(:,indice));
        else
            moy_b(k) = +epsilon + Y(indice) - prodwr(X,alpha,X(:,indice));
        end
        k=k+1;
    end
end

b=mean(moy_b);   %  on moyenne les b des points supports

%--------------------------------------------------------------

save Modele X Y alpha b  % On stocke le materiel qui permetra de reconstruire
                         % la fonction apprise à des fins d'optim.
                         % Construire cette fonction
%
% Dessin du "modèle" (fonction de regression). 
% On peut aussi l'appeler en dehors de ce programme
Dessin_Modele

% dessin des performances de la convergence 
figure(2)  % dessin des performances de convergence
plot(log(vit_conv),'linewidth',2)
xlabel(' nombre itérations')
ylabel(' log (ecart)')
title('comportement convergence')
grid
fprintf(' nb iterations: %d critère: %e \r',npas,crit_arret)

%% --------------------------------------------------------------------------

load antennes_2d_test  C  % Dans le fichier antenne, C contient X' et 

X=C';
npoints=size(X,2);

S_test = zeros(1,npoints);

for i=1:npoints
    S_test(i) = pred(X(:,i));
end

save Benechehab_Gueddari_Adref_test_2d S_test
