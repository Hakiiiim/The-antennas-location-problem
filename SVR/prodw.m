function prodw = prodw(X,lab,alph,u)
% on calcule ici le produit scalaire de w avec le vecteur u qui doit etre 
% l'image d'un vecteur de Rn par la retranscription
% X est le tableau des points d'apprentissage

prodw=0;
na=length(lab); % nombre de points 
        for k=1:na
            prodw=prodw+ lab(k)*alph(k)*kernel(X(:,k),u);
        end