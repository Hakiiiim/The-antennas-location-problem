function prodw = prodwr(X,alph,u)
% on calcule ici le produit scalaire de w avec le vecteur u qui doit etre 
% l'image d'un vecteur de Rn par la retranscription
% X est le tableau des points d'apprentissage
prodw=0;
n=length(X);
        for k=1:n
            prodw=prodw + (alph(k)-alph(k+n))*kernel(X(:,k),u);
        end