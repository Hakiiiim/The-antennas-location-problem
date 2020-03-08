function S = pred(C)
    load Modele X Y alpha b
    S = prodwr(X,alpha,C')+b;
end

