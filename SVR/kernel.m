function kern=kernel(u,v)
%  fonction renvoyant la valeur du noyau K(u,v)
% u et v sont deux vecteurs (colonnes) de Rn
global ikernel % passé en global pour alleger

if ikernel==1
    kern=u'*v;
elseif ikernel==2
    kern=exp(-0.1*norm(u-v)^2);
end

