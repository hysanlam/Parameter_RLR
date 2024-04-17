function [W,V]=GN_factors(S,U,Theta)
% output W,V such that G=W*V is the Generlized nystorm apporixmation.
    [Q,R]=qr(Theta'*S,0);
     W=S*pinv(R,10*eps(1));
     V=Q'*U; 
    
end