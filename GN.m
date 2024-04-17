function G=GN(S,U,Theta)
% output W, V such that X=W*V where X is the Generlized Nystorm approximation, tol=10*eps(1)
    [Q,R]=qr(Theta'*S,0);
     W=S*pinv(R,10*eps(1));
     V=Q'*U; 
     G=W*V;
end