function [A,B] = RSVD_linear_factors(Q,Y,Z,U)
    % return RSVD approximation based on Online phase discussed in the
    % paper.
        [Q_tilde,~]=qr(ttv(Y,U',3).data,0);
        A=(Q_tilde'*ttv(Z,U',3).data);
        B=Q*Q_tilde;
        
end