function [output] = RSVD_linear(Q,Y,Z,U)
    % return RSVD approximation based on Online phase discussed in the
    % paper.

    [Q_tilde,~]=qr(ttv(Y,U',3).data,0);   
    output=Q*Q_tilde*(Q_tilde'*ttv(Z,U',3).data);

end