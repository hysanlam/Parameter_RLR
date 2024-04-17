function [output] = RSVD(A,S)

%% Return RSVD approximation. By Randomized methods for matrix computations p.9 -Per-Gunnar Martinsson 

[Q,~] = qr(S,0);
output=Q*(Q'*A);

end