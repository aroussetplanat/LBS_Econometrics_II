function [L, x_t] = kf_lik(gamma,Y,a1,p1,fct)

R = 0;
F = [gamma(2), 1 ; 0, 0];
H = [1,0];
Q = [gamma(4)^2, gamma(3)*gamma(4); gamma(3)*gamma(4), gamma(3)^2*gamma(4)^2];
A = gamma(1);

[x_t, x_t_1, sigma_t, sigma_t_1, L] = fct(Y, F,Q,A,H,R,a1,p1);

end