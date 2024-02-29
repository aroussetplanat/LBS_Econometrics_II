function [x_t, x_t_1, sigma_t, sigma_t_1, L] = KalmanFilter(Y, F,Q,A,H,R,a1,p1)
%3.1
[T,n]= size(Y);
k = size(a1,1);

%3.2 
x_t =zeros(k,T);
x_t_1 =zeros(k,T);
sigma_t =zeros(k,k,T);
sigma_t_1 =zeros(k,k,T);
L = zeros(T,1);
lik = zeros(T,1);

%3.3
tilde_x_t = 0;
tilde_x_t_1 = 0;
tilde_sigma_t = 0;
tilde_sigma_t_1 = 0;

%3.4
tilde_x_t = a1;
tilde_sigma_t = p1;

    for i=1:T
        %i. 
        tilde_x_t_1 = F*tilde_x_t;
        tilde_sigma_t_1 = F*tilde_sigma_t * F' + Q;

        %ii. 
        tilde_y = Y(i,:)'-A-H *tilde_x_t_1;
        tilde_omega = H * tilde_sigma_t_1 * H' + R;

        %iii.
        K = tilde_sigma_t_1 *H' *inv(tilde_omega);
        tilde_x_t = tilde_x_t_1 + K *tilde_y;
        tilde_sigma_t = tilde_sigma_t_1 - K*H*tilde_sigma_t_1;

        %iv.
        lik(i) = -.5*log(det(tilde_omega)) - .5*tilde_y' *inv(tilde_omega)*tilde_y;

        %v. 
        x_t(:,i) =tilde_x_t;
        x_t_1(:,i) =tilde_x_t_1;
        sigma_t(:,:,i) =tilde_sigma_t;
        sigma_t_1(:,:,i) =tilde_sigma_t_1;
    end

%3.5 
L = -(sum(lik,1) - (T/2)*log(2)*pi);
x_t = x_t';
x_t_1 = x_t_1';

end