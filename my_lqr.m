function K = my_lqr(A, B, Q, R)
    % solve ARE:
    % find the positive definite solution of the Riccati equation
    syms P;
    eqn_ARE = P*A + transpose(A)*P - P*B*inv(R)*transpose(B)*P + Q == 0;
    tau = [A -B*inv(R)*transpose(B); -Q -transpose(A)];
    % get eigenvalues
    [V, D] = eig(tau);
    [row, col] = size(D);
    % store stable eignvalues' index
    index = 1;
    for i = 1: row
        if D(i,i) < 0
            d(1,index) = i;
            index = index + 1;
        end
    end
    [row_ele, col_ele] = size(d);
    index = 1;
    for i = 1: col_ele
        ele(:, index) = V(:, d(1, i));
        index = index + 1;
    end
    [r, c] = size(ele);
    v = ele(1:0.5*r, :);
    mu = ele(0.5*r+1:r, :);
    P = mu/v;
    K = real(inv(R)*transpose(B)*P);
end