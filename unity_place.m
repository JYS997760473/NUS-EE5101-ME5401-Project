% place use Ackermann's formula unity rank method
% Author Jia Yansong
function K = unity_place(A, B, P)
    [num_state, num_input] = size(B);
    syms s;
    q = create_q(num_input);
    B = B*q;
    W_c = [B A*B A^2*B A^3*B A^4*B A^5*B];
    [row, col] = size(W_c);
    r = min(row, col);
    if rank(W_c) < r
        error('it is uncontrollable!');
    end
    % select num_state independent vectors out of num_state*num_input
    % vectors from the controllability matrix in the strict order form left
    % to right:
    C = W_c;
    % inverse C:
    C_inv = C^(-1);
    co = 1;
    for i = 1: num_state
        co = co * (s - P(1,i));
    end
    co = expand(co);
    coco = sym2poly(co);
    phid = polyvalm(coco, A);
    K = create_one(num_state)*C_inv*phid;
    K = q*K;
end
function res = create_q(num_input)
    res = zeros(num_input, 1);
    res(num_input, 1) = 1;
    res = [0.1;1];
end
function res = create_one(num_state)
    res = zeros(1, num_state);
    res(1, num_state) = 1;
end