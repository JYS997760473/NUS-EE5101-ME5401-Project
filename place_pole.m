% use cannocial form and pole placement method to place poles 
% guide by nus ME5401 Chapter 7 Multiply input pole placement
% Author---Jia Yansong---2022.11.6-----
function K = place_pole(A, B, P)
    % this function is used to place poles by pole placement 
    % write this for nus ME5401 mini project's six order system. 
    % first three poles in P belong to first third order system,
    % second three poles in P belong to second third order system.
    % use controllable canonical form to calculate poles.
    [num_state, num_input] = size(B);
    syms s;
    %W_c = [B A*B A^2*B A^3*B A^4*B A^5*B];
    %[row, col] = size(W_c);
    index = 1;
    W_c = [];
    for i = 0: num_state - 1
        middle = A^i*B;
        W_c = [W_c middle];
        index = index + 1;
    end
    [row, col] = size(W_c);
    r = min(row, col);
    if rank(W_c) < r
        error('it is uncontrollable!');
    end
    % select num_state independent vectors out of num_state*num_input
    % vectors from the controllability matrix in the strict order form left
    % to right:
    d = zeros(1, num_input); % di imply the number of vectors in C related to the ith input, ui.
    CC = W_c(:, 1);
    j = 1;
    d(1,1) = 1;
    for i = 2: num_state*num_input
        j = j + 1;
        CC(:, j) = W_c(:, i);
        % judge this vector whether independent.
        if rank(CC) < j
            CC(:, j) = [];
            j = j - 1;
            continue;
        end
        % update d matrix in order to help next step to group them
        index = rem(i, num_input);
        if index == 0
            d(1, num_input) = i / num_input;
        else
            d(1, rem(i, num_input)) = ceil(i / num_input); 
        end
        if j == num_state
            break;
        end
    end
    % group them in a square matrix:
    C = CC;
    index = 1;
    for i = 1: num_input
        for j = 0: d(1, i)-1
            C(:, index) = A^j*B(:, i);
            index = index + 1;
            if index > num_state
                break;
            end
        end
    end   
    % inverse C:
    C_inv = C^(-1);
    % form T according the process:
    T = zeros(num_state, num_state);
    index = 1;
    d_index = 0;
    for i = 1: num_input
        d_index = d_index + d(1, i);
        for j = 0: d(1, i)-1
            T(index, :) = C_inv(d_index, :)*A^j;
            index = index + 1;
        end
    end
    % calculate A_bar and B_bar
    A_bar = T*A/T; 
    B_bar = T*B;
    K_bar = sym('k', [num_input num_state]); 
    closed_loop_matrix= A_bar - B_bar*K_bar;
    % desired close loop matrix in canonical form create num_input sub-canonical form in desired close loop matrix:
    % create this matrix from top left to right down:
    A_cl = zeros(num_state, num_state);
    d_index = 0;
    for i = 1: num_input
        pre_index = d_index + 1;
        d_index = d_index + d(1,i);
        sub_poles = P(1, pre_index : d_index);
        sub_co = 1;
        [~, num_po] = size(sub_poles);
        for j = 1: num_po
            sub_co = sub_co * (s - sub_poles(1, j));
        end
        co = expand(sub_co); 
        co1 = sym2poly(co);
        canonical_co = -1*flip(co1);
        caco = canonical_co(1,1:num_state/num_input);
        sub_matrix = create_canonical(caco);
        A_cl((i-1)*num_state/num_input + 1 : i*num_state/num_input, (i-1)*num_state/num_input + 1 : i*num_state/num_input) = sub_matrix;
    end
    A_cl;
    d;
    % solve k
    d_index = 0;
    for i = 1: num_input
        d_index = d_index + d(1,i);
        eqn(i, :) = closed_loop_matrix(d_index, :) == A_cl(d_index, :);
    end
    sol = solve(eqn, K_bar);
    sol = struct2cell(sol);
    index = 1;
    for i = 1: num_state
        for j = 1: num_input
            K_bar(j, i) = vpa(str2num(string(sol(index))));

            index = index + 1;
        end
    end
    K = double(K_bar*T);
end
function matrix = create_canonical(co)
% co is a matrix contains non-trail row of desried canonical form matrix
    [~, num] = size(co);
    if num == 1
        matrix = co;
    else 
        right_top = eye(num - 1);
        top = [zeros(num-1, 1) right_top];
        whole = [top; co];
        matrix = whole;
    end
end