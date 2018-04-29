 %Tanner: Decodes a bit vector using soft LDPC
    % Exchanges bit probabilities instead of binary values between v-nodes
    % and c-nodes.
H = [0 1 0 1 1 0 0 1;
    1 1 1 0 0 1 0 0;
    0 0 1 0 0 1 1 1;
    1 0 0 1 1 0 1 0];
    transmitted = [1 0 0 1 0 1 0 1];
    sent = [1 1 -1 1 -1 1 -1 1];

    sigma = 0.2;
    rcvd = sent + sigma*randn(1,length(sent));

    [NbrCheckn,NbrVarn] = size(H);
    % Compute variance
    v = [rcvd(rcvd>0)-1, rcvd(rcvd<0)+1];
    sigma2 = var(v);    % estimate of variance of noise
    %sigma2 = sigma^2;
    % P0 = P(ci=0|ri)
    % P1 = P(ci=1|ri)
    P0 = 1./(1+exp(2*rcvd/sigma2));
    P1 = 1 - P0;
    
    q0 = repmat(P0',1,NbrCheckn);
    q1 = repmat(P1',1,NbrCheckn);
    
    Q0 = zeros(1,NbrVarn);
    Q1 = zeros(1,NbrVarn);
    
    decision = zeros(1,NbrVarn);
    % begin iterative process
    nbrIterations = 4;
     % r0(j,i): response of fj to ci, probability that ci is a zero
    r0 = zeros(NbrCheckn,NbrVarn); r1 = zeros(NbrCheckn,NbrVarn);
    for cntr = 1:nbrIterations
        for i = 1:NbrVarn
            for j = 1:NbrCheckn
                % set of row locs of the ones in the jth column of H
                tmp = H(j,:);  % jth row of H
                % R: the set of column locations of the 1's in the jth row of H, excluding column i
                tmp(i) = 0;  % negate ci
                Ri = find(tmp); % gives the indices of nonzero elements in the jth row 
                r0(j,i) = 1/2 + 1/2 * prod(1-2*q1(Ri,j));
                r1(j,i) = 1 - r0(j,i);
            end
        end
        K   = zeros(NbrVarn,NbrCheckn);
        for i = 1:NbrVarn
            for j = 1:NbrCheckn
                % set of row locs of the ones in the jth column of H
                tmp = H(:,i);  % ith column of H
                tmp(j) = 0; % negate j
                Cj = find(tmp); % gives the indices of nonzero elements in the jth row
                q0(i,j) = P0(i)*prod(r0(Cj,i));
                q1(i,j) = P1(i)*prod(r1(Cj,i));
                K(i,j) = 1/(q0(i,j) + q1(i,j));
                q0(i,j) = K(i,j) * q0(i,j);
                q1(i,j) = K(i,j) * q1(i,j);
            end
            K = zeros(1,NbrVarn);
            C = H(:,i); C = find(C);
            Q0(i) = P0(i) * prod(r0(C,i));
            Q1(i) = P1(i) * prod(r1(C,i));
            K(i) = 1/(Q0(i)+Q1(i));
            Q0(i) = K(i) * Q0(i);
            Q1(i) = K(i) * Q1(i);
            if Q1(i)>0.5
                decision(i) = 1;
            else
                decision(i) = 0;
            end
        end
    end