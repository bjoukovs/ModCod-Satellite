function [ decision ] = Tanner( H, rcvdBits, nbrIterations )
    %Tanner: Decodes a bit vector using LDPC
    [M,N] = size(H);

    % begin iterative process
    for cntr = 1:nbrIterations

        % Send values of vnodes to cnodes 
        % rcvd{i}: vector with received bits at cnode i
        rcvd = {};
        for i = 1:M   %go over all check nodes
           tmp = [];
           connectedNodes = find(H(i, :));
           for k = connectedNodes
               tmp = [tmp rcvdBits(k)];
           end
           rcvd{i} = tmp;
        end

        % Compute what cnode will send back using xor of rcvd by excluding an
        % element
        % sent{i}: vector with the bits cnode i will send back to the corresponding
        % vnodes
        sent = {};
        for i = 1:M
            tmp = [];
            for j = 1:numel(rcvd{i})
                v = rcvd{i}; v(j) = [];          % exclude jth element from rcvd vector
                tmp(j) = mod(sum(v),2);         % compute xor
            end
            sent{i} = tmp;
        end

        nonZeroIndices = {};
        for j = 1:M
            nonZeroIndices{j} = find(H(j, :));
        end
                

        % make decision (majority voting)
        decision = zeros(1, N);
        for i = 1:N
            % 3) messages from checknodes
            tmp =[];
            for j = 1:M             % go over all check nodes
                nzi = nonZeroIndices{j};
                kMax = numel(nzi); 
                for k = 1:kMax
                    if (nzi(k) == i)                         % if variable node i connected to check node j and it's the K'th
                        tmp = [tmp sent{j}(k)];               % sent bit by cnode j to vnode i
                    end
                end
            end

            % 5) decision: majority voting
            vec = [rcvdBits(i) tmp];
            decision(i) = round(sum(vec) / numel(vec));             % most common element in vector
        end
        rcvdBits = decision;
    end
end

