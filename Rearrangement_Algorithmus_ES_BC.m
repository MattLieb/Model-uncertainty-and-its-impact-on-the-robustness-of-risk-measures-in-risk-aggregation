function [ES_WC, Zeilensumme2] = Rearrangement_Algorithmus_ES_BC(X, tol, N, alpha)
    % Random column permutations and basic setup
    X = arrayfun(@(x) X(randperm(size(X,1)),x), 1:size(X,2), 'UniformOutput', false);
    X = cat(2, X{:});
    d = size(X, 2);
    Y = X;
    Y_rs = sum(Y, 2); % Y row sums (only computed once)
    m_rs_old = max(Y_rs); % initial minimal row sums (to compare against)
    
    % Main
    while true
        % Loop over all columns and oppositely reorder the jth with respect to
        % the sum of all others
        for j = 1:d
            Y_j = Y(:,j); % jth column
            Y_rs_mj = Y_rs - Y_j; % row sums without jth = total row sums - jth column
            
            % Oppositely order the jth column with respect to the sum of all others
            [~, ranks] = sort(Y_rs_mj, 'descend');
            [~, inv_ranks] = sort(ranks, 'descend');
            [~, sort_ind] = sort(Y_j, 'ascend');
            Y(:,j) = Y_j(sort_ind(inv_ranks));
            Y_rs = Y_rs_mj + Y(:,j); % update total row sum
        end
        m_rs_new = max(Y_rs); % update minimal row sum
        
        % Check stopping criterion
        tol_ = abs((m_rs_new - m_rs_old) / m_rs_old); % relative change of minimal row sum
        if isempty(tol)
            tol_reached = isequal(Y, X); % stop only if matrix did not change after d iterations
        else
            tol_reached = tol_ <= tol; % reached tolerance if tol_ <= tol
        end
        
        % If fulfilled, stop, otherwise update and continue
        if tol_reached % if tolerance was reached
            break % break while()
        else % update (and continue)
            m_rs_old = m_rs_new;
            X = Y;
        end
    end
    
    % Return
  
     % Return
    %lower bound
    Zeilensummen = sum(Y, 2);
    Zeilensumme2=sort(Zeilensummen);
    ES_WC=sum(Zeilensumme2(floor((N*alpha))+1:N))/N/(1-alpha);
end
