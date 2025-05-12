function [best_VaR] = Rearrangement_Algorithmus_VaR_BC(X, tol)
    % Zufällige Spaltenpermutationen und Grundaufbau
    % Hier werden die Spalten der Matrix X zufällig permutiert
    % Jede Spalte von X wird unabhängig zufällig durchmischt.
    X = arrayfun(@(x) X(randperm(size(X,1)),x), 1:size(X,2), 'UniformOutput', false);
    X = cat(2, X{:}); % Zusammenfügen der permutierten Spalten zu einer neuen Matrix
    d = size(X, 2);% Anzahl der Spalten (Dimension)
    Z = X;% Initialisierung der Matrix Z, die umgeordnet wird
    Z_rs = sum(Z, 2); % Berechnung der Zeilensummen von Z wird nur einmal am Anfang berechnet
    m_rs_old = max(Z_rs); % Anfangsmaximum der Zeilensummen für Vergleichszwecke
    
     % Hauptschleife
    while true
        % Schleife über alle Spalten und Neuanordnung der j-ten Spalte in
        % Bezug auf die Summe aller anderen Spalten
        for j = 1:d
            Z_j = Z(:,j); % Extrahiere die j-te Spalte
            Z_rs_mj = Z_rs - Z_j; % Berechne die Zeilensummen ohne die j-te Spalte (Z_rs - Z_j)
            
             % Sortiere die j-te Spalte in Bezug auf die Summe der anderen Spalten
            [~, ranks] = sort(Z_rs_mj, 'descend');% Sortiere Zeilensummen ohne j-te Spalte in absteigender Reihenfolge
            [~, inv_ranks] = sort(ranks, 'descend');% Bestimme die inverse Rangordnung der sortierten Werte
            [~, sort_ind] = sort(Z_j, 'ascend'); % Sortiere die j-te Spalte aufsteigend
            Z(:,j) = Z_j(sort_ind(inv_ranks));% Neuanordnung der j-ten Spalte basierend auf der inversen Rangfolge
            Z_rs = Z_rs_mj + Z(:,j); % Aktualisiere die Zeilensummen, nachdem die j-te Spalte neu geordnet wurde
        end
        m_rs_new = max(Z_rs); % Aktualisiere das Maximum der Zeilensummen
        
        % Überprüfung des Abbruchkriteriums
        tol_ = abs((m_rs_new - m_rs_old) / m_rs_old); % Berechne relative Änderung der maximalen Zeilensummen
        if isempty(tol)
            tol_reached = isequal(Z, X); % Beende, wenn sich die Matrix nach d Iterationen nicht geändert hat
        else
            tol_reached = tol_ <= tol; % Beende, wenn die relative Änderung kleiner oder gleich der Toleranz ist
        end
        
         % Wenn das Abbruchkriterium erfüllt ist, beende die Schleife, andernfalls
        % aktualisiere die Werte und fahre fort
        if tol_reached % Wenn die Toleranz erreicht wurde
            break % Beende die while-Schleife
        else  % Aktualisiere die Werte und fahre fort
            m_rs_old = m_rs_new; % Aktualisiere das alte Maximum der Zeilensummen
            X = Z; % Setze die permutierte Matrix als neue Vergleichsmatrix
        end
    end
    
    % Rückgabewert
    best_VaR = max(Z_rs); % Gib das maximale Value-at-Risk (VaR) zurück, d.h. die größte Zeilensumme
      % Optional: Die letzte aktualisierte Matrix mit den Zeilensummen kann auch verwendet werden

    output.X_rearranged = Z; % Optional: Die umgeordnete Matrix kann ausgegeben werden
end
