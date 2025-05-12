function TCopula(Hauptpfad, Dim, alpha, N, Size)
    % TCopula function
    % Inputs:
    % Hauptpfad: Main path for file operations
    % Dim: Array of dimensions
    % alpha: Array of confidence levels
    % N: Number of paths
    % Size: Array of sample sizes

    rng(2307); % Reproducibility
    eps = 0;

    % Load parameter file
    Parameter = fullfile(Hauptpfad, 'RV_Parameter.xls');
    RV_Parameter = xlsread(Parameter);

    % Initialize variables
    VaR_BC = zeros(N, 1);
    VaR_WC = zeros(N, 1);
    ES_BC = zeros(N, 1);
    ES_WC = zeros(N, 1);
    VaR_Como = zeros(N, 1);
    ErgebnisMatrix_t = zeros(28, 48);
    Kennzahlen_all = [];
    Konvergenz_all = [];
    Nuhat_all2 = [];

    % Main loop over dimensions
    for u = 1:length(Dim)
        D = Dim(u);
        Rhohat_all = [];
        Nuhat_all = [];

        for r = 1:length(Size)
            % Load raw data for current dimension and sample size
            filename = sprintf('Basis%02d_%02d.xlsx', Dim(u), Size(r));
            RV_t = xlsread(fullfile(Hauptpfad, filename));

            % Fit t-Copula
            [Rhohat, Nuhat] = copulafit('t', RV_t);
            Nuhat_all(r, 1) = Nuhat;
            Rhohat_all = [Rhohat_all, Rhohat, NaN(length(Rhohat(1, :)), 1)];

            % Initialize matrices for random variables
            tR1 = zeros(Size(r), D/5);
            tR2 = zeros(Size(r), D/5);
            tR3 = zeros(Size(r), D/5);
            tR4 = zeros(Size(r), D/5);
            tR5 = zeros(Size(r), D/5);

            % Initialize result matrices
            RM_t = zeros(length(alpha), 4);
            VaR_Comontonic = zeros(length(alpha), 2);
            Mean = zeros(length(alpha), 5);

            % Loop over confidence levels
            for l = 1:length(alpha)
                for i = 1:N
                    % Simulate t-Copula
                    tcopula = copularnd('t', Rhohat, Nuhat, Size(r));
                    Ma_t = sort(tcopula);

                    % Transform marginals using distribution parameters
                    for j = 1:D/5
                        tR1(:, j) = icdf('gp', Ma_t(:, D - (D - j)), RV_Parameter(1, j), RV_Parameter(2, j), RV_Parameter(3, j));
                        tR2(:, j) = icdf('LogNormal', Ma_t(:, D - (D - (j + D/5))), RV_Parameter(4, j), RV_Parameter(5, j));
                        tR3(:, j) = icdf('Exponential', Ma_t(:, D - (D - (j + 2*D/5))), RV_Parameter(6, j));
                        tR4(:, j) = icdf('wbl', Ma_t(:, D - (D/5*2 - j)), RV_Parameter(7, j), RV_Parameter(8, j));
                        tR5(:, j) = icdf('Gamma', Ma_t(:, D - ((D/5) - j)), RV_Parameter(9, j), RV_Parameter(10, j));
                    end

                    t_Matrix = [tR1 tR2 tR3 tR4 tR5];
                    alpha_RM = alpha(l);

                    % Calculate percentiles for VaR and ES
                    idx = ceil(alpha_RM * Size(r));
                    Matrix_BC = t_Matrix(1:idx, :);
                    Matrix_WC = t_Matrix(idx+1:end, :);

                    % Calculate risk measures
                    VaR_Como(i, :) = sum(t_Matrix(idx+1, :));
                    VaR_BC(i, :) = Rearrangement_Algorithmus_VaR_BC(Matrix_BC, eps);
                    VaR_WC(i, :) = Rearrangement_Algorithmus_VaR_WC(Matrix_WC, eps);
                    ES_BC(i, :) = Rearrangement_Algorithmus_ES_BC(t_Matrix, eps, Size(r), alpha_RM);
                    Zeilensumme = sum(t_Matrix, 2);
                    ES_WC(i, :) = sum(Zeilensumme(floor((Size(r) * alpha_RM)) + 1:end)) / Size(r) / (1 - alpha_RM);
                end

                % Aggregate results for current alpha
                RM_t(l, :) = [min(VaR_BC) max(VaR_WC) min(ES_BC) max(ES_WC)];
                Mean(l, :) = [mean(VaR_BC) mean(VaR_WC) mean(ES_BC) mean(ES_WC) mean(VaR_Como)];
                VaR_Comontonic(l, :) = [min(VaR_Como) max(VaR_Como)];
            end

            % Calculate additional metrics
            Delta_WC = RM_t(:, 2) ./ Mean(:, 5);
            Spread_VaR = RM_t(:, 2) - RM_t(:, 1);
            Spread_ES = RM_t(:, 4) - RM_t(:, 3);
            Kennzahlen = [Spread_VaR Spread_ES Delta_WC];

            % Combine data
            DATEN_C = [RM_t Spread_VaR Spread_ES];

            % Save results for current dimension and sample size
            ErgebnisMatrix_t(u, r) = DATEN_C;
            Kennzahlen_all = [Kennzahlen_all; Kennzahlen];
        end

        % Save t-Copula parameters for the current dimension
        Nuhat_all2 = [Nuhat_all2; Nuhat_all];
    end

    % Write results to Excel
    xlswrite(fullfile(Hauptpfad, 'Rohdaten_t.xlsx'), Rhohat_all, 'Theta', 'B1');
    xlswrite(fullfile(Hauptpfad, 'Rohdaten_t.xlsx'), ErgebnisMatrix_t, 'Tabellen', 'B2');
    xlswrite(fullfile(Hauptpfad, 'Rohdaten_t.xlsx'), Kennzahlen_all, 'Kennzahlen', 'B2');
    xlswrite(fullfile(Hauptpfad, 'Rohdaten_t.xlsx'), Nuhat_all2, 'Nu', 'A1');
end