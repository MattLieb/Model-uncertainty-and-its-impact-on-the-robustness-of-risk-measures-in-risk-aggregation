function GumbelCopula(Hauptpfad, Dim, alpha, N, Size)
    % GumbelCopula function
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
    ErgebnisMatrix_Gumbel = zeros(28, 48);
    Kennzahlen_all = [];
    Konvergenz_all = [];

    % Main loop over dimensions
    for u = 1:length(Dim)
        D = Dim(u);

        for r = 1:length(Size)
            % Load raw data for current dimension and sample size
            filename = sprintf('Basis%02d_%02d.xlsx', Dim(u), Size(r));
            RV_C = xlsread(fullfile(Hauptpfad, filename));

            % Fit Gumbel Copula
            families1 = {'G'};
            [HACObjectC, fitLogC] = HACopulafit(RV_C, families1);

            % Initialize matrices for random variables
            GumbelR1 = zeros(Size(r), D/5);
            GumbelR2 = zeros(Size(r), D/5);
            GumbelR3 = zeros(Size(r), D/5);
            GumbelR4 = zeros(Size(r), D/5);
            GumbelR5 = zeros(Size(r), D/5);

            % Initialize result matrices
            RM_Gumbel = zeros(length(alpha), 4);
            VaR_Comontonic = zeros(length(alpha), 2);
            Mean = zeros(length(alpha), 5);

            % Loop over confidence levels
            for l = 1:length(alpha)
                for i = 1:N
                    % Simulate Gumbel Copula
                    GumbelCopula = HACopularnd(HACObjectC, Size(r));
                    Ma_Gumbel = sort(GumbelCopula);

                    % Transform marginals using distribution parameters
                    for j = 1:D/5
                        GumbelR1(:, j) = icdf('gp', Ma_Gumbel(:, D - (D - j)), RV_Parameter(1, j), RV_Parameter(2, j), RV_Parameter(3, j));
                        GumbelR2(:, j) = icdf('LogNormal', Ma_Gumbel(:, D - (D - (j + D/5))), RV_Parameter(4, j), RV_Parameter(5, j));
                        GumbelR3(:, j) = icdf('Exponential', Ma_Gumbel(:, D - (D - (j + 2*D/5))), RV_Parameter(6, j));
                        GumbelR4(:, j) = icdf('wbl', Ma_Gumbel(:, D - (D/5*2 - j)), RV_Parameter(7, j), RV_Parameter(8, j));
                        GumbelR5(:, j) = icdf('Gamma', Ma_Gumbel(:, D - ((D/5) - j)), RV_Parameter(9, j), RV_Parameter(10, j));
                    end

                    Matrix = [GumbelR1 GumbelR2 GumbelR3 GumbelR4 GumbelR5];
                    alpha_RM = alpha(l);

                    % Calculate percentiles for VaR and ES
                    idx = ceil(alpha_RM * Size(r));
                    Matrix_BC = Matrix(1:idx, :);
                    Matrix_WC = Matrix(idx+1:end, :);

                    % Calculate risk measures
                    VaR_Como(i, :) = sum(Matrix(idx+1, :));
                    VaR_BC(i, :) = Rearrangement_Algorithmus_VaR_BC(Matrix_BC, eps);
                    VaR_WC(i, :) = Rearrangement_Algorithmus_VaR_WC(Matrix_WC, eps);
                    ES_BC(i, :) = Rearrangement_Algorithmus_ES_BC(Matrix, eps, Size(r), alpha_RM);
                    Zeilensumme = sum(Matrix, 2);
                    ES_WC(i, :) = sum(Zeilensumme(floor((Size(r) * alpha_RM)) + 1:end)) / Size(r) / (1 - alpha_RM);
                end

                % Aggregate results for current alpha
                RM_Gumbel(l, :) = [min(VaR_BC) max(VaR_WC) min(ES_BC) max(ES_WC)];
                Mean(l, :) = [mean(VaR_BC) mean(VaR_WC) mean(ES_BC) mean(ES_WC) mean(VaR_Como)];
                VaR_Comontonic(l, :) = [min(VaR_Como) max(VaR_Como)];
            end

            % Calculate additional metrics
            Delta_WC = RM_Gumbel(:, 2) ./ Mean(:, 5);
            Spread_VaR = RM_Gumbel(:, 2) - RM_Gumbel(:, 1);
            Spread_ES = RM_Gumbel(:, 4) - RM_Gumbel(:, 3);
            Kennzahlen = [Spread_VaR Spread_ES Delta_WC];

            % Combine data
            DATEN_C = [RM_Gumbel Spread_VaR Spread_ES];

            % Save results for current dimension and sample size
            ErgebnisMatrix_Gumbel(u, r) = DATEN_C;
            Kennzahlen_all = [Kennzahlen_all; Kennzahlen];
        end

        % Save Gumbel Copula parameters for the current dimension
        CopulaParameter(u, 1) = HACObjectC.Parameter;
    end

    % Write results to Excel
    xlswrite(fullfile(Hauptpfad, 'RohdatenGumbel.xlsx'), CopulaParameter, 'Theta', 'B1');
    xlswrite(fullfile(Hauptpfad, 'RohdatenGumbel.xlsx'), ErgebnisMatrix_Gumbel, 'Tabellen', 'B2');
    xlswrite(fullfile(Hauptpfad, 'RohdatenGumbel.xlsx'), Kennzahlen_all, 'Kennzahlen', 'B2');
end