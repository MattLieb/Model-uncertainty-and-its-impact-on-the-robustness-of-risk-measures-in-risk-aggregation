function GaussianCopula(Hauptpfad, Dim, alpha, N, Size)
    % GaussianCopula function
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
    ErgebnisMatrix_Gauss = zeros(28, 48);
    Kennzahlen_all = [];
    Konvergenz_all = [];

    % Main loop over dimensions
    for u = 1:length(Dim)
        D = Dim(u);
        Rhohat_all = [];

        for r = 1:length(Size)
            % Load raw data for current dimension and sample size
            filename = sprintf('Basis%02d_%02d.xlsx', Dim(u), Size(r));
            RV_Gauss = xlsread(fullfile(Hauptpfad, filename));

            % Fit Gaussian Copula
            RHOHATg = copulafit('Gaussian', RV_Gauss);
            Rhohat_all = [Rhohat_all, RHOHATg, NaN(length(RHOHATg(1, :)), 1)];

            % Initialize matrices for random variables
            GaussR1 = zeros(Size(r), D/5);
            GaussR2 = zeros(Size(r), D/5);
            GaussR3 = zeros(Size(r), D/5);
            GaussR4 = zeros(Size(r), D/5);
            GaussR5 = zeros(Size(r), D/5);

            % Initialize result matrices
            RM_Gauss = zeros(length(alpha), 4);
            VaR_Comontonic = zeros(length(alpha), 2);
            Mean = zeros(length(alpha), 5);

            % Loop over confidence levels
            for l = 1:length(alpha)
                for i = 1:N
                    % Simulate Gaussian Copula
                    Gauss_copula = copularnd('Gaussian', RHOHATg, Size(r));
                    Ma_Gauss = sort(Gauss_copula);

                    % Transform marginals using distribution parameters
                    for j = 1:D/5
                        GaussR1(:, j) = icdf('gp', Ma_Gauss(:, D - (D - j)), RV_Parameter(1, j), RV_Parameter(2, j), RV_Parameter(3, j));
                        GaussR2(:, j) = icdf('LogNormal', Ma_Gauss(:, D - (D - (j + D/5))), RV_Parameter(4, j), RV_Parameter(5, j));
                        GaussR3(:, j) = icdf('Exponential', Ma_Gauss(:, D - (D - (j + 2*D/5))), RV_Parameter(6, j));
                        GaussR4(:, j) = icdf('wbl', Ma_Gauss(:, D - (D/5*2 - j)), RV_Parameter(7, j), RV_Parameter(8, j));
                        GaussR5(:, j) = icdf('Gamma', Ma_Gauss(:, D - ((D/5) - j)), RV_Parameter(9, j), RV_Parameter(10, j));
                    end

                    Gauss_Matrix = [GaussR1 GaussR2 GaussR3 GaussR4 GaussR5];
                    alpha_RM = alpha(l);

                    % Calculate percentiles for VaR and ES
                    idx = ceil(alpha_RM * Size(r));
                    Gauss_Matrix_BC = Gauss_Matrix(1:idx, :);
                    Gauss_Matrix_WC = Gauss_Matrix(idx+1:end, :);

                    % Calculate risk measures
                    VaR_Como(i, :) = sum(Gauss_Matrix(idx+1, :));
                    VaR_BC(i, :) = Rearrangement_Algorithmus_VaR_BC(Gauss_Matrix_BC, eps);
                    VaR_WC(i, :) = Rearrangement_Algorithmus_VaR_WC(Gauss_Matrix_WC, eps);
                    ES_BC(i, :) = Rearrangement_Algorithmus_ES_BC(Gauss_Matrix, eps, Size(r), alpha_RM);
                    Zeilensumme_Clayton = sum(Gauss_Matrix, 2);
                    ES_WC(i, :) = sum(Zeilensumme_Clayton(floor((Size(r) * alpha_RM)) + 1:end)) / Size(r) / (1 - alpha_RM);
                end

                % Aggregate results for current alpha
                RM_Gauss(l, :) = [min(VaR_BC) max(VaR_WC) min(ES_BC) max(ES_WC)];
                Mean(l, :) = [mean(VaR_BC) mean(VaR_WC) mean(ES_BC) mean(ES_WC) mean(VaR_Como)];
                VaR_Comontonic(l, :) = [min(VaR_Como) max(VaR_Como)];
            end

            % Calculate additional metrics
            Delta_WC = RM_Gauss(:, 2) ./ Mean(:, 5);
            Spread_VaR = RM_Gauss(:, 2) - RM_Gauss(:, 1);
            Spread_ES = RM_Gauss(:, 4) - RM_Gauss(:, 3);
            Kennzahlen = [Spread_VaR Spread_ES Delta_WC];

            % Combine data
            DATEN_C = [RM_Gauss Spread_VaR Spread_ES];

            % Save results for current dimension and sample size
            ErgebnisMatrix_Gauss(u, r) = DATEN_C;
            Kennzahlen_all = [Kennzahlen_all; Kennzahlen];
        end

        % Save Gaussian Copula parameters
        xlswrite(fullfile(Hauptpfad, 'RohdatenGauss'), Rhohat_all, 'Theta', sprintf('B%d', u));
        toc
    end

    % Write results to Excel
    xlswrite(fullfile(Hauptpfad, 'RohdatenGauss'), Kennzahlen_all, 'Kennzahlen', 'B2');
    xlswrite(fullfile(Hauptpfad, 'RohdatenGauss'), ErgebnisMatrix_Gauss, 'Tabellen', 'B2');
end