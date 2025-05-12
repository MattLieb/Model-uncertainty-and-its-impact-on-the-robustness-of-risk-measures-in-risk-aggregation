function FrankCopula(Hauptpfad, Dim, alpha, N, Size)
    % FrankCopula function
    % Inputs:
    % Hauptpfad: Main path for file operations
    % Dim: Array of dimensions
    % alpha: Array of confidence levels
    % N: Number of paths
    % Size: Array of sample sizes

    rng(2307); % Reproducibility
    eps = 0;
    tol = 0;

    % Load parameter file
    Parameter = fullfile(Hauptpfad, 'RV_Parameter.xls');
    RV_Parameter = xlsread(Parameter);

    % Initialize variables
    VaR_BC = zeros(N, 1);
    VaR_WC = zeros(N, 1);
    ES_BC = zeros(N, 1);
    ES_WC = zeros(N, 1);
    VaR_Como = zeros(N, 1);
    ErgebnisMatrix_Frank = zeros(28, 48);
    Kennzahlen_all = [];
    Konvergenz_all = [];

    % Main loop over dimensions
    for u = 1:length(Dim)
        D = Dim(u);

        for r = 1:length(Size)
            % Load raw data for current dimension and sample size
            filename = sprintf('Basis%02d_%02d.xlsx', Dim(u), Size(r));
            RV_C = xlsread(fullfile(Hauptpfad, filename));

            % Fit Frank Copula
            families1 = {'C'};
            [HACObjectC, fitLogC] = HACopulafit(RV_C, families1);

            % Initialize matrices for random variables
            FrankR1 = zeros(Size(r), D/5);
            FrankR2 = zeros(Size(r), D/5);
            FrankR3 = zeros(Size(r), D/5);
            FrankR4 = zeros(Size(r), D/5);
            FrankR5 = zeros(Size(r), D/5);

            % Initialize result matrices
            RM_Frank = zeros(length(alpha), 4);
            VaR_Comontonic = zeros(length(alpha), 2);
            Mean = zeros(length(alpha), 5);

            % Loop over confidence levels
            for l = 1:length(alpha)
                for i = 1:N
                    % Simulate Frank Copula
                    Frank = HACopularnd(HACObjectC, Size(r));
                    Ma_Frank = sort(Frank);

                    % Transform marginals using distribution parameters
                    for j = 1:D/5
                        FrankR1(:, j) = icdf('gp', Ma_Frank(:, D - (D - j)), RV_Parameter(1, j), RV_Parameter(2, j), RV_Parameter(3, j));
                        FrankR2(:, j) = icdf('LogNormal', Ma_Frank(:, D - (D - (j + D/5))), RV_Parameter(4, j), RV_Parameter(5, j));
                        FrankR3(:, j) = icdf('Exponential', Ma_Frank(:, D - (D - (j + 2*D/5))), RV_Parameter(6, j));
                        FrankR4(:, j) = icdf('wbl', Ma_Frank(:, D - (D/5*2 - j)), RV_Parameter(7, j), RV_Parameter(8, j));
                        FrankR5(:, j) = icdf('Gamma', Ma_Frank(:, D - ((D/5) - j)), RV_Parameter(9, j), RV_Parameter(10, j));
                    end

                    Frank_Matrix = [FrankR1 FrankR2 FrankR3 FrankR4 FrankR5];
                    alpha_RM = alpha(l);

                    % Calculate percentiles for VaR and ES
                    idx = ceil(alpha_RM * Size(r));
                    Frank_Matrix_BC = Frank_Matrix(1:idx, :);
                    Frank_Matrix_WC = Frank_Matrix(idx+1:end, :);

                    % Calculate risk measures
                    VaR_Como(i, :) = sum(Frank_Matrix(idx+1, :));
                    VaR_BC(i, :) = Rearrangement_Algorithmus_VaR_BC(Frank_Matrix_BC, eps);
                    VaR_WC(i, :) = Rearrangement_Algorithmus_VaR_WC(Frank_Matrix_WC, eps);
                    ES_BC(i, :) = Rearrangement_Algorithmus_ES_BC(Frank_Matrix, eps, Size(r), alpha_RM);
                    Zeilensumme_Clayton = sum(Frank_Matrix, 2);
                    ES_WC(i, :) = sum(Zeilensumme_Clayton(floor((Size(r) * alpha_RM)) + 1:end)) / Size(r) / (1 - alpha_RM);
                end

                % Aggregate results for current alpha
                RM_Frank(l, :) = [min(VaR_BC) max(VaR_WC) min(ES_BC) max(ES_WC)];
                Mean(l, :) = [mean(VaR_BC) mean(VaR_WC) mean(ES_BC) mean(ES_WC) mean(VaR_Como)];
                VaR_Comontonic(l, :) = [min(VaR_Como) max(VaR_Como)];
            end

            % Calculate additional metrics
            Delta_WC = RM_Frank(:, 2) ./ Mean(:, 5);
            Delta_Quo_max = RM_Frank(:, 4) ./ VaR_Comontonic(:, 2);
            Spread_VaR = RM_Frank(:, 2) - RM_Frank(:, 1);
            Spread_ES = RM_Frank(:, 4) - RM_Frank(:, 3);
            Kennzahlen = [Spread_VaR Spread_ES Delta_WC Delta_Quo_max];

            % Combine data
            Haupttabelle = [RM_Frank Spread_VaR Spread_ES];

            % Save results for current dimension and sample size
            ErgebnisMatrix_Frank(u, r) = Haupttabelle;
            Kennzahlen_all = [Kennzahlen_all; Kennzahlen];
        end

        % Save Copula parameters for the current dimension
        CopulaParameter(u, 1) = HACObjectC.Parameter;
    end

    % Write results to Excel
    xlswrite(fullfile(Hauptpfad, 'RohDaten_Frank.xlsx'), CopulaParameter, 'Theta', 'B1');
    xlswrite(fullfile(Hauptpfad, 'RohDaten_Frank.xlsx'), ErgebnisMatrix_Frank, 'Tabellen', 'B2');
    xlswrite(fullfile(Hauptpfad, 'RohDaten_Frank.xlsx'), Kennzahlen_all, 'Kennzahlen', 'B2');
end