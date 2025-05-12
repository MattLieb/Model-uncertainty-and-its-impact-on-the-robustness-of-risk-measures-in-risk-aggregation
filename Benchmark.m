function Benchmark(HP, Dim, alpha, N, Size)
    % Benchmark function
    % Inputs:
    % HP (Hauptpfad): Main path for file operations
    % Dim: Array of dimensions
    % alpha: Array of confidence levels
    % N: Number of paths
    % Sizes: Sample sizes

    rng(2307); % Reproducibility

    % Initialize results matrices
    ErgebnisMatrix = zeros(28, 48);
    KENN_ZAHLEN = [];
    Konvergenz_all = [];

    % Generate distribution parameters
    Pareto_k = 0.1 + (0.5-0.1).*rand(1, length(Dim));
    Pareto_sigma = 2 + (4-2).*rand(1, length(Dim));
    Pareto_mu = 0.00001 + (0.8-0.00001).*rand(1, length(Dim));
    LogNor_mu = (2.2).*rand(1, length(Dim));
    LogNor_sigma = 1.5.*rand(1, length(Dim));
    Expo_mu = 15 + (20-15).*rand(1, length(Dim));
    Weibull_alpha = 20 + (20-15).*rand(1, length(Dim));
    Weibull_beta = 0.7 + (2-0.7).*rand(1, length(Dim));
    Gamma_alpha = 0.8 + (2-0.8).*rand(1, length(Dim));
    Gamma_beta = 15 + (20-15).*rand(1, length(Dim));

    % Save distribution parameters to Excel
    rowname1 = {'Dim05', 'Dim10', 'Dim15', 'Dim20'};
    column_Names = {'Pareto_k', 'Pareto_sigma', 'Pareto_mu', 'LogNor_mu', 'LogNor_sigma', ...
                    'Expo_mu', 'Weibull_alpha', 'Weibull_beta', 'Gamma_alpha', 'Gamma_beta'};
    xlswrite([HP 'RV_Parameter.xlsx'], rowname1', 'B1');
    xlswrite([HP 'RV_Parameter.xlsx'], column_Names', 'A2:A11');
    xlswrite([HP 'RV_Parameter.xlsx'], [Pareto_k; Pareto_sigma; Pareto_mu; LogNor_mu; ...
        LogNor_sigma; Expo_mu; Weibull_alpha; Weibull_beta; Gamma_alpha; Gamma_beta], 'B2:E11');

    % Copula initialization
    C = cell(1, length(Dim));
    for u = 1:length(Dim)
        elements = 1:Dim(u);
        C{u} = {{'C', 5}, elements};
    end

    % Loop over dimensions and sample sizes
    for u = 1:length(Dim)
        for r = 1:length(Size)
            sz = [1 Size(r)];
            D = Dim(u);
            myHAC = HACopula();
            cell2tree(myHAC, C{u});
            BenchMark = HACopularnd(myHAC, Size(r));

            % Save benchmark data to Excel
            xlswrite([HP sprintf('Basis_%02d_%02d.xlsx', Dim(u), Size(r))], BenchMark);

            % Transform data using distribution parameters
            BM = zeros(Size(r), D);
            for j = 1:D/5
                BM(:, j) = icdf('gp', BenchMark(:, D - (D - j)), Pareto_k(j), Pareto_sigma(j), Pareto_mu(j));
                BM(:, j+D/5) = icdf('LogNormal', BenchMark(:, D - (D - (j + D/5))), LogNor_mu(j), LogNor_sigma(j));
                BM(:, j+2*D/5) = icdf('Exponential', BenchMark(:, D - (D - (j + 2*D/5))), Expo_mu(j));
                BM(:, j+3*D/5) = icdf('wbl', BenchMark(:, D - (D/5*2-j)), Weibull_alpha(j), Weibull_beta(j));
                BM(:, j+4*D/5) = icdf('Gamma', BenchMark(:, D-((D/5)-j)), Gamma_alpha(j), Gamma_beta(j));
            end
            BM = sort(BM);

            % Confidence level loop
            RM_BM = zeros(length(alpha), 4);
            MEAN_BM = zeros(length(alpha), 5);
            VaR_Comontonic_BM = zeros(length(alpha), 2);
            for l = 1:length(alpha)
                VaR_BC = zeros(N, 1);
                VaR_WC = zeros(N, 1);
                ES_BC = zeros(N, 1);
                ES_WC = zeros(N, 1);
                VaR_Como_BM = zeros(N, 1);

                % Path simulation loop
                for i = 1:N
                    BM_copula = HACopularnd(myHAC, Size(r));
                    for j = 1:D/5
                        BM(:, j) = icdf('gp', BM_copula(:, D - (D - j)), Pareto_k(j), Pareto_sigma(j), Pareto_mu(j));
                        BM(:, j+D/5) = icdf('LogNormal', BM_copula(:, D - (D - (j + D/5))), LogNor_mu(j), LogNor_sigma(j));
                        BM(:, j+2*D/5) = icdf('Exponential', BM_copula(:, D - (D - (j + 2*D/5))), Expo_mu(j));
                        BM(:, j+3*D/5) = icdf('wbl', BM_copula(:, D - (D/5*2-j)), Weibull_alpha(j), Weibull_beta(j));
                        BM(:, j+4*D/5) = icdf('Gamma', BM_copula(:, D-((D/5)-j)), Gamma_alpha(j), Gamma_beta(j));
                    end

                    BM = sort(BM);
                    alpha_RM = alpha(l);
                    idx = ceil(alpha_RM * Size(r));

                    VaR_Como_BM(i, :) = sum(BM(idx+1, :));
                    VaR_BC(i, :) = Rearrangement_Algorithmus_VaR_BC(BM(1:idx, :), 0);
                    VaR_WC(i, :) = Rearrangement_Algorithmus_VaR_WC(BM(idx+1:end, :), 0);
                    ES_BC(i, :) = Rearrangement_Algorithmus_ES_BC(BM, 0, Size(r), alpha_RM);
                    ES_WC(i, :) = sum(sum(BM(idx+1:end, :))) / Size(r) / (1 - alpha_RM);
                end

                RM_BM(l, :) = [min(VaR_BC), max(VaR_WC), min(ES_BC), max(ES_WC)];
                MEAN_BM(l, :) = [mean(VaR_BC), mean(VaR_WC), mean(ES_BC), mean(ES_WC), mean(VaR_Como_BM)];
                VaR_Comontonic_BM(l, :) = [min(VaR_Como_BM), max(VaR_Como_BM)];
            end

            % Save results in main matrices
            ErgebnisMatrix(u, r) = RM_BM;

            % Save Excel files
            xlswrite([HP 'Rohdaten_BM.xlsx'], RM_BM, 'Tabellen', 'B2');
        end
    end
end