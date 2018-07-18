function [mw, rho, sol] = loadOdorSolub()
fileName = fullfile('data', 'odor_density_mw_solubility.csv');

dataT = readtable(fileName);    % load Excel file

mw  = dataT.MW_g_mol_;
rho = dataT.Density_g_mL_;
sol = dataT.SolubilityInWater_g_L_;

end