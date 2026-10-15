function [doses_annual, doses_distributed] = brazil_dose_allocation(cfg)
%BRAZIL_DOSE_ALLOCATION State shares of the 2024/25 dengue vaccine distribution.
%   doses_distributed: [S x 1] doses distributed by the Ministry of Health
%                      (shipments R1-R8 and 2025, column "sum")
%   doses_annual     : [S x 1] the same shares scaled to cfg.annual_doses
%
%   The header row is skipped and columns are read by position
%   (Var1 = IBGE state code, Var12 = sum), as in the original analysis.

tbl = readtable(cfg.file_doses, 'NumHeaderLines', 1);
ibge_in_file = tbl.Var1;
doses_sum    = tbl.Var12;

doses_distributed = zeros(cfg.S, 1);
for si = 1:cfg.S
    row = find(ibge_in_file == cfg.ibge_codes(si));
    assert(numel(row) == 1, 'State %s not found in %s', cfg.state_codes{si}, cfg.file_doses);
    doses_distributed(si) = doses_sum(row);
end
doses_annual = cfg.annual_doses * doses_distributed / sum(doses_distributed);
end
