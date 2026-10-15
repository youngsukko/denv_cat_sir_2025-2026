function [doses_annual, doses_distributed] = brazil_dose_allocation(cfg)
%BRAZIL_DOSE_ALLOCATION State shares of the 2024/25 dengue vaccine distribution.
%   doses_distributed: [S x 1] doses distributed by the Ministry of Health
%                      (shipments R1-R8 and 2025)
%   doses_annual     : [S x 1] the same shares scaled to cfg.annual_doses

tbl = readtable(cfg.file_doses, 'TextType', 'char');
[found, row] = ismember(cfg.state_codes, tbl.UF);
assert(all(found), 'Missing states in %s', cfg.file_doses);
doses_distributed = tbl.sum(row);
doses_annual = cfg.annual_doses * doses_distributed / sum(doses_distributed);
end
