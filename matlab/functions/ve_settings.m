function ve = ve_settings(name)
%VE_SETTINGS Vaccine effectiveness for the three vaccine-derived immunity states.
%   ve.inf, ve.dis, ve.hsp, ve.mort: [3 x 2] fractions; rows are <1 year,
%   1-5 years and >5 years after vaccination, columns are seronegative and
%   seropositive. Daniels estimates are averaged over serotypes. Effectiveness
%   against death equals that against hospitalization (Table S9).

src = ve_source_data();
yr  = [1 2 10];   % representative year for each immunity state

switch lower(name)
    case 'shen'
        s = src.shen;
        ve.inf = [s.neg_inf(yr);  s.pos_inf(yr)]' / 100;
        ve.dis = [s.neg_symp(yr); s.pos_symp(yr)]' / 100;
        ve.hsp = [s.hosp(yr);     s.hosp(yr)]' / 100;
        ve.label = 'Shen';
    case 'daniels'
        d = src.daniels;
        ve.inf = [mean(d.neg_inf(:, yr), 1);  mean(d.pos_inf(:, yr), 1)]' / 100;
        ve.dis = [mean(d.neg_symp(:, yr), 1); mean(d.pos_symp(:, yr), 1)]' / 100;
        ve.hsp = [mean(d.neg_hosp(:, yr), 1); mean(d.pos_hosp(:, yr), 1)]' / 100;
        ve.label = 'Daniels';
    otherwise
        error('Unknown VE setting "%s" (use shen or daniels)', name);
end
ve.mort = ve.hsp;
ve.name = lower(name);
end
