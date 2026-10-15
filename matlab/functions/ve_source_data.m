function src = ve_source_data()
%VE_SOURCE_DATA Published TAK-003 vaccine effectiveness (%) by year since vaccination.
%   Daniels et al. (Nat Med 2025): serotype-specific estimates, rows DENV1-4,
%   columns years 1..11 after vaccination.
%   Shen et al. (PLoS Med 2025): serotype-averaged estimates against
%   symptomatic disease, years 1..11; effectiveness against infection is
%   half of that against disease and against hospitalization is constant.

src.daniels.neg_symp = [ 84  55  41  29   16    2  -10  -28  -37 -45 -51;
                        100 100 100 100  100  100  100  100   97  96  93;
                         76  33  18   2  -10  -22  -34  -45  -49 -59 -65;
                         69  24  10  -2  -16  -29  -40  -47  -55 -59 -67];
src.daniels.pos_symp = [ 95  74  60  50   42   35   29   23   18  14  12;
                        100  91  81  75   68   60   53   45   38  32  25;
                         82  48  32  25   19   15   12   10    7   6   5;
                         94  76  65  57   50   43   36   30   24  20  15];
src.daniels.neg_hosp = [100  87  83  70   62   49   40   23    6 -10 -23;
                        100 100 100 100  100  100  100  100  100 100 100;
                         92  62  49  32   15   -6  -27  -45  -70 -91 -100;
                        100  70  61  49   36   27    2  -10  -19 -32 -40];
src.daniels.pos_hosp = [ 97  91  83  78   72   64   56   49   40  33  29;
                        100  95  94  90   87   84   79   74   66  58  49;
                         95  75  60  51   43   35   30   24   19  16  13;
                         98  91  84  80   75   68   62   56   49  43  36];
src.daniels.neg_inf  = [  0   0   0   0    0    0    0    0    0   0   0;
                         66  40  32  26   20   17   13   10    8   6   5;
                          0   0   0   0    0    0    0    0    0   0   0;
                          0   0   0   0    0    0    0    0    0   0   0];
src.daniels.pos_inf  = [ 48  15   8   7    5    3    2    2    2   1   1;
                         71  32  21  16   13   11    7    7    5   5   2;
                         20   5   3   1    1    0    0    0    0   0   0;
                         53  20  12  10    7    5    5    5    2   1   1];

src.shen.neg_symp = [56.2 40.5 34.8 31.4 29.2 27.5 26.1 25.0 24.1 23.3 22.6];
src.shen.pos_symp = [72.0 53.1 46.2 42.1 39.2 37.1 35.5 34.1 32.9 31.9 31.1];
src.shen.neg_inf  = src.shen.neg_symp / 2;
src.shen.pos_inf  = src.shen.pos_symp / 2;
src.shen.hosp     = 84.4 * ones(1, 11);
end
