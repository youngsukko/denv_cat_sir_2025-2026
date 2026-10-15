function dydt = func_ode_alt_v6(t,y,inputs)
% func_ode_alt_v5과 완전히 동일하나, 백신접종 파트만 다름:
%   v5: inputs.vacc(스칼라) + inputs.vacorder(1~9 정수 스위치, 대상연령 선택) — 그룹 1개만 가능
%   v6: inputs.vacc_targets — struct array, 각 원소 {idx: 대상연령 인덱스들, dose: 그 그룹 총도즈}
%       그룹별로 v5와 동일한 "그 시점 감수성 pool 비례 rate" 메커니즘을 독립 적용.
%       여러 그룹(예: children 연령군 + elderly 연령군)을 동시에, 서로 다른 도즈로 넣을 수 있음.

y(y<0) = 0;
agenum     = inputs.agenum;
beta0      = inputs.beta0;
beta1      = inputs.beta1(floor(t)+1);

gamma      = inputs.gamma;
ohm        = inputs.ohm;
i1         = inputs.i1;
i2         = inputs.i2;
i3         = inputs.i3;
mu         = inputs.mu;
mu_final   = inputs.mu_final;
q          = inputs.q;
birthrate  = inputs.birthrate;

vacc_targets = inputs.vacc_targets;   % struct array: .idx(연령인덱스들), .dose(스칼라 총도즈)

e_p_i  = inputs.e_p_i;
e_p_d  = inputs.e_p_d;
e_p_h  = inputs.e_p_h;
e_p_m  = inputs.e_p_m;
e_n_i  = inputs.e_n_i;
e_n_d  = inputs.e_n_d;
e_n_h  = inputs.e_n_h;
e_n_m  = inputs.e_n_m;

ew_p_i = inputs.ew_p_i;
ew_p_d = inputs.ew_p_d;
ew_p_h = inputs.ew_p_h;
ew_p_m = inputs.ew_p_m;
ew_n_i = inputs.ew_n_i;
ew_n_d = inputs.ew_n_d;
ew_n_h = inputs.ew_n_h;
ew_n_m = inputs.ew_n_m;

eww_p_i = inputs.eww_p_i;
eww_p_d = inputs.eww_p_d;
eww_p_h = inputs.eww_p_h;
eww_p_m = inputs.eww_p_m;
eww_n_i = inputs.eww_n_i;
eww_n_d = inputs.eww_n_d;
eww_n_h = inputs.eww_n_h;
eww_n_m = inputs.eww_n_m;

ohm_vac    = inputs.ohm_vac;
ohm_vac2   = inputs.ohm_vac2;
ohm_vac3   = inputs.ohm_vac3;

h1         = inputs.h1;
h2         = inputs.h2;
h3         = inputs.h3;
h4         = inputs.h4;
f1         = inputs.f1;
f2         = inputs.f2;
f3         = inputs.f3;
f4         = inputs.f4;
cprop_1    = inputs.cprop_1;
cprop_2    = inputs.cprop_2;
cprop_3    = inputs.cprop_3;
cprop_4    = inputs.cprop_4;

cageadj = inputs.cageadj;

S0  = y(agenum* 0 + (1:agenum));
I1  = y(agenum* 1 + (1:agenum));
R1  = y(agenum* 2 + (1:agenum));
S1  = y(agenum* 3 + (1:agenum));
I2  = y(agenum* 4 + (1:agenum));
R2  = y(agenum* 5 + (1:agenum));
S2  = y(agenum* 6 + (1:agenum));
I3  = y(agenum* 7 + (1:agenum));
R3  = y(agenum* 8 + (1:agenum));
S3  = y(agenum* 9 + (1:agenum));
I4  = y(agenum*10 + (1:agenum));
R4  = y(agenum*11 + (1:agenum));

Sv0 = y(agenum*12 + (1:agenum));
Iv1 = y(agenum*13 + (1:agenum));
Rv1 = y(agenum*14 + (1:agenum));
Sv1 = y(agenum*15 + (1:agenum));
Iv2 = y(agenum*16 + (1:agenum));
Rv2 = y(agenum*17 + (1:agenum));
Sv2 = y(agenum*18 + (1:agenum));
Iv3 = y(agenum*19 + (1:agenum));
Rv3 = y(agenum*20 + (1:agenum));
Sv3 = y(agenum*21 + (1:agenum));
Iv4 = y(agenum*22 + (1:agenum));
Rv4 = y(agenum*23 + (1:agenum));

Sw0 = y(agenum*24 + (1:agenum));
Iw1 = y(agenum*25 + (1:agenum));
Rw1 = y(agenum*26 + (1:agenum));
Sw1 = y(agenum*27 + (1:agenum));
Iw2 = y(agenum*28 + (1:agenum));
Rw2 = y(agenum*29 + (1:agenum));
Sw2 = y(agenum*30 + (1:agenum));
Iw3 = y(agenum*31 + (1:agenum));
Rw3 = y(agenum*32 + (1:agenum));
Sw3 = y(agenum*33 + (1:agenum));
Iw4 = y(agenum*34 + (1:agenum));
Rw4 = y(agenum*35 + (1:agenum));

Sww0 = y(agenum*36 + (1:agenum));
Iww1 = y(agenum*37 + (1:agenum));
Rww1 = y(agenum*38 + (1:agenum));
Sww1 = y(agenum*39 + (1:agenum));
Iww2 = y(agenum*40 + (1:agenum));
Rww2 = y(agenum*41 + (1:agenum));
Sww2 = y(agenum*42 + (1:agenum));
Iww3 = y(agenum*43 + (1:agenum));
Rww3 = y(agenum*44 + (1:agenum));
Sww3 = y(agenum*45 + (1:agenum));
Iww4 = y(agenum*46 + (1:agenum));
Rww4 = y(agenum*47 + (1:agenum));

N = sum(S0 + I1 + R1 + S1 + I2 + R2 + S2 + I3 + R3 + S3 + I4 + R4 + ...
    Sv0 + Iv1 + Rv1 + Sv1 + Iv2 + Rv2 + Sv2 + Iv3 + Rv3 + Sv3 + Iv4 + Rv4 + ...
    Sw0 + Iw1 + Rw1 + Sw1 + Iw2 + Rw2 + Sw2 + Iw3 + Rw3 + Sw3 + Iw4 + Rw4 + ...
    Sww0 + Iww1 + Rww1 + Sww1 + Iww2 + Rww2 + Sww2 + Iww3 + Rww3 + Sww3 + Iww4 + Rww4);

FOI = sum(beta0.*(I1 + I2 + I3 + I4 + Iv1 + Iv2 + Iv3 + Iv4 + Iw1 + Iw2 + Iw3 + Iw4 + Iww1 + Iww2 + Iww3 + Iww4))/N;
dS0 = -beta1 * beta0 .* FOI .* S0;
dI1 = beta1 * beta0 .* FOI .* S0 - gamma * I1;
dR1 = gamma * I1 - ohm * R1;
dS1 = ohm * R1 - i1 * beta1 * beta0 .* FOI .* S1;
dI2 = i1 * beta1 * beta0 .* FOI .* S1 - gamma * I2;
dR2 = gamma * I2 - ohm * R2;
dS2 = ohm * R2 - i2 * beta1 * beta0 .* FOI .* S2;
dI3 = i2 * beta1 * beta0 .* FOI .* S2 - gamma * I3;
dR3 = gamma * I3 - ohm * R3;
dS3 = ohm * R3 - i3 * beta1 * beta0 .* FOI .* S3;
dI4 = i3 * beta1 * beta0 .* FOI .* S3 - gamma * I4;
dR4 = gamma * I4;
dSv0 = - (1 - e_n_i) * beta1 * beta0 .* FOI .* Sv0;
dIv1 = (1 - e_n_i) * beta1 * beta0 .* FOI .* Sv0 - gamma * Iv1;
dRv1 = gamma * Iv1 - ohm * Rv1;
dSv1 = ohm * Rww1 + ohm * Rw1 + ohm * Rv1 - i1 * (1 - e_p_i) * beta1 * beta0 .* FOI .* Sv1;
dIv2 = i1 * (1 - e_p_i) * beta1 * beta0 .* FOI .* Sv1 - gamma * Iv2;
dRv2 = gamma * Iv2 - ohm * Rv2;
dSv2 = ohm * Rww2 + ohm * Rw2 + ohm * Rv2 - i2 * (1 - e_p_i) * beta1 * beta0 .* FOI .* Sv2;
dIv3 = i2 * (1 - e_p_i) * beta1 * beta0 .* FOI .* Sv2 - gamma * Iv3;
dRv3 = gamma * Iv3 - ohm * Rv3;
dSv3 =  ohm * Rww3  + ohm * Rw3 + ohm * Rv3 - i3 * (1 - e_p_i) * beta1 * beta0 .* FOI .* Sv3;
dIv4 = i3 * (1 - e_p_i) * beta1 * beta0 .* FOI .* Sv3 - gamma * Iv4;
dRv4 = gamma * Iv4;
dSw0 = - (1 - ew_p_i) * beta1 * beta0 .* FOI .* Sw0;
dIw1 = (1 - ew_p_i) * beta1 * beta0 .* FOI .* Sw0 - gamma * Iw1;
dRw1 = gamma * Iw1 - ohm * Rw1;
dSw1 = - i1 * (1 - ew_p_i) * beta1 * beta0 .* FOI .* Sw1;
dIw2 = i1 * (1 - ew_p_i) * beta1 * beta0 .* FOI .* Sw1 - gamma * Iw2;
dRw2 = gamma * Iw2 - ohm * Rw2;
dSw2 = - i2 * (1 - ew_p_i) * beta1 * beta0 .* FOI .* Sw2;
dIw3 = i2 * (1 - ew_p_i) * beta1 * beta0 .* FOI .* Sw2 - gamma * Iw3;
dRw3 = gamma * Iw3 - ohm * Rw3;
dSw3 =  - i3 * (1 - ew_p_i) * beta1 * beta0 .* FOI .* Sw3;
dIw4 = i3 * (1 - ew_p_i) * beta1 * beta0 .* FOI .* Sw3 - gamma * Iw4;
dRw4 = gamma * Iw4;
dSww0 = - (1 - eww_p_i) * beta1 * beta0 .* FOI .* Sww0;
dIww1 = (1 - eww_p_i) * beta1 * beta0 .* FOI .* Sww0 - gamma * Iww1;
dRww1 = gamma * Iww1 - ohm * Rww1;
dSww1 = - i1 * (1 - eww_p_i) * beta1 * beta0 .* FOI .* Sww1;
dIww2 = i1 * (1 - eww_p_i) * beta1 * beta0 .* FOI .* Sww1 - gamma * Iww2;
dRww2 = gamma * Iww2 - ohm * Rww2;
dSww2 = - i2 * (1 - eww_p_i) * beta1 * beta0 .* FOI .* Sww2;
dIww3 = i2 * (1 - eww_p_i) * beta1 * beta0 .* FOI .* Sww2 - gamma * Iww3;
dRww3 = gamma * Iww3 - ohm * Rww3;
dSww3 = - i3 * (1 - eww_p_i) * beta1 * beta0 .* FOI .* Sww3;
dIww4 = i3 * (1 - eww_p_i) * beta1 * beta0 .* FOI .* Sww3 - gamma * Iww4;
dRww4 = gamma * Iww4;

% ── 백신접종 (v5와 다른 부분): 그룹별로 v5와 동일한 "그 시점 pool 비례
%    rate" 방식을 독립 적용 (여러 그룹을 동시에 다른 도즈로 넣을 수 있음)
for gi = 1:numel(vacc_targets)
    idx  = vacc_targets(gi).idx;
    dose = vacc_targets(gi).dose;

    pool = sum(S0(idx) + S1(idx) + R1(idx) + S2(idx) + R2(idx) + S3(idx) + R3(idx) + R4(idx));
    pool_safe = max(pool, eps);
    grate = min(1, dose/pool_safe);   % v5의 rate = vacc/vactargetpop 와 동일한 정의, 그룹별로

    dS0(idx) = dS0(idx) - grate.*S0(idx);   dSv0(idx) = dSv0(idx) + grate.*S0(idx);
    dS1(idx) = dS1(idx) - grate.*S1(idx);   dSv1(idx) = dSv1(idx) + grate.*S1(idx);
    dS2(idx) = dS2(idx) - grate.*S2(idx);   dSv2(idx) = dSv2(idx) + grate.*S2(idx);
    dS3(idx) = dS3(idx) - grate.*S3(idx);   dSv3(idx) = dSv3(idx) + grate.*S3(idx);
    dR1(idx) = dR1(idx) - grate.*R1(idx);   dRv1(idx) = dRv1(idx) + grate.*R1(idx);
    dR2(idx) = dR2(idx) - grate.*R2(idx);   dRv2(idx) = dRv2(idx) + grate.*R2(idx);
    dR3(idx) = dR3(idx) - grate.*R3(idx);   dRv3(idx) = dRv3(idx) + grate.*R3(idx);
    dR4(idx) = dR4(idx) - grate.*R4(idx);   dRv4(idx) = dRv4(idx) + grate.*R4(idx);
end

% vaccine wanning
dSv0 = dSv0 - ohm_vac * Sv0 ;
dSw0 = dSw0 + ohm_vac * Sv0 - ohm_vac2 * Sw0;
dSww0 = dSww0 + ohm_vac2 * Sw0 - ohm_vac3 * Sww0;
dS0 = dS0 + ohm_vac3 * Sww0;
dSv1   = dSv1   - ohm_vac  * Sv1;
dSw1   = dSw1   + ohm_vac  * Sv1   - ohm_vac2 * Sw1;
dSww1  = dSww1  + ohm_vac2 * Sw1   - ohm_vac3 * Sww1;
dS1    = dS1    + ohm_vac3 * Sww1;
dSv2   = dSv2   - ohm_vac  * Sv2;
dSw2   = dSw2   + ohm_vac  * Sv2   - ohm_vac2 * Sw2;
dSww2  = dSww2  + ohm_vac2 * Sw2   - ohm_vac3 * Sww2;
dS2    = dS2    + ohm_vac3 * Sww2;
dSv3   = dSv3   - ohm_vac  * Sv3;
dSw3   = dSw3   + ohm_vac  * Sv3   - ohm_vac2 * Sw3;
dSww3  = dSww3  + ohm_vac2 * Sw3   - ohm_vac3 * Sww3;
dS3    = dS3    + ohm_vac3 * Sww3;

% demographic
dS0(1) = dS0(1) + birthrate * N/1000/1;
dS0(2:agenum) = dS0(2:agenum) + mu .* q .* S0(1:agenum-1);
dS0 = dS0 - mu .* S0;    dS0(agenum) = dS0(agenum) - mu_final * S0(agenum);
dS1(2:agenum) = dS1(2:agenum) + mu .* q .* S1(1:agenum-1);
dS1 = dS1 - mu .* S1;    dS1(agenum) = dS1(agenum) - mu_final * S1(agenum);
dS2(2:agenum) = dS2(2:agenum) + mu .* q .* S2(1:agenum-1);
dS2 = dS2 - mu .* S2;    dS2(agenum) = dS2(agenum) - mu_final * S2(agenum);
dS3(2:agenum) = dS3(2:agenum) + mu .* q .* S3(1:agenum-1);
dS3 = dS3 - mu .* S3;    dS3(agenum) = dS3(agenum) - mu_final * S3(agenum);
dR1(2:agenum) = dR1(2:agenum) + mu .* q .* R1(1:agenum-1);
dR1 = dR1 - mu .* R1;    dR1(agenum) = dR1(agenum) - mu_final * R1(agenum);
dR2(2:agenum) = dR2(2:agenum) + mu .* q .* R2(1:agenum-1);
dR2 = dR2 - mu .* R2;    dR2(agenum) = dR2(agenum) - mu_final * R2(agenum);
dR3(2:agenum) = dR3(2:agenum) + mu .* q .* R3(1:agenum-1);
dR3 = dR3 - mu .* R3;    dR3(agenum) = dR3(agenum) - mu_final * R3(agenum);
dR4(2:agenum) = dR4(2:agenum) + mu .* q .* R4(1:agenum-1);
dR4 = dR4 - mu .* R4;    dR4(agenum) = dR4(agenum) - mu_final * R4(agenum);
dSv0(2:agenum) = dSv0(2:agenum) + mu .* q .* Sv0(1:agenum-1);
dSv0 = dSv0 - mu .* Sv0;    dSv0(agenum) = dSv0(agenum) - mu_final * Sv0(agenum);
dSv1(2:agenum) = dSv1(2:agenum) + mu .* q .* Sv1(1:agenum-1);
dSv1 = dSv1 - mu .* Sv1;    dSv1(agenum) = dSv1(agenum) - mu_final * Sv1(agenum);
dSv2(2:agenum) = dSv2(2:agenum) + mu .* q .* Sv2(1:agenum-1);
dSv2 = dSv2 - mu .* Sv2;    dSv2(agenum) = dSv2(agenum) - mu_final * Sv2(agenum);
dSv3(2:agenum) = dSv3(2:agenum) + mu .* q .* Sv3(1:agenum-1);
dSv3 = dSv3 - mu .* Sv3;    dSv3(agenum) = dSv3(agenum) - mu_final * Sv3(agenum);
dRv1(2:agenum) = dRv1(2:agenum) + mu .* q .* Rv1(1:agenum-1);
dRv1 = dRv1 - mu .* Rv1;    dRv1(agenum) = dRv1(agenum) - mu_final * Rv1(agenum);
dRv2(2:agenum) = dRv2(2:agenum) + mu .* q .* Rv2(1:agenum-1);
dRv2 = dRv2 - mu .* Rv2;    dRv2(agenum) = dRv2(agenum) - mu_final * Rv2(agenum);
dRv3(2:agenum) = dRv3(2:agenum) + mu .* q .* Rv3(1:agenum-1);
dRv3 = dRv3 - mu .* Rv3;    dRv3(agenum) = dRv3(agenum) - mu_final * Rv3(agenum);
dRv4(2:agenum) = dRv4(2:agenum) + mu .* q .* Rv4(1:agenum-1);
dRv4 = dRv4 - mu .* Rv4;    dRv4(agenum) = dRv4(agenum) - mu_final * Rv4(agenum);
dSw0(2:agenum) = dSw0(2:agenum) + mu .* q .* Sw0(1:agenum-1);
dSw0 = dSw0 - mu .* Sw0;    dSw0(agenum) = dSw0(agenum) - mu_final * Sw0(agenum);
dSw1(2:agenum) = dSw1(2:agenum) + mu .* q .* Sw1(1:agenum-1);
dSw1 = dSw1 - mu .* Sw1;    dSw1(agenum) = dSw1(agenum) - mu_final * Sw1(agenum);
dSw2(2:agenum) = dSw2(2:agenum) + mu .* q .* Sw2(1:agenum-1);
dSw2 = dSw2 - mu .* Sw2;    dSw2(agenum) = dSw2(agenum) - mu_final * Sw2(agenum);
dSw3(2:agenum) = dSw3(2:agenum) + mu .* q .* Sw3(1:agenum-1);
dSw3 = dSw3 - mu .* Sw3;    dSw3(agenum) = dSw3(agenum) - mu_final * Sw3(agenum);
dRw1(2:agenum) = dRw1(2:agenum) + mu .* q .* Rw1(1:agenum-1);
dRw1 = dRw1 - mu .* Rw1;    dRw1(agenum) = dRw1(agenum) - mu_final * Rw1(agenum);
dRw2(2:agenum) = dRw2(2:agenum) + mu .* q .* Rw2(1:agenum-1);
dRw2 = dRw2 - mu .* Rw2;    dRw2(agenum) = dRw2(agenum) - mu_final * Rw2(agenum);
dRw3(2:agenum) = dRw3(2:agenum) + mu .* q .* Rw3(1:agenum-1);
dRw3 = dRw3 - mu .* Rw3;    dRw3(agenum) = dRw3(agenum) - mu_final * Rw3(agenum);
dRw4(2:agenum) = dRw4(2:agenum) + mu .* q .* Rw4(1:agenum-1);
dRw4 = dRw4 - mu .* Rw4;    dRw4(agenum) = dRw4(agenum) - mu_final * Rw4(agenum);
dSww0(2:agenum) = dSww0(2:agenum) + mu .* q .* Sww0(1:agenum-1);
dSww0 = dSww0 - mu .* Sww0;    dSww0(agenum) = dSww0(agenum) - mu_final * Sww0(agenum);
dSww1(2:agenum) = dSww1(2:agenum) + mu .* q .* Sww1(1:agenum-1);
dSww1 = dSww1 - mu .* Sww1;    dSww1(agenum) = dSww1(agenum) - mu_final * Sww1(agenum);
dSww2(2:agenum) = dSww2(2:agenum) + mu .* q .* Sww2(1:agenum-1);
dSww2 = dSww2 - mu .* Sww2;    dSww2(agenum) = dSww2(agenum) - mu_final * Sww2(agenum);
dSww3(2:agenum) = dSww3(2:agenum) + mu .* q .* Sww3(1:agenum-1);
dSww3 = dSww3 - mu .* Sww3;    dSww3(agenum) = dSww3(agenum) - mu_final * Sww3(agenum);
dRww1(2:agenum) = dRww1(2:agenum) + mu .* q .* Rww1(1:agenum-1);
dRww1 = dRww1 - mu .* Rww1;    dRww1(agenum) = dRww1(agenum) - mu_final * Rww1(agenum);
dRww2(2:agenum) = dRww2(2:agenum) + mu .* q .* Rww2(1:agenum-1);
dRww2 = dRww2 - mu .* Rww2;    dRww2(agenum) = dRww2(agenum) - mu_final * Rww2(agenum);
dRww3(2:agenum) = dRww3(2:agenum) + mu .* q .* Rww3(1:agenum-1);
dRww3 = dRww3 - mu .* Rww3;    dRww3(agenum) = dRww3(agenum) - mu_final * Rww3(agenum);
dRww4(2:agenum) = dRww4(2:agenum) + mu .* q .* Rww4(1:agenum-1);
dRww4 = dRww4 - mu .* Rww4;    dRww4(agenum) = dRww4(agenum) - mu_final * Rww4(agenum);

% hosp cumm and death
cum_hosp = zeros(agenum,1);
cum_death = zeros(agenum,1);
cum_case = zeros(agenum,1);

cum_case = cum_case + cageadj.*cprop_1.*(beta1 * beta0 .* FOI .* S0);
cum_case = cum_case + cageadj.*cprop_2.*(i1 * beta1 * beta0 .* FOI .* S1);
cum_case = cum_case + cageadj.*cprop_3.*(i2 * beta1 * beta0 .* FOI .* S2);
cum_case = cum_case + cageadj.*cprop_4.*(i3 * beta1 * beta0 .* FOI .* S3);
cum_case = cum_case + cageadj.*(1 - e_n_d)*cprop_1.*((1 - e_n_i) * beta1 * beta0 .* FOI .* Sv0);
cum_case = cum_case + cageadj.*(1 - e_p_d)*cprop_2.*(i1 * (1 - e_p_i) * beta1 * beta0 .* FOI .* Sv1);
cum_case = cum_case + cageadj.*(1 - e_p_d)*cprop_3.*(i2 * (1 - e_p_i) * beta1 * beta0 .* FOI .* Sv2);
cum_case = cum_case + cageadj.*(1 - e_p_d)*cprop_4.*(i3 * (1 - e_p_i) * beta1 * beta0 .* FOI .* Sv3);
cum_case = cum_case + cageadj.*(1 - ew_n_d)*cprop_1.*((1 - ew_n_i) * beta1 * beta0 .* FOI .* Sw0);
cum_case = cum_case + cageadj.*(1 - ew_p_d)*cprop_2.*(i1 * (1 - ew_p_i) * beta1 * beta0 .* FOI .* Sw1);
cum_case = cum_case + cageadj.*(1 - ew_p_d)*cprop_3.*(i2 * (1 - ew_p_i) * beta1 * beta0 .* FOI .* Sw2);
cum_case = cum_case + cageadj.*(1 - ew_p_d)*cprop_4.*(i3 * (1 - ew_p_i) * beta1 * beta0 .* FOI .* Sw3);
cum_case = cum_case + cageadj.*(1 - eww_n_d)*cprop_1.*((1 - eww_n_i) * beta1 * beta0 .* FOI .* Sww0);
cum_case = cum_case + cageadj.*(1 - eww_p_d)*cprop_2.*(i1 * (1 - eww_p_i) * beta1 * beta0 .* FOI .* Sww1);
cum_case = cum_case + cageadj.*(1 - eww_p_d)*cprop_3.*(i2 * (1 - eww_p_i) * beta1 * beta0 .* FOI .* Sww2);
cum_case = cum_case + cageadj.*(1 - eww_p_d)*cprop_4.*(i3 * (1 - eww_p_i) * beta1 * beta0 .* FOI .* Sww3);

cum_hosp = cum_hosp + cageadj.*cprop_1.*h1 .* (beta1 * beta0 .* FOI .* S0);
cum_hosp = cum_hosp + cageadj.*cprop_2.*h2 .* (i1 * beta1 * beta0 .* FOI .* S1);
cum_hosp = cum_hosp + cageadj.*cprop_3.*h3 .* (i2 * beta1 * beta0 .* FOI .* S2);
cum_hosp = cum_hosp + cageadj.*cprop_4.*h4 .* (i3 * beta1 * beta0 .* FOI .* S3);
cum_hosp = cum_hosp + cageadj.*(1 - 0)*(1 - e_n_h)*cprop_1.*h1 .* ((1 - e_n_i) * beta1 * beta0 .* FOI .* Sv0);
cum_hosp = cum_hosp + cageadj.*(1 - 0)*(1 - e_p_h)*cprop_2.*h2 .* (i1 * (1 - e_p_i) * beta1 * beta0 .* FOI .* Sv1);
cum_hosp = cum_hosp + cageadj.*(1 - 0)*(1 - e_p_h)*cprop_3.*h3 .* (i2 * (1 - e_p_i) * beta1 * beta0 .* FOI .* Sv2);
cum_hosp = cum_hosp + cageadj.*(1 - 0)*(1 - e_p_h)*cprop_4.*h4 .* (i3 * (1 - e_p_i) * beta1 * beta0 .* FOI .* Sv3);
cum_hosp = cum_hosp + cageadj.*(1 - 0)*(1 - ew_n_h)*cprop_1.*h1 .* ((1 - ew_n_i) * beta1 * beta0 .* FOI .* Sw0);
cum_hosp = cum_hosp + cageadj.*(1 - 0)*(1 - ew_p_h)*cprop_2.*h2 .* (i1 * (1 - ew_p_i) * beta1 * beta0 .* FOI .* Sw1);
cum_hosp = cum_hosp + cageadj.*(1 - 0)*(1 - ew_p_h)*cprop_3.*h3 .* (i2 * (1 - ew_p_i) * beta1 * beta0 .* FOI .* Sw2);
cum_hosp = cum_hosp + cageadj.*(1 - 0)*(1 - ew_p_h)*cprop_4.*h4 .* (i3 * (1 - ew_p_i) * beta1 * beta0 .* FOI .* Sw3);
cum_hosp = cum_hosp + cageadj.*(1 - 0)*(1 - eww_n_h)*cprop_1.*h1 .* ((1 - eww_n_i) * beta1 * beta0 .* FOI .* Sww0);
cum_hosp = cum_hosp + cageadj.*(1 - 0)*(1 - eww_p_h)*cprop_2.*h2 .* (i1 * (1 - eww_p_i) * beta1 * beta0 .* FOI .* Sww1);
cum_hosp = cum_hosp + cageadj.*(1 - 0)*(1 - eww_p_h)*cprop_3.*h3 .* (i2 * (1 - eww_p_i) * beta1 * beta0 .* FOI .* Sww2);
cum_hosp = cum_hosp + cageadj.*(1 - 0)*(1 - eww_p_h)*cprop_4.*h4 .* (i3 * (1 - eww_p_i) * beta1 * beta0 .* FOI .* Sww3);

cum_death = cum_death + cageadj.*cprop_1.*f1 .* (beta1 * beta0 .* FOI .* S0);
cum_death = cum_death + cageadj.*cprop_2.*f2 .* (i1 * beta1 * beta0 .* FOI .* S1);
cum_death = cum_death + cageadj.*cprop_3.*f3 .* (i2 * beta1 * beta0 .* FOI .* S2);
cum_death = cum_death + cageadj.*cprop_4.*f4 .* (i3 * beta1 * beta0 .* FOI .* S3);
cum_death = cum_death + cageadj.*(1 - 0)*(1 - e_n_m)*cprop_1.*f1 .* ((1 - e_n_i) * beta1 * beta0 .* FOI .* Sv0);
cum_death = cum_death + cageadj.*(1 - 0)*(1 - e_p_m)*cprop_2.*f2 .* (i1 * (1 - e_p_i) * beta1 * beta0 .* FOI .* Sv1);
cum_death = cum_death + cageadj.*(1 - 0)*(1 - e_p_m)*cprop_3.*f3 .* (i2 * (1 - e_p_i) * beta1 * beta0 .* FOI .* Sv2);
cum_death = cum_death + cageadj.*(1 - 0)*(1 - e_p_m)*cprop_4.*f4 .* (i3 * (1 - e_p_i) * beta1 * beta0 .* FOI .* Sv3);
cum_death = cum_death + cageadj.*(1 - 0)*(1 - ew_n_m)*cprop_1.*f1 .* ((1 - ew_n_i) * beta1 * beta0 .* FOI .* Sw0);
cum_death = cum_death + cageadj.*(1 - 0)*(1 - ew_p_m)*cprop_2.*f2 .* (i1 * (1 - ew_p_i) * beta1 * beta0 .* FOI .* Sw1);
cum_death = cum_death + cageadj.*(1 - 0)*(1 - ew_p_m)*cprop_3.*f3 .* (i2 * (1 - ew_p_i) * beta1 * beta0 .* FOI .* Sw2);
cum_death = cum_death + cageadj.*(1 - 0)*(1 - ew_p_m)*cprop_4.*f4 .* (i3 * (1 - ew_p_i) * beta1 * beta0 .* FOI .* Sw3);
cum_death = cum_death + cageadj.*(1 - 0)*(1 - eww_n_m)*cprop_1.*f1 .* ((1 - eww_n_i) * beta1 * beta0 .* FOI .* Sww0);
cum_death = cum_death + cageadj.*(1 - 0)*(1 - eww_p_m)*cprop_2.*f2 .* (i1 * (1 - eww_p_i) * beta1 * beta0 .* FOI .* Sww1);
cum_death = cum_death + cageadj.*(1 - 0)*(1 - eww_p_m)*cprop_3.*f3 .* (i2 * (1 - eww_p_i) * beta1 * beta0 .* FOI .* Sww2);
cum_death = cum_death + cageadj.*(1 - 0)*(1 - eww_p_m)*cprop_4.*f4 .* (i3 * (1 - eww_p_i) * beta1 * beta0 .* FOI .* Sww3);

dydt = [
    dS0; dI1; dR1; dS1; dI2; dR2; dS2; dI3; dR3; dS3; dI4; dR4;
    dSv0; dIv1; dRv1; dSv1; dIv2; dRv2; dSv2; dIv3; dRv3; dSv3; dIv4; dRv4;
    dSw0; dIw1; dRw1; dSw1; dIw2; dRw2; dSw2; dIw3; dRw3; dSw3; dIw4; dRw4;
    dSww0; dIww1; dRww1; dSww1; dIww2; dRww2; dSww2; dIww3; dRww3; dSww3; dIww4; dRww4;
    cum_case;
    cum_hosp;
    cum_death
];