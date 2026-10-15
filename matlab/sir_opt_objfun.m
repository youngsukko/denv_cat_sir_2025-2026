function daly_total = sir_opt_objfun(x, pre, ode_inp, tspan, A5, S, annual_national)
%% sir_opt_objfun
%  x         : [S×1] 상대 배분 비율 (비음수)
%  pre       : 사전계산된 초기조건 struct (A0~A4, R0, popbystate)
%  ode_inp   : inputs struct (beta0, cageadj, ... vacorder 포함)
%  returns   : 전국 5년 누적 DALY 합계

x = max(x, 0);
doses_annual = annual_national * x / sum(x);   % 비율 → 실제 접종수

daly_total = 0;
for si = 1:S
    ode_inp.beta1 = pre.R0pred(si, pre.T+1:pre.T+5) * 52;
    ode_inp.vacc  = doses_annual(si);

    y0 = sir_make_y0(pre.A0{si}, pre.A1{si}, pre.A2{si}, ...
                     pre.A3{si}, pre.A4{si}, pre.R0si(si), A5);

    try
        [~,y] = ode23(@(t,yy) func_ode_alt_v5(t,yy,ode_inp), tspan, y0);
        Casen = y(:, A5*48+(1:A5));
        Hosp  = y(:, A5*49+(1:A5));
        Dead  = y(:, A5*50+(1:A5));
        YLD   = Casen*0.00084 + Hosp*0.0051;
        YLL   = Dead .* pre.lifeexp_age;
        daly_total = daly_total + sum(YLD(end,:)) + sum(YLL(end,:));
    catch
        daly_total = daly_total + 1e12;   % 실패시 큰 값
    end
end
end
