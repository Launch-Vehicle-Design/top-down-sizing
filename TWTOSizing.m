%% FUNCTION - Two Stage Iteration
function sizing_vec = TWTOSizing(dvrq_stg_2,dvrq_stg_1,mPL,sig1,sig2,ve1,ve2)

    sizing_vec = nan(10,1);

    % second stage
    mPL_stg_2 = mPL;
    sigma_stg_2 = sig2;
    ve_stg_2 = ve2;
    [m0_stg_2,ms_stg_2,mp_stg_2,mu_stg_2,PI_stg_2] = sizing(dvrq_stg_2,mPL_stg_2,sigma_stg_2,ve_stg_2);
    
    % first stage
    mPL_stg_1 = m0_stg_2;
    sigma_stg_1 = sig1;
    ve_stg_1 = ve1;
    [m0_stg_1,ms_stg_1,mp_stg_1,mu_stg_1,PI_stg_1] = sizing(dvrq_stg_1,mPL_stg_1,sigma_stg_1,ve_stg_1);

    sizing_vec(1:5) = [m0_stg_2 ms_stg_2 mp_stg_2 mu_stg_2 PI_stg_2]';
    sizing_vec(6:10) = [m0_stg_1 ms_stg_1 mp_stg_1 mu_stg_1 PI_stg_1]';
    if sum(sizing_vec < 0) ~= 0
        sizing_vec = nan(10,1);
    end