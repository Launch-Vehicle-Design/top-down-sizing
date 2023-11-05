function [optimal,dvdisb] = iterativeConceptDesign(param, is_3stg, plot_imperial)
    
    set(0,'DefaultTextInterpreter','latex')
    set(0,'DefaultFigureColor',[1,1,1])
    set(groot,'defaultAxesFontSize',16)
    
    % parameter conditioning
    if ~exist("param","var")
        param = sysParam();
    end
    if ~exist("plot_imperial","var")
        plot_imperial = false;
    end

    max_dv_stg_1_perc = 0.99;
    if param.is_scram
        lambda = 1.401; R = 287;
        gloss = @(h1,h2) sqrt(2*param.Mearth*param.G*(1/(h1+param.Rearth)-1/(h2+param.Rearth)));
        max_dv_stg_1 = param.scram_mach*sqrt(lambda*R*param.scram_ceiling_temp) + ...
            1*gloss(param.release_alt,param.scram_ceiling) + 0.9*(param.dragloss+param.proploss+param.steeloss);
        max_dv_stg_1_perc = max_dv_stg_1/param.vrq;

        delta_v_types = categorical({'\Delta v_{orbital}','\Delta v_{grav loss}','\Delta v_{drag loss}', ...
            '\Delta v_{prop loss}','\Delta v_{steering loss}'});
        delta_v_types = reordercats(delta_v_types,{'\Delta v_{orbital}','\Delta v_{grav loss}','\Delta v_{drag loss}', ...
            '\Delta v_{prop loss}','\Delta v_{steering loss}'});
        total_delta_v = [param.orb_vel 0.8*gloss(param.release_alt,param.orb_alt) param.dragloss param.proploss param.steeloss];
        step_1_delta_v = [param.scram_mach*sqrt(lambda*R*param.scram_ceiling_temp) 0.8*gloss(param.release_alt,param.scram_ceiling) ...
            0.9*param.dragloss 0.9*param.proploss 0.9*param.steeloss];
        delta_v = [step_1_delta_v; total_delta_v-step_1_delta_v];

        if plot_imperial
            plot_delta_v = delta_v*3.28084;
            label_delta_v = "$\Delta$ v ($\frac{ft}{s}$)";
        else
            plot_delta_v = delta_v;
            label_delta_v = "$\Delta$ v ($\frac{m}{s}$)";
        end
        figure; b = bar(delta_v_types,plot_delta_v,"stacked",'EdgeColor','none','BarWidth',0.4); grid on
        b(1).FaceColor = [107 175 226]/255; % light blue
        b(2).FaceColor = [0 104 56]/255; % dark green
        title("Non Smeared Air Breathing Propulsion $\Delta v$ Allocation")
        ylabel(label_delta_v); legend("Step 1 Dual Ram/Scram","Step 2 - RDRE","interpreter","latex")
    end
    
    solid_boost_dv_stg_1_perc = 0;
    if param.is_scram && param.is_scram_solid_boost
        lambda = 1.4; R = 287;
        solid_boost_dv = param.scram_solid_boost_mach*sqrt(lambda*R*param.release_temp);
        solid_boost_dv_stg_1_perc = solid_boost_dv/param.vrq;

        solid_boost_delta_v = [param.scram_solid_boost_mach*sqrt(lambda*R*param.release_temp) 0 0.3*param.dragloss 0 0];
        delta_v = [solid_boost_delta_v; step_1_delta_v-solid_boost_delta_v; total_delta_v-step_1_delta_v];

        if plot_imperial
            plot_delta_v = delta_v*3.28084;
            label_delta_v = "$\Delta$ v ($\frac{ft}{s}$)";
        else
            plot_delta_v = delta_v;
            label_delta_v = "$\Delta$ v ($\frac{m}{s}$)";
        end
        figure; b = bar(delta_v_types,plot_delta_v,"stacked",'EdgeColor','none','BarWidth',0.4); grid on
        b(1).FaceColor = [247 129 52]/255; % orange
        b(2).FaceColor = [107 175 226]/255; % light blue
        b(3).FaceColor = [0 104 56]/255; % dark green
        title("Non Smeared Air Breathing Propulsion $\Delta v$ Allocation")
        ylabel(label_delta_v); legend("Step 1 Solid Boost","Step 2 - Scramjet","Step 3 - RDRE","interpreter","latex")
    end
    
    %% sizing iteration - two stage rocket
    vrq_ratio_stg_1 = 0.01:0.001:max_dv_stg_1_perc;
    vrq_ratio_stg_2 = 1-vrq_ratio_stg_1;
    
    dvrq_stg_1 = vrq_ratio_stg_1*param.vrq;
    dvrq_stg_2 = vrq_ratio_stg_2*param.vrq;
    
    sizing_table = nan([10 size(vrq_ratio_stg_1,2)]);
    for i = 1:size(vrq_ratio_stg_1,2)
        sizing_table(:,i) = TSTOSizing(dvrq_stg_2(i),dvrq_stg_1(i),...
            param.mPL,param.sigma_stg1,param.sigma_stg2,param.Isp_stg1*param.g0,param.Isp_stg2*param.g0);
    end
    fuel_vol_2stg = sizing_table(3,:)/param.density_stg2 + sizing_table(8,:)/param.density_stg1;

    [min_mass_2stg,ind_2stg] = min(sizing_table(6,:));
    fuel_vol_min_mass_2stg = fuel_vol_2stg(ind_2stg);
    
    if plot_imperial
        plot_m0 = sizing_table(6,:)*2.20462;
        plot_vp = fuel_vol_2stg*61023.7;
        limit_bounding_mass = param.bounding_mass*2.20462;
        limit_bounding_volu = param.bounding_box_volu*61023.7;
        min_identified_m0 = min_mass_2stg*2.20462;
        min_identified_vp = fuel_vol_min_mass_2stg*61023.7;
        label_m0 = "Vehicle Mass $m_0$ (lb)";
        label_vp = "Propellant Volume $V_p (inch^3)$";
        text_placement_y = 0.2*61023.7;
    else
        plot_m0 = sizing_table(6,:);
        plot_vp = fuel_vol_2stg;
        limit_bounding_mass = param.bounding_mass;
        limit_bounding_volu = param.bounding_box_volu;
        min_identified_m0 = min_mass_2stg;
        min_identified_vp = fuel_vol_min_mass_2stg;
        label_m0 = "Vehicle Mass $m_0$ (kg)";
        label_vp = "Propellant Volume $V_p (m^3)$";
        text_placement_y = 0.2;
    end
    figure; subplot(2,1,1)
    plot(vrq_ratio_stg_1*100,plot_m0,'k','LineWidth',1.2); hold on
    scatter(vrq_ratio_stg_1(ind_2stg)*100, min_identified_m0, 'ko','LineWidth',1.2);
    xline(vrq_ratio_stg_1(ind_2stg)*100,"k--",'LineWidth',1);
    title("TSTO Brute Force Sizing - Mass");
    xlabel("Step 1 ${\Delta v\%}$"); ylabel(label_m0);
    legend("Sizing Iteration","Optimal Point Design","interpreter","latex");
    xlim([0 100]); ylim([0 limit_bounding_mass]); grid on
    
    subplot(2,1,2)
    plot(vrq_ratio_stg_1*100,plot_vp,'k','LineWidth',1.2); hold on
    scatter(vrq_ratio_stg_1(ind_2stg)*100, min_identified_vp, 'ko','LineWidth',1.2);
    xline(vrq_ratio_stg_1(ind_2stg)*100,"k--",'LineWidth',1);
    text(vrq_ratio_stg_1(ind_2stg)*100+1,text_placement_y,num2str(vrq_ratio_stg_1(ind_2stg)*100)+"\%");
    title("TSTO Brute Force Sizing - Volume");
    xlabel("Step 1 ${\Delta v\%}$"); ylabel(label_vp);
    legend("Sizing Iteration","Optimal Point Design","interpreter","latex");
    xlim([0 100]); ylim([0 limit_bounding_volu]); grid on
    
    %% sizing iteration - three stage rocket
    vrq_ratio_stg_1_3stg = 0.01:0.002:0.99;
    if param.is_scram_solid_boost
        vrq_ratio_stg_1_3stg = solid_boost_dv_stg_1_perc;
    end
    vrq_ratio_stg_2_3stg = 0.01:0.002:max_dv_stg_1_perc-solid_boost_dv_stg_1_perc;
    [VRQ1,VRQ2] = meshgrid(vrq_ratio_stg_1_3stg,vrq_ratio_stg_2_3stg);
    
    VRQ3 = 1-VRQ1-VRQ2; VRQ3(VRQ3<0) = nan;
    
    dVRQ1 = VRQ1*param.vrq; dVRQ2 = VRQ2*param.vrq;
    dVRQ3 = VRQ3*param.vrq;
    
    sizing_table_3stg = nan([size(dVRQ1,1) size(dVRQ1,2) 15]);
    for i = 1:size(dVRQ1,1)
        for j = 1:size(dVRQ1,2)
            sizing_vec = TSTOSizing(dVRQ3(i,j),dVRQ2(i,j), ...
                param.mPL,param.sigma_stg2,param.sigma_stg3,param.Isp_stg2*param.g0,param.Isp_stg3*param.g0);
            % first stage
            mPL_stg_1 = sizing_vec(6); sigma_stg_1 = param.sigma_stg1; ve_stg_1 = param.Isp_stg1*param.g0;
            [m0_stg_1,ms_stg_1,mp_stg_1,mu_stg_1,PI_stg_1] = sizing(dVRQ1(i,j),mPL_stg_1,sigma_stg_1,ve_stg_1);
            sizing_vec(11:15) = [m0_stg_1 ms_stg_1 mp_stg_1 mu_stg_1 PI_stg_1]';
            if sum(sizing_vec < 0) ~= 0
                sizing_vec = nan(15,1);
            end
            sizing_table_3stg(i,j,:) = sizing_vec;
        end
    end
    
    fuel_vol_3stg = sizing_table_3stg(:,:,3)/param.density_stg3 + ...
        sizing_table_3stg(:,:,8)/param.density_stg2 + sizing_table_3stg(:,:,13)/param.density_stg1;
    [min_mass_3stg,ind_3stg] = min(sizing_table_3stg(:,:,11),[],"all");
    fuel_vol_min_mass_3stg = fuel_vol_3stg(ind_3stg);

    if is_3stg
        % SI
        disp("3 Stage - Min Mass = "+min_mass_3stg+" kg Fuel Vol = "+fuel_vol_min_mass_3stg+" m^3 as "+ ...
            fuel_vol_min_mass_3stg/param.bounding_box_volu*100+"% of the bounding box")
        optimal_3stg = nan([size(sizing_table_3stg,3),1]);
        for i = 1:size(sizing_table_3stg,3)
            temp = sizing_table_3stg(:,:,i);
            optimal_3stg(i) = temp(ind_3stg);
        end
        disp("3 Stage (kg) - m03="+optimal_3stg(1)+" ms3="+optimal_3stg(2)+" mp3="+optimal_3stg(3)+ ...
            " m02="+optimal_3stg(6)+" ms2="+optimal_3stg(7)+" mp2="+optimal_3stg(8)+ ...
            " m01="+optimal_3stg(11)+" ms1="+optimal_3stg(12)+" mp1="+optimal_3stg(13))

        % Imperial
        disp("3 Stage - Min Mass = "+min_mass_3stg*2.20462+" lb Fuel Vol = "+fuel_vol_min_mass_3stg*61023.7+" inch^3 as "+ ...
            fuel_vol_min_mass_3stg/param.bounding_box_volu*100+"% of the bounding box")
        optimal_3stg = nan([size(sizing_table_3stg,3),1]);
        for i = 1:size(sizing_table_3stg,3)
            temp = sizing_table_3stg(:,:,i);
            optimal_3stg(i) = temp(ind_3stg);
        end
        disp("3 Stage (lb) - m03="+optimal_3stg(1)*2.20462+" ms3="+optimal_3stg(2)*2.20462+" mp3="+optimal_3stg(3)*2.20462+ ...
            " m02="+optimal_3stg(6)*2.20462+" ms2="+optimal_3stg(7)*2.20462+" mp2="+optimal_3stg(8)*2.20462+ ...
            " m01="+optimal_3stg(11)*2.20462+" ms1="+optimal_3stg(12)*2.20462+" mp1="+optimal_3stg(13)*2.20462)

        disp("3 Stage DV - 1st stg: "+VRQ1(ind_3stg)+", 2nd stg: "+VRQ2(ind_3stg)+", 3rd stg: "+VRQ3(ind_3stg))
        disp("################## 3 stage concepts ######################")

        optimal = optimal_3stg;
        dvdisb = [VRQ1(ind_3stg) VRQ2(ind_3stg) VRQ3(ind_3stg)];
    else
        % SI
        disp("2 Stage - Min Mass = "+min_mass_2stg+" kg Fuel Vol = "+fuel_vol_min_mass_2stg+" m^3 as "+ ...
            fuel_vol_min_mass_2stg/param.bounding_box_volu*100+"% of the bounding box")
        optimal_2stg = sizing_table(:,ind_2stg);
        disp("2 Stage (kg) - m02="+optimal_2stg(1)+" ms2="+optimal_2stg(2)+" mp2="+optimal_2stg(3)+ ...
            " m01="+optimal_2stg(6)+" ms1="+optimal_2stg(7)+" mp1="+optimal_2stg(8))
        % imperial
        disp("2 Stage - Min Mass = "+min_mass_2stg*2.20462+" lb Fuel Vol = "+fuel_vol_min_mass_2stg*61023.7+" inch^3 as "+ ...
            fuel_vol_min_mass_2stg/param.bounding_box_volu*100+"% of the bounding box")
        optimal_2stg = sizing_table(:,ind_2stg);
        disp("2 Stage (lb) - m02="+optimal_2stg(1)*2.20462+" ms2="+optimal_2stg(2)*2.20462+" mp2="+optimal_2stg(3)*2.20462+ ...
            " m01="+optimal_2stg(6)*2.20462+" ms1="+optimal_2stg(7)*2.20462+" mp1="+optimal_2stg(8)*2.20462)

        disp("2 Stage DV - 1st stg: "+vrq_ratio_stg_1(ind_2stg)+", 2nd stg: "+vrq_ratio_stg_2(ind_2stg))
        disp("################## 2 stage concepts ######################")

        optimal = optimal_2stg;
        dvdisb = [vrq_ratio_stg_1(ind_2stg) vrq_ratio_stg_2(ind_2stg)];
    end
    disp("")

    
    if ~(param.is_scram && param.is_scram_solid_boost)
        [X,Y,Z] = sphere; C(:,:,1) = zeros(size(X)); C(:,:,2) = zeros(size(X)); C(:,:,3) = zeros(size(X));

        if plot_imperial
            plot_m0 = sizing_table_3stg(:,:,11)*2.20462;
            plot_vp = fuel_vol_3stg*61023.7;
            limit_bounding_mass = param.bounding_mass*2.20462;
            limit_bounding_volu = param.bounding_box_volu*61023.7;
            min_identified_m0 = min_mass_3stg*2.20462;
            min_identified_vp = fuel_vol_min_mass_3stg*61023.7;
            label_m0 = "Vehicle Mass (lb)";
            label_vp = "Propellant Volume $V_p (inch^3)$";
        else
            plot_m0 = sizing_table_3stg(:,:,11);
            plot_vp = fuel_vol_3stg;
            limit_bounding_mass = param.bounding_mass;
            limit_bounding_volu = param.bounding_box_volu;
            min_identified_m0 = min_mass_3stg;
            min_identified_vp = fuel_vol_min_mass_3stg;
            label_m0 = "Vehicle Mass (kg)";
            label_vp = "Propellant Volume $V_p (m^3)$";
        end
    
        figure; surf(VRQ1*100,VRQ2*100,plot_m0,"EdgeColor","interp",'FaceAlpha',0.5); hold on
        clim([0,limit_bounding_mass]); zlim([0,limit_bounding_mass]); xlim([0,100]); ylim([0,100])
        surf(X+VRQ1(ind_3stg)*100,Y+VRQ2(ind_3stg)*100,Z*20+min_identified_m0,C,"EdgeColor","none");
        title("ThreeSTO Brute Force Sizing - Mass"); colorbar
        zlabel(label_m0); xlabel("Step 1 ${\Delta v\%}$"); ylabel("Step 2 ${\Delta v\%}$");
        legend("Sizing Iteration","Optimal Point Design","interpreter","latex");

        figure; contourf(VRQ1*100,VRQ2*100,plot_m0,...
            min_identified_m0:(limit_bounding_mass-min_identified_m0)/41:limit_bounding_mass,"EdgeColor","none"); hold on
        clim([min_identified_m0,limit_bounding_mass]); xlim([0,100]); ylim([0,100])
        scatter(VRQ1(ind_3stg)*100,VRQ2(ind_3stg)*100,"ro","LineWidth",1.2);
        xline(VRQ1(ind_3stg)*100,"k--",'LineWidth',1); yline(VRQ2(ind_3stg)*100,"k--",'LineWidth',1);
        text(VRQ1(ind_3stg)*100+0.3,30.2,""+VRQ1(ind_3stg)*100+"\%","interpreter","latex");
        text(38,VRQ2(ind_3stg)*100+0.3,num2str(VRQ2(ind_3stg)*100)+"\%","interpreter","latex");
        title("ThreeSTO Brute Force Sizing - Mass"); c = colorbar; c.Label.String = label_m0;
        xlabel("Step 1 ${\Delta v\%}$"); ylabel("Step 2 ${\Delta v\%}$");
        legend("Sizing Iteration","Optimal Point Design","interpreter","latex");
        
        figure; surf(VRQ1*100,VRQ2*100,plot_vp,"EdgeColor","interp",'FaceAlpha',0.5); hold on
        surf(X+VRQ1(ind_3stg)*100,Y+VRQ2(ind_3stg)*100,Z*0.02+min_identified_vp,C,"EdgeColor","none");
        clim([0,limit_bounding_volu]); zlim([0,limit_bounding_volu]); xlim([0,100]); ylim([0,100])
        title("ThreeSTO Brute Force Sizing - Volume"); colorbar
        zlabel(label_vp); xlabel("Step 1 ${\Delta v\%}$"); ylabel("Step 2 ${\Delta v\%}$");
        legend("Sizing Iteration","Optimal Point Design","interpreter","latex");
        
        % figure; surf(VRQ1*100,VRQ2*100,100*sizing_table_3stg(:,:,11)/param.bounding_mass,"EdgeColor","interp"); hold on
        % surf(VRQ1*100,VRQ2*100,100*fuel_vol_3stg/param.bounding_box_volu,"EdgeColor","interp"); 
        % clim([0,100]); zlim([0,100]); xlim([0,100]); ylim([0,100])
        % xlabel("Step 1 ${\Delta v\%}$"); ylabel("Step 2 ${\Delta v\%}$");
    else
        if plot_imperial
            plot_m0 = sizing_table_3stg(:,:,11)*2.20462;
            plot_vp = fuel_vol_3stg*61023.7;
            limit_bounding_mass = param.bounding_mass*2.20462;
            limit_bounding_volu = param.bounding_box_volu*61023.7;
            min_identified_m0 = min_mass_3stg*2.20462;
            min_identified_vp = fuel_vol_min_mass_3stg*61023.7;
            label_m0 = "Vehicle Mass $m_0$ (lb)";
            label_vp = "Propellant Volume $V_p (inch^3)$";
            text_placement_y = 0.2*61023.7;
        else
            plot_m0 = sizing_table_3stg(:,:,11);
            plot_vp = fuel_vol_3stg;
            limit_bounding_mass = param.bounding_mass;
            limit_bounding_volu = param.bounding_box_volu;
            min_identified_m0 = min_mass_3stg;
            min_identified_vp = fuel_vol_min_mass_3stg;
            label_m0 = "Vehicle Mass $m_0$ (kg)";
            label_vp = "Propellant Volume $V_p (m^3)$";
            text_placement_y = 0.2;
        end

        figure; subplot(2,1,1)
        comb_VRQ12 = vrq_ratio_stg_1_3stg+VRQ2;
        plot(comb_VRQ12*100,plot_m0,'k','LineWidth',1.2); hold on
        scatter(comb_VRQ12(ind_3stg)*100, min_identified_m0, 'ko','LineWidth',1.2);
        xline(comb_VRQ12(ind_3stg)*100,"k--",'LineWidth',1);
        xline(vrq_ratio_stg_1_3stg*100,"k:",'LineWidth',1);
        title("ThreeSTO Scramjet with Solid Booster Brute Force Sizing - Mass");
        xlabel("Combined Step 1 and 2 ${\Delta v\%}$"); ylabel(label_m0);
        legend("Sizing Iteration","Optimal Point Design","interpreter","latex");
        xlim([0 100]); ylim([0 limit_bounding_mass]); grid on
        
        subplot(2,1,2)
        plot(comb_VRQ12*100,plot_vp,'k','LineWidth',1.2); hold on
        scatter(comb_VRQ12(ind_3stg)*100, min_identified_vp, 'ko','LineWidth',1.2);
        xline(comb_VRQ12(ind_3stg)*100,"k--",'LineWidth',1);
        text(comb_VRQ12(ind_3stg)*100+1,text_placement_y,"2nd Step "+num2str(VRQ2(ind_3stg)*100)+"\%");
        xline(vrq_ratio_stg_1_3stg*100,"k:",'LineWidth',1);
        text(vrq_ratio_stg_1_3stg*100+1,text_placement_y,"1st Step "+num2str(vrq_ratio_stg_1_3stg*100,3)+"\%");
        title("ThreeSTO Scramjet with Solid Booster Brute Force Sizing - Volume");
        xlabel("Combined Step 1 and 2 ${\Delta v\%}$"); ylabel(label_vp);
        legend("Sizing Iteration","Optimal Point Design","interpreter","latex");
        xlim([0 100]); ylim([0 limit_bounding_volu]); grid on
    end
    
    % close all
end