function [coincidence] = solar_load_factor(available_solar,loadprof)


loadprof.lighting_all = loadprof.cumulative.non_critical_light + loadprof.cumulative.critical_light;
loadprof.cell_all = loadprof.cumulative.non_critical_phone + loadprof.cumulative.critical_phone;
loadprof.fan_all = loadprof.cumulative.non_critical_fan + loadprof.cumulative.critical_fan;

coincidence.lights = zeros(length(available_solar),1);
coincidence.cell = zeros(length(available_solar),1);
coincidence.fan = zeros(length(available_solar),1);


for i=1:length(available_solar)
    
    if loadprof.lighting_all(i,end) > 0
        coincidence.lights(i) = available_solar(i) / loadprof.lighting_all(i,end);
        
        if coincidence.lights(i) > 1
            coincidence.lights(i) = 1;
        end
    end
    
    
    if loadprof.cell_all(i,end) > 0
        coincidence.cell(i) = available_solar(i) / loadprof.cell_all(i,end);
        
        if coincidence.cell(i) > 1
            coincidence.cell(i) = 1;
        end
    end
    
    if loadprof.fan_all(i,end) > 0
        coincidence.fan(i) = available_solar(i) / loadprof.fan_all(i,end);
        
        if coincidence.fan(i) > 1
            coincidence.fan(i) = 1;
        end
        
    end
    
    
end



coincidence.fan_avg = mean(coincidence.fan(coincidence.fan>0));
coincidence.lights_avg = mean(coincidence.lights(coincidence.lights>0));
coincidence.cell_avg = mean(coincidence.cell(coincidence.cell>0));



end