function [ vth_graph, vth_sd] = vth_calculation( result_address_vgs_vth,type )
    
%%%%%%%%%%% Vth değerinin eğri çizilerek bulunması %%%%%%%%%%%

A = importdata(result_address_vgs_vth);
hedef = A.data;
V_aim_vth = hedef(:,1);
I_aim_vth = hedef(:,2);


if type =='n'
[~,loc]=max(diff(I_aim_vth));
elseif type == 'p'
[~,loc]=max(abs(diff(I_aim_vth)));   
end

if loc>=3
m= (I_aim_vth(loc)-I_aim_vth(loc-2))/(V_aim_vth(loc)-V_aim_vth(loc-2));
else
m= (I_aim_vth(loc)-I_aim_vth(loc+2))/(V_aim_vth(loc)-V_aim_vth(loc+2)); 
end
b = I_aim_vth(loc) - m*V_aim_vth(loc);


vth_graph = -b/m;

[~,loc]=max(smooth(V_aim_vth(1:size(V_aim_vth)-2,1),diff(I_aim_vth,2)));
vth_sd = V_aim_vth(loc);
if type == 'p'
   vth_sd = 1.8-vth_sd; 
   vth_graph = 1.8-vth_graph;
end   

end

