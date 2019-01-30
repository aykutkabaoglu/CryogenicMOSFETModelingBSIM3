function [ V, I ] = simulation( ocean_address, result_address )
% ocean dosyasi cagirilarak simulasyon baslatilir
ocean_command = sprintf('ocean -restore %s',ocean_address); 
unix(ocean_command)

result = importdata(result_address);
new_data = result.data;
V = new_data(:,1);
I = new_data(:,2);

    
end

