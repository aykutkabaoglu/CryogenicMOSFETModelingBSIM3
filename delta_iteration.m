function [ delta_n, Vds, Ids ] = delta_iteration( model_address, baslangic, son, command, step, Ids, Vds, Ids_aim, Vds_aim, ocean_address_vds, result_address_vds, delta_n_max, error_vds, vdsat_difference )
% delta parametresi grafikteki hata oranini azaltmak icin ve Vdsat'i daha iyi yakinsamak icin iterasyonlarla degistirilir
file = fopen(model_address);
formatSpec = '%s';
C = textscan(file,formatSpec);
[boyut, ~] = size(C{:});

for i=baslangic:son
     k=strfind(C{1}{i},'delta=');
     if(k==1)
        delta_n = textscan(C{1}{i},'delta=%f');
        delta_n = delta_n{1}
        delta_n_char = C{1}{i};
     end

end
    delta_n_def = delta_n;
    %%%%%%%%%%%%%% DELTA %%%%%%%%%%%%%%%%%%%

%if max(abs(diff(Ids,2))) < max(abs(diff(Ids_aim,2))) && delta_n ~= 0
if vdsat_difference<0 && delta_n ~= 0
    
    error_diff=0;
    best_diff=0;

best_diff = max(abs(diff(Ids_aim,2)))-max(abs(diff(Ids,2)))^2; % maksimum error

while error_diff>=0 
    
   delta_n = delta_n-(delta_n_def/step); 
   
    % Değişkenlerin tekrardan karaktere dönüştürülmesi
    delta_n_char_new = sprintf('delta=%d', delta_n);
    delta_n_char_old = delta_n_char;
    
    % Parametrelerin yeni değerleri ile aynı dosya üzerine yazdırılması  
    execute = sprintf(command,delta_n_char,delta_n_char_new,model_address);
    system(execute);
    
%% Parametrelerin dosyadan tekrar okunması (hata almamak için)

file = fopen(model_address);
formatSpec = '%s';
C = textscan(file,formatSpec);

for i=baslangic:son
    
     k=strfind(C{1}{i},'delta=');
     if(k==1)
        delta_n = textscan(C{1}{i},'delta=%f');
        delta_n = delta_n{1};
        delta_n_char = C{1}{i};
     end
end  
fclose(file);
delta_n

%%

    %Yeni parametreler ile analiz
    ocean_command = sprintf('ocean -restore %s',ocean_address_vds); 
    unix(ocean_command)
       
    result = importdata(result_address_vds);
    new_data = result.data;
    Vds = new_data(:,1);
    Ids = new_data(:,2);
    %%

    %% Bir önceki error ile karşılaştırma
        error_ds=0;
        for i=1:size(Ids)
             error_ds=(Ids_aim(i)-Ids(i))^2+error_ds; % değiştirilmiş parametrelerle hesaplanan error
        end 
                   
        if(error_ds<error_vds) % en düşük error değeri ile karşılaştırma
           error_vds = error_ds;

       else
           % Parametrelerin eski değerleri ile aynı dosya üzerine yazdırılması
           execute = sprintf(command,delta_n_char_new,delta_n_char_old,model_address);
           system(execute);
           break
        end
        
       if delta_n>=delta_n_max
            break
       end     
end  
    
    
%elseif max(abs(diff(Ids,2))) > max(abs(diff(Ids_aim,2)))
elseif vdsat_difference > 0
  
    error_diff=0;
    best_diff=0;
    
best_diff = max(abs(diff(Ids_aim,2)))-max(abs(diff(Ids,2)))^2; % maksimum error

while error_diff>=0 
    
   delta_n = delta_n+(delta_n_max-delta_n_def)/step; 
   
    % Değişkenlerin tekrardan karaktere dönüştürülmesi
    delta_n_char_new = sprintf('delta=%d', delta_n);
    delta_n_char_old = delta_n_char;
    
    % Parametrelerin yeni değerleri ile aynı dosya üzerine yazdırılması  
    execute = sprintf(command,delta_n_char,delta_n_char_new,model_address);
    system(execute);
    
%% Parametrelerin dosyadan tekrar okunması (hata almamak için)

file = fopen(model_address);
formatSpec = '%s';
C = textscan(file,formatSpec);

for i=baslangic:son
    
     k=strfind(C{1}{i},'delta=');
     if(k==1)
        delta_n = textscan(C{1}{i},'delta=%f');
        delta_n = delta_n{1};
        delta_n_char = C{1}{i};
     end
end  
fclose(file);
delta_n

%%

    %Yeni parametreler ile analiz
    ocean_command = sprintf('ocean -restore %s',ocean_address_vds); 
    unix(ocean_command)
       
    result = importdata(result_address_vds);
    new_data = result.data;
    Vds = new_data(:,1);
    Ids = new_data(:,2);
    %%

    %% Bir önceki error ile karşılaştırma       
        error_ds=0;
        for i=1:size(Ids)
             error_ds=(Ids_aim(i)-Ids(i))+error_ds; % değiştirilmiş parametrelerle hesaplanan error
        end 

        toplam=0;
        for i=1:size(Ids)
            toplam = abs(Ids_aim(i))+toplam;
        end  
        
        error_ds = error_ds/toplam*100;
        if(error_ds<error_vds) % en düşük error değeri ile karşılaştırma
           error_vds = error_ds;

       else
           % Parametrelerin eski değerleri ile aynı dosya üzerine yazdırılması
           execute = sprintf(command,delta_n_char_new,delta_n_char_old,model_address);
           system(execute);
           break
        end
        
       if delta_n>=delta_n_max
            break
       end     
end
        
end
    
end

