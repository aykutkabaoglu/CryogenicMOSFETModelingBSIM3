function [ ute_n, Ids, I ] = ute_iteration_algorithm(command, model_address,step, I, I_aim, Ids, Ids_aim, ocean_address_vds, ocean_address_vgs, result_address_vds, result_address_vgs, baslangic, son, type )
% UTE parametresi denklem ile cozdurulemediginde MATLAB hata donugunde iterasyonlarla UTE degeri cozumu bulunamayan degerinden kaydilir ve tekrar denklemlere sokulur boylece koku bulunamayan denklem degistirilmis olur
%%% Dosya Okuma    
    file = fopen(model_address);
formatSpec = '%s';
C = textscan(file,formatSpec);
[boyut, ~] = size(C{:});

for i=baslangic:son
    
     k=strfind(C{1}{i},'kt1=');
     if(k==1)
        kt1_n = textscan(C{1}{i},'kt1=%f');
        kt1_n = kt1_n{1}
        kt1_n_char = C{1}{i};
     end
     
     k=strfind(C{1}{i},'kt2=');
     if(k==1)
        kt2_n = textscan(C{1}{i},'kt2=%f');
        kt2_n = kt2_n{1}
        kt2_n_char = C{1}{i};
     end
     
     k=strfind(C{1}{i},'ute=');
     if(k==1)
        ute_n = textscan(C{1}{i},'ute=%f');
        ute_n = ute_n{1}
        ute_n_char = C{1}{i};        
     end 

     k=strfind(C{1}{i},'ua1=');
     if(k==1)
        ua1_n = textscan(C{1}{i},'ua1=%f');
        ua1_n = ua1_n{1}
        ua1_n_char = C{1}{i};
     end

     k=strfind(C{1}{i},'ub1=');
     if(k==1)
        ub1_n = textscan(C{1}{i},'ub1=%f');
        ub1_n = ub1_n{1}
        ub1_n_char = C{1}{i};        
     end

     k=strfind(C{1}{i},'uc1=');
     if(k==1)
        uc1_n = textscan(C{1}{i},'uc1=%f');
        uc1_n = uc1_n{1}
        uc1_n_char = C{1}{i};
     end
     
     k=strfind(C{1}{i},'at=');
     if(k==1)
        at_n = textscan(C{1}{i},'at=%f');
        at_n = at_n{1}
        at_n_char = C{1}{i};
     end
     
     k=strfind(C{1}{i},'prt=');
     if(k==1)
        prt_n = textscan(C{1}{i},'prt=%f');
        prt_n = prt_n{1}
        prt_n_char = C{1}{i};
     end
     
     k=strfind(C{1}{i},'delta=');
     if(k==1)
        delta_n = textscan(C{1}{i},'delta=%f');
        delta_n = delta_n{1}
        delta_n_char = C{1}{i};
     end

end
    
fclose(file);

if type == 'n'
ute_n_max = -0.280;
ute_n_min = -1.82;
elseif type == 'p'
ute_n_max = -0.020;
ute_n_min = -1.84;
end

ute_n_def = ute_n;
konum = 1;
%%%
    
    
% --------------------------ALGORİTMA--------------------------------------
if step~=0

%     % Limit değerler ile dosyadan okunan değerlerin karşılaştırılıp ekrana uyarı basılması
% if abs(kt1_n)>abs(kt1_n_max)
%     disp('UYARI: KT1 degeri limit degerinden daha buyuk')
% end
% if abs(ute_n)>abs(ute_n_max)
%     disp('UYARI: UTE degeri limit degerinden daha buyuk')    
% end
% if abs(ua1_n)>abs(ua1_n_max)
%     disp('UYARI: UA1 degeri limit degerinden daha buyuk')    
% end

% %%%%%%%%%%%% UTE %%%%%%%%%%%%%%%%%%%%%%%

    error=0; 
    best=0;
    for i=konum:size(I)
best = (I_aim(i)-I(i))^2+best; % maksimum error
   end  
    error_ds=0; 
    best_ds=0;
    for i=1:size(Ids)
best_ds = (Ids_aim(i)-Ids(i))^2+best_ds; % maksimum error
   end  
   
%    if abs(ute_n)<abs(ute_n_max)  
    
while error_ds>=0 
   ute_n = ute_n+(ute_n_max-ute_n_def)/step; 
   
           
        if ute_n<=ute_n_min||ute_n>=ute_n_max
            break
        end         
   
    % Değişkenlerin tekrardan karaktere dönüştürülmesi
    ute_n_char_new = sprintf('ute=%d', ute_n);
    ute_n_char_old = ute_n_char;
    % Parametrelerin yeni değerleri ile aynı dosya üzerine yazdırılması 
    execute = sprintf(command,ute_n_char,ute_n_char_new, model_address);
    system(execute);
    
%% Parametrelerin dosyadan tekrar okunması (hata almamak için)

file = fopen(model_address);
formatSpec = '%s';
C = textscan(file,formatSpec);

for i=baslangic:son
    
     k=strfind(C{1}{i},'ute=');
     if(k==1)
        ute_n = textscan(C{1}{i},'ute=%f');
        ute_n = ute_n{1};
        ute_n_char = C{1}{i};
     end
end  
fclose(file);
ute_n

%%
%%%%%%% yeni değerler ile analiz %%%%%%%%
[ Vds, Ids ] = simulation( ocean_address_vds, result_address_vds );
[ V, I] = simulation(ocean_address_vgs,result_address_vgs);


    %%

    %% Bir önceki error ile karşılaştırma
%         error=0;
%         for i=konum:size(I)
%              error=(I_aim(i)-I(i))^2+error; % değiştirilmiş parametrelerle hesaplanan error
%         end
        error_ds=0;
        for i=1:size(Ids)
             error_ds=(Ids_aim(i)-Ids(i))^2+error_ds; % değiştirilmiş parametrelerle hesaplanan error
        end
    %%    
%         if find(diff(I,2)<-6.0e-09)
%            % Parametrelerin eski değerleri ile aynı dosya üzerine yazdırılması
%            execute = sprintf(command,ute_n_char_new,ute_n_char_old, model_address);
%            system(execute);
%            break            
            
        if(error_ds<best_ds) % en düşük error değeri ile karşılaştırma
           best = error;
           best_ds = error_ds;
            % grafiğin yeni değerlerle güncellenmesi
           subplot(2,2,1)
           plot(V,I,'-')

           subplot(2,2,2)
           plot(Vds,Ids,'-')
           drawnow
           best
        else        
           % Parametrelerin eski değerleri ile aynı dosya üzerine yazdırılması
           execute = sprintf(command,ute_n_char_new,ute_n_char_old, model_address);
           system(execute);
           
               %% Parametrelerin dosyadan tekrar okunması (hata almamak için)
                file = fopen(model_address);
                formatSpec = '%s';
                C = textscan(file,formatSpec);
                    for i=baslangic:son

                 k=strfind(C{1}{i},'ute=');
                 if(k==1)
                    ute_n = textscan(C{1}{i},'ute=%f');
                    ute_n = ute_n{1};
                    ute_n_char = C{1}{i};
                 end
                    end  
                    fclose(file);
               %----
           
           %%%%%%% yeni değerler ile analiz %%%%%%%%
                [ Vds, Ids ] = simulation( ocean_address_vds, result_address_vds );
                [ V, I] = simulation(ocean_address_vgs,result_address_vgs);

                           subplot(2,2,1)
                           plot(V,I,'-')

                           subplot(2,2,2)
                           plot(Vds,Ids,'-')
                           drawnow
           
                
                            error=0; 
                            best=0;
                            for i=konum:size(I)
                        best = (I_aim(i)-I(i))^2+best; % maksimum error
                           end  
                            error_ds=0; 
                            best_ds=0;
                            for i=1:size(Ids)
                        best_ds = (Ids_aim(i)-Ids(i))^2+best_ds; % maksimum error
                           end  

                        %    if abs(ute_n)<abs(ute_n_max)  

                        while error_ds>=0 
                           ute_n = ute_n+(ute_n_min-ute_n_def)/step; 
                           
                                  
                                if ute_n<=ute_n_min||ute_n>=ute_n_max
                                    break
                                end         

                            % Değişkenlerin tekrardan karaktere dönüştürülmesi
                            ute_n_char_new = sprintf('ute=%d', ute_n);
                            ute_n_char_old = ute_n_char;
                            % Parametrelerin yeni değerleri ile aynı dosya üzerine yazdırılması 
                            execute = sprintf(command,ute_n_char,ute_n_char_new, model_address);
                            system(execute);

                        %% Parametrelerin dosyadan tekrar okunması (hata almamak için)

                        file = fopen(model_address);
                        formatSpec = '%s';
                        C = textscan(file,formatSpec);

                        for i=baslangic:son

                             k=strfind(C{1}{i},'ute=');
                             if(k==1)
                                ute_n = textscan(C{1}{i},'ute=%f');
                                ute_n = ute_n{1};
                                ute_n_char = C{1}{i};
                             end
                        end  
                        fclose(file);
                        ute_n

                        %%
                        %%%%%%% yeni değerler ile analiz %%%%%%%%
                        [ Vds, Ids ] = simulation( ocean_address_vds, result_address_vds );
                        [ V, I] = simulation(ocean_address_vgs,result_address_vgs);


                            %%

                            %% Bir önceki error ile karşılaştırma
                                error=0;
                                for i=konum:size(I)
                                     error=(I_aim(i)-I(i))^2+error; % değiştirilmiş parametrelerle hesaplanan error
                                end
                                error_ds=0;
                                for i=1:size(Ids)
                                     error_ds=(Ids_aim(i)-Ids(i))^2+error_ds; % değiştirilmiş parametrelerle hesaplanan error
                                end
                            %%    
                        %         if find(diff(I,2)<-6.0e-09)
                        %            % Parametrelerin eski değerleri ile aynı dosya üzerine yazdırılması
                        %            execute = sprintf(command,ute_n_char_new,ute_n_char_old, model_address);
                        %            system(execute);
                        %            break            

                                if(error_ds<best_ds) % en düşük error değeri ile karşılaştırma
                                   best = error;
                                   best_ds = error_ds;
                                    % grafiğin yeni değerlerle güncellenmesi
                                   subplot(2,2,1)
                                   plot(V,I,'-')

                                   subplot(2,2,2)
                                   plot(Vds,Ids,'-')
                                   drawnow
                                   best
                                else        
                                   % Parametrelerin eski değerleri ile aynı dosya üzerine yazdırılması
                                   execute = sprintf(command,ute_n_char_new,ute_n_char_old, model_address);
                                   system(execute);
                                   
                                       %% Parametrelerin dosyadan tekrar okunması (hata almamak için)
                                        file = fopen(model_address);
                                        formatSpec = '%s';
                                        C = textscan(file,formatSpec);
                                            for i=baslangic:son

                                         k=strfind(C{1}{i},'ute=');
                                         if(k==1)
                                            ute_n = textscan(C{1}{i},'ute=%f');
                                            ute_n = ute_n{1};
                                            ute_n_char = C{1}{i};
                                         end
                                            end  
                                            fclose(file);
                                        %----
                                   
                                   %%%%%%% yeni değerler ile analiz %%%%%%%%
                                        [ Vds, Ids ] = simulation( ocean_address_vds, result_address_vds );
                                        [ V, I] = simulation(ocean_address_vgs,result_address_vgs);
                                           subplot(2,2,1)
                                           plot(V,I,'-')

                                           subplot(2,2,2)
                                           plot(Vds,Ids,'-')
                                           drawnow

                                   break
                                end
         
                        end

           
           break
        end
   
end



end
% -------------------------------------------------------------------------

    
end

