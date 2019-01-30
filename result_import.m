% tum olcumler dosyaya kaydedilirken ayni formatta kaydedilmistir ve bu .csv dosyalari wafer ve chip isimleri ile klasorlenmistir.
% asagidaki klasor isimlerini yani wafer_name ve chip_name'leri dogru yazdiginizda o klasor altindaki tum .csv dosyalari otomatik olarak okunarak transistor adlari ve olcum gerilimlerine gore .mat dosyasina cevirilir
% dosya adi formati '%d_%3s_%3s%d-%d_%3s-%3s','.csv'
% ornek1: 01_VDS_VGS0.6-1.8_VBS0.csv
% ornek2: 01_VGS_VBS0-0.4_VDS1.8.csv
% daha fazla ornek icin olcum sonuclari klasoru incelenebilir
clc
clear
format short
wafer_name = 'W06';
chip_name = 'C8';
fileDir = sprintf('./meas_results_backup/%s/%s/',wafer_name,chip_name); % klasor yolu
j=1;
formatSpec = '%f%f%f%f%f%f%*s%*s%*s%*s%*s%*s%[^\n\r]';

F = dir(fullfile(fileDir,'*.csv'));
for ii = 1:length(F)
    filename = F(ii).name;
    fileID = fopen(fullfile(fileDir,F(ii).name));
    dataArray = textscan(fileID, formatSpec, 'Delimiter', ',', 'HeaderLines' , 6 , 'ReturnOnError', false);
    
    format2 = strcat('%d_%3s_%3s%d-%d_%3s-%3s','.csv');
    name_array=textscan(filename,format2);
    no = name_array{:,1};
    sweep_type = name_array{:,2,1};
    sweep_type = sweep_type{1};
    second_sweep = name_array{:,3};
    second_sweep = second_sweep{1};
    bias_type = name_array{:,6};
    bias_type = bias_type{1};
    if bias_type == 'VBS'
        bias = name_array{:,7};
        if isempty(bias)
            bias = 0;
        else
            bias = bias{1};
            if bias == '0.1'
                bias = 1;
            elseif bias == '0.2'
                bias = 2;
            elseif bias == '0.3'
                bias = 3;
            elseif bias == '0.4'
                bias = 4;
            end
        end
    elseif bias_type == 'VDS'
        format2 = strcat('%d_%3s_%3s%d-%d_%3s%f','.csv');
        name_array=textscan(filename,format2);
        bias = name_array{:,7};
        if bias == 1.8
            bias = 1;
        elseif bias == 0.05
            bias = 2;
        elseif bias == 0.01
            bias = 3;
        end
    end
    
    %% Close the text file.
    fclose(fileID);
    
    if sweep_type =='VDS'
        VD = dataArray{:, 1};
        ID = dataArray{:, 2};
        VG = dataArray{:, 3};
        VB = bias;
        
        %        eval([sprintf('%s%sT%d_VD',wafer_name,chip_name,no),'= VD;']);
        %        eval([sprintf('%s%sT%d_ID_VB%d',wafer_name,chip_name,no,bias),'= ID;']); %ID-VD values for all bias conditions
        if bias == 0
            i = 1;
            [boyut, ~] = size(VG);
            while VG(boyut)== 1.8
                VD0(i,1) = VD(boyut);
                ID0(i,1) = ID(boyut);
                boyut = boyut-1;
                i = i+1;
            end
            VD0 = flipud(VD0);
            ID0 = flipud(ID0);
            eval([sprintf('%s%sT%d_VD0',wafer_name,chip_name,no),'= VD0;']);
            eval([sprintf('%s%sT%d_ID0',wafer_name,chip_name,no),'= ID0;']); %ID-VD values for VGS=1.8 VBS=0
        elseif bias == 4
            i = 1;
            [boyut, ~] = size(VG);
            while VG(boyut)== 1.8
                VD0(i,1) = VD(boyut);
                ID0(i,1) = ID(boyut);
                boyut = boyut-1;
                i = i+1;
            end
            VD0 = flipud(VD0);
            ID0 = flipud(ID0);
            eval([sprintf('%s%sT%d_VD0',wafer_name,chip_name,no),'= VD0;']);
            eval([sprintf('%s%sT%d_ID0_VBlow',wafer_name,chip_name,no),'= ID0;']); %ID-VD values for VGS=1.8 VBS=-0.4
        end
        
    elseif sweep_type =='VGS'
        VG = dataArray{:, 1};
        ID = dataArray{:, 2};
        VB = dataArray{:, 6};
        VD = bias;
        
        %        eval([sprintf('%s%sT%d_VG',wafer_name,chip_name,no),'= VG;']);
        %        eval([sprintf('%s%sT%d_ID_VD%d',wafer_name,chip_name,no,bias),'= ID;']); %ID-VG values for all bias conditions
        
        if bias == 1
            i = 1;
            while VB(i)== 0
                VG1(i,1) = VG(i);
                ID1(i,1) = ID(i);
                i = i+1;
            end
            eval([sprintf('%s%sT%d_VG0',wafer_name,chip_name,no),'= VG1;']);
            eval([sprintf('%s%sT%d_ID0_VG',wafer_name,chip_name,no),'= ID1;']); %ID-VG values for VDS=1.8 VBS=0
            
        %% VTH CALCULATION %%
        elseif bias == 2
            i = 1;
            while VB(i)== 0
                VG1(i,1) = VG(i);
                ID1(i,1) = ID(i);
                i = i+1;
            end
            eval([sprintf('%s%sT%d_VG0',wafer_name,chip_name,no),'= VG1;']);
            eval([sprintf('%s%sT%d_ID0_VG_VDlow',wafer_name,chip_name,no),'= ID1;']); %ID-VG values for VDS=0.05 VBS=0
            
            if no >=9 && no <=13
                type = 'p';
            else
                type = 'n';
            end
            % ID1 = smooth(ID1);
            [vth_graph, vth_sd] = vth_calculation2(VG1, ID1, type);
            vth_app = vth_sd;
            if abs(abs(vth_graph)-abs(vth_sd))>0.2
                vth_app = vth_graph;
            end
            %            eval([sprintf('%s%sT%d_Vth',wafer_name,chip_name,no),'= vth_app;']); % VTH values for VDS = 0.05
            %            VTH_array{j,1} = sprintf('%s%sT%d',wafer_name,chip_name,no); % VTH array for all size of devices
            %            VTH_array{j,2} = vth_app;
            j = j+1;
        end
    end
    
end
%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;





