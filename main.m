clc
clear
tic

%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%
manual_input = false;
% olcum verileri manual olarak girildiginde buradaki dosya yollari guncellenecek
if manual_input
    Id_Vg_Vdmin_path = './meas_results_manual/T1/Id_Vg_Vdmin.csv';
    Id_Vg_Vdmin1_path = './meas_results_manual/T1/Id_Vg_Vdmin1.csv'; %isminde 1 olanlar farkli vbs gerilimindeki olcumler
    Id_Vg_Vdmax_path = './meas_results_manual/T1/Id_Vg_Vdmax.csv';
    Id_Vg_Vdmax1_path = './meas_results_manual/T1/Id_Vg_Vdmax1.csv';
    Id_Vd_Vgmax_path = './meas_results_manual/T1/Id_Vd_Vgmax.csv';
    Id_Vd_Vgmax1_path = './meas_results_manual/T1/Id_Vd_Vgmax1.csv';
    %    Id_Vg_Vdmin_s_path = './meas_results_manual/T1/Id_Vg_Vdmin.csv'; % farkli kanal boyutundaki olcum sonucu dahil edilerek kt1l hesaplanmak istendiginde burasi aktif edilir
else
    transistor_no = 2; % ölçüm dosyalarından sonuçları otomatik çekebilmek için kullanilir
    % No - W(um), L(um) , type
    %  1 - 0.24 , 0.18  , n
    %  2 - 10   , 0.18  , n
    %  3 - 0.24 , 10    , n
    %  4 - 10   , 10    , n
    %  5 - 35   , 0.5   , n
    %  9 - 0.24 , 0.18  , p
    % 10 - 10   , 0.18  , p
    % 11 - 0.24 , 10    , p
    % 12 - 10   , 10    , p
    % 13 - 35   , 0.5   , p
    % 15 - 35   , 0.18  , rf-n
end
temp = -196; % olcumlerin yapildigi sicaklik girilir model bu sicaklikta en az hatayi verecek sekilde sonuc verir

Wdrawn=10e-06; % transistor width
Ldrawn=0.18e-06; % transistor length

short_ch = 0; %%% threshold'a short channel etkisinin katılıp katılmayacağını belirler 1 yapildiginda manual inputta Id_Vg_Vdmin_s_path girilir
Ldrawn_s = 0.18e-06; % short_ch aktif ise kullanilir ve kt1l'nin hesaplanacagi transistorun kanal boyu girilir, Ldrawn'un, Ldrawn_s'den buyuk olmasi gerekir

%rf transistorlerde girilen width length degerlerinin bir onemi yok degisiklikler cadence'da simulasyonlarin yapildigi devre uzerinde yapilmali
transistor_name = 'n_18'; % n_18 / p_18 / n_l18w500_18_rf / p_l18w500_18_rf (p_18 ile n_18 ayni)
type = 'n'; % n / p - transistor_name ne olursa olsun burada n veya p girilir


vds_second_bias = 8; % Id-Vds için VGS degerleri 8 9 10 11 12 (VGS = 1.8 1.5 1.2 0.9 0.6)

step = 10; % algoritma içerisindeki iterasyon sayısı UB1 için kullaniliyor
maksimum_simulasyon_sayisi = 10; % yapilacak maksimum simulasyon limiti
target_min_error = 0.5; % hedeflenen minimum hata degeri, bu degere ulasana kadar parametreler guncelleniyor (0.5%)

% sweeplerde kullanilan max gerilim ocean dosyasindan geliyor. burada
% ikinci kutuplama degerleri yer aliyor, daha fazla detay icin ocean dosyalari incelenebilir
if type == 'n'
    vd = 1.8;
    vg = 1.8;
    vs = 0;
    vb = 0;
elseif type == 'p'
    vd = 0;
    vg = 0;
    vs = 1.8;
    vb = 1.8;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ax1 = subplot(2,2,1);
ax2 = subplot(2,2,2);
ax3 = subplot(2,2,3);
ax4 = subplot(2,2,4);

hold (ax1,'on')
hold (ax2,'on')
hold (ax3,'on')
hold (ax4,'on')
hold all

% girilen transistor_name'e gore hangi ocean dosyalarinin cagirilicagi secilir farkli devreler icin ocean dosyalari da degistirilmeli
[~,transistor_size] = size(transistor_name);
if transistor_size == 4 & (transistor_name =='n_18' | transistor_name =='p_18')
    %-------- N_18,P_18 --------%
    ocean_address_n_vds = './ocean_files/ocean_n_vds.ocn'; % id-vd sweep
    ocean_address_n_vgs = './ocean_files/ocean_n_vgs.ocn'; % id-vg sweep
    ocean_address_n_vgs_low_vds = './ocean_files/ocean_n_vgs_low_vds.ocn';
    ocean_address_p_vds = './ocean_files/ocean_p_vds.ocn';
    ocean_address_p_vgs = './ocean_files/ocean_p_vgs.ocn';
    ocean_address_p_vgs_low_vds = './ocean_files/ocean_p_vgs_low_vds.ocn';
    model_address = '/vlsi/projects/sergiomodeling2/UMC_180nm_02_16/Cadence_IC6_RF/UMC_18_CMOS/../Models/Spectre/mm180_reg18_v124.mdl.scs';
    
    file = fopen(model_address);
    formatSpec = '%s';
    C = textscan(file, formatSpec);
    [boyut, ~] = size(C{:});
    
    % Normal model dosyasi icerisinde hangi bolgenin degistirilecegi bulunur
    for i=1:boyut
        k=strfind(C{1}{i},'type=');
        if k==1
            model_mosfet_type = textscan(C{1}{i},'type=%s');
            model_mosfet_type = model_mosfet_type{1}{1};
            if type == 'n'
                if model_mosfet_type=='n'
                    baslangic = i; %nmos konum
                elseif model_mosfet_type=='p'
                    son = i; %pmos konum
                end
            end
            
            if type == 'p'
                if model_mosfet_type=='p'
                    baslangic = i; %pmos konum
                    son = boyut;
                end
            end
        end
    end
    fclose(file);
    
else % eger transistor_name rf'lerden ise
    %-------- RF --------%
    ocean_address_n_vds = './ocean_files/ocean_n_vds_rf.ocn';
    ocean_address_n_vgs = './ocean_files/ocean_n_vgs_rf.ocn';
    ocean_address_n_vgs_low_vds = './ocean_files/ocean_n_vgs_low_vds_rf.ocn';
    ocean_address_p_vds = './ocean_files/ocean_p_vds_rf.ocn';
    ocean_address_p_vgs = './ocean_files/ocean_p_vgs_rf.ocn';
    ocean_address_p_vgs_low_vds = './ocean_files/ocean_p_vgs_low_vds_rf.ocn';
    model_address = '/vlsi/projects/sergiomodeling2/UMC_180nm_02_16/Cadence_IC6_RF/UMC_18_CMOS/../Models/Spectre/core_rf_v2d4.mdl.scs'; % rf model adresi
    
    file = fopen(model_address);
    formatSpec = '%s';
    C = textscan(file, formatSpec);
    [boyut, ~] = size(C{:});
    
    % RF Transistor model dosyası
    n = 1;
    for i=1:boyut
        k = strfind(C{1}{i}, transistor_name);
        if k==1
            if n==1
                baslangic = i;
                n = 2;
            elseif n == 2
                son = i;
            end
        end
    end
end

% yapilan otomatik simulasyonlarin sonuclarinin kaydedilecegi dizin
result_address_vgs = './sonuclar/n_vgs.dat';
result_address_vgs_low_vds = './sonuclar/n_vgs_low_vds.dat';
result_address_vds = './sonuclar/n_vds.dat';

% iterasyon yapilabilecek parametrelerin max değerleri
ub1_n_max=-10.6730e-18; % vgs'deki bozulmayı düzeltmek için mutlak değeri artar
%ute_n_max=-0.280;    % normal kosullarda kullanilmaz, aktif olmayan iteration fonksiyonuna gonderilen bir parametre
%at_n_max= 50e+03; % normal kosullarda kullanilmaz
delta_n_max = 0.2;% vds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

command = 'sed -i -e ''s/''"%s"''/''"%s"''/g'' %s';
if type =='n'
    ocean_address_vds = ocean_address_n_vds;
    ocean_address_vgs = ocean_address_n_vgs;
    ocean_address_low_vds = ocean_address_n_vgs_low_vds;
    
elseif type =='p'
    ocean_address_vds = ocean_address_p_vds;
    ocean_address_vgs = ocean_address_p_vgs;
    ocean_address_low_vds = ocean_address_p_vgs_low_vds;
end

% grafiklerde soguk olcumlerin yaninda oda sicakliginin da gosterilmesi istendiginde bu kisim kullanilir onun disinda gerekli degil
%%%%%%%%%% 24 derece için Input Import / Ocean File %%%%%%%%%%
%
% ocean_file_edit(24, vg, vd, vs, vb, Wdrawn, Ldrawn, ocean_address_vds);
% ocean_file_edit(24, vg, vd, vs, vb, Wdrawn, Ldrawn, ocean_address_vgs);
% if type == 'n'
% vd = 0.05;
% ocean_file_edit(24, vg, vd, vs, vb, Wdrawn, Ldrawn, ocean_address_low_vds);
% vd = 1.8;
% elseif type == 'p'
% vd = 1.75;
% ocean_file_edit(24, vg, vd, vs, vb, Wdrawn, Ldrawn, ocean_address_low_vds);
% vd = 0;
% end
% %%%%%%%%% 24 derecede default değerlerde analiz %%%%%%%%%%
% drawnow
%
% [ Vds, Ids ] = simulation( ocean_address_vds, result_address_vds );
% [ V, I ] = simulation( ocean_address_vgs, result_address_vgs );
% [ V_low_vds, I_low_vds ] = simulation( ocean_address_low_vds, result_address_vgs_low_vds );
% vth_def = importdata('./sonuclar/vth.dat');
% Ids_cad_def = Ids;
% I_cad_def = I;
% I_cad_low_vds_def = I_low_vds;
%
% if type == 'n'
%     vb = -0.4;
% elseif type == 'p'
%     vb = 2.2;
% end
% ocean_file_edit(24, vg, vd, vs, vb, Wdrawn, Ldrawn, ocean_address_vds);
% [ Vds, Ids_cad_vbs_def ] = simulation( ocean_address_vds, result_address_vds );
% if type=='n'
%     vb = 0;
% elseif type=='p'
%     vb = 1.8;
% end
%
%     subplot(2,2,1)
%         ax1_I_cad_def = plot(V,I_cad_def,'g-');
%         % Create xlabel
%         xlabel('\bf Vgs');
%
%         % Create ylabel
%         ylabel('\bf Id');
%
%         % Create title
%         title('\fontsize{12} \bf Id-Vgs Curve Vbs=0, Vds=1.8V');
%
%     subplot(2,2,3)
%         ax2_I_cad_def = plot(Vds,I_cad_low_vds_def,'g-');
%         % Create xlabel
%         xlabel('\bf Vgs');
%
%         % Create ylabel
%         ylabel('\bf Id');
%
%         % Create title
%         title('\fontsize{12} \bf Id-Vgs Curve Vbs=0, Vds=0.05V');
%
%
%     subplot(2,2,2)
%         ax3_I_cad_def = plot(Vds,Ids_cad_def,'g-');
%         % Create xlabel
%         xlabel('\bf Vds');
%
%         % Create ylabel
%         ylabel('Id');
%
%         % Create title
%         title('\fontsize{12} \bf Id - Vds Curve Vbs=0, Vgs=1.8V');
%
%
%     subplot(2,2,4)
%         ax4_I_cad_def = plot(Vds,Ids_cad_vbs_def,'g-');
%         % Create xlabel
%         xlabel('\bf  Vds');
%
%         % Create ylabel
%         ylabel('\bf Id');
%
%         % Create title
%         title('\fontsize{12} \bf Id - Vds Curve Vbs=-04, Vgs=1.8V');
%     drawnow
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%% Algoritma dongusu buradan baslar %%%%%%%%%
simulasyon_sayisi=1;
level = 1; % algoritma en bastan calistirilmak istenmez ise level 2 ya da 3 olacak secilebilir
% level 1 - vbs=0 icin gerekli parametreleri degistirir
% level 2 - bulk effect icin parametreleri degistirir
% level 3 - vdsat ve delta degerleri icin kullanilir
while level <= 3
    if level == 2
        if type == 'n'
            vb = -0.4; % bulk effect icin kullanılacak vbs degerini belirler (vbs = -0.4)
        elseif type == 'p'
            vb = 2.2; % bulk effect icin kullanılacak vbs degerini belirler (vbs = 0.4)
        end
    elseif level == 3 % level 3'te body voltage eski haline getirilir
        if type == 'n'
            vb = 0;
        elseif type == 'p'
            vb = 1.8;
        end
    end
    
    
    
    %%%%%% Ölçüm Verileri %%%%%%
    
    if manual_input
        %%%%% Ölçüm verisi el ile girilecek ise %%%%%
        % en basta girilen dosya yollarindan olcum sonuclari alinir
        if type == 'n'
            if (level == 1 || level == 3)
                FF = importdata(Id_Vg_Vdmin_path);
                V_aim_low_vds = FF(:,1);
                I_aim_low_vds = FF(:,2);
                FF = importdata(Id_Vg_Vdmax_path);
                V_aim = FF(:,1);
                I_aim = FF(:,2);
                FF = importdata(Id_Vd_Vgmax_path);
                Vds_aim = FF(:,1);
                Ids_aim = FF(:,2);
            elseif level == 2
                FF = importdata(Id_Vg_Vdmin1_path);
                V_aim_low_vds = FF(:,1);
                I_aim_low_vds = FF(:,2);
                I_aim_low_vds_vbs = I_aim_low_vds;
                %         FF = importdata(Id_Vg_Vdmax1_path);
                %         V_aim = FF(:,1);
                %         I_aim = FF(:,2);
                FF = importdata(Id_Vd_Vgmax1_path);
                Vds_aim = FF(:,1);
                Ids_aim = FF(:,2);
                Ids_aim_vbs = Ids_aim;
            end
        elseif type=='p'
            if (level == 1 || level == 3)
                FF = importdata(Id_Vg_Vdmin_path);
                V_aim_low_vds = FF(:,1)+1.8;
                I_aim_low_vds = FF(:,2)*-1;
                FF = importdata(Id_Vg_Vdmax_path);
                V_aim = FF(:,1)+1.8;
                I_aim = FF(:,2)*-1;
                FF = importdata(Id_Vd_Vgmax_path);
                Vds_aim = FF(:,1)+1.8;
                Ids_aim = FF(:,2)*-1;
            elseif level == 2
                FF = importdata(Id_Vg_Vdmin1_path);
                V_aim_low_vds = FF(:,1)+1.8;
                I_aim_low_vds = FF(:,2)*-1;
                I_aim_low_vds_vbs = I_aim_low_vds;
                %         FF = importdata(Id_Vg_Vdmax1_path);
                %         V_aim = FF(:,1)+1.8;
                %         I_aim = FF(:,2)*-1;
                FF = importdata(Id_Vd_Vgmax1_path);
                Vds_aim = FF(:,1)+1.8;
                Ids_aim = FF(:,2)*-1;
                Ids_aim_vbs = Ids_aim;
            end
        end
    end
    
    
    if ~manual_input
        %%%%% Ölçüm verileri otomatik kullanılacağı zaman - Tüm ölçümlerin ortalamasını çıkartarak simülasyonlarda kullanır %%%%%%
        
        if type == 'n'
            V = 0:0.01:1.8;
            if (level == 1 || level == 3)
                [ Id_mean_ln2,  Id_3sigma, Id_e3sigma] = target_meas_data( transistor_no, 'VGS',  2, 8 );
                V_aim_low_vds = V;
                I_aim_low_vds = Id_mean_ln2;
                
                [ Id_mean_ln2,  Id_3sigma, Id_e3sigma] = target_meas_data( transistor_no, 'VGS',  1, 8 );
                V_aim = V;
                I_aim = Id_mean_ln2;
                
                [ Id_mean_ln2,  Id_3sigma, Id_e3sigma] = target_meas_data( transistor_no, 'VDS',  0, vds_second_bias );
                Vds_aim = V;
                Ids_aim = Id_mean_ln2;
            elseif level == 2
                [Id_mean_ln2,  Id_3sigma, Id_e3sigma] = target_meas_data( transistor_no, 'VGS',  2, 12 );
                V_aim_low_vds = V;
                I_aim_low_vds = Id_mean_ln2;
                I_aim_low_vds_vbs = I_aim_low_vds;
                %         FF = importdata(Id_Vg_Vdmax1_path);
                %         V_aim = FF(:,1);
                %         I_aim = FF(:,2);
                [ Id_mean_ln2,  Id_3sigma, Id_e3sigma] = target_meas_data( transistor_no, 'VDS',  4, vds_second_bias );
                Vds_aim = V;
                Ids_aim = Id_mean_ln2;
                Ids_aim_vbs = Ids_aim;
            end
        elseif type=='p'
            V = 1.8:-0.01:0;
            if (level == 1 || level == 3)
                [Id_mean_ln2,  Id_3sigma, Id_e3sigma] = target_meas_data( transistor_no, 'VGS',  2, 8 );
                V_aim_low_vds = V;
                I_aim_low_vds = Id_mean_ln2*-1;
                
                [ Id_mean_ln2,  Id_3sigma, Id_e3sigma] = target_meas_data( transistor_no, 'VGS',  1, 8 );
                V_aim = V;
                I_aim = Id_mean_ln2*-1;
                
                [ Id_mean_ln2,  Id_3sigma, Id_e3sigma] = target_meas_data( transistor_no, 'VDS',  0, vds_second_bias );
                Vds_aim = V;
                Ids_aim = Id_mean_ln2*-1;
            elseif level == 2
                [Id_mean_ln2,  Id_3sigma, Id_e3sigma] = target_meas_data( transistor_no, 'VGS',  2, 12 );
                V_aim_low_vds = V;
                I_aim_low_vds = Id_mean_ln2*-1;
                I_aim_low_vds_vbs = I_aim_low_vds;
                %         FF = importdata(Id_Vg_Vdmax1_path);
                %         V_aim = FF(:,1)+1.8;
                %         I_aim = FF(:,2)*-1;
                [ Id_mean_ln2,  Id_3sigma, Id_e3sigma] = target_meas_data( transistor_no, 'VDS',  4, vds_second_bias );
                Vds_aim = V;
                Ids_aim = Id_mean_ln2*-1;
                Ids_aim_vbs = Ids_aim;
            end
        end
        
    end
    
    
    
    if level == 1  %%olcum verilerinden gelen grafikleri çizmek için
        subplot(2,2,1)
        ax1_I_aim = plot(V_aim,I_aim,'linewidth',2,'color','r');
        %            ax1_I_def = plot(V,I_def,'linewidth',2,'color','r')  ;
        
        subplot(2,2,3)
        ax2_I_aim = plot(V_aim_low_vds,I_aim_low_vds,'linewidth',2,'color','r');
        %           ax2_I_def = plot(V_low_vds,I_low_vds_def,'linewidth',2,'color','r');
        
        subplot(2,2,2)
        ax3_I_aim = plot(Vds_aim,Ids_aim,'linewidth',2,'color','r');
        %           ax3_I_def = plot(Vds,Ids_def,'linewidth',2,'color','r');
        
    elseif level == 2
        subplot(2,2,4)
        ax4_I_aim = plot(Vds_aim,Ids_aim_vbs,'linewidth',2,'color','r');
        %       ax4_I_def = plot(Vds,Ids_def_vbs,'linewidth',2,'color','r');
    end
    drawnow
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%% Input Import / Ocean File %%%%%%%%%%
    % ocean dosyalari uzerinde gerilim degerleri ve transistor boyutlari degistirilir
    ocean_file_edit(temp, vg, vd, vs, vb, Wdrawn, Ldrawn, ocean_address_vds);
    ocean_file_edit(temp, vg, vd, vs, vb, Wdrawn, Ldrawn, ocean_address_vgs);
    if type == 'n'
        vd = 0.05;
        ocean_file_edit(temp, vg, vd, vs, vb, Wdrawn, Ldrawn, ocean_address_low_vds);
        vd = 1.8;
    elseif type == 'p'
        vd = 1.75;
        ocean_file_edit(temp, vg, vd, vs, vb, Wdrawn, Ldrawn, ocean_address_low_vds);
        vd = 0;
    end
    %%%%%%%%% default değerlerde analiz %%%%%%%%%%
    
    drawnow
    
    [ Vds, Ids ] = simulation( ocean_address_vds, result_address_vds ); % cadence'ta simulasyon kosturulur sonuclari alinir (ID-VDS)
    [ V, I ] = simulation( ocean_address_vgs, result_address_vgs ); % (ID-VGS)
    [ V_low_vds, I_low_vds ] = simulation( ocean_address_low_vds, result_address_vgs_low_vds ); % (ID-VGS @ low VDS)
    vth = importdata('./sonuclar/vth.dat'); % cadence'in hesapladigi vth
    vdsat = importdata('./sonuclar/vdsat.dat'); % cadence'in hesapladigi vdsat
    if level ==1
        Ids_cad_ln2 = Ids;
        I_cad_ln2 = I;
        
        % plot(V,Ids) % dongu icerisinde gerceklesen her simulasyon sonucunun grafikte gorulmesi istenirse buradaki yorum kaldirilir
        % plot(V,Ids_aim)
        
        subplot(2,2,1)
        ax1_I_cad_ln2 = plot(V,I_cad_ln2,'-');
        
        subplot(2,2,2)
        ax3_I_cad_ln2 = plot(Vds,Ids_cad_ln2,'-');
        
        drawnow
    end
    
    
    %%%%%%%%%%%%%%%%%%%% VTH %%%%%%%%%%%%%%%%%%%%%%
    % olcum sonuclarinin grafikleri uzerinden vth elde edilir
    I_aim_low_vds = smooth(I_aim_low_vds);
    
    if type =='n'
        [~,loc]=max(diff(I_aim_low_vds));
    elseif type == 'p'
        [~,loc]=max(abs(diff(I_aim_low_vds)));
    end
    
    if loc>=3
        m= (I_aim_low_vds(loc)-I_aim_low_vds(loc-2))/(V_aim_low_vds(loc)-V_aim_low_vds(loc-2));
    else
        m= (I_aim_low_vds(loc)-I_aim_low_vds(loc+2))/(V_aim_low_vds(loc)-V_aim_low_vds(loc+2));
    end
    b = I_aim_low_vds(loc) - m*V_aim_low_vds(loc);
    
    vth_graph = -b/m; % lineer metot ile vth hesaplama
    
    if ~manual_input
        [~,loc]=max(smooth(transpose(V_aim_low_vds(1,1:size(V_aim_low_vds,2)-2)),diff(I_aim_low_vds,2)));
    else
        [~,loc]=max(smooth(transpose(V_aim_low_vds(1:size(V_aim_low_vds,1)-2)),diff(I_aim_low_vds,2)));
    end
    
    vth_sd = V_aim_low_vds(loc); % ikinci derece turevden vth hesaplama
    if type == 'p'
        vth_sd = 1.8-vth_sd;
        vth_graph = 1.8-vth_graph;
    end
    
    % hangi yontemin bundan sonra vth olarak kullanilcagina karar verilir (burasi direkt olarak vth_app = vth_graph olarak kullanilabilir boylece hep ayni yontem kullanilmis olur)
    vth_app = vth_sd;
    if abs(abs(vth_graph)-abs(vth_sd))>0.2
        vth_app = vth_graph;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    error_vds=100; % her dongude hesaplanan hata sifirlanir
    [Ids_boyut,~] = size(Ids);
    
    while error_vds>target_min_error || (Ids_aim(Ids_boyut,1)-Ids(Ids_boyut,1))/Ids(Ids_boyut,1)>0.01 % hedeflenen hata oranina ulasana kadar ya da maksimum simulasyon sayisina erisene kadar algoritma calismaya devam eder
    % hedef hata orani degistirilebilir    
        %%%%%%%%%%%%%%%%%%% VDSAT %%%%%%%%%%%%%%%%%%%%
        
        [ vdsat_difference ] = vdsat_finding( Vds, Ids, Vds, Ids_aim );
              
        %%%%%%%%%%%%% Hesaplamalar %%%%%%%%%%%%%%%%%%
        % BSIM esitlikleri kullanilarak parametreler hesaplanir
        if type =='n'
            [level, hata, Ids_app] = parametric_calculations_n ( temp, vg, vd, vs, vb, vth, vth_app, Wdrawn, Ldrawn, model_address, command, Ids_aim, Ids, I, I_aim, vdsat, vdsat_difference, ocean_address_vgs, ocean_address_low_vds, ocean_address_vds, result_address_vgs, result_address_vgs_low_vds, result_address_vds, ub1_n_max, step, level, baslangic, son, short_ch, type );
            
        elseif type =='p'
            [level, hata, Ids_app] = parametric_calculations_p ( temp, vg, vd, vs, vb, vth, vth_app, Wdrawn, Ldrawn, model_address, command, Ids_aim, Ids, I, I_aim, vdsat, vdsat_difference, ocean_address_vgs,  ocean_address_low_vds, ocean_address_vds, result_address_vgs, result_address_vgs_low_vds, result_address_vds, ub1_n_max, step, level, baslangic, son, short_ch, type );
        end
        
        % BSIM esitliklerini MATLAB cozemezse hata doner, genellikle akim denkleminden UTE bulunurken hata olusuyor; bu yüzden ute icin ayri bir hata parametresi var
        if hata == 1
            if level == 1
                hata_ute = 1;
            end
            %        break % hata olustugun programin calismayi durdurmasi isteniyorsa 
        end
        
        
        
        % kisa kanal etkisi ayrica dahil edilmek istendiginde level 1'e ek bir level atanir ve KT1L parametresi hesaplanir
        if short_ch == 1 && level == 1 && simulasyon_sayisi > 2
            level = 4;
            
            if type =='n'
                FF = importdata(Id_Vg_Vdmin_s_path);
                V_aim_low_vds_s = FF(:,1);
                I_aim_low_vds_s = FF(:,2);
                
            elseif type == 'p'
                FF = importdata(Id_Vg_Vdmin_s_path);
                V_aim_low_vds_s = FF(:,1)+1.8;
                I_aim_low_vds_s = FF(:,2)*-1;
            end          
            
            %%%%%%%%%%%%%%%%%%%% VTH %%%%%%%%%%%%%%%%%%%%%%
            I_aim_low_vds_s = smooth(I_aim_low_vds_s);
            
            if type =='n'
                [~,loc]=max(diff(I_aim_low_vds_s));
            elseif type == 'p'
                [~,loc]=max(abs(diff(I_aim_low_vds_s)));
            end
            
            if loc>=3
                m= (I_aim_low_vds_s(loc)-I_aim_low_vds_s(loc-2))/(V_aim_low_vds_s(loc)-V_aim_low_vds_s(loc-2));
            else
                m= (I_aim_low_vds_s(loc)-I_aim_low_vds_s(loc+2))/(V_aim_low_vds_s(loc)-V_aim_low_vds_s(loc+2));
            end
            b = I_aim_low_vds_s(loc) - m*V_aim_low_vds_s(loc);
            
            vth_graph = -b/m;
            
            [~,loc]=max(smooth(V_aim_low_vds_s(1:size(V_aim_low_vds_s)-2,1),diff(I_aim_low_vds_s,2)));
            vth_sd = V_aim_low_vds_s(loc);
            if type == 'p'
                vth_sd = 1.8-vth_sd;
                vth_graph = 1.8-vth_graph;
            end
            
            vth_app = vth_sd;
            if abs(abs(vth_graph)-abs(vth_sd))>0.2
                vth_app = vth_graph;
            end
            
            Ldrawn_old = Ldrawn;
            Ldrawn = Ldrawn_s;
            
            %%%%%%%%%%% Input Import / Ocean File %%%%%%%%%%
            
            ocean_file_edit(temp, vg, vd, vs, vb, Wdrawn, Ldrawn, ocean_address_vgs);
            if type == 'n'
                vd = 0.05;
                ocean_file_edit(temp, vg, vd, vs, vb, Wdrawn, Ldrawn, ocean_address_low_vds);
                vd = 1.8;
            elseif type == 'p'
                vd = 1.75;
                ocean_file_edit(temp, vg, vd, vs, vb, Wdrawn, Ldrawn, ocean_address_low_vds);
                vd = 0;
            end
            
            %%%%%%% analiz %%%%%
            
            [ V_low_vds, I_low_vds ] = simulation( ocean_address_low_vds, result_address_vgs_low_vds );
            [ V, I ] = simulation( ocean_address_vgs, result_address_vgs );
            A = importdata('./sonuclar/vth.dat');
            vth = A;
            
            
            %%%%%%%%%%%%% Hesaplamalar %%%%%%%%%%%%%%%%%%
            
            if type =='n'
                [level, hata, Ids_app] = parametric_calculations_n ( temp, vg, vd, vs, vb, vth, vth_app, Wdrawn, Ldrawn, model_address, command, Ids_aim, Ids, I, vdsat, vdsat_difference, ocean_address_vgs, ocean_address_low_vds, ocean_address_vds, result_address_vgs, result_address_vgs_low_vds, result_address_vds, ub1_n_max, step, level, baslangic, son, short_ch );
                
            elseif type =='p'
                [level, hata, Ids_app] = parametric_calculations_p ( temp, vg, vd, vs, vb, vth, vth_app, Wdrawn, Ldrawn, model_address, command, Ids_aim, Ids, I, vdsat, vdsat_difference, ocean_address_vgs,  ocean_address_low_vds, ocean_address_vds, result_address_vgs, result_address_vgs_low_vds, result_address_vds, ub1_n_max, step, level, baslangic, son, short_ch );
            end
            
            
            Ldrawn = Ldrawn_old;
            level = 1; % kisa kanal icin KT1L bulunduktan sonra program Ldrawn degeri icin level=1'den normal calismasina devam eder
        end
        
       
        %%%%%%% bulunun yeni parametreler ile analiz %%%%%%%%
        [ Vds, Ids ] = simulation( ocean_address_vds, result_address_vds );
        [ V_low_vds, I_low_vds ] = simulation( ocean_address_low_vds, result_address_vgs_low_vds );
        [ V, I] = simulation(ocean_address_vgs,result_address_vgs);
        A = importdata('./sonuclar/vth.dat');
        vth = A;
        A = importdata('./sonuclar/vdsat.dat');
        vdsat = A;
        if level == 2
            Ids_vbs = Ids;
            I_low_vds_vbs = I_low_vds;
            I_vbs = I;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % olcum egrisi ile yeni parametrelerin eklendigi cadence simulasyonu grafikleri arasindaki hata orani hesaplanir
        %%%%%% ERROR %%%%%%
        %         error=0;
        %     for i=1:size(I)
        % error = abs(I_aim_low_vds(i)-I_low_vds(i))+error;
        %     end
        %         toplam=0;
        %     for i=1:size(I)
        % toplam = abs(I_aim_low_vds(i))+toplam;
        %     end
        %
        % error_vgs = error/toplam*100 % maksimum error
        %
        error=0;
        for i=1:size(Ids)
            error = abs(Ids_aim(i)-Ids(i))+error;
        end
        toplam=0;
        for i=1:size(Ids)
            toplam = abs(Ids_aim(i))+toplam;
        end
        
        error_vds_new = error/toplam*100 % maksimum error
        
        if abs(error_vds_new-error_vds) < 0.1
            break
        else
            error_vds = error_vds_new
        end
        %%%%%%%%%%%%%%%%%%%
        
        % maksimum simulasyon sayisi asildiginda program durdurulur
        simulasyon_sayisi=simulasyon_sayisi+1
        if level == 1 && simulasyon_sayisi==maksimum_simulasyon_sayisi
            break
        end
        
        if level == 2 && simulasyon_sayisi==2*maksimum_simulasyon_sayisi
            break
        end
        
        if level == 3 && simulasyon_sayisi==2*maksimum_simulasyon_sayisi + 2
            break
        end
        
        if hata == 3
            break
        end
        
        
    end
    
    %%%%%%%%%%% İterasyonlar %%%%%%%%%%%%%%%%%
    if level == 3
        % sadece delta parametresi icin iterasyon yapilir ancak diger parametrelerde de iterasyon yapilmasi gerekirse iteration_algorithm fonksiyonu aktif hale getirilir. bu eski bir fonksiyon ve hata olusabilir
        %[ kt1_n, ute_n, at_n, delta_n, ua1_n ] = iteration_algorithm(command, model_address,step, ute_n_max, I, I_aim, Ids, Ids_aim, ocean_address_vds, ocean_address_vgs, result_address_vds, result_address_vgs );
        [ delta_n, Vds, Ids ] = delta_iteration( model_address, baslangic, son, command, step, Ids, Vds, Ids_aim, Vds_aim, ocean_address_vds, result_address_vds, delta_n_max, error_vds, vdsat_difference );
    end
    
    level = level + 1
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% algoritma dongusu sonu, tum parametreler hesaplandi ve guncellendi

%%%%%%%%%%% Input Import / Ocean File %%%%%%%%%%

ocean_file_edit(temp, vg, vd, vs, vb, Wdrawn, Ldrawn, ocean_address_vds);
ocean_file_edit(temp, vg, vd, vs, vb, Wdrawn, Ldrawn, ocean_address_vgs);
if type == 'n'
    vd = 0.05;
    ocean_file_edit(temp, vg, vd, vs, vb, Wdrawn, Ldrawn, ocean_address_low_vds);
    vd = 1.8;
elseif type == 'p'
    vd = 1.75;
    ocean_file_edit(temp, vg, vd, vs, vb, Wdrawn, Ldrawn, ocean_address_low_vds);
    vd = 0;
end

%%%%%%% Belirlenen parametrelerle son analiz %%%%%

[ Vds, Ids ] = simulation( ocean_address_vds, result_address_vds );
[ V_low_vds, I_low_vds ] = simulation( ocean_address_low_vds, result_address_vgs_low_vds );
[ V, I ] = simulation( ocean_address_vgs, result_address_vgs );
A = importdata('./sonuclar/vth.dat');
vth = A;
% A = importdata('./sonuclar/vdsat.dat');
% vdsat = A;

if type =='n'
    vb = -0.4;
elseif type =='p'
    vb = 2.2;
end
ocean_file_edit(temp, vg, vd, vs, vb, Wdrawn, Ldrawn, ocean_address_vds);
ocean_file_edit(temp, vg, vd, vs, vb, Wdrawn, Ldrawn, ocean_address_vgs);
[ Vds, Ids_vbs ] = simulation( ocean_address_vds, result_address_vds );
[ V, I_vbs ] = simulation( ocean_address_vgs, result_address_vgs );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parametrelerin dosyadan tekrar okunması (dosyadaki son değerleri görebilmek için)

file = fopen(model_address);
formatSpec = '%s';
C = textscan(file,formatSpec);
[boyut, ~] = size(C{:});
for i=baslangic:son
    
    k=strfind(C{1}{i},'kt1=');
    if(k==1)
        kt1_n = textscan(C{1}{i},'kt1=%f');
        kt1_n = kt1_n{1};
    end
    
    k=strfind(C{1}{i},'kt1l=');
    if(k==1)
        kt1l_n = textscan(C{1}{i},'kt1l=%f');
        kt1l_n = kt1l_n{1};
    end
    
    k=strfind(C{1}{i},'kt2=');
    if(k==1)
        kt2_n = textscan(C{1}{i},'kt2=%f');
        kt2_n = kt2_n{1};
    end
    
    k=strfind(C{1}{i},'ute=');
    if(k==1)
        ute_n = textscan(C{1}{i},'ute=%f');
        ute_n = ute_n{1};
    end
    
    k=strfind(C{1}{i},'at=');
    if(k==1)
        at_n = textscan(C{1}{i},'at=%f');
        at_n = at_n{1};
    end
    
    k=strfind(C{1}{i},'ua1=');
    if(k==1)
        ua1_n = textscan(C{1}{i},'ua1=%f');
        ua1_n = ua1_n{1};
    end
    
    k=strfind(C{1}{i},'ub1=');
    if(k==1)
        ub1_n = textscan(C{1}{i},'ub1=%f');
        ub1_n = ub1_n{1};
    end
    
    k=strfind(C{1}{i},'uc1=');
    if(k==1)
        uc1_n = textscan(C{1}{i},'uc1=%f');
        uc1_n = uc1_n{1};
    end
    
    k=strfind(C{1}{i},'delta=');
    if(k==1)
        delta_n = textscan(C{1}{i},'delta=%f');
        delta_n = delta_n{1};
    end
end
fclose(file);

if exist('delta_n','var')== 0
    delta_n = 0;
end


subplot(2,2,1)
ax1_I = plot(V,I,'--+');

subplot(2,2,3)
ax2_I = plot(V_low_vds,I_low_vds,'--+');

subplot(2,2,2)
ax3_I = plot(Vds,Ids,'--+');

subplot(2,2,4)
ax4_I = plot(Vds,Ids_vbs,'--+');


%%%%%% ERROR %%%%%%
error=0;
for i=1:size(I)
    error = abs(I_aim(i)-I(i))^2+error;
end
toplam=0;
for i=1:size(I)
    toplam = abs(I_aim(i))^2+toplam;
end
% ID-VGS error
error_vgs = sqrt(error/toplam)*100; % maksimum error 

error=0;
for i=1:size(Ids)
    error = abs(Ids_aim(i)-Ids(i))^2+error;
end
toplam=0;
for i=1:size(Ids)
    toplam = abs(Ids_aim(i))^2+toplam;
end
% ID-VDS error
error_vds = sqrt(error/toplam)*100; % maksimum error


error=0;
for i=1:size(Ids_vbs)
    error = abs(Ids_aim_vbs(i)-Ids_vbs(i))^2+error;
end
toplam=0;
for i=1:size(Ids_vbs)
    toplam = abs(Ids_aim_vbs(i))^2+toplam;
end
% ID-VDS @ |VBS|>0
error_vds_vbs = sqrt(error/toplam)*100; % maksimum error


error=0;
for i=1:size(I)
    error = abs(I_aim(i)-I_cad_ln2(i))^2+error;
end
toplam=0;
for i=1:size(I)
    toplam = abs(I_aim(i))^2+toplam;
end
% ID-VGS error (cadence'in default parametreler ile dusuk sicaklikta)
error_vgs_cad = sqrt(error/toplam)*100; % maksimum error

error=0;
for i=1:size(Ids)
    error = abs(Ids_aim(i)-Ids_cad_ln2(i))^2+error;
end
toplam=0;
for i=1:size(Ids)
    toplam = abs(Ids_aim(i))^2+toplam;
end
% ID-VDS error (cadence'in default parametreler ile dusuk sicaklikta)
error_vds_cad = sqrt(error/toplam)*100; % maksimum error

toc
%%%%%%%%%%%%%%%%%%%
kt1_n
kt1l_n
kt2_n
ute_n
ua1_n
ub1_n
uc1_n
at_n
delta_n

%vth_def
vth_app
vth

error_vgs_cad
error_vgs
error_vds_cad
error_vds
error_vds_vbs

