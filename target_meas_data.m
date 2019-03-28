function [ Id_mean_ln2,  Id_3sigma, Id_e3sigma] = target_meas_data( transistor_no, sweep,  first_bias, second_bias )
% asagidaki formatta yazilmis .mat dosyasi icerisindeki tum olcumleri alarak istenen boyuttaki transistorun tum olcum sonuclari ya da tum sonuclarin ortalamalai, +-3sigma degerleri grafik olarak cizdirilebilir
% ID_cell{1,1} = 'wafer';
% ID_cell{2,1} = 'chip';
% ID_cell{3,1} = 'transistor';
% ID_cell{4,1} = 'ln2';
% ID_cell{5,1} = 'sweep';
% ID_cell{6,1} = 'bias';
% ID_cell{7,1} = 'V';
% ID_cell{8,1} = 'ID';

    % format icin bu dosya incelenebilir bu dosya da yine otomatik olarak result_import programi sayesinde olusturulur
    load('id_cell_cleared.mat') % id_cell_cleared yapilmis olan tum duzgun olcum sonuclarini barindirir

% grafigi cizdirilmesi istenen transistor(ler), wafer(lar), chip(ler) array icerisine dogru isimlendirme ile yazildiginda grafikleri ekrana bastirilir
wafer = [3 5 6]; % 0 3 5 6 wafer numbers
chip = [1:11]; % 1:10 chip numbers
transistor = [transistor_no]; % 1:22 transistor numbers
ln2 = [1];  % 0 (24C), 1 (under LN2) 
bias = [first_bias]; % For VDS, VBS: 0 1 2 3 4 (0 -0.1 -0.2 -0.3 -0.4) / For VGS, VDS: 1 2 3 (1.8 0.05 0.01)

if sweep == 'VDS'
    curve = second_bias;  % 8 9 10 11 12 (VGS = 1.8 1.5 1.2 0.9 0.6)
elseif sweep == 'VGS'
    curve = second_bias;  % 8 9 10 11 12 (VBS = 0 -0.1 -0.2 -0.3 -0.4)
end

hold all
counter = 0;
counter_ln2 = 0;
toplam = 0;
toplam_ln2 = 0;
toplam_low = 0;
toplam_ln2_low = 0;
for i=2:size(ID_cell,2)
    for n=1:size(wafer,2)
        if wafer(n) == ID_cell{1,i}
            for m=1:size(chip,2)
                if chip(m) == ID_cell{2,i}
                    for l=1:size(transistor,2)
                        if transistor(l) == ID_cell{3,i}
                            for k=1:size(ln2,2)
                                if ln2(k) == ID_cell{4,i}
                                    if sweep == char(ID_cell{5,i})
                                        for j=1:size(bias,2)
                                            if bias(j) == ID_cell{6,i}
                                                if ln2(k) == 0
                                                    toplam = toplam + ID_cell{curve,i};
                                                    counter = counter + 1;
                                                elseif ln2(k) == 1
                                                    if ~isempty(ID_cell{curve,i})
                                                    toplam_ln2 = toplam_ln2 + ID_cell{curve,i};
                                                    counter_ln2 = counter_ln2 + 1;
                                                    end
                                                end
                                                name = sprintf('W%dC%dT%d', wafer(n),chip(m),transistor(l));
                                                name_array{1,counter+counter_ln2} = name;
                                                name_array{2,counter+counter_ln2} = ID_cell{curve,i};
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end


Id_mean = toplam / counter;
Id_mean_ln2 = toplam_ln2 / counter_ln2;

cell_boyut = size(name_array);
cell_boyut = cell_boyut(2);
for n=1:cell_boyut
    if ~isempty(name_array{2,n})
       sigma_array(:,n) = name_array{2,n};  
    end
end
sigma_array = transpose(sigma_array);
sigma = std(sigma_array);

Id_3sigma = Id_mean_ln2 + 3*transpose(sigma);
Id_e3sigma = Id_mean_ln2 - 3*transpose(sigma);

clearvars sigma_array sigma name_array ID_cell

end

