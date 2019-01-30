 function ocean_file_edit(temp, vg, vd, vs, vb, Wdrawn, Ldrawn, ocean_address)
% ocean dosyasi icerisindeki parametreleri degistirmek icin kullanilir
% shell komutlari ile text olarak degisiklik yapilir
     
command = 'sed -i -e ''s/%s/%s/g'' %s';
    
    file = fopen(ocean_address);
formatSpec = '%s';
C = textscan(file,formatSpec);
[boyut,~] = size(C{:});

for i=1:boyut
    
     k=strfind(C{1}{i},'temp(');
     if(k==1)
        temp_old = textscan(C{1}{i},'temp(%d)');
        temp_old = temp_old{1};
        temp_old_char = C{1}{i};
     end
     
     k=strfind(C{1}{i},'desVar("vg",');
     if(k==1)
        vg_old = textscan(C{1}{i},'desVar("vg",%d)');
        vg_old = vg_old{1};
        vg_old_char = C{1}{i};
     end
     
     k=strfind(C{1}{i},'desVar("vd",');
     if(k==1)
        vd_old = textscan(C{1}{i},'desVar("vd",%d)');
        vd_old = vd_old{1};
        vd_old_char = C{1}{i};
     end
     
     k=strfind(C{1}{i},'desVar("vs",');
     if(k==1)
        vs_old = textscan(C{1}{i},'desVar("vs",%d)');
        vs_old = vs_old{1};
        vs_old_char = C{1}{i};
     end
     
     k=strfind(C{1}{i},'desVar("vb",');
     if(k==1)
        vb_old = textscan(C{1}{i},'desVar("vb",%d)');
        vb_old = vb_old{1};
        vb_old_char = C{1}{i};
     end
     
     k=strfind(C{1}{i},'desVar("width",');
     if(k==1)
        width_old = textscan(C{1}{i},'desVar("width",%d)');
        width_old = width_old{1};
        width_old_char = C{1}{i};
     end
     
     k=strfind(C{1}{i},'desVar("length",');
     if(k==1)
        length_old = textscan(C{1}{i},'desVar("length",%d)');
        length_old = length_old{1};
        length_old_char = C{1}{i};
     end
     
     k=strfind(C{1}{i},'?start');
     if(k==1)
        start_old = textscan(C{1}{i+1},'"%d"');
        start_old = start_old{1};
        start_old_char = C{1}{i+1};
     end
     
     k=strfind(C{1}{i},'?stop');
     if(k==1)
        stop_old = textscan(C{1}{i+1},'"%.2f"');
        stop_old = stop_old{1};
        stop_old_char = C{1}{i+1};
     end
     
     k=strfind(C{1}{i},'?step');
     if(k==1)
        step_old = textscan(C{1}{i+1},'"%.2f"');
        step_old = step_old{1};
        step_old_char = C{1}{i+1};
        step_old_char = sprintf('?step "%.2f"',step_old_char);
     end
end

fclose(file);

temp_char = sprintf('temp(%d)',temp);
execute = sprintf(command,temp_old_char,temp_char,ocean_address);
unix(execute);

vg_char = sprintf('desVar("vg",%d)',vg);
execute = sprintf(command,vg_old_char,vg_char,ocean_address);
unix(execute);

vd_char = sprintf('desVar("vd",%d)',vd);
execute = sprintf(command,vd_old_char,vd_char,ocean_address);
unix(execute);

vs_char = sprintf('desVar("vs",%d)',vs);
execute = sprintf(command,vs_old_char,vs_char,ocean_address);
unix(execute);

vb_char = sprintf('desVar("vb",%d)',vb);
execute = sprintf(command,vb_old_char,vb_char,ocean_address);
unix(execute);
% 
% step_char = sprintf('?step "%.2f"',step);
% execute = sprintf(command,step_old_char,step_char,ocean_address);
% unix(execute);

if exist('width_old','var')
width_char = sprintf('desVar("width",%d)',Wdrawn);
execute = sprintf(command,width_old_char,width_char,ocean_address);
unix(execute);

length_char = sprintf('desVar("length",%d)',Ldrawn);
execute = sprintf(command,length_old_char,length_char,ocean_address);
unix(execute);
end

end

