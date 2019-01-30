function [ vdsat_difference ] = vdsat_finding( Vds, Ids, Vds_aim, Ids_aim )
% minimum ve maksimum turev degerlerine sahip noktalardan dogrusal egriler cizdirilir ve bu egrilerin kesisim noktasi hesaplanir bu degerin X ekseni Vdsat olarak alinir    
Vds_draw = Vds;
Ids_draw = Ids;
    
[~,loc1] = max(diff(Ids_draw));
if loc1<178
m1 = (Ids_draw(loc1) - Ids_draw(loc1+2))/(Vds_draw(loc1)-Vds_draw(loc1+2));
else
m1 = (Ids_draw(loc1) - Ids_draw(loc1-2))/(Vds_draw(loc1)-Vds_draw(loc1-2));
end
b1 = Ids_draw(loc1) - m1*Vds_draw(loc1);

x1 = 0: 0.01: 1.8;
y1 = m1*x1 + b1;
% subplot(2,2,2)
% plot(x1,y1)

[~,loc2] = min(diff(Ids_draw));
if loc2>=3
m2 = (Ids_draw(loc2) - Ids_draw(loc2-2))/(Vds_draw(loc2)-Vds_draw(loc2-2));
else
m2 = (Ids_draw(loc2) - Ids_draw(loc2+2))/(Vds_draw(loc2)-Vds_draw(loc2+2));  
end
b2 = Ids_draw(loc2) - m2*Vds_draw(loc2);

x2 = 0: 0.01: 1.8;
y2 = m2*x2 + b2;
% subplot(2,2,2)
% plot(x2,y2)

[~,loc3]=min(abs((y1-y2)));

Vdsat_def = Vds_draw(loc3);

   
Vds_draw = Vds_aim;
Ids_draw = Ids_aim;
    
[~,loc1] = max(diff(Ids_draw));
if loc1>=3
m1 = (Ids_draw(loc1) - Ids_draw(loc1-2))/(Vds_draw(loc1)-Vds_draw(loc1-2));
else
m1 = (Ids_draw(loc1) - Ids_draw(loc1+2))/(Vds_draw(loc1)-Vds_draw(loc1+2));    
end
b1 = Ids_draw(loc1) - m1*Vds_draw(loc1);

x1 = 0: 0.01: 1.8;
y1 = m1*x1 + b1;
% subplot(2,2,2)
% plot(x1,y1)


[~,loc2] = min(diff(Ids_draw));
if loc2>=3
m2 = (Ids_draw(loc2) - Ids_draw(loc2-2))/(Vds_draw(loc2)-Vds_draw(loc2-2));
else
m2 = (Ids_draw(loc2) - Ids_draw(loc2+2))/(Vds_draw(loc2)-Vds_draw(loc2+2));   
end
b2 = Ids_draw(loc2) - m2*Vds_draw(loc2);

x2 = 0: 0.01: 1.8;
y2 = m2*x2 + b2;
% subplot(2,2,2)
% plot(x2,y2)

[~,loc3]=min(abs((y1-y2)));

Vdsat_aim = Vds_draw(loc3);

vdsat_difference = Vdsat_aim - Vdsat_def;
    
end

