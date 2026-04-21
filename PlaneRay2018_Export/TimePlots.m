% TimePlots
clear

load trans;load para; load env;
load signal;load range; load t_real
load para;
scale=50000; threshold=0;
figure(88);clf;

M_range=length(para.eigrange);
for m=1:M_range;
    range=para.eigrange(m)
    newpara.eigrange=range;
    
    [ eigresult] = validate( newpara, env, SUFBOT,eigen, threshold);
    range=eigresult.range;
    
    plot(t_real(m,:),range+scale*signal(m,:))
    
    hold on
end