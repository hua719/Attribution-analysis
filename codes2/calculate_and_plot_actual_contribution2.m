%
clear;
nameset={'Jimai','Tuotuohe','Zhimenda','Shigu','Xiangda','Daojieba','Nuxia','Shigatse','Lasha','Gengzhang','Lazi','Jiayuqiao'};
nameset1={'Nuxia','Tuotuohe','Jiayuqiao','Shigatse','Lasha','Gengzhang','Lazi'};
nameset2={'Lazi','Shigatse','Lasha','Gengzhang','Nuxia','Jiayuqiao','Tuotuohe'};
namesetx={'S1','S2','S3','S4','S5','S6','S7'};
for i=1:length(nameset1)
    for j=1:length(nameset)
a=strcmp(nameset1{i},nameset{j});
if a==1
    cindex(i,1)=j;
end
    end
end

for i=1:length(nameset2)
    for j=1:length(nameset1)
a=strcmp(nameset2{i},nameset1{j});
if a==1
    cindex2(i,1)=j;
end
    end
end
%% read sensitivity
sensitivity_set=textread('G:\variation_and_attribution_for_runoff_and_sediment_flux\codes2\input_data\ridge_sensitivity_plus_ndvi.txt');
varset={'pre','ta','growing_season_noaa_cdr_ndvi','gleam_es','q'};
startyear=1950;
endyear=2022;
turning_point=1995;
xset1=ones(endyear-startyear+1,length(varset))*NaN;
num_vars=length(varset);
num_stations=length(nameset1);

%% read data and calculate relative changes in climate variables and runoff
for s=1:num_stations
for i=1:length(varset)
if i<3
    xset=textread(['G:\variation_and_attribution_for_runoff_and_sediment_flux\codes2\input_data\obs_annual_' varset{i} '.txt']);
   xset1((1980-1950+1):(2020-1950+1),i)=xset(:,s);
elseif i==3
        xset=textread(['G:\variation_and_attribution_for_runoff_and_sediment_flux\codes2\input_data\annual_' varset{i} '.txt']);
           xset1((1982-1950+1):(2020-1950+1),i)=xset(:,s);

elseif i==4
        xset=textread(['G:\variation_and_attribution_for_runoff_and_sediment_flux\codes2\input_data\annual_' varset{i} '.txt']);
   xset1((1980-1950+1):(2020-1950+1),i)=xset(:,s);
else
    
    
    xset=textread(['G:\variation_and_attribution_for_runoff_and_sediment_flux\codes2\input_data\annual_' varset{i} '.txt']);
   xset1(:,i)=xset(:,cindex(s)); 
end   

end

index=find(~isnan(xset1(:,end)) & ~isnan(xset1(:,3)));
xset2=xset1(index,:);
start_endset(s,1)=index(1)+startyear-1;
start_endset(s,2)=index(end)+startyear-1;
x=(index(1)+startyear-1):1:(index(end)+startyear-1);
xx=[ones(length(x),1),x'];

index1=turning_point-x(1)+1;
for i=1:length(varset)
    ref=xset2(1:index1,i);
    var=xset2(index1+1:end,i);
    ave_ref_set(s,i)=mean(ref);
    if i==2
        change_set(s,i)=mean(var)-mean(ref);   % compared to prior-1995 period
        
    else
   change_set(s,i)=(mean(var)-mean(ref))/mean(ref);
   
    end
end

xset2_all{s,1}=xset2;
end
sensitivity_set1=sensitivity_set;
sensitivity_set2=sensitivity_set1';
change_set1=change_set(:,[1,3,4,5]);
change_set2=change_set1*100;

%% calculate contribution
for i=1:num_vars-2
    for j=1:size(sensitivity_set2,1)
        if i==2 || i==3
            a0=sensitivity_set2(j,i)*change_set1(j,i)*100;
            if a0>0
     actual_con(j,i)=-1*a0;
            else
               actual_con(j,i)=a0; 
            end
        else
        actual_con(j,i)=sensitivity_set2(j,i)*change_set1(j,i)*100;
        end
    end
end
other_con=change_set2(:,end)-sum(actual_con,2);
actual_con_set=[actual_con,other_con];
actual_con_set1=actual_con_set./abs(change_set2(:,end))*100;
%% actual contribution in km3
ave_ref_setx=ave_ref_set(:,5)';
for s=1:7
actual_con_setx(s,:)=actual_con_set(s,:)*ave_ref_setx(s)/100;
end
actual_con_setx2=actual_con_setx(cindex2,:);

%% adjust the order of stations in the grah.
actual_con_set2=actual_con_set1(cindex2,:);
change_set3=change_set2(cindex2,:);
index_x_y=find(abs(actual_con_set(:,1))>abs(actual_con_set(:,2)));
[~,index_x_y1]=max(actual_con_set,[],2);
index_x_y2=find(actual_con_set(:,2)>actual_con_set(:,3));
%% relative contribution
actual_con_set3=abs(actual_con_set2);
sum_con=sum(actual_con_set3,2);
for s=1:num_stations
relative_con_set3(s,:)=actual_con_set3(s,:)./sum_con(s,1)*100;
end

%% set color
   str = '#a6cee3';
color_pre = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

str = '#b2df8a';
color_ndvi = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

str = '#fdbf6f';
color_et = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

str = '#1f78b4';
color_q = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

str = '#fb9a99';
color_others = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
%% plot relative changes and contributions
rows=1;
cols=2;
nplots=2;

af=1; 
width=6.8*af;
height=4.5*af; 
dhor=2.4*af;
dver=1*af;
white_bottom=0.7*af;
white_left=1.6*af;
white_right=0.4*af;
white_top=1*af;
FontS=15*af; 
qq=rem(nplots,cols);
titleset={'h','i','c','d'};
figure
set(gcf,'units','centimeters','position',[3 3 (white_left+(cols-1)*dhor+cols*width+white_right) (white_bottom+white_top+(rows-1)*dver+height*rows)]);
xx=[1:1:size(actual_con_set2,1)]';
for order=1:nplots
        ax=subplot(rows,cols,order);
    
    q=rem(order,cols);
    if q==0
        row=floor(order/cols);
        col=cols;
    else
        row=floor(order/cols)+1;
        col=q;
    end
if order==1
xx0=1:1:length(nameset1);
y1=change_set3;
b1=bar(ax,xx0,y1,1);
b1(1).LineWidth=0.5;
b1(2).LineWidth=0.5;
b1(3).LineWidth=0.5;
b1(4).LineWidth=0.5;

b1(1).FaceColor=color_pre;
b1(2).FaceColor=color_ndvi;
b1(3).FaceColor=color_et;
b1(4).FaceColor=color_q;

set(ax,'xlim',[xx0(1)-0.5,xx0(end)+0.5]);
set(ax,'ytick',[-20:20:100],'ylim',[-10,100]);
ylabel(ax,'Relative change/%');
h1=legend(ax,[b1(1),b1(2),b1(3),b1(4)],{'P','NDVI','Es','Q'},'units','centimeters','location','north',...
    'box','off','numcolumns',4,'fontsize',FontS);
h1_pos=get(h1,'position');
set(h1,'position',[h1_pos(1)+0.6 h1_pos(2)+0.8 h1_pos(3) h1_pos(4)]);
set(ax,'units','centimeters','position',[white_left+(col-1)*(dhor+width) white_bottom+(rows-row)*(dver+height) width height]);
set(ax,'fontname','Times New Roman','fontweight','bold','fontsize',FontS,'box','off','linewidth',0.8);
hold on
xlim=get(ax,'xlim');
ylim=get(ax,'ylim');
set(ax,'xticklabels',namesetx);
else
  b1=bar(ax,xx,actual_con_set2,0.5,'stack');
b1(1).FaceColor=color_pre;
b1(2).FaceColor=color_ndvi;
b1(3).FaceColor=color_et;
b1(4).FaceColor=color_others;

b1(1).LineWidth=0.5;
b1(2).LineWidth=0.5;
b1(3).LineWidth=0.5;
b1(4).LineWidth=0.5;

set(ax,'xlim',[0.5,size(actual_con_set2,1)+0.5]);
l=legend({'P','NDVI','Es','C'},...
    'fontsize',FontS);
set(l,'box','off','location','north','numcolumns',5,'units','centimeters','fontsize',FontS);
pos=get(l,'position');
l.Position=[pos(1)+1.8,pos(2)+0.8,pos(3),pos(4)];
ylabel('Contribution/%');
set(ax,'ylim',[-200,300]);
set(ax,'units','centimeters','position',[white_left+(col-1)*(dhor+width) white_bottom+(rows-row)*(dver+height) width height]);
set(ax,'fontname','Times New Roman');
set(ax,'FontSize',FontS,'FontWeight','bold','linewidth',0.8,'box','off','Fontname','Times New Roman');
hold on
xlim=get(ax,'xlim');
ylim=get(ax,'ylim');
a1=0.1;
text(5-a1,170,'*','fontsize',FontS);   % nuxia
text(6.8,200,'**','fontsize',FontS);   %  tuotuohe
text(6-a1,180,'*','fontsize',FontS);    %   jiayuqiao
text(1-a1,170,'*','fontsize',FontS);    %   lazi
set(ax,'xticklabels',namesetx);


end
ax1=axes;
pos=get(ax,'position');
set(ax1,'units','centimeters','position',pos);
set(ax1,'XaxisLocation','top','YaxisLocation','Right','xlim',[xlim(1),xlim(2)],'ylim',[ylim(1),ylim(2)],'Color','none','LineWidth',0.8,...
    'Xtick',[],'Ytick',[]);
end
%% export figures
    delete('G:\variation_and_attribution_for_runoff_and_sediment_flux\codes2\output\ridge_contribution_based_on_relative_change_and_obs6.tif');
print('-dtiff','-r300','G:\variation_and_attribution_for_runoff_and_sediment_flux\codes2\output\ridge_contribution_based_on_relative_change_and_obs6.tif');

