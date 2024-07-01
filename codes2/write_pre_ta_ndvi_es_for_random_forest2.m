%This section of code is used to generate input data for random forest modeling.
clear;
nameset={'Jimai','Tuotuohe','Zhimenda','Shigu','Xiangda','Daojieba','Nuxia','Rikaze','Lasha','Gengzhang','Lazi','Jiayuqiao'};
            nameset1={'Nuxia','Tuotuohe','Jiayuqiao','Rikaze','Lasha','Gengzhang','Lazi'};
num_stations=length(nameset1);
for i=1:length(nameset1)
    for j=1:length(nameset)
        a=strcmp(nameset1{i},nameset{j});
        if a==1
            cindex(i,1)=j;
        end
    end
end
%% read data
varset={'pre','ta','growing_season_noaa_cdr_ndvi','gleam_es','q'};
startyear=1950;
endyear=2022;
turning_point=1995;
num_vars=length(varset);

% 站点序号
sset=1:1:7;
% sset=[sset(1:11),sset(end)];  %去掉了昌都站，该站为第12个站

for s=1:length(nameset1)
    xset1=ones(endyear-startyear+1,length(varset))*NaN;
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

    for i=1:length(varset)
        cvset(s,i)=std(xset2(:,i))/mean(xset2(:,i));
        [b,~,~,~,stat]=regress(xset2(:,i),xx);
        if i==1
            bset(s,i)=b(2)/mean(xset2(:,i))*100;
        else
        bset(s,i)=b(2);
        end
        if stat(3)<0.05 && stat(3)>0.01
            pset(s,i)=1;
        elseif stat(3)<0.01
            pset(s,i)=11;
        else
            pset(s,i)=0;
        end
    end

      xset2_all{s,1}=xset2;  
    
    x_before=(index(1)+startyear-1):1:turning_point;
    xx_before=[ones(length(x_before),1),x_before'];
    y1=length(x_before);
        for i=1:length(varset)
       cvset_before(s,i)=std(xset2(1:y1,i))/mean(xset2(1:y1,i));
         [b_before,~,~,~,stat_before]=regress(xset2(1:y1,i),xx_before);
         if i==1
             bset_before(s,i)=b_before(2)/mean(xset2(:,i))*100;
         else
         bset_before(s,i)=b_before(2);
         end
        if stat_before(3)<0.05 && stat_before(3)>0.01
            pset_before(s,i)=1;
        elseif stat_before(3)<0.01
            pset_before(s,i)=11;
        else
            pset_before(s,i)=0;
        end
        end
    
        
            x_after=(turning_point+1):1:(index(end)+startyear-1);
    xx_after=[ones(length(x_after),1),x_after'];
    y2=length(x_after);
        for i=1:length(varset)
                    d1=mean(xset2((y1+1):end,i))-mean(xset2(1:y1,i));
                    if i==2
                       rset(s,i)=d1;  
                    else
      rset(s,i)=d1/mean(xset2(1:y1,i))*100;
                    end
                          cvset_after(s,i)=std(xset2((y1+1):end,i))/mean(xset2((y1+1):end,i));    
   [b_after,~,~,~,stat_after]=regress(xset2((y1+1):end,i),xx_after);
   if i==1
       bset_after(s,i)=b_after(2)/mean(xset2(:,i))*100;
   else
         bset_after(s,i)=b_after(2);
   end
   
        if stat_after(3)<0.05 && stat_after(3)>0.01
            pset_after(s,i)=1;
        elseif stat_after(3)<0.01
            pset_after(s,i)=11;
        else
            pset_after(s,i)=0;
        end
        end
        

       
          
end

bset1=bset;
pset1=pset;
xset3_all=xset2_all;
start_endsetx=start_endset;
v_num=length(varset);

        nameset2={'Lazi','Rikaze','Lasha','Gengzhang','Nuxia','Jiayuqiao','Tuotuohe'};
for i=1:length(nameset2)
    for j=1:length(nameset1)
a=strcmp(nameset2{i},nameset1{j});
if a==1
    cindex2(i,1)=j;
end
    end
end
xset4_all=xset3_all(cindex2,1);
start_endsetx1=start_endsetx(cindex2,:);
varset1={'pre','ta','ndvi','es','q'};
name1={'Three-year moving average serial number'};
%% write output
   output_filename='G:\variation_and_attribution_for_runoff_and_sediment_flux\codes2\input_data\relative_changes_in_pre_ta_ndvi_es_q_for_random_forest.xlsx';
delete(output_filename);
num_moving=3;
for s=1:num_stations
   xsetx1=xset4_all{s,1};
   xsetx2=xsetx1(:,[1,2,3,4,5]);
   num_years=size(xsetx2,1);
    r1=num_years-(num_moving-1);  % How many data sets are there after the moving average
      index0=1995-start_endsetx1(s,1)+1;
   ave_xset0=mean(xsetx2(1:index0,:));  % for the prior-1995 period
   
   
for i=1:v_num
    for j=1:num_years
        delta_x=xsetx2(j,i)-ave_xset0(i);
        if i==2
             delta_x_set(j,i)=delta_x;
        else
      delta_x_set(j,i)=delta_x/ave_xset0(i);
        end
    end 
    
for i1=1:r1
    delta_x_set_moving(i1,i)=mean(delta_x_set(i1:(i1+2),i));
end
end


%Performing z-score normalization on data to mitigate the influence of magnitude differences between different variables on the results
zset=[];
for i=1:v_num
    ave_std_delta(i,1)=mean(delta_x_set_moving(:,i));
ave_std_delta(i,2)=std(delta_x_set_moving(:,i));
       zset(:,i)=(delta_x_set_moving(:,i)-ave_std_delta(i,1))./ave_std_delta(i,2);
end
   
    
serial_number=[1:1:size(zset,1)]';

   
   xlswrite(output_filename,name1,['s' num2str(s)],'a1');
xlswrite(output_filename,serial_number,['s' num2str(s)],'a2');
xlswrite(output_filename,varset1,['s' num2str(s)],'b1');
xlswrite(output_filename,zset,['s' num2str(s)],'b2');
   
end






