close all 
clc
clear 
%-----------------------------------------------------------------------
%Input data
%-----------------------------------------------------------------------
country='United States of America';
gauss = 'gauss1'; 
%-----------------------------------------------------------------------
%import data
data=importdata('WHO-COVID-19-global-data.csv');
CountryName=data.textdata(2:end,3);
day1=data.textdata(2:end,1);
CumulativeConfirmed=data.data(1:end,4);
%-----------------------------------------------------------------------
flagFirst=true;
for i=1:length(day1)-1
    if strcmp(CountryName(i+1),country)&&CumulativeConfirmed(i+1)>0 && flagFirst
            first =i+1;
            flagFirst=false;
    end
    if strcmp(CountryName(i),country)&& ~strcmp(CountryName(i+1),country)
            ending=i;
    end
end
ending =ending-1; %removing today, maybe do not actualized yet
diffwofilter = diff(CumulativeConfirmed(first:ending)); %derivative without filter
ma_order=10;                                            %moving average order
h=ones(1,ma_order)/ma_order;                            %h of filter
tcf =filtfilt(h,1,CumulativeConfirmed(first:ending));           
dif = diff(tcf);                                        %filtered derivative
t=0:ending-first-1;
model = fit(t',dif,gauss);                              %fit gaussian model
%--------------------------------------------------------------------------
v=[];
predictionsDays=0;
while length(v) <1
    gf=model(0:(ending-first+predictionsDays));         %gaussian prediction
    [m,I]=max(gf);             
    v=find(gf(I:end)<1);                                %find first day with cases less than 1
    predictionsDays=predictionsDays+1;
end
dt=datetime(day1(first:ending));
dd=datetime((dt(end)+ days(1:predictionsDays)))';
dd.Format='dd-MMM-yyyy';
df=[dt;dd];
pos=v(1)+I-1;
dataending=df(pos);
disp(['first day with cases <1: ',datestr(dataending)]); 
gfacc=cumtrapz(gf);                                     %sigmoide (gaussian integration)
%-----------------------------------------------------------------------
%Graphics plot
%-----------------------------------------------------------------------
subplot(2,1,1);
plot(CumulativeConfirmed(first:ending),'-o')
hold on
title(strcat(country,' - Confirmed Cases'));
ylabel('cases');
plot(gfacc)
grid
legend('total cases',strcat(gauss,' fit acc.'),'location','northwest');
%----------------------------------------------------------------------
subplot(2,1,2);
plot(diff(CumulativeConfirmed(first:ending)))
hold on
plot(gf)
plot(pos,2,'+r');
title('Derivative Analysis');
xlabel('days after first case');
ylabel('daily cases');
legend('derivative',strcat(gauss,' fit'),'location','northwest');
grid


