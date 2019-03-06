function climate=climterp_linear(deltaT)

climate=zeros(2200000,25);

%weather/radiation data
climdat=csvread('SampleWeather.csv');

[l,w]=size(climdat);
mult=linspace(1,l,(1+(l-1)*3600)/deltaT);

for col=2:w
    new=interp1(climdat(:,1),climdat(:,col),mult);
    if (col==2)
        climate=new';
    else
        climate=[climate new'];
    end
end

check_nonneg=climate(:,3:20)>=0; 
climate(:,3:20)=check_nonneg.*climate(:,3:20);

save('climate.mat','climate');

end