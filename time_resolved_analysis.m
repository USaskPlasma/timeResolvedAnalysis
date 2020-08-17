clear;
clc;
close all;
filename = dir('TEK*.CSV'); %find the scope files
Numelement=length(filename); %get the number of scope files
delimiterIn = ','; %define the delimiter of the file
headerlinesIn = 16; %define the number of header lines
%
R_Sensor=150; %sensor resistor
Z_probe=4.1; %probe impedance
Z_total=R_Sensor+Z_probe;
%length probe
Lp=5e-3;
%probe radius
Rp=0.5e-3; %1mm diameter
Ap=2*pi*Rp*Lp;
%
%read the first file to extract the Time base it should be the same for all
%file if the trigger setting wasn't changed
A = importdata(filename(1).name,delimiterIn,headerlinesIn); %read the data
TIME = A.data(:,1);  %time array
%preallocate the memory
TIME_check = zeros(length(TIME),Numelement);%check that the time vector is really a constant can be done manually at the end of analysis
I_target = zeros(length(TIME),Numelement);%current throught target
V_pulse = zeros(length(TIME),Numelement);%high voltage pulse
V_powersupply = zeros(length(TIME),Numelement);%Voltage given by the power supply connected to probe (not the tip voltage)
I_probe = zeros(length(TIME),Numelement);%probe current
V_probe = zeros(length(TIME),Numelement);%take into account the voltage drop due to Rsensor and Z_probe;
%creat the vector of qunatities of interest
Te=zeros(length(TIME),1); %electron temperature
V_plasma=zeros(length(TIME),1); %plasma potential
Ie_sat=zeros(length(TIME),1); % electron saturation current
%Ne=zeros(length(TIME),1); %electron density (OML theory)
errTe_min=zeros(length(TIME),1); %errors
errTe_max=zeros(length(TIME),1);
errV_plasma_min=zeros(length(TIME),1);
errV_plasma_max=zeros(length(TIME),1);
errIe_sat_min=zeros(length(TIME),1);
errIe_sat_max=zeros(length(TIME),1);
%errNe_max=zeros(length(TIME),1);
%errNe_min=zeros(length(TIME),1);
%loop through all the data
f = waitbar(0,'1','Name','Analysing data...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
for j=1:Numelement 
    A = importdata(filename(j).name,delimiterIn,headerlinesIn); %read the data
    TIME_check(:,j) = A.data(1:length(TIME),1);
    I_target(:,j) =  A.data(1:length(TIME),3); %current throught target
    V_pulse(:,j) = A.data(1:length(TIME),4); %high voltage pulse
    V_powersupply(:,j) = A.data(1:length(TIME),5); %Voltage given by the power supply connected to probe (not the tip voltage)
    I_probe(:,j) = (A.data(1:length(TIME),2)-0.05)/R_Sensor; %probe current it seems like there is an offset in the measurement of +0.05V so i remove it
    V_probe(:,j) = V_powersupply(:,j)-Z_total.*I_probe(:,j); %take into account the voltage drop due to Rsensor and Z_probe;
end
%let's order the value
[V_probe_ordered,Indexes]=sort(V_probe,2); %ascending oreder of probe voltage at all times
I_probe_ordered = zeros(length(TIME),Numelement);%ordered probe current according to voltage ordering
for k=1:length(TIME)-1
    % Check for clicked Cancel button
    if getappdata(f,'canceling')
        break
    end
    
    % Update waitbar and message
    waitbar(k/length(TIME),f,sprintf('Wait for it...%5.5f',k/length(TIME)))
    
    I_probe_ordered(k,:)=I_probe(k,Indexes(k,:)); %current ordered with respect to probe voltage
    %lest analyse the data while we are in the loop
    V_analysis=V_probe_ordered(k,:);
    I_analysis=I_probe_ordered(k,:);
%     %lets remove the ion current (assuming OML sqrt(V) dependance)
%     I2=I_analysis.^2;
%     %prepare the data
%     [V_analysis, I2] = prepareCurveData( V_analysis, I2);
%     % Set up fittype and options.
%     ft = fittype( 'poly1' );
%     excludedPoints = (V_analysis > -20);
%     opts = fitoptions( 'Method', 'LinearLeastSquares' );
%     opts.Exclude = excludedPoints;
%     % Fit model to data.
%     [ion_current, ~] = fit( V_analysis, I2, ft, opts );
%     I_ion=(ion_current(V_analysis));
%     I_ion(I_ion<0)=0;
%     I_ion=sqrt(I_ion);
%     I_analysis=I_analysis+I_ion'; %+ sign because ion current is negative
    %let's do the langmuir analysis
    %we remove all negative current that remain
    indices=I_analysis>0;
    logI=log(I_analysis(indices));
    V_log_analysis=V_analysis(indices);
    [V_log_analysis, logI] = prepareCurveData( V_log_analysis, logI );
    % Set up fittype and options.
    ft = fittype( 'poly1' );
    excludedPoints = (logI < -6.9) | (logI > max(logI)-1.3);
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    opts.Exclude = excludedPoints;
    % Fit model to data.
    [electron_current, ~] = fit( V_log_analysis, logI, ft, opts );
    %check that the fit is not crazy otherwise redo with slighttly extended
    %boundaries
    ci = confint(electron_current);
    if ci(1,1)<0.067
        excludedPoints = (logI < -7.0) | (logI > max(logI)-1.2);
        opts.Exclude = excludedPoints;
        [electron_current, ~] = fit( V_log_analysis, logI, ft, opts );
        ci = confint(electron_current);
    end
    %extract the electron temperature in eV
    Te(k)=1/(electron_current.p1);
    errTe_min(k)=Te(k)-1/ci(2,1);
    errTe_max(k)=1/ci(1,1)-Te(k);
    %fit with imposed max Te
    opts.Upper = [ci(1,1) Inf];
    opts.Lower = [ci(1,1) -Inf];
    [electron_current_Temax, ~] = fit( V_log_analysis, logI, ft, opts );
    %fit with imposed min Te
    opts.Upper = [ci(2,1) Inf];
    opts.Lower = [ci(2,1) -Inf];
    [electron_current_Temin, ~] = fit( V_log_analysis, logI, ft, opts );
    %fit the electron saturation current
    % Set up fittype and options.
    ft = fittype( 'poly1' );
    excludedPoints =  (logI < max(logI)-0.7);
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    opts.Exclude = excludedPoints;
    % Fit model to data.
    [electron_saturation, ~] = fit( V_log_analysis, logI, ft, opts );
    ci_es = confint(electron_saturation);
    %fit saturation imposed max slop
    opts.Upper = [ci_es(2,1) Inf];
    opts.Lower = [ci_es(2,1) -Inf];
    [electron_saturation_maxslope, ~] = fit( V_log_analysis, logI, ft, opts );
    %fit saturation imposed min slop
    opts.Upper = [ci_es(1,1) Inf];
    opts.Lower = [ci_es(1,1) -Inf];
    [electron_saturation_minslope, ~] = fit( V_log_analysis, logI, ft, opts );
    %extract the quantities
    V_plasma(k)=roots(coeffvalues(electron_current)-coeffvalues(electron_saturation));
    errV_plasma_max(k)=roots(coeffvalues(electron_current_Temax)-coeffvalues(electron_saturation_minslope))-V_plasma(k);
    errV_plasma_min(k)=V_plasma(k)-roots(coeffvalues(electron_current_Temin)-coeffvalues(electron_saturation_maxslope));
    Ie_sat(k)=exp(electron_current(V_plasma(k)));
    errIe_sat_max(k)=exp(electron_current_Temax(V_plasma(k)+errV_plasma_max(k)))-Ie_sat(k);
    errIe_sat_min(k)=Ie_sat(k)-exp(electron_current_Temin(V_plasma(k)-errV_plasma_min(k)));
    %
end
delete(f)
%convert saturation current to electron densities
%Merlino
%Am. J. Phys. 75 (12), December 2007
elm=1.602e-19;
me=9.11e-31;
vthe=sqrt(8*elm.*Te./(pi*me));
vthemin=sqrt(8*elm.*(Te-errTe_min)./(pi*me));
vthemax=sqrt(8*elm.*(Te+errTe_max)./(pi*me));
Ne=4.*Ie_sat./(elm.*vthe.*Ap);
errNe_min=Ne-4.*(Ie_sat-errIe_sat_min)./(elm.*vthemax.*Ap);
errNe_max=4.*(Ie_sat+errIe_sat_max)./(elm.*vthemin.*Ap)-Ne;
%
figure
errTe=vertcat(smooth(errTe_max,7)',smooth(errTe_min,7)');
shadedErrorBar((TIME.*1e6)',smooth(Te,7)',errTe);
set(gca,'xlim',[-10 35]);
box on
set(gca,'linewidth',2);
set(gca,'FontSize',16,'FontName','Times');
xlabel('Time ({\mu}s)','FontSize',16,'FontName','Times');
ylabel('T_e (eV)','FontSize',16,'FontName','Times');
%
figure
errNe=vertcat(smooth(errNe_max,7)',smooth(errNe_min,7)');
errNe(isnan(errNe))=0.;
shadedErrorBar((TIME.*1e6)',1e-15.*smooth(Ne,7)',1e-15.*errNe);
set(gca,'xlim',[-10 35]);
box on
set(gca,'linewidth',2);
set(gca,'FontSize',16,'FontName','Times');
xlabel('Time ({\mu}s)','FontSize',16,'FontName','Times');
ylabel('n_e ({\times}10^9 cm^{-3})','FontSize',16,'FontName','Times');
%
figure
errVp=vertcat(smooth(errV_plasma_max,7)',smooth(errV_plasma_min,7)');
shadedErrorBar((TIME.*1e6)',smooth(V_plasma,7)',errVp);
set(gca,'xlim',[-10 35]);
box on
set(gca,'linewidth',2);
set(gca,'FontSize',16,'FontName','Times');
xlabel('Time ({\mu}s)','FontSize',16,'FontName','Times');
ylabel('V_p (V)','FontSize',16,'FontName','Times');
%
%ion saturation current
Iisat=mean(I_probe_ordered(:,1:5),2);
errIisat=std(I_probe_ordered(:,1:5),0,2);
%calulate the ion density
%Bohm velocity
mi=40.*1.66e-27;
ub=sqrt(elm.*Te./mi);
ub_min=sqrt(elm.*(Te-errTe_min)./mi);
ub_max=sqrt(elm.*(Te+errTe_max)./mi);
Ni=-Iisat./(0.6.*elm.*ub.*Ap);
errNi_min=Ni-(-Iisat+errIisat)./(0.6.*elm.*ub_max.*Ap);
errNi_max=-Ni+(-Iisat-errIisat)./(0.6.*elm.*ub_min.*Ap);
figure
errNi=vertcat(smooth(errNi_max,7)',smooth(errNi_min,7)');
errNi(isnan(errNi))=0.;
shadedErrorBar((TIME.*1e6)',1e-15.*smooth(Ni,7)',1e-15.*errNi);
set(gca,'xlim',[-10 35]);
box on
set(gca,'linewidth',2);
set(gca,'FontSize',16,'FontName','Times');
xlabel('Time ({\mu}s)','FontSize',16,'FontName','Times');
ylabel('n_i ({\times}10^9 cm^{-3})','FontSize',16,'FontName','Times');
