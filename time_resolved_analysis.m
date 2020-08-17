clear;
clc;
cd 'C:\Users\alex_\Desktop\magnetron\2020 probe\7-8-2020 Si ICP time resolved';
filename = dir('*.CSV'); %find the scope files
Numelement=length(filename); %get the number of scope files
delimiterIn = ','; %define the delimiter of the file
headerlinesIn = 12; %12 for Al %define the number of header lines
R1_Sensor=100; %sense resistor 1
Z_probe=0; %probe impedance = 0 just a piece of wire
Z_totalR1=R1_Sensor+Z_probe;
Lp=6E-3; %probe length 
Rp=(.31e-3)/2; %probe radius
Ap=2*pi*Rp*Lp;
A = importdata(filename(1).name,delimiterIn,headerlinesIn); %read the data
TIME = A.data(:,1); %time array
%% preallocate the memory
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
errTe_min=zeros(length(TIME),1); %errors
errTe_max=zeros(length(TIME),1);
errV_plasma_min=zeros(length(TIME),1);
errV_plasma_max=zeros(length(TIME),1);
errIe_sat_min=zeros(length(TIME),1);
errIe_sat_max=zeros(length(TIME),1);
cantfit = zeros(length(TIME),1);
%% Load in data loop
for j=1:Numelement 
    A = importdata(filename(j).name,delimiterIn,headerlinesIn); %read the data
    TIME_check(:,j) = A.data(1:length(TIME),1);
    I_target(:,j) =  A.data(1:length(TIME),3); %current through target
    V_pulse(:,j) = A.data(1:length(TIME),2); %high voltage pulse
    V_powersupply(:,j) = A.data(1:length(TIME),5); %Voltage given by the power supply connected to probe (not the tip voltage)     
    I_probe(:,j) = (A.data(1:length(TIME),4))/R1_Sensor;
    V_probe(:,j) = V_powersupply(:,j)-Z_totalR1.*I_probe(:,j); %take into account the voltage drop due to Rsensor and Z_probe;
end
%% let's order the value
[V_probe_ordered,Indexes]=sort(V_probe,2); %ascending oreder of probe voltage at all times
I_probe_ordered = zeros(length(TIME),Numelement);%ordered probe current according to voltage ordering
%% Plots
plotcurves2d = 1;
plotcurves3d = 0;
if plotcurves3d == 1
    figure
    hold on   
    view(3)
    xlabel('Voltage(V)')
    ylabel('Time(\mus)')
    zlabel('Log Current(log(A))')
    set(gca, 'YDir','reverse','fontname','times','fontsize',14)
    zlim([-10 -0])
    ylim([TIME(1)*1e6 TIME(end)*1e6])
    grid on
    box on
end
% Anaylsis loop
for k = NearVal(10e-6,TIME)%1:5000:size(TIME,1)
    I_probe_ordered(k,:)=I_probe(k,Indexes(k,:)); %current ordered with respect to probe voltage
    V_analysis=V_probe_ordered(k,:);
    I_analysis=I_probe_ordered(k,:);
    indices=(I_analysis>0);
    logI=log(I_analysis(indices));
    V_log_analysis=V_analysis(indices);
    [V_log_analysis, logI] = prepareCurveData( V_log_analysis, logI );
    %% Transition region fit
    try
        [xData, yData] = prepareCurveData( V_log_analysis, logI );
        I_outlie = isoutlier(logI,'ThresholdFactor',8.5); %Threshold to reject outliers
        V_outlie = isoutlier(logI,'ThresholdFactor',8.5); 
        outlie = V_outlie | I_outlie;
        if plotcurves2d == 1 
            figure; hold on; box on; grid on
            scatter(V_log_analysis,logI,'filled')
            scatter(V_log_analysis(I_outlie),logI(I_outlie),'r','filled')
            scatter(V_log_analysis(V_outlie),logI(V_outlie),'g','filled')
            legend('Original Data','I outliers','V outliers','location','southeast');
        end
        ft = fittype( 'a*atan(x*b+c)+d', 'independent', 'x', 'dependent', 'y' ); %Fit arctan to I-V curve
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.Upper = [10,10,inf,100];
        opts.Lower = [0,0.1,-inf,-100];
        excludedPoints0 = V_log_analysis > .6 * max(V_log_analysis) | outlie; 
        opts.Exclude = excludedPoints0;
        [fitresult, gof] = fit( xData, yData, ft, opts );
        coef = coeffvalues(fitresult);
        a = coef(1); b = coef(2); c = coef(3); d =coef(4);
        ddfit =@(x) (-(2*a*b^2*(c + b.*x))./(1 +(c + b.*x).^2).^2); %2nd derivative of arctan
        [~,corner] = min(ddfit(V_log_analysis));
        corner = corner + 13; % Change value to move claculated plasma potiential
        %% Electron transition
        ft = fittype( 'poly1' );
        excludedPoints1 = V_log_analysis >=  V_log_analysis(corner) | outlie;
        opts = fitoptions( 'Method', 'LinearLeastSquares' );
        opts.Exclude = excludedPoints1;
        [electron_current, ~] = fit( V_log_analysis, logI, ft, opts );
        ci = confint(electron_current);
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
        %% Electron saturation fit
        % Set up fittype and options.
        ft = fittype( 'poly1' );
        excludedPoints2 =   V_log_analysis < V_log_analysis(corner+5)| outlie; 
        opts = fitoptions( 'Method', 'LinearLeastSquares' );
        opts.Exclude = excludedPoints2;
        % Fit model to data.
        [electron_saturation, ~] = fit( V_log_analysis, logI, ft, opts );
        %% Plot IV curves
        if plotcurves3d == 1
            plot3(V_log_analysis,TIME(k)*ones(size(V_log_analysis))*1e6,logI,'color','k')
            electron_currentfunction = electron_current.p1.*V_log_analysis+electron_current.p2;
            e1plt = plot3(V_log_analysis,TIME(k)*ones(size(V_log_analysis))*1e6,electron_currentfunction,'-b','linewidth',.01);
            e1plt.Color=[0,0,0,0.4];
            electron_saturationfunction = electron_saturation.p1.*V_log_analysis+electron_saturation.p2;
            e2plt = plot3(V_log_analysis,TIME(k)*ones(size(V_log_analysis))*1e6,electron_saturationfunction,'-b','linewidth',.01); 
            e2plt.Color=[0,0,0,0.4];
        end
        if plotcurves2d == 1
            figure; hold on;
            title(sprintf("%.3f \\mus",k*2e-9*1e6));
            scatter(V_log_analysis,logI,'k','filled')
            scatter(V_log_analysis(~excludedPoints2),logI(~excludedPoints2),'r','filled')
            scatter(V_log_analysis(~excludedPoints1),logI(~excludedPoints1),'b','filled')
            xline(V_log_analysis(corner),'--b');
            electron_currentfunction = electron_current.p1.*V_log_analysis+electron_current.p2;
            e1plt = plot(V_log_analysis,electron_currentfunction,'-b','linewidth',.01);
            e1plt.Color=[0,0,0,0.4];
            electron_saturationfunction = electron_saturation.p1.*V_log_analysis+electron_saturation.p2;
            e2plt = plot(V_log_analysis,electron_saturationfunction,'-b','linewidth',.01); %TIME(k)*ones(size(V_log_analysis))*1e6
            e2plt.Color=[0,0,0,0.4];
        end
        %%
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

    catch e
        disp(e.identifier)
        disp(e.message)
        cantfit(k) = 1;
        V_plasma(k)=NaN;
        errV_plasma_max(k)=NaN;
        errV_plasma_min(k)=NaN;
        Ie_sat(k)=NaN;
        errIe_sat_max(k)=NaN;
        errIe_sat_min(k)=NaN; 
    end
end

return
%% Set start of pulse as T = 0
Von = V_pulse(:,1) < mean(V_pulse(:,1));
PulseOnIndex = find(Von == 1,1);
TIME = TIME - TIME(PulseOnIndex);
%% convert saturation current to electron densities
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
%% Electron Temperature plot
figure; grid on; box on
errTe=vertcat(errTe_max',errTe_min');
Teplot = Te;
errTeplot = errTe;
shadedErrorBar((TIME.*1e6)',Teplot',errTeplot);
set(gca,'FontSize',14,'FontName','Times');
xlabel('Time (\mus)','FontSize',14,'FontName','Times');
ylabel('T_e (eV)','FontSize',14,'FontName','Times');
%% Electron Density plot
figure; grid on; box on
errNe=vertcat(errNe_max'*2/3,errNe_min'*2/3);
Neplot = Ne;
errNeplot = errNe;
shadedErrorBar((TIME.*1e6)',abs(1e-15.*Neplot'),abs(1e-15.*errNeplot));
set(gca,'FontSize',14,'FontName','Times');
xlabel('Time (\mus)','FontSize',14,'FontName','Times');
ylabel('n_e (\times 10^9 cm^{-3})','FontSize',14,'FontName','Times');
%% Plasma Potiential plot
figure; grid on; box on
errVp=vertcat(errV_plasma_max',errV_plasma_min');
errVpplot = errVp;
V_plasmaplot = V_plasma;
shadedErrorBar((TIME.*1e6)',V_plasmaplot',abs(errVpplot));
set(gca,'FontSize',14,'FontName','Times');
xlabel('Time (\mus)','FontSize',14,'FontName','Times');
ylabel('V_p (V)','FontSize',14,'FontName','Times');
%% Ion Density
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
figure; grid on; box on
errNi=vertcat(errNi_max',errNi_min');
errNi(isnan(errNi))=0.;
errNiplot = errNi;
Niplot = Ni;
shadedErrorBar((TIME.*1e6)',1e-15.*Niplot',1e-15.*errNiplot);
set(gca,'FontSize',14,'FontName','Times');
xlabel('Time (\mus)','FontSize',14,'FontName','Times');
ylabel('n_i (\times10^{9} cm^{-3})','FontSize',14,'FontName','Times');
%% Plot all in one figure
figure('position',[300 50 600 820])
ha = tight_subplot(5,1,[.02 .03],[.04 .01],[.1 .08]);
axes(ha(1)); hold on; grid on ; box on;
yyaxis left
plot(TIME,mean(V_pulse,2)*1000,'k');
ylabel('Voltage (V)','FontSize',11);
set(gca,'YTickLabel',cellstr(num2str(get(gca,'YTick')')));
yyaxis right
plot(TIME,mean(I_target,2)*2,'k','linestyle','--');
ylabel('Current (A)','FontSize',11);
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
legend('Voltage','Current','location','east');
set(gca,'YTickLabel',cellstr(num2str(get(gca,'YTick')')),'FontName','Times')

axes(ha(2)); grid on; box on;
shadedErrorBar((TIME.*1e6)',Teplot',errTeplot);
xlim([0 140]); ylim([0 20]);
ylabel('T_e (eV)','FontSize',11);
set(gca,'YTickLabel',cellstr(num2str(get(gca,'YTick')')),'FontName','Times')

axes(ha(3)); grid on; box on;
shadedErrorBar((TIME.*1e6)',abs(1e-15.*Neplot'),abs(1e-15.*errNeplot));
xlim([0 140]); ylim([0 200]);
ylabel('n_e (\times 10^9 cm^{-3})','FontSize',11);
set(gca,'YTickLabel',cellstr(num2str(get(gca,'YTick')')),'FontName','Times')

axes(ha(4)); grid on; box on;
shadedErrorBar((TIME.*1e6)',V_plasmaplot',abs(errVpplot));
xlim([0 140]); ylim([0 60]);
ylabel('V_p (V)','FontSize',11);
set(gca,'YTickLabel',cellstr(num2str(get(gca,'YTick')')),'FontName','Times')

axes(ha(5)); grid on; box on;
shadedErrorBar((TIME.*1e6)',1e-15.*Niplot',1e-15.*errNiplot);
xlim([0 140]); ylim([0 800]);
xlabel('Time (\mus)','FontSize',11);
ylabel('n_i (\times10^{9} cm^{-3})','FontSize',11);
set(gca,'YTickLabel',cellstr(num2str(get(gca,'YTick')')), ...
    'xTickLabel',cellstr(num2str(get(gca,'xTick')')) ,'FontName','Times')

%%
%Takes time value x and time vector t
%Returns element of t that is closest to x
function nearest_value = NearVal(x,t) 
[~,nearest_value] = min(abs(x-t));
end
