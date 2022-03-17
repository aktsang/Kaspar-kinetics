function fitOnKinetics (datasets)

%window of data to fit
g = fittype('a-b*exp(-c*(x-d))');
hF = figure; 

for Dix = 1:length(datasets)
    D = datasets(Dix);
    
    color = 'k';
    switch D.construct
        case 'WT'
            color = [0 0 1];
        case 'v857'
            color = [0.8 0.8 0];
        case 'v82'
            color = [1 0.25 0];
    end
    
    kon = zeros(1, length(D.conc));
    for Cix = 1:length(D.conc) %for each concentration
        
        avg = mean(D.rawTraces{Cix},2);
        
        tEnd = 4*find(avg>(0.1*avg(1) + 0.9*avg(200)), 1, 'first');
        tStart = min(find(D.t>0.002, 1, 'first'), find(avg>(0.66*avg(1) + 0.33*avg(200)), 1, 'first'));
        selT = tStart:tEnd;
        
        x =D.t(selT)';
        y = avg(selT);
        
        gOpts = fitoptions(g);
        gOpts.StartPoint = [mean(avg) (mean(avg)-y(1)) 1000 0];
        gOpts.Upper = [1.5*max(y) max(y) 10000 0.0005];
        gOpts.Lower = [0.5*max(y) 0.01 10 0];
        f0 = fit(x,y,g,gOpts);
        
        kon(Cix) = f0.c;
    end
    
    %fit KM
    selKM = kon>max(kon)/4; %data at very low ligand concentrations don't obey assumptions about excess ligand
    h = fittype('a*x/(b+x)'); %initial rate; a = k12, b = KM
    hOpts = fitoptions(h);
    hOpts.StartPoint = [10000 10000];
    hOpts.Upper = [1e6 50000];
    hOpts.Lower = [1 1];
    f1 = fit(D.conc(selKM)',kon(selKM)',h,hOpts);
    
    scatter(D.conc,kon, 'markerfacecolor', color, 'markeredgecolor', color); hold on, plot(0:1e4, f1(0:1e4), 'color', color, 'linewidth', 2);
    KM(Dix) = f1.b;
    tmp = confint(f1,normcdf(1) - normcdf(-1));
    KME(Dix,:) = tmp(:,2)';
    disp([D.construct ' K_M = ' int2str(f1.b) 'uM'])
    disp([D.construct ' v_max = ' int2str(f1.a) 'uM'])
    
end

set(gca, 'tickdir', 'out', 'linewidth', 2, 'xlim', [0 5000])
xlabel('[glutamate] (\muM)');
ylabel('k_{obs} (s^{-1})');

figure, bar(KM, 'linestyle', 'none'); hold on, errorbar(1:length(KM), KM, abs(KM-KME(:,1)'), abs(KM-KME(:,2)'), 'k', 'linewidth', 2, 'linestyle', 'none');
set(gca, 'tickdir', 'out', 'xticklabel',{datasets.construct}, 'xlim', [0.5 4.5], 'box', 'off', 'ytick', (0:2:30)*1000)
end
