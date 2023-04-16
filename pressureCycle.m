%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                      NAME: TRAITEMENT_PRESSIONS                         %
%                      AUTHOR: PabDawan                                   %
%                      DATE: 06/2022                                      %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Description: Traitement repartition des pressions (Intellisoles , 9 capteurs de pressions)
%% Traitement marche et statique

tic
clc                                                                         %Clear Command Window
clear                                                                       %Clear Workspace
close all                                                                   %Close figure


%% Chemin d'accès
start_path=pwd;
allFiles = genpath(pwd);
allFiles = split(allFiles,';');

rawDataFile = ~cellfun(@isempty,regexp(allFiles,'rawData'));
initialPath=allFiles{rawDataFile};

poids = [60 53];
cd(initialPath)

list=dir(initialPath);
list={list.name}'; %Liste noms des fichiers
%% Import Options and import the DATA
opts = delimitedTextImportOptions("NumVariables", 10);
opts.DataLines = [31, Inf];
opts.Delimiter = ";";
opts.VariableNames = ["mode", "VarName1", "VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

rawdata=cellfun(@(x) readtable(x,opts),list(3:end),'uniformOutput',false);  %Stockage des tableaux dans des cellules
names=extractBefore(list(3:end),'_STATIQUE');                               %Noms des sujets pour chaque fichier
condition=extractBefore(extractAfter(list(3:end),'INTELLI_'),'.csv');       %Nom des conditions pour chaque fichier
[e,r,t]=unique(names,'stable');                                             %Groupe et occurence des noms de sujets
[e1,r1,t1]=unique(condition,'stable');
noms=lower(e);

%% Attribution des régions plantaires
% 4 régions plantaires différentes: talon, midfoot distal, midfoot proximal et forefoot
for k = 1: numel(rawdata)
    currentData = rawdata{k};
    
    talon_moy=mean([currentData.VarName1, currentData.VarName6],2);
    midfoot_moy_dist=mean([currentData.VarName3,currentData.VarName9],2);
    midfoot_moy_prox=mean([currentData.VarName2,currentData.VarName4],2);
    forefoot_moy=mean([currentData.VarName5,currentData.VarName7,currentData.VarName8],2);
    
    % Vecteur temps
    t=linspace(0,length(talon_moy)/100,length(talon_moy));
    t=t';
    
    signaux_brut=[talon_moy,midfoot_moy_dist,midfoot_moy_prox,forefoot_moy];
    % plot(t,signaux_brut)
    
    
    %% Application filtre passe bas (fc=6Hz)
    fs=100;                                                                 %Fréquence d'acquisition
    fc = 6;
    [b,a] = butter(3,fc/(fs/2),'low');                                      %Butterworth, lowpass, fc = 6Hz
    S_filt_lp=filtfilt(b,a,signaux_brut);                                   %Zero-phase digital filtering
    S_filt=arrayfun(@(x) S_filt_lp(:,x),1:4,'uni',0);
    
    %% Découpe du signal: selection auto de 6 pics
    seuil=0.4*max(S_filt{4});                                                   %seuil à 40% du max
    seuil2=0.4*max(S_filt{1});
    [~,o]=findpeaks(S_filt{4},t,'MinPeakHeight',seuil,'MinPeakDistance',0.4,'minpeakprominence',10);
    [~,p]=findpeaks(S_filt{1},t,'MinPeakHeight',seuil2,'MinPeakDistance',0.4,'minpeakprominence',10);
    
    aller=o(4:6);
    retour=o(end-4:end-2);
    
    aller(1)=aller(1)-0.8;
    retour(1)=retour(1)-0.8;
    
    aller(end)=aller(end)+0.4;                                                  %je prends 0.4 sec avant le premier pic et 0.4 sec apres le dernier pics
    retour(end)=retour(end)+0.4;
    
    %% Découper le signal pour ne traiter que la fenetre entre les deux points
    
    mark1=aller(1);                                                             % marqueur auquel je veux commencer mon traitement
    mark2=aller(end);                                                           % marqueur auquel je veux finir mon traitement
    index1=find(abs(t-mark1)==min(abs(t-mark1)));                               % index1 sert a stocker la position de la valeur la plus proche de 80sec
    index2=find(abs(t-mark2)==min(abs(t-mark2)));                               % pareil: je ne suis pas sur de trouver la valeur exacte 80 ou 100 donc je cherche la plus proche
    
    t_aller=t((index1:index2),:);                                               % apres avoir pris les indexes, on crée une sous matrice avec les 20 secondes d'interet
    S_filt_aller{1}=S_filt{1}((index1:index2),1);                               % apres avoir pris les indexes, on crée une sous matrice avec les 20 secondes d'interet
    S_filt_aller{2}=S_filt{2}((index1:index2),1);                               % apres avoir pris les indexes, on crée une sous matrice avec les 20 secondes d'interet
    S_filt_aller{3}=S_filt{3}((index1:index2),1);                               % apres avoir pris les indexes, on crée une sous matrice avec les 20 secondes d'interet
    S_filt_aller{4}=S_filt{4}((index1:index2),1);                               % apres avoir pris les indexes, on crée une sous matrice avec les 20 secondes d'interet
    
    mark3=retour(1);                                                            % marqueur auquel je veux commencer mon traitement
    mark4=retour(end);                                                          % marqueur auquel je veux finir mon traitement
    index1=find(abs(t-mark3)==min(abs(t-mark3)));                               % index1 sert a stocker la position de la valeur la plus proche de 80sec
    index2=find(abs(t-mark4)==min(abs(t-mark4)));                               % pareil: je ne suis pas sur de trouver la valeur exacte 80 ou 100 donc je cherche la plus proche
    
    t_retour=t((index1:index2),:);                                              % apres avoir pris les indexes, on crée une sous matrice avec les 20 secondes d'interet
    S_filt_retour{1}=S_filt{1}((index1:index2),1);                              % apres avoir pris les indexes, on crée une sous matrice avec les 20 secondes d'interet
    S_filt_retour{2}=S_filt{2}((index1:index2),1);                              % apres avoir pris les indexes, on crée une sous matrice avec les 20 secondes d'interet
    S_filt_retour{3}=S_filt{3}((index1:index2),1);                              % apres avoir pris les indexes, on crée une sous matrice avec les 20 secondes d'interet
    S_filt_retour{4}=S_filt{4}((index1:index2),1);
    
    %% Visualisation et detection des debut de cycle et fin de cycle
    S_filt=cellfun(@(x,y) [x;y],S_filt_aller,S_filt_retour,'UniformOutput',false); %j'assemble les deux signaux pour n'en faire qu'un
    t_n=0:1/100:numel(S_filt{1})/100-1/100;
    t_n=t_n';
    % Pour la fin du cycle
    figure
    plot(t_n,[S_filt{1} S_filt{2} S_filt{3} S_filt{4}]);hold on
    deriv_fin=diff(S_filt{4})./diff(t_n);                                       % Je dérive le signal 4 (forefoot) car c'est la derniere zone en contact avec le sol
    %plot(t_n(1:end-1),deriv_fin);hold on
    find(deriv_fin(1:end-1)>0 & deriv_fin(2:end) < 0);
    
    crossZero=deriv_fin(1:end-1)>0 & deriv_fin(2:end) < 0;                      % Je trouve les points lorsque la dérivée traverse 0
    xPics=find(crossZero);
    locspic=find(S_filt{4}(crossZero)>seuil)+1;
    locspic_fin=xPics(locspic);
    scatter(t_n(locspic_fin),S_filt{4}(locspic_fin),'bo','filled')
    xline(t_n(locspic_fin),"Color","r")
    
    % Pour le debut du cycle
    deriv_deb=diff(S_filt{1})./diff(t_n);                                       % Je dérive le signal 1 (talon) car c'est la première zone en contact avec le sol (ATTENTION: pas forcément vrai en course à pied)
    %plot(t_aller(1:end-1),deriv);hold on
    %scatter(t_aller(find(deriv_deb(1:end-1)>0 & deriv_deb(2:end) < 0)),S_filt_aller{1}(find(deriv_deb(1:end-1)>0 & deriv_deb(2:end) < 0)))
    crossZero=deriv_deb(1:end-1)>0 & deriv_deb(2:end) < 0;
    xPics=find(crossZero);
    locspic=find(S_filt{1}(crossZero)>seuil2)-1;
    locspic_deb=xPics(locspic);
    scatter(t_n(locspic_deb),S_filt{1}(locspic_deb),'magenta','filled')
    xline(t_n(locspic_deb),"Color","r")
    %% rescale surface et poids
    S_filt=cellfun(@(x) x/((pi*0.6^2)*poids(k)),S_filt,'uni',false);               %divisé par la surface d'un capteur de rayon 6 mm et par le poids du sujet (kg)
    
    %% Normalisation du temps (cycle 0-100%)
    localMax_FIN=locspic_fin;
    localMax_DEB=locspic_deb;
    Contact_Norm=cell(numel(S_filt),length(localMax_FIN));
    
    for j = 1:4
        for i=1:length(localMax_FIN)
            Signal     = S_filt{j};
            HeelStrike = localMax_DEB(i);
            ToeOff     = localMax_FIN(i);
            Contact    = Signal(HeelStrike:ToeOff);
            Contact_Norm{j,i}  = interp1(1:numel(Contact), Contact, linspace(1, numel(Contact), 101))';
            
        end
    end
    
    Contact_Norm = arrayfun(@(x) cell2mat(Contact_Norm(x,:)),1:4,'uni',0);
    
    %% Je trace le résultat du sujet avec son cycle moyen sur 6 appuis (0-100%)
    % TALON
    figure
    moyenne_1=mean(Contact_Norm{1},2);
    ET_1=std(Contact_Norm{1},0,2);
    borne_sup=moyenne_1+ET_1;
    borne_inf=moyenne_1-ET_1;
    
    max_1=max(moyenne_1);                                                       % Pic de pression
    intergale_1=max(cumtrapz(moyenne_1));                                       %Pressure Time Integral (PTI) (kPa.s)
    atteinte_max1=find(moyenne_1==max(moyenne_1));                              % (%)
    
    %plot([Contact_Norm_1{1},Contact_Norm_1{2},Contact_Norm_1{3},Contact_Norm_1{4},Contact_Norm_1{5},Contact_Norm_1{6},Contact_Norm_1{7},Contact_Norm_1{8}],':')
    hold on
    plot(moyenne_1,'Color',[57 106 177]./255,'LineWidth',2)
    hold on
    plot([borne_inf,borne_sup],'Color',[57 106 177]./255,'LineWidth',0.5,'LineStyle','--','HandleVisibility','off')
    
    % Visualisation cycle moyen +- ecartype
    y = moyenne_1'; % your mean vector;
    x = 1:numel(y);
    std_dev = ET_1';
    curve1 = y + std_dev;
    curve2 = y - std_dev;
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    
    fill(x2, inBetween,[114 147 203]./255,'LineStyle','-.','FaceAlpha',.2,'HandleVisibility','off');
    %plot(x, y, 'b', 'LineWidth', 1);
    hold on
    %% Midfoot distal
    moyenne_2=mean(Contact_Norm{2},2);
    ET_2=std(Contact_Norm{2},0,2);
    borne_sup=moyenne_2+ET_2;
    borne_inf=moyenne_2-ET_2;
    
    max_2=max(moyenne_2);
    intergale_2=max(cumtrapz(moyenne_2)); %Pressure Time Integral (PTI) (kPa.s)
    atteinte_max2=find(moyenne_2*((pi*0.6^2)*poids)==max(moyenne_2*((pi*0.6^2)*poids)));
    
    %plot([Contact_Norm_1{1},Contact_Norm_1{2},Contact_Norm_1{3},Contact_Norm_1{4},Contact_Norm_1{5},Contact_Norm_1{6},Contact_Norm_1{7},Contact_Norm_1{8}],':')
    hold on
    plot(moyenne_2,'Color',[62 150 81]./255,'LineWidth',2)
    hold on
    plot([borne_inf,borne_sup],'Color',[132 186 91]./255,'LineWidth',0.5,'LineStyle','--','HandleVisibility','off')
    
    % Visualisation cycle moyen +- ecartype
    y = moyenne_2'; % your mean vector;
    x = 1:numel(y);
    std_dev = ET_2';
    curve1 = y + std_dev;
    curve2 = y - std_dev;
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    
    %Remplir et mettre transparent
    fill(x2, inBetween,[132 186 91]./255,'LineStyle','-.','FaceAlpha',.2,'HandleVisibility','off');
    hold on
    %% Midfoot proximal
    
    moyenne_3=mean(Contact_Norm{3},2);
    ET_3=std(Contact_Norm{3},0,2);
    borne_sup=moyenne_3+ET_3;
    borne_inf=moyenne_3-ET_3;
    
    max_3=max(moyenne_3);
    intergale_3=max(cumtrapz(moyenne_3));%Pressure Time Integral (PTI) (kPa.s)
    atteinte_max3=find(moyenne_3*((pi*0.6^2)*poids)==max(moyenne_3*((pi*0.6^2)*poids)));
    
    %plot([Contact_Norm_1{1},Contact_Norm_1{2},Contact_Norm_1{3},Contact_Norm_1{4},Contact_Norm_1{5},Contact_Norm_1{6},Contact_Norm_1{7},Contact_Norm_1{8}],':')
    plot(moyenne_3,'Color',[204 37 41]./255,'LineWidth',2)
    hold on
    plot([borne_inf,borne_sup],'Color',[211 94 96]./255,'LineWidth',0.5,'LineStyle','--','HandleVisibility','off')
    
    % Visualisation cycle moyen +- ecartype
    y = moyenne_3'; % your mean vector;
    x = 1:numel(y);
    std_dev = ET_3';
    curve1 = y + std_dev;
    curve2 = y - std_dev;
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    
    fill(x2, inBetween, [211 94 96]./255,'LineStyle','-.','FaceAlpha',.2,'HandleVisibility','off');
    
    hold on
    %% Forefoot
    %Stats descriptives
    moyenne_4=mean(Contact_Norm{4},2);
    ET_4=std(Contact_Norm{4},0,2);
    borne_sup=moyenne_4+ET_4;
    borne_inf=moyenne_4-ET_4;
    
    max_4=max(moyenne_4);
    intergale_4=max(cumtrapz(moyenne_4));%Pressure Time Integral (PTI) (kPa.s)
    atteinte_max4=find(moyenne_4*((pi*0.6^2)*poids)==max(moyenne_4*((pi*0.6^2)*poids)));
    
    %plot([Contact_Norm_1{1},Contact_Norm_1{2},Contact_Norm_1{3},Contact_Norm_1{4},Contact_Norm_1{5},Contact_Norm_1{6},Contact_Norm_1{7},Contact_Norm_1{8}],':')
    plot(moyenne_4,'Color',[107 76 154]./255,'LineWidth',2)
    hold on
    plot([borne_inf,borne_sup],'Color',[144 103 167]./255,'LineWidth',0.5,'LineStyle','-.','HandleVisibility','off')
    
    % Visualisation cycle moyen +- ecartype
    y = moyenne_4'; % your mean vector;
    x = 1:numel(y);
    std_dev = ET_4';    %your SD
    curve1 = y + std_dev;
    curve2 = y - std_dev;
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    
    fill(x2, inBetween, [144 103 167]./255,'LineStyle','-.','FaceAlpha',.2,'HandleVisibility','off');
    
    
    
    %% Esthetique du graphe
    hfig= gcf;  % save the figure handle in a variable
    title('Pressure evolution during stance phase','Interpreter','latex')
    subtitle(lower(names{k}),'Interpreter','latex')
    legend('Talon','midfoot distal','midfoot proximal','forefoot')
    xlabel('Stance phase (\%)')
    ylabel('Pressure ($kPa$.$kg^{-1}$)')
    picturewidth = 20; % set this parameter and keep it forever
    hw_ratio = 0.65; % feel free to play with this ratio
    set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document
    
    %set(findall(hfig,'-property','Box'),'Box','off') % optional
    set(findall(hfig,'-property','Interpreter'),'Interpreter','latex')
    set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
    set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
    pos = get(hfig,'Position');
    set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
    xlim([0 100])
    ylim([0 50])
    
    
end
