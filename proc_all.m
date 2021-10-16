% Import cases spreadsheet
clear
cas=importdata('./InputFiles/US_COVID_Cases_Deaths_Sorted.csv');
cas_states=unique(cas.textdata(2:end,2));
cas_dates=unique(cas.textdata(2:end,1));
[Y,I]=sort(datenum(cas_dates));
cas_dates=cas_dates(I);

% Import age-stratified state population spreadsheet
pop=importdata('./InputFiles/sc-est2019-agesex-civ_Sorted.csv');
pop_states=unique(pop.textdata(2:end,2));

% Import vaccine coverage data
vax=importdata('./InputFiles/COVID-19_US_Vaccinations_Sorted.csv');
vax_dates=unique(vax.textdata(2:end,1));
vax_states=unique(vax.textdata(2:end,2));
[Y,I]=sort(datenum(vax_dates));
vax_dates=vax_dates(I);

% import deaths by age and sex
% 21 months, 3 sex groups (all, male, female) and 17 age groups
% and 53 states
dth=importdata('./InputFiles/Deaths_by_Sex_and_Age_Sorted.csv');
dth_dates=unique(dth.textdata(2:end,1));
dth_states=unique(dth.textdata(2:end,4));
[Y,I]=sort(datenum(dth_dates));
dth_dates=dth_dates(I);
dth_dates=dth_dates(2:21); % excluded '1/1/2020 and 9/1/2021' since it's not a complete month

% get unique states, including NYC as separate state
states_wnyc=intersect(cas_states,dth_states);

% add values from NYC to the values for NY in the cas and dth datasets
[LIA,LOCB] = ismember(cas.textdata(:,2),'NYC');
cas_NYC_data=cas.data(find(LIA==1)-1,:);
[LIA,LOCB] = ismember(cas.textdata(:,2),'NY');
cas_NY_data=cas.data(find(LIA==1)-1,:);
cas_NY_NYC_data=cas_NY_data+cas_NYC_data;
cas.data(find(LIA==1)-1,:)=cas_NY_NYC_data;

[LIA,LOCB] = ismember(dth.textdata(:,4),'NYC');
dth_NYC_data=dth.data(find(LIA==1)-1,:);
[LIA,LOCB] = ismember(dth.textdata(:,4),'NY');
dth_NY_data=dth.data(find(LIA==1)-1,:);
dth_NY_NYC_data=dth_NY_data+dth_NYC_data;
dth.data(find(LIA==1)-1,:)=dth_NY_NYC_data;

states=intersect(states_wnyc,vax_states);
states=intersect(states,pop_states);

% Filter out states not in all spreadsheets
[LIA,LOCB] = ismember(cas.textdata(:,2),states);
cas.textdata=[cas.textdata(1,:); cas.textdata(find(LIA==1),:)];
cas.data=cas.data(find(LIA==1)-1,:);

[LIA,LOCB] = ismember(vax.textdata(:,2),states);
vax.textdata=[vax.textdata(1,:); vax.textdata(find(LIA==1),:)];
vax.data=vax.data(find(LIA==1)-1,:);

[LIA,LOCB] = ismember(dth.textdata(:,4),states);
dth.textdata=[dth.textdata(1,:); dth.textdata(find(LIA==1),:)];
dth.data=dth.data(find(LIA==1)-1,:);

[LIA,LOCB] = ismember(pop.textdata(:,2),states);
pop.textdata=[pop.textdata(1,:); pop.textdata(find(LIA==1),:)];
pop.data=pop.data(find(LIA==1)-1,:);

% Filter out dates. Should only have 19 dates (for each month starting 2/1/2020)
% since will be extracting cumulative cases, cumulative vaccinations at beginning of
% each month, and also number of deaths for the total month
[LIA,LOCB] = ismember(vax.textdata(:,1),dth_dates);
vax.textdata=[vax.textdata(1,:); vax.textdata(find(LIA==1),:)];
vax.data=vax.data(find(LIA==1)-1,:);

[LIA,LOCB] = ismember(cas.textdata(:,1),dth_dates);
cas.textdata=[cas.textdata(1,:); cas.textdata(find(LIA==1),:)];
cas.data=cas.data(find(LIA==1)-1,:);

[LIA,LOCB] = ismember(dth.textdata(:,1),dth_dates);
dth.textdata=[dth.textdata(1,:); dth.textdata(find(LIA==1),:)];
dth.data=dth.data(find(LIA==1)-1,:);

% dth will have 53040 rows (20 months * 3 sex groups * 17 age groups * 52 states)
% vax will have 468 rows (9 months * 52 states) 
% cas will have 1040 rows (20 months * 52 states)
% pop will have 4524 rows (87 age groups * 52)

%%
cas_labels={'date','state','tot_cases','conf_cases','prob_cases','new_case',...
'pnew_case','tot_death','conf_death','prob_death','new_death','pnew_death'};
close_figures 

make_figs=0; % Boolean to indicate with want to save out all plots or just write out results tables
months_idx=[2:8]; % Start February, end with August

%atypes={'vax-admin-log-rob'};
atypes={'vax-admin-log-rob'};
%atypes={'vax-admin-ages-gt65-log-rob'};
%atypes={'vax-adj-cas-log-rob'};
%atypes={'vax-adj-cas-log-rob-lag0'};
%atypes={'vax-adj-cas-ages-gt65-log-rob'};

for at=1:length(atypes)
    clear vax_id2 deaths
    if ~isempty(strfind(atypes{at},'lt18-gt18'))
        vax_id=[2,4]; vax_id2=[4];
    elseif ~isempty(strfind(atypes{at},'ages-18-65'))
        vax_id=[2,4,5]; % subtracts gt18 from total, then subtracts gt65
    elseif ~isempty(strfind(atypes{at},'ages-lt12'))
        vax_id=[2,3];
    elseif ~isempty(strfind(atypes{at},'ages-lt18'))
        vax_id=[2,4]; 
    elseif ~isempty(strfind(atypes{at},'ages-lt65'))
        vax_id=[2,5]; 
    elseif ~isempty(strfind(atypes{at},'ages-gt12'))
        vax_id=3; 
    elseif ~isempty(strfind(atypes{at},'ages-gt18'))
        vax_id=4; 
    elseif ~isempty(strfind(atypes{at},'ages-gt65'))
        vax_id=5;
    elseif ~isempty(strfind(atypes{at},'janssen'))
        vax_id=6; 
    elseif ~isempty(strfind(atypes{at},'moderna'))
        vax_id=7; 
    elseif ~isempty(strfind(atypes{at},'pfizer'))
        vax_id=8; % 2 is Administered, 
    else
        vax_id=2; % 2 is Administered, 
    end
    
    % Define column indices to extract relevant variables from the data
    % matrices
    cas_id=1; % For total cases (cumulative)
    dth_cell={4,3,[4,8]}; % 3 for COVID deaths, 4 for total deaths, [4,8] for non covid, pneumonia, flue deaths
    pop_id=2; % For 2019 census populations
    
    % Initial indices mapping age groups to respective rows in the death tables 
    % 3  : 0-17 years, all sexes, AK
    % 7  : 18-29 years, all sexes, AK
    % 9  : 30-39 years, all sexes, AK
    % 11 : 40-49 years, all sexes, AK
    % 13 : 50-64 years, all sexes, AK
    % 15 : 65-74 years, all sexes, AK
    % 16 : 75-84 years, all sexes, AK
    % 17 : >85 years, all sexes, AK
    
    % 1-18   : 0-17 years, all sexes, AK
    % 19-30  : 18-29 years, all sexes, AK
    % 31-40  : 30-39 years, all sexes, AK
    % 41-50  : 40-49 years, all sexes, AK
    % 51-65  : 50-64 years, all sexes, AK
    % 66-75  : 65-74 years, all sexes, AK
    % 76-85  : 75-84 years, all sexes, AK
    % 86    : 85 years, all sexes, AK (DON'T USE)
    
    age_groups=[3, 7, 9, 11, 13, 15, 16, 17]; % start indeces for age groups in the dth spreadsheet
    pop_age_idx={[1:18],[19:30],[31:40],[41:50],[51:65],[66:75],[76:85],[86]}; % age groups idxs in the pop sheet
    months={'Jan','Feb','March','April','May','June','July','Aug'};
    age_labels={'0-17','18-29','30-39','40-49','50-64','65-74','75-84','85-plus'};

    for death_type=1:length(dth_cell)
        clear state_vax state_dth state_cas out state_pop
        dth_id=dth_cell{death_type};

        for a=1:length(age_groups) % iterate over all age groups

            % initialize variables that define spreadsheet row structure
            month_gap=51; state_gap=20; 
            vax_starti=1; vax_month_gap=1; vax_state_gap=9;
            cas_starti=1; cas_month_gap=1; cas_state_gap=20;
            pop_starti=1; pop_month_gap=87;
            
            dth_starti=age_groups(a);

            for s=1:length(states)      
                % extract deaths
                state_idx=[dth_starti:month_gap:(state_gap-1)*month_gap+dth_starti]; % Should be length 20
                if length(dth_id)>1
                    state_dth(:,s)=dth.data(state_idx,dth_id(1))-dth.data(state_idx(2),dth_id(2)); % Months in rows, states in columns
                    dth_type='Non-COVFLUPNU-deaths';
                else
                    state_dth(:,s)=dth.data(state_idx,dth_id); % Months in rows, states in columns
                    if dth_id==3
                        dth_type='COVID-deaths';
                    elseif dth_id==4
                        dth_type='Total-deaths';
                    end
                end

                dth_starti=state_idx(end)+month_gap;

                % extract vaccination rates
                vax_idx=[vax_starti:(vax_state_gap-1)+vax_starti];
                if length(vax_id)==3
                    state_vax(:,s)=vax.data(vax_idx,vax_id(1))-vax.data(vax_idx,vax_id(2))-...
                        vax.data(vax_idx,vax_id(3));
                    i=findstr(atypes{at},'-ages-');
                    prepend=[dth_type atypes{at}(i:i+11)];
                elseif length(vax_id)==2
                    state_vax(:,s)=vax.data(vax_idx,vax_id(1))-vax.data(vax_idx,vax_id(2));
                    i=findstr(atypes{at},'-ages-');
                    prepend=[dth_type atypes{at}(i:i+10)];
                else
                    state_vax(:,s)=vax.data(vax_idx,vax_id);
                    prepend=[dth_type '-vax-pop-'];
                end
                if exist('vax_id2')
                    state_vax2(:,s)=vax.data(vax_idx,vax_id2);
                    prepend=[dth_type '-vax-lt18-gt18-'];
                end
                
                vax_starti=vax_idx(end)+vax_month_gap;

                % extract cas
                cas_idx=[cas_starti:(cas_state_gap-1)+cas_starti];
                state_cas(:,s)=cas.data(cas_idx,cas_id);
                cas_starti=cas_idx(end)+cas_month_gap;
                
                % extract age stratified populations                    
                state_pop(1,s)=sum(pop.data(pop_age_idx{a}+(pop_starti-1)*pop_month_gap,pop_id));
                pop_starti=pop_starti+1;
            end
            
            close_figures

            for m=months_idx, % Compute for months April through August
            %for m=7:7
                state_labels=states; state19_pop=state_pop;
            
                % predict deaths in month after the vaccination rate
                dth_vec=(state_dth(m+11,:)-state_dth(m-1,:))./state_dth(m-1,:);

                dth21_abs=(state_dth(m+11,:)); 
                dth20_abs=(state_dth(m-1,:));
                vax_cum=state_vax(m,:);
                
                % Extract vaccination administred for current month and
                % previous month
                vax_cur_vec=state_vax(m+1,:)-state_vax(m,:);
                vax_vec=state_vax(m,:)-state_vax(m-1,:); 
                if ~isempty(findstr(atypes{at},'lag0'))
                    vax_vec=vax_cur_vec;
                end
                if exist('vax_id2') % modify to accommodate same month lag if requested
                    vax_vec2=state_vax2(m,:)-state_vax2(m-1,:);
                end

                % define vectors for the COVID cases
                cas_cum=state_cas(m+11,:); % Should be total cases at start of the month
                cas_abs=state_cas(m+11,:)-state_cas(m+10,:); % Should be new cases previous month 

                % write out table of vax rates and age-stratified death counts for Y20 and Y21       
                fid = fopen(['ages-' age_labels{a} '-month-' months{m}  '.csv'],'w');
                if a==1
                    fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s\n','State','Cases Cum','Cases Change','Vax Cum',...
                        'Vax Change',['Y19 pop ages:' age_labels{a}],['Y20: ages ' age_labels{a}],['Y21: ages ' age_labels{a}]);
                    for s=1:length(state_labels)
                        fprintf(fid,'%s,%i,%i,%0.2f,%0.2f,%i,%i,%i\n',state_labels{s},cas_cum(s),cas_abs(s),...
                            vax_cum(s),vax_vec(s),state_pop(s),dth20_abs(s),dth21_abs(s));
                    end
                elseif a==length(age_groups)
                    fprintf(fid,'%s,%s,%s\n',['Y19 pop ages:' age_labels{a}],['Y20: ages ' age_labels{a}],['Y21: ages ' age_labels{a}]);
                    for s=1:length(state_labels)
                        fprintf(fid,'%i,%i,%i\n',state_pop(s),dth20_abs(s),dth21_abs(s));
                    end
                else
                    fprintf(fid,'%s,%s,%s\n',['Y19 pop ages:' age_labels{a}],['Y20: ages ' age_labels{a}],['Y21: ages ' age_labels{a}]);
                    for s=1:length(state_labels)
                        fprintf(fid,'%i,%i,%i\n',state_pop(s),dth20_abs(s),dth21_abs(s));
                    end
                end
                fclose(fid);

                % take log
                dth20_log=log(dth20_abs); dth21_log=log(dth21_abs); cas_log=log(cas_abs); 
                vax_log=log(vax_vec); pop_log=log(state19_pop); vax_cur_log=log(vax_cur_vec);
                % Store vectors for estimating deaths later on
                vax_full=vax_vec; dth20_full=dth20_abs;
                
                if exist('vax_id2')
                    vax_log2=log(vax_vec2);
                end
                
                % remove negative values where not allowed, NaN and inf for
                % log transform
                rm_ind=find(dth20_abs<0 | dth21_abs<0 | cas_abs<0 | vax_vec<0 |...
                            isnan(dth20_abs) | isnan(dth21_abs) | isnan(cas_abs) | isnan(vax_vec) | isnan(dth_vec) |...
                            isinf(dth20_log) | isinf(dth21_log) | isinf(cas_log) | isinf(vax_log));
                dth20_abs(rm_ind)=[]; dth21_abs(rm_ind)=[]; cas_abs(rm_ind)=[]; 
                vax_vec(rm_ind)=[]; dth_vec(rm_ind)=[]; state_labels(rm_ind)=[];
                dth20_log(rm_ind)=[]; dth21_log(rm_ind)=[]; cas_log(rm_ind)=[]; vax_log(rm_ind)=[];
                vax_cur_log(rm_ind)=[];
                state19_pop(rm_ind)=[]; pop_log(rm_ind)=[]; 
                if exist('vax_id2')
                    vax_vec2(rm_ind)=[]; vax_log2(rm_ind)=[];
                end
      
                fname  = [dth_type '-' atypes{at} '-' months{m} '-ages-' age_labels{a}]; 
                % plot histograms
                if make_figs
                    figure
                    %subplot(2,1,1); hist(dth_vec); subplot(2,1,2); hist(vax_vec)
                    subplot(3,1,1); hist(dth20_log); 
                    subplot(3,1,2); hist(dth21_log); 
                    subplot(3,1,3); hist(vax_log)
                    saveas(gca,[fname '-HIST.png']);
                end
                
                ant=atypes{at};
                            
                if ~isempty(dth21_abs)
                    switch(ant)
                        
                        case {'vax-adj-cas-log-rob','vax-adj-cas-log-rob-lag0','vax-adj-cas-ages-gt65-log-rob'}
                            disp(atypes{at})
                            try 
                                [b,stats]=robustfit([dth20_log',cas_log',vax_log'],[dth21_log']);
                                [bres,resid]=robustfit([dth20_log',cas_log'],[dth21_log']);
                            catch 
                                [b,dev,stats]=glmfit([dth20_log',cas_log',vax_log'],[dth21_log']);
                                [~,~,resid]=glmfit([dth20_log',cas_log'],[dth21_log']);
                            end
                            xvec=vax_log'; yvec=resid.resid; robust=1;
                            if ~isempty(strfind(atypes{at},'lag0')) 
                                xlab   = [months{m} ': vaccination rate (log(#) administered)'];
                            else
                                xlab   = [months{m-1} ': vaccination rate (log(#) administered)'];
                            end
                            ylab   = [months{m} ': residual log(Y21) deaths (vs. Y20), ages ' age_labels{a}];
                            ftitle = ['Vax rates vs. ' dth_type ': tval=' num2str(stats.t(end),'%0.2f') ', beta=' num2str(b(end),'%0.4f') ', p=' num2str(stats.p(end),'%0.5f')];
                        case {'vax-adj-cas-log-glm','vax-adj-cas-log-glm-lag0','vax-adj-cas-gt65-log-glm'}
                            disp(atypes{at})

                            [b,dev,stats]=glmfit([dth20_log',cas_log',vax_log'],[dth21_log']);
                            [~,~,resid]=glmfit([dth20_log',cas_log'],[dth21_log']);
                        
                            xvec=vax_log'; yvec=resid.resid; robust=0;
                            if ~isempty(strfind(atypes{at},'lag0')) 
                                xlab   = [months{m} ': vaccination rate (log(#) administered)'];
                            else
                                xlab   = [months{m-1} ': vaccination rate (log(#) administered)'];
                            end
                            ylab   = [months{m} ': residual log(Y21) deaths (vs. Y20), ages ' age_labels{a}];
                            ftitle = ['Vax rates vs. ' dth_type ': tval=' num2str(stats.t(end),'%0.2f') ', beta=' num2str(b(end),'%0.4f') ', p=' num2str(stats.p(end),'%0.5f')];
                        case 'cas-adj-vax-log-rob'
                            disp(atypes{at})
                            try
                                [b,stats]=robustfit([dth20_log',vax_log',cas_log'],[dth21_log']);
                                [bres,resid]=robustfit([dth20_log',cas_log'],[dth21_log']);
                            catch
                                [b,~,stats]=glmfit([dth20_log',vax_log',cas_log'],[dth21_log']);
                                [bres,~,resid]=glmfit([dth20_log',cas_log'],[dth21_log']);
                            end
                            if isstruct(resid)
                                xvec=cas_log'; yvec=resid.resid; 
                            else
                                xvec=NaN; yvec=NaN;
                            end
                            robust=1;
                            xlab   = [months{m-1} ': log(#) COVID cases'];
                            ylab   = [months{m} ': residual log(#) deaths (vs. Y20), ages ' age_labels{a}];
                            ftitle = ['Cases vs. ' dth_type ': tval=' num2str(stats.t(end),'%0.2f') ', beta=' num2str(b(end),'%0.4f') ', p=' num2str(stats.p(end),'%0.5f')];
                        case 'cas-abs-log-rob'
                            disp(atypes{at})
                            try 
                                [b,stats]=robustfit([dth20_log',cas_log'],[dth21_log']);
                                [bres,resid]=robustfit([dth20_log'],[dth21_log']);
                            catch 
                                [b,dev,stats]=glmfit([dth20_log',cas_log'],[dth21_log']);
                                [~,~,resid]=glmfit([dth20_log'],[dth21_log']);
                            end
                            xvec=cas_log'; yvec=resid.resid; robust=1;
                            xlab   = [months{m-1} ': log(#) COVID cases'];
                            ylab   = [months{m} ': residual log(#) deaths (vs. Y20), ages ' age_labels{a}];
                            ftitle = ['Cases vs. ' dth_type ': tval=' num2str(stats.t(end),'%0.2f') ', beta=' num2str(b(end),'%0.4f') ', p=' num2str(stats.p(end),'%0.5f')];
                        case {'vax-admin-log-rob','vax-admin-log-rob-lag0','vax-admin-ages-gt65-log-rob'}
                            disp(atypes{at})
                            try
                                [b,stats]=robustfit([dth20_log',vax_log'],[dth21_log']);
                                [bres,resid]=robustfit([dth20_log'],[dth21_log']);
                            catch
                                [b,~,stats]=glmfit([dth20_log',vax_log'],[dth21_log']);
                                [bres,~,resid]=glmfit([dth20_log'],[dth21_log']);
                            end
                            if isstruct(resid)
                                xvec=vax_log'; yvec=resid.resid; 
                            else
                                xvec=NaN; yvec=NaN;
                            end
                            robust=1;
                            if ~isempty(strfind(atypes{at},'lag0')) 
                                xlab   = [months{m} ': vaccination rate (log(#) administered)'];
                            else
                                xlab   = [months{m-1} ': vaccination rate (log(#) administered)'];
                            end
                            ylab   = [months{m} ': residual log(Y21) deaths (vs. Y20), ages ' age_labels{a}];
                            ftitle = ['Vax rates vs. ' dth_type ': tval=' num2str(stats.t(end),'%0.2f') ', beta=' num2str(b(end),'%0.4f') ', p=' num2str(stats.p(end),'%0.5f')];        
                        case {'vax-admin-log-glm','vax-admin-log-glm-lag0','vax-admin-ages-gt65-log-glm'}
                            disp(atypes{at})
                            [b,~,stats]=glmfit([dth20_log',vax_log'],[dth21_log']);
                            [bres,~,resid]=glmfit([dth20_log'],[dth21_log']);
                            if isstruct(resid)
                                xvec=vax_log'; yvec=resid.resid; 
                            else
                                xvec=NaN; yvec=NaN;
                            end
                            robust=0;
                            if ~isempty(strfind(atypes{at},'lag0'))
                                xlab   = [months{m} ': vaccination rate (log(#) administered)'];
                            else
                                xlab   = [months{m-1} ': vaccination rate (log(#) administered)'];
                            end
                            ylab   = [months{m} ': residual log(Y21) deaths (vs. Y20), ages ' age_labels{a}];
                            ftitle = ['Vax rates vs. ' dth_type ': tval=' num2str(stats.t(end),'%0.2f') ', beta=' num2str(b(end),'%0.4f') ', p=' num2str(stats.p(end),'%0.5f')];
                        otherwise
                            disp('Invalid case')
                    end

                    % plot and save out the graphs
                    try 
                        if make_figs
                            figure
                            plot_fit(xvec,yvec,robust,'k.','r');
                            xlabel(xlab,'fontSize',15); ylabel(ylab,'fontSize',15); title(ftitle,'fontSize',15); 
                            saveas(gca,[fname '.png']);
                        end
                        
                        out.corr(a,m)=b(end); 
                        out.pval(a,m)=stats.p(end);
                        
                        % Compute and store correlation between age
                        % stratified state population and 
                        %[r,p]=corr(state19_pop',yvec,'type','spearman');
                        %pop_res.corr(a,m)=r; pop_res.pval(a,m)=p;    
                    catch
                        
                        out.corr(a,m)=NaN; out.pval(a,m)=NaN;
                    end
                    % Estimate deaths attributed to vaccine
                    % yfit=b(1)+b(2)*dth20_log+b(3)*vax_log;
                    clear numD numVax
                    % remove NaN or missing values or zero
                    rm_full=find(isnan(dth20_full) | dth20_full==0);
                    dth20_full(rm_full)=[]; vax_full(rm_full)=[];
                    FF=0.1; % Sample at 10% increases
                    for i=1:length(dth20_full)
                        Y1_log=b(1)+b(2)*log(dth20_full(i))+b(3)*log(vax_full(i));
                        Y1=exp(Y1_log);
                        Y2=Y1*exp(b(3)*log(1+FF));
                        numD(i)=Y2-Y1; numVax(i)=FF*vax_full(i);   
                    end
                    
                    rate=sum(numD)/sum(numVax);
                    out.sum_vax(1,m)=sum(vax_full);
                    out.deaths{death_type}(a,m)=rate*sum(vax_full);
                    out.perc_rate{death_type}(a,m)=rate*100;
                else        
                    out.corr(a,m)=NaN; out.pval(a,m)=NaN;
                end
                
            end 
        end

        % Write out the result table
        resultsFileName=[prepend 'results.csv'];
        fid = fopen(resultsFileName,'w');
        fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n','Ages',...
            [months{months_idx(1)} ' (es)'],[months{months_idx(1)} ' (pval)'],...
            [months{months_idx(2)} ' (es)'],[months{months_idx(2)} ' (pval)'],...
            [months{months_idx(3)} ' (es)'],[months{months_idx(3)} ' (pval)'],...
            [months{months_idx(4)} ' (es)'],[months{months_idx(4)} ' (pval)'],...
            [months{months_idx(5)} ' (es)'],[months{months_idx(5)} ' (pval)'],...
            [months{months_idx(6)} ' (es)'],[months{months_idx(6)} ' (pval)'],...
            [months{months_idx(7)} ' (es)'],[months{months_idx(7)} ' (pval)']);

        % Write row to table
        for a=1:length(age_labels)
            %fprintf(fid,'%s,%0.2f,%0.4f,%0.2f,%0.4f,%0.2f,%0.4f,%0.2f,%0.4f,%0.2f,%0.4f\n',...
            fprintf(fid,'%s,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f\n',...
                        age_labels{a},...
                        out.corr(a,months_idx(1)),out.pval(a,months_idx(1)),...
                        out.corr(a,months_idx(2)),out.pval(a,months_idx(2)),...
                        out.corr(a,months_idx(3)),out.pval(a,months_idx(3)),...
                        out.corr(a,months_idx(4)),out.pval(a,months_idx(4)),...
                        out.corr(a,months_idx(5)),out.pval(a,months_idx(5)),...
                        out.corr(a,months_idx(6)),out.pval(a,months_idx(6)),...
                        out.corr(a,months_idx(7)),out.pval(a,months_idx(7)));
        end
        fclose(fid);
        
        % Write out the result table with adjusted p-values
        % run FDR correction 
        [h crit_p adj_p]=fdr_bh([out.pval(:,months_idx(1):months_idx(end))]);
        out.adj_pval=[zeros(size(out.pval,1),months_idx(1)-1) adj_p];
        
        % Mask out deaths if less than p<0.05 FDR corrected. Can also change to use p<0.05
        % uncorrected
        out.deaths{death_type}(out.adj_pval>0.05)=NaN;
        out.perc_rate{death_type}(out.adj_pval>0.05)=NaN;
        
        resultsFileName=[prepend 'adj_pval_results.csv'];
        fid = fopen(resultsFileName,'w');
        fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n','Ages',...
            [months{months_idx(1)} ' (es)'],[months{months_idx(1)} ' (pval)'],...
            [months{months_idx(2)} ' (es)'],[months{months_idx(2)} ' (pval)'],...
            [months{months_idx(3)} ' (es)'],[months{months_idx(3)} ' (pval)'],...
            [months{months_idx(4)} ' (es)'],[months{months_idx(4)} ' (pval)'],...
            [months{months_idx(5)} ' (es)'],[months{months_idx(5)} ' (pval)'],...
            [months{months_idx(6)} ' (es)'],[months{months_idx(6)} ' (pval)'],...
            [months{months_idx(7)} ' (es)'],[months{months_idx(7)} ' (pval)']);

        % Write row to table
        for a=1:length(age_labels)
            %fprintf(fid,'%s,%0.2f,%0.4f,%0.2f,%0.4f,%0.2f,%0.4f,%0.2f,%0.4f,%0.2f,%0.4f\n',...
            fprintf(fid,'%s,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f,%0.4f\n',...
                        age_labels{a},...
                        out.corr(a,months_idx(1)),out.adj_pval(a,months_idx(1)),...
                        out.corr(a,months_idx(2)),out.adj_pval(a,months_idx(2)),...
                        out.corr(a,months_idx(3)),out.adj_pval(a,months_idx(3)),...
                        out.corr(a,months_idx(4)),out.adj_pval(a,months_idx(4)),...
                        out.corr(a,months_idx(5)),out.adj_pval(a,months_idx(5)),...
                        out.corr(a,months_idx(6)),out.adj_pval(a,months_idx(6)),...
                        out.corr(a,months_idx(7)),out.adj_pval(a,months_idx(7)));
        end
        fclose(fid);
        
        % Write out the death estimate tables
        resultsFileName=[prepend 'deaths.csv'];
        fid = fopen(resultsFileName,'w');
        fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s\n','Ages',months{months_idx(1)},months{months_idx(2)},...
            months{months_idx(3)},months{months_idx(4)},months{months_idx(5)},months{months_idx(6)},...
            months{months_idx(7)});

        % Write row to top table
        for a=1:length(age_labels)
            %fprintf(fid,'%s,%0.2f,%0.4f,%0.2f,%0.4f,%0.2f,%0.4f,%0.2f,%0.4f,%0.2f,%0.4f\n',...
            fprintf(fid,'%s,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f\n',...
                        age_labels{a},out.deaths{death_type}(a,months_idx(1)),out.deaths{death_type}(a,months_idx(2)),out.deaths{death_type}(a,months_idx(3)),...
                        out.deaths{death_type}(a,months_idx(4)),out.deaths{death_type}(a,months_idx(5)),out.deaths{death_type}(a,months_idx(6)),out.deaths{death_type}(a,months_idx(7)))
                        
        end
        % Write percentages to bottom table
        for a=1:length(age_labels)
            %fprintf(fid,'%s,%0.2f,%0.4f,%0.2f,%0.4f,%0.2f,%0.4f,%0.2f,%0.4f,%0.2f,%0.4f\n',...
            fprintf(fid,'%s,%0.5f,%0.5f,%0.5f,%0.5f,%0.5f,%0.5f,%0.5f\n',...
                        age_labels{a},out.perc_rate{death_type}(a,months_idx(1)),out.perc_rate{death_type}(a,months_idx(2)),out.perc_rate{death_type}(a,months_idx(3)),...
                        out.perc_rate{death_type}(a,months_idx(4)),out.perc_rate{death_type}(a,months_idx(5)),out.perc_rate{death_type}(a,months_idx(6)),out.perc_rate{death_type}(a,months_idx(7)))
                        
        end
        % Add the total vax (from vax_full) 
        fprintf(fid,'%s,%i,%i,%i,%i,%i,%i,%i\n','Total vax',out.sum_vax(months_idx(1)),out.sum_vax(months_idx(2)),...
            out.sum_vax(months_idx(3)),out.sum_vax(months_idx(4)),out.sum_vax(months_idx(5)),out.sum_vax(months_idx(6)),out.sum_vax(months_idx(7)))
        fclose(fid);

        % combine the csv files into one table for each month
        for m=months_idx %size(state_vax,1)

            outputFileName=[prepend 'month-' months{m} '-full-table.csv'];
            allCsv=[];
            for a=1:length(age_groups)
                if a==1 
                    startc=1;
                else 
                    startc=0;
                end
                tmp = csvread(['ages-' age_labels{a} '-month-' months{m}  '.csv'],1,startc);
                allCsv = [allCsv tmp];
            end
            
            if size(allCsv,2)==10
                tmp=nan(size(allCsv,1),10);
                allCsv=[tmp allCsv];
            end
            
            fid = fopen(outputFileName,'w');
            fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n','State','Cases Cum','Cases Change','Vax Cum','Vax Change',...
                ['Y19 pop: ages ' age_labels{1}],['Y20: ages ' age_labels{1}],['Y21: ages ' age_labels{1}],...
                ['Y19 pop: ages ' age_labels{2}],['Y20: ages ' age_labels{2}],['Y21: ages ' age_labels{2}],...
                ['Y19 pop: ages ' age_labels{3}],['Y20: ages ' age_labels{3}],['Y21: ages ' age_labels{3}],...
                ['Y19 pop: ages ' age_labels{4}],['Y20: ages ' age_labels{4}],['Y21: ages ' age_labels{4}],...
                ['Y19 pop: ages ' age_labels{5}],['Y20: ages ' age_labels{5}],['Y21: ages ' age_labels{5}],....
                ['Y19 pop: ages ' age_labels{6}],['Y20: ages ' age_labels{6}],['Y21: ages ' age_labels{6}],...
                ['Y19 pop: ages ' age_labels{7}],['Y20: ages ' age_labels{7}],['Y21: ages ' age_labels{7}],...
                ['Y19 pop: ages ' age_labels{8}],['Y20: ages ' age_labels{8}],['Y21: ages ' age_labels{8}]);

            for s=1:length(states)
                fprintf(fid,'%s,%i,%i,%0.2f,%0.2f,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i\n',states{s},...
                    allCsv(s,1),allCsv(s,2),allCsv(s,3),allCsv(s,4),allCsv(s,5),allCsv(s,6),allCsv(s,7),...
                    allCsv(s,8),allCsv(s,9),allCsv(s,10),allCsv(s,11),allCsv(s,12),allCsv(s,13),allCsv(s,14),...
                    allCsv(s,15),allCsv(s,16),allCsv(s,17),allCsv(s,18),allCsv(s,19),allCsv(s,20),...
                    allCsv(s,21),allCsv(s,22),allCsv(s,23),allCsv(s,24),allCsv(s,25),allCsv(s,26),allCsv(s,27),allCsv(s,28));
            end
            fclose(fid);
        end

    delete('ages-*.csv');
    end % end for loop for death types
    
    table_fol=['./Tables/' atypes{at}]; mkdir(table_fol)
    pause(5)
    movefile('*full-table.csv', table_fol); 
    results_fol=['./Results/' atypes{at}]; mkdir(results_fol)
    pause(5)
    movefile('*results.csv', results_fol);
    pause(5)
    movefile('*deaths.csv', results_fol);
    if make_figs
        fig_fol=['./Figures/' atypes{at}]; mkdir(fig_fol)
        pause(5)
        movefile('*ages-*.png', fig_fol);
    end
    
end

%% Here do sign tests on the EU lagged correlations and apply FDR correction
clear
dat=importdata('lagr_EU_SupTable1.csv') % Same table as Supplementary Table 3

dat(:,2:7)=dat(:,2:7)./100;
g={[0,1,2],[3,4,5],[6,7,8],[9,10,11],[12,13,14],[15,16,17],[18,19,20],...
    [21,22,23,24],[25,26,27,28]};

%g={[0,1,2,3],[4,5,6,7],[8,9,10,11],[12,13,14,15],[16,17,18,19],[20,21,22,23],[24,25,26,27]}

for lag=1:length(g)
    %tmp=dat(find(dat(:,1)==lag),:);
    tmp=dat(ismember(dat(:,1),g{lag}),:);
    
    for a=1:6
        [pval,~,out.stats(a,lag).stats]=signtest(tmp(:,1+a));
        
        out.signp(a,lag)=pval;
    end
    
end
[h crit_p out.adj_p]=fdr_bh(out.signp,.05,'pdep','yes');

%[h crit_p out.adj_p]=fdr_bh(out.signp);
    

   


    