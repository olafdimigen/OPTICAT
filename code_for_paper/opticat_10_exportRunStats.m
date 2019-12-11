% generate output for statistics in "R"
% olaf.dimigen@hu-berlin.de, 2018-06-10

clear, clc
load Z:/OPTICA/results/opticaresults_20180604.mat % used only for trial count information

% metric_sac(1,SUBJECTS,2,:,1,THRESHOLD2PLOT))
line1 = 1; line2 = 1; line = 1;
THRESHOLD = 7; % corrected to "1.1" = 7 on August 7

% % go tru factors
% for exp = 1:2
%     for s = 1:12
%         for hc_level = 1:2
%             for lc_level = 1:20
%                 for ow_level = 1:2
% 
%                     all_results(line,1) = exp;
%                     all_results(line,2) = (exp-1)*12 + s; % subject ID (between-subject factor)
%                     all_results(line,3) = hc_level;
%                     all_results(line,4) = lc_level;
%                     all_results(line,5) = ow_level;
%                     all_results(line,6) = opticaresults(exp).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).thresh(THRESHOLD).metric_sac;   % CR
%                     %all_results(line,7) = opticaresults(exp).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).thresh(THRESHOLD).metric_sp2;   % SP
%                     all_results(line,7) = opticaresults(exp).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).thresh(THRESHOLD).metric_sp;   % SP
%                     all_results(line,8) = opticaresults(exp).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).thresh(THRESHOLD).metric_stim;  % Overcorr.
%                     line = line+1;
%                     
%                 end % ow
%             end % lc
%         end % hc
%     end % s
% end % exp
% 
% 
% %% Control plot
% for exp = 1:2
%     for lc = 1:20
%         for hc = 1:2
%             for ow = 1:2
%                 ix = find( all_results(:,1) == exp & all_results(:,3) == hc & all_results(:,4) == lc & all_results(:,5) == ow);
%                 controldata(exp,lc,hc,ow) = mean(all_results(ix,7));
%                 
%             end
%         end
%     end
% end
% 
% %% control plot
% figure 
% subplot(1,2,1); hold on; title('SP: Scenes')
% plot(1:20, controldata(1,:,1,1),'ro-')
% plot(1:20, controldata(1,:,1,2),'go-')
% plot(1:20, controldata(1,:,2,1),'bo-')
% plot(1:20, controldata(1,:,2,2),'co-')
% ylim([0 6]) % for metric_sp     --> GFP in window 
% % ylim([0 10]) % for metric_sp2 --> GFP at SP peak
% 
% subplot(1,2,2);  hold on; title('SP: Reading')
% plot(1:20, controldata(2,:,1,1),'ro-')
% plot(1:20, controldata(2,:,1,2),'go-')
% plot(1:20, controldata(2,:,2,1),'bo-')
% plot(1:20, controldata(2,:,2,2),'co-')
% ylim([0 6])
% 
% dlmwrite('statistic/Statistics_20181004_asInPaper.txt', all_results);
% 
% 
% %% now also include the used threshold
% line = 1;
% % go tru factors
% for exp = 1:2
%     for s = 1:12
%         for hc_level = 1:2
%             for lc_level = 1:20
%                 for ow_level = 1:2
%                     for thresh = 2:11 % 0.6 to 1.5
%                         
%                         all_results_withThreshold(line, 1) = exp;
%                         all_results_withThreshold(line, 2) = (exp-1)*12 + s;
%                         all_results_withThreshold(line, 3) = hc_level;
%                         all_results_withThreshold(line, 4) = lc_level;
%                         all_results_withThreshold(line, 5) = ow_level;
%                         all_results_withThreshold(line, 6) = thresh;
%                         all_results_withThreshold(line, 7) = opticaresults(exp).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).thresh(THRESHOLD).metric_sac;
%                         all_results_withThreshold(line, 8) = opticaresults(exp).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).thresh(THRESHOLD).metric_sp;
%                         all_results_withThreshold(line, 9) = opticaresults(exp).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).thresh(THRESHOLD).metric_stim;
%                         all_results_withThreshold(line,10) = opticaresults(exp).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).thresh(THRESHOLD).nbadcomps;
%                         line = line+1;
%                         
%                     end % ow
%                 end % lc
%             end % hc
%         end % s
%     end % exp
% end
% dlmwrite('statistic/Statistics_20181004_asInPaper_Threshold.txt', all_results_withThreshold);
% 
%  
% 
% % consider threshold, but only at 2 Hz, 100 Hz, with overweighting (i.e. settings from figure)
% line = 1;
% 
% hc_level = 2;
% lc_level = 8;
% ow_level = 2;
% 
% % go tru factors
% for exp = 1:2
%     for s = 1:12
%         for thresh = 2:11 % 0.6 to 1.5
%             
%             all_results_withThreshold_minimal(line, 1) = exp;
%             all_results_withThreshold_minimal(line, 2) = (exp-1)*12 + s;
%             all_results_withThreshold_minimal(line, 3) = hc_level;
%             all_results_withThreshold_minimal(line, 4) = lc_level;
%             all_results_withThreshold_minimal(line, 5) = ow_level;
%             all_results_withThreshold_minimal(line, 6) = thresh;
%             all_results_withThreshold_minimal(line, 7) = opticaresults(exp).subj(s).hc(hc_level).lc(lc_level).ow(ow_level).thresh(THRESHOLD).metric_stim;
%             line = line+1;
%         end % thresh
%     end % s
% end % exp
%  
% dlmwrite('statistic/Statistics_20181004_asInPaper_Threshold_minimal.txt', all_results_withThreshold_minimal);
% fprintf('\nStats exported as txt.')


%% compare best ICA solutions to MSEC solution
% do this for Scenes and for Reading, separately
% use simple post-hoc paired t-tests
LC_LEVEL = 8 % 2 Hz

for exp = 1:2
    for s = 1:12
        cr_bestica(s,exp)   = opticaresults(exp).subj(s).hc(2).lc(LC_LEVEL).ow(2).thresh(7).metric_sac;
        sp_bestica(s,exp)   = opticaresults(exp).subj(s).hc(2).lc(LC_LEVEL).ow(2).thresh(7).metric_sp;
        stim_bestica(s,exp) = opticaresults(exp).subj(s).hc(2).lc(LC_LEVEL).ow(2).thresh(7).metric_stim;
                
        cr_msec(s,exp)      = opticaresults(exp).subj(s).msec.metric_sac;  % cr
        sp_msec(s,exp)      = opticaresults(exp).subj(s).msec.metric_sp;   % sp
        stim_msec(s,exp)    = opticaresults(exp).subj(s).msec.metric_stim; % overcorr
    end
end


%%

% scenes
[h, p1, ci1, stats1] = ttest(cr_msec(:,1),cr_bestica(:,1));   % cr
[mean(cr_msec(:,1)) mean(cr_bestica(:,1)) mean(cr_msec(:,1)-cr_bestica(:,1)) ci1']
[p1, stats1.tstat]
% ci1a = get95CIs(cr_bestica(:,1))  % CIs (unsubtracted)
% ci1b = get95CIs(cr_msec(:,1))      % CIs (unsubtracted)

[h, p2, ci2, stats2] = ttest(sp_msec(:,1),sp_bestica(:,1));   % sp
[mean(sp_msec(:,1)) mean(sp_bestica(:,1)) mean(sp_msec(:,1)-sp_bestica(:,1)) ci2']
[p2, stats2.tstat]
% ci2a = get95CIs(sp_bestica(:,1))  % CIs (unsubtracted)
% ci2b = get95CIs(sp_msec(:,1))      % CIs (unsubtracted)

[h, p3, ci3, stats3] = ttest(stim_msec(:,1),stim_bestica(:,1));  % overcorr
[mean(stim_msec(:,1)) mean(stim_bestica(:,1)) mean(stim_msec(:,1)-stim_bestica(:,1)) ci3']
[p3, stats3.tstat]
% ci3a = get95CIs(stim_bestica(:,1)) % CIs (unsubtracted)
% ci3b = get95CIs(stim_msec(:,1))     % CIs (unsubtracted)

% reading
[h, p4, ci4, stats4] = ttest(cr_msec(:,2),cr_bestica(:,2));   
[mean(cr_msec(:,2)) mean(cr_bestica(:,2)) mean(cr_msec(:,2)-cr_bestica(:,2)) ci4']
[p4, stats4.tstat] % cr
% ci4a = get95CIs(cr_bestica(:,2))  % CIs (unsubtracted)
% ci4b = get95CIs(cr_msec(:,2))     % CIs (unsubtracted)

[h, p5, ci5, stats5] = ttest(sp_msec(:,2),sp_bestica(:,2)); % sp
[mean(sp_msec(:,2)) mean(sp_bestica(:,2)) mean(sp_msec(:,2)-sp_bestica(:,2)) ci5']
[p5, stats5.tstat]
% ci5a = get95CIs(sp_bestica(:,2))  % CIs (unsubtracted)
% ci5b = get95CIs(sp_msec(:,2))     % CIs (unsubtracted)

[h, p6, ci6, stats6] = ttest(stim_msec(:,2),stim_bestica(:,2)); % overcorr
[mean(stim_msec(:,2)) mean(stim_bestica(:,2)) mean(stim_msec(:,2)-stim_bestica(:,2)) ci6']
[p6, stats6.tstat]
% ci6a = get95CIs(stim_bestica(:,2)) % CIs (unsubtracted)
% ci6b = get95CIs(stim_msec(:,2))    % CIs (unsubtracted)
 

 
%% compare threshold of 1.5 to threshold of 0.9 or 1.0. for "best" ICA solution
% do this for Scenes and for Reading, separately; use simple t-test
LC_LEVEL = 8 % is this the "best" LC setting?

for exp = 1:2 % exp loop
    for s = 1:12 % subject loop
        % get values for near-optimal ICA (2-100 Hz, with oveweighting) at...
        
        % ...threshold of 0.9
        overcorr_0point9(s,exp) = opticaresults(exp).subj(s).hc(2).lc(LC_LEVEL).ow(2).thresh(5).metric_stim;  % threshold level 5 = 0.9
        
        % ...threshold of 1.0
        overcorr_1point0(s,exp) = opticaresults(exp).subj(s).hc(2).lc(LC_LEVEL).ow(2).thresh(6).metric_stim;  % threshold level 6 = 1.0
        
        % ...threshold of 1.5 (lenient reference threshold)
        overcorr_1point5(s,exp) = opticaresults(exp).subj(s).hc(2).lc(LC_LEVEL).ow(2).thresh(11).metric_stim; % threshold level 11 = 1.5
    end
end

% show correponding t-test results to quantify pattern already seen in 
% the 95% confidence interval plotted in Figure 8. (Side note: Looking for 
% the first threshold that does not produce "significant" overcorrection of 
% course does not imply that the null hypotheses is correct at this threshold)


% exp 1: scenes: 0.9 vs. 1.5
[h1 p1 ci1 stats1] = ttest(overcorr_0point9(:,1),overcorr_1point5(:,1)); 
p1, stats1.tstat

% exp 1: scenes: 1.0 vs. 1.5
[h2 p2 ci2 stats2] = ttest(overcorr_1point0(:,1),overcorr_1point5(:,1)); 
p2, stats2.tstat % *significant

% exp 2: reading: 0.9 vs. 1.5
[h3 p3 ci3 stats3] = ttest(overcorr_0point9(:,2),overcorr_1point5(:,2)); 
p3, stats3.tstat % *significant

% exp 2: reading: 1.0 vs. 1.5
[h4 p4 ci4 stats4] = ttest(overcorr_1point0(:,2),overcorr_1point5(:,2)); 
p4, stats4.tstat



%% to get 95% confidence interval around threshold of 1.5 (level 11)
for exp = 1:2
    for s = 1:12
        overcorr_1point0(s,exp) = opticaresults(exp).subj(s).hc(2).lc(LC_LEVEL).ow(2).thresh(11).metric_stim;
    end
end
ALPHA = 0.05;
% exp 1: scenes
[~,~,ci_crit1,stats] = ttest(overcorr_1point0(:,1),0,ALPHA); ci_crit1
% exp 2: reading
[~,~,ci_crit2,stats] = ttest(overcorr_1point0(:,2),0,ALPHA); ci_crit2


%% Add this to top line of output statistics text files for import to "R" 
% (quick fix, this should be fprintf()'ed properly into the text file as a header line...)
% Exp,Subject,HC,LC,OW,DV_CR,DV_SP,DV_STIM
% Exp,Subject,HC,LC,OW,THRESH,DV_CR,DV_SP,DV_STIM,DV_NUMIC
% Exp,Subject,HC,LC,OW,THRESH,DV_STIM

