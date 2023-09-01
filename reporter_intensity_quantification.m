%% Koleske Lab
% Last updated September 1, 2023
% reporter_intensity_quantification.m
%
% DESCRIPTION: Quantify the intensity of tubulin reporter
% on 23C and 37C storage treated biotinylated rhodamine
% taxol-stabilized MTs to determine whether Abl2 assists in recruiting tubulin to damaged sites
%
% INPUT: 
% 1) filaments: total number of MTs
% 2) mt_tot: array of MT lengths and respective MT intensity per unit of length (left column = unit of length measurement; right column = intensity of MT at that particular locus of filament)
% 3) tub_tot: array of MT lengths and respective tubulin reporter intensity per unit of length (left column = unit of MT length measurement; right column = intensity of Alexa647tub at that particular locus of filament)
% 4) threshold: the normalized intensity value that user denotes whether a segment is damaged or not (default value: 0.45)
% 5) rep_threshold: the normalized intensity value of Alexa647 tubulin reporter that user denotes whether reporter is present or not (above background noise; default value: 0.4)
% --- NOTE: the number of filaments in both the mt_tot and tub_tot (every
% other 2 columns in the mt_tot and tub_tot arrays) should match ---
% 
% OUTPUT: 
% 1) reporter_fraction: length of tubulin reporter divided by overall length of damaged MT
% 2) avg_reporter_length: the length (microns) of reporter per damaged MT
% 3) count: number of reporter events within a 23C/37C O/N stored damaged MT
% 4) avg_reporter_intensities: mean intensities of tubulin reporter per damaged MT

function [reporter_fraction,avg_reporter_length,count,avg_reporter_intensities] = ...
reporter_intensity_quantification(filaments,mt_tot,tub_tot,threshold,rep_threshold)

reporter_fraction = [];
avg_reporter_length = [];
count = [];
avg_reporter_intensities = [];


for i=1:2:(filaments-1)
    reporter_length = []; 
    mt=rmmissing(mt_tot(:,i:i+1)); 
    tub=rmmissing(tub_tot(:,i:i+1)); 
   
    try
    mt_norm = [mt(:,1) (mt(:,2)-min(mt(:,2)))./(max(mt(:,2))-min(mt(:,2)))];
    tub_norm = [tub(:,1) (tub(:,2)-min(tub(:,2)))./(max(tub(:,2))-min(tub(:,2)))];

    ind_dam = find(mt_norm(:,2) < threshold);

    if isempty(ind_dam) == 0
        a=diff(ind_dam);
        b=find(a ~= 1);

        if isempty(b) == 1 
            reporter_intensity = mean(tub_norm(ind_dam(1:end),2));
            if reporter_intensity >= rep_threshold
                avg_reporter_intensities = [avg_reporter_intensities; reporter_intensity];
                reporter_length = tub_norm(ind_dam(end),1)-tub_norm(ind_dam(1),1); % reporting this twice for reporter_fraction quantificaton (see end of loop)
                avg_reporter_length = [avg_reporter_length; tub_norm(ind_dam(end),1)-tub_norm(ind_ok(1),1)]; % still need to retain for outside the loop (output of the function)
                count = [count; 1./(mt(end,1)-mt(1,1))];
            end
        else
            if length(b) == 1 
                reporter_intensity = mean([tub_norm(ind_dam(b):a(b)+ind_dam(b),2) tub_norm(ind_dam(b+1):ind_dam(end),2)]);

                if reporter_intensity >= rep_threshold
                    total = tub_norm(a(b)+ind_dam(b),1)-tub_norm(ind_dam(b),1);
                    total = total + (tub_norm(ind_dam(end),1) - tub_norm(ind_dam(b+1),1))
                    reporter_length = total;
                    if reporter_length ~= 0
                        avg_reporter_length = [avg_reporter_length; reporter_length];
                    end
                    count = [count; 2./(mt(end,1)-mt(1,1))];
                    avg_reporter_intensities = [avg_reporter_intensities; reporter_intensity];

                end
            else 
                
                index=1; 
                num = 0;
                reporter_length=0;
                for j=1:length(b)-1
                    old=b(j);
                    if a(b(j)+1) == 1 && a(b(j+1)-1) == 1
                        c = tub_norm(ind_dam(b(j),1)) - tub_norm(ind_dam(index),1);                        
                        reporter_intensity = mean(tub_norm(ind_dam(index:b(j)),2));

                        if reporter_intensity >= rep_threshold
                            reporter_length = reporter_length+c; 
                            avg_reporter_length = [avg_reporter_length; reporter_length];
                            avg_reporter_intensities = [avg_reporter_intensities; reporter_intensity];
                        end
                        num=num+1;
                    end
                    index=old;
                end
                count = [count; num./(mt(end,1)-mt(1,1))];
            end
        end
    end
    if reporter_length ~= 0
        reporter_fraction = [reporter_fraction; reporter_length./(mt(end,1)-mt(1,1))];
    end
    catch
    end
end
end