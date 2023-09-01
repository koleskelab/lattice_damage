%% Koleske Lab
% Last updated September 1, 2023
% abl2_intensity_quantification.m
%
% DESCRIPTION: Quantify the intensity of Abl2
% on 23C and 37C storage treated biotinylated rhodamine
% taxol-stabilized MTs to determine where Abl2 resides (boundary analysis)
%
% INPUT: 
% 1) filaments: total number of MTs
% 2) mt_tot: array of MT lengths and respective MT intensity per unit of length (left column = unit of length measurement; right column = intensity of MT at that particular locus of filament)
% 3) abl2_tot: array of MT lengths and respective Abl2 intensity per unit of length (left column = unit of MT length measurement; right column = intensity of Abl2 at that particular locus of filament)
% 4) mt_threshold: the normalized intensity value that user denotes whether a segment is damaged or not (default value: 0.5)
% 5) increment: the # of pieces for user to define border from lattice or break the healthy MT segment up into (default value: 7)
% --- NOTE: the number of filaments in both the mt_tot and abl2_tot (every
% other 2 columns in the mt_tot and abl2_tot arrays) should match ---
% 
% OUTPUTt: 
% 1) abl2_intensity: mean global intensities of Abl2 per MT
% 2) mt_intesnity: mean MT intensities
% 3) avg_mt_length: 'healthy' MT segments (adjacent to damaged segments)
% 4) abl2_loci_beg: Abl2 intensity at the 'beg' portions of protected/healthy MT
% segments (end #1) that borders damaged sites
% 5) abl2_loci_middle: Abl2 intensity at the 'middle' portions of protected/healthy MT
% segments
% 6) abl2_loci_end: Abl2 intensity at the 'end' portions of protected/healthy MT
% segments (end #2) that borders damaged sites


function [abl2_intensity,mt_intensity,...
    avg_mt_length, abl2_loci_beg, abl2_loci_middle, ...
    abl2_loci_end] = abl2_reporter_intensity_quantification(filaments,mt_tot,abl2_tot,mt_threshold,increment)

abl2_intensity = [];
mt_intensity = [];

avg_mt_length = []; % the protected segments/boundaries
abl2_loci_beg = []; % abl2 intensity on the beginning portions of protected segments/boundaries
abl2_loci_middle = []; % abl2 intensity on middle of protected segments/boundaries
abl2_loci_end = []; % abl2 intensity on the ends of protected segments/boundaries


for i=1:2:(filaments-1)
    mt=rmmissing(mt_tot(:,i:i+1)); 
    abl2=rmmissing(abl2_tot(:,i:i+1)); 
   
    %need to internally normalize
    mt_norm = [mt(:,1) (mt(:,2)-min(mt(:,2)))./(max(mt(:,2))-min(mt(:,2)))];
    abl2_norm = [abl2(:,1) (abl2(:,2)-min(abl2(:,2)))./(max(abl2(:,2))-min(abl2(:,2)))];
    
    ratio_mt = abl2_norm(:,2)./mt_norm(:,2);  
    mt_norm = mt_norm(find(ratio_mt ~= Inf),:);
    abl2_norm = abl2_norm(find(ratio_mt ~= Inf),:);

    ind_ok = find(abl2_norm(:,2) >= mt_threshold); % find areas along the single MT that are above user-defined threshold of 'healthy' versus 'damaged'
    tot_intensity = sum(abl2_norm(ind_ok,2)); 
    abl2_intensity = [abl2_intensity; tot_intensity./mt(end,1)];

    ind_condensation = find(abl2_norm(:,2) > 0.3);
    mt_intensity = [mt_intensity; mt_norm(ind_condensation,2)];
    
    % The following section will focus on determining where Abl2 binds and
    % if it binds at boundaries between 'healthy' and 'damaged' 

    pos_ind = find(mt_norm(:,2) >= mt_threshold); % only selecting the regions of the MT that are there to differentiate the segments that are 'damaged'
    
        a=diff(pos_ind);
        b=find(a ~= 1);
        if isempty(b) == 1 
            mt_length = mt_norm(pos_ind(end),1)-mt_norm(pos_ind(1),1); % segment of the MT that isn't too damaged (right next to the damaged areas)
            mt_length_tot = mt_length/increment; % divide into 3 distinct sections: first border, middle, end border
            
            inc = round(length(pos_ind)/increment,0);
            first_portion = pos_ind(1):pos_ind(1)+inc; 
            end_portion = pos_ind(end)-inc:pos_ind(end);
            middle_portion = pos_ind(1)+inc+1:pos_ind(end)-inc-1; 

            % now need to determine the locus of the Abl2 condensates
            avg_mt_length = [avg_mt_length; mt_norm(pos_ind(end),1)-mt_norm(pos_ind(1),1)]; % still need to retain for outside the loop (output of the function)
            abl2_loci_middle = [abl2_loci_middle; mean(abl2_norm(middle_portion,2))]; 
            
        else
            if length(b) == 1 
                seg1_len = mt_norm(a(b)+pos_ind(b),1)-mt_norm(pos_ind(b),1); % first healthy segment
                seg2_len = (mt_norm(pos_ind(end),1) - mt_norm(pos_ind(b+1),1)); % second healthy segment
                total = seg1_len + seg2_len;
                mt_length = total;
                avg_mt_length = [avg_mt_length; mt_length];
                seg1 = mt_norm(pos_ind(b): pos_ind(b)+a(b), 1);
                seg2 = mt_norm(pos_ind(b+1): pos_ind(end), 1);
                    if mt_norm(end,1) > 5 
                        inc1 = round(length(pos_ind(b):pos_ind(b)+a(b))/increment,0)

                        end_portion1 = (pos_ind(b)+a(b)-inc1):pos_ind(b)+a(b)
                        middle_portion1 = (pos_ind(b)+inc1+1):pos_ind(end)-inc1-1                
                        abl2_loci_middle = [abl2_loci_middle; mean(abl2_norm(middle_portion1,2))]; 
                        abl2_loci_end = [abl2_loci_end; mean(abl2_norm(end_portion1,2))]; 

                        inc2 = round(length(pos_ind(b+1):pos_ind(end))/increment,0)
                        
                        first_portion2 = pos_ind(b+1):pos_ind(b+1)+inc2
                        middle_portion2 = (pos_ind(b+1)+inc2+1):pos_ind(end)-inc2-1
                        
                        % now need to determine the locus of the Abl2 condensates
                        abl2_loci_beg = [abl2_loci_beg; mean(abl2_norm(first_portion2,2))]; 
                        abl2_loci_middle = [abl2_loci_middle; mean(abl2_norm(middle_portion2,2))]; 
                    end
            else % this is for length(b) >= 1
                    index=1; % this variable will change within the loop to reflect the 'old' b index
                    num = 0;
                    mt_length=0;
                    for j=1:length(b)-1
                        old=b(j);
                        if a(b(j)+1) == 1 && a(b(j+1)-1) == 1
                            c = mt_norm(pos_ind(b(j+1),1)) - mt_norm(pos_ind(b(j)),1);
                            mt_length = mt_length+c;                        
                            seg = mt_norm(pos_ind(b(j)):pos_ind(b(j+1)),1);
                            len = mt(pos_ind(b(j)),1) - mt(pos_ind(b(j+1)),1); 

                            if length(pos_ind(b(j):b(j+1))) >= 5
                                inc = round(length(pos_ind(b(j):b(j+1)))/increment,0) 
                                if isempty(find(mt_norm(pos_ind(b(j)):pos_ind(b(j+1)),1) == mt_norm(1,1))) == 1 &&  isempty(find(mt_norm(pos_ind(b(j)):pos_ind(b(j+1)),1) == mt_norm(end,1))) == 1
                                    pos_ind(b(j))
                                    pos_ind(b(j+1))
                                    first_portion = pos_ind(b(j)):pos_ind(b(j))+inc
                                    end_portion = pos_ind(b(j+1))-inc:pos_ind(b(j+1))
                                    middle_portion = pos_ind(b(j))+inc+1:pos_ind(b(j+1))-inc-1
                  
                                    % now need to determine the locus of the Abl2 condensates
                                    abl2_loci_beg = [abl2_loci_beg; mean(abl2_norm(first_portion,2))]; 
                                    abl2_loci_middle = [abl2_loci_middle; mean(abl2_norm(middle_portion,2))]; 
                                    abl2_loci_end = [abl2_loci_end; mean(abl2_norm(end_portion,2))]; 
                                end
                            elseif length(pos_ind(b(j):b(j+1))) < 5 && length(pos_ind(b(j):b(j+1))) >= 3
                                inc = round(length(pos_ind(b(j):b(j+1)))/3,0) 
                                if isempty(find(mt_norm(pos_ind(b(j)):pos_ind(b(j+1)),1) == mt_norm(1,1))) == 1 &&  isempty(find(mt_norm(pos_ind(b(j)):pos_ind(b(j+1)),1) == mt_norm(end,1))) == 1
                                    pos_ind(b(j))
                                    pos_ind(b(j+1))
                                    first_portion = pos_ind(b(j))
                                    end_portion = pos_ind(b(j+1))
                                    middle_portion = pos_ind(b(j))+inc:pos_ind(b(j+1))-inc

                                    abl2_loci_beg = [abl2_loci_beg; mean(abl2_norm(first_portion,2))]; 
                                    abl2_loci_middle = [abl2_loci_middle; mean(abl2_norm(middle_portion,2))]; 
                                    abl2_loci_end = [abl2_loci_end; mean(abl2_norm(end_portion,2))]; 
                                end

                        end
                    end
                    avg_mt_length = [avg_mt_length; mt_length];             
                    end
            end
        end
end
end