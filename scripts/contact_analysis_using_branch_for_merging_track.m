% Contact Analysis Using Branch Segments for Live-Cell Tracking
% --------------------------------------------------------------
% This script computes contact statistics for branches before merging events, using spot-to-mitochondria distances.

% Clear existing variables to avoid conflicts
clearvars;

% Load spot distance data (output from Python script)
%    'result_whole_spot_distance.xlsx' contains [SpotID, Distance]
all_spots_D = readmatrix('result_whole_spot_distance.xlsx', 'Range',[2,2]);

% Load the TrackMate graph structure (nodes & edges)
%    Nodes: merging_g.Nodes.ID contains spot IDs
%    Edges: merging_g.Edges.SPOT_SOURCE_ID / SPOT_TARGET_ID for links
merging_g = trackmateGraph('merging_G.xml',[],[],true);

% Read branch information for merging events
%    Columns: predecessor flag (col2), successor flag (col3), ..., startID (col8), endID (col9)merging_b=readmatrix('merging_b.csv', 'Range',[5,2]);
% Initialize containers for branch node index lists
before_merge_branch_ID = {};
after_merge_branch_ID  = {};
before_merge_branch_sID = []; % succersor's ID of before-merging-branch
whole_merge_branch_ID  = {};

% Extract node index paths for each branch segment
for i = 1 : size(merging_b,1)
    start_s = find(merging_g.Nodes.ID==merging_b(i,8));
    end_s = find(merging_g.Nodes.ID==merging_b(i,9));
    % Compute shortest path (node indices) between start and end
    whole_merge_branch_ID{end+1} = [shortestpath(merging_g,start_s,end_s)];

    % Pre-merge branch: successor flag == 1
    if merging_b(i,3)==1
        before_merge_branch_ID{end+1}=[shortestpath(merging_g,start_s,end_s)];
        % Record the spot ID at merge point for pairing
        end_edge= find(merging_g.Edges.SPOT_SOURCE_ID==merging_b(i,9));
        before_merge_branch_sID(end+1)= merging_g.Edges(end_edge,3).SPOT_TARGET_ID;
    end
    % Post-merge branch: predecessor flag == 2
    if merging_b(i,2)==2
        after_merge_branch_ID{end+1}=[shortestpath(merging_g,start_s,end_s)];
    end
end


% Organize distances for each branch by mapping node indices to distances
before_merge_branch_Ds={};
after_merge_branch_Ds={};
whole_merge_branch_Ds={};

% Loop through pre-merge branches
for i = 1:size(before_merge_branch_ID,2)
    ds = arrayfun(@(x) get_distance(merging_g.Nodes(x,1).ID,all_spots_D), before_merge_branch_ID{1,i});
    before_merge_branch_Ds{end+1}=ds;
end
% Loop through post-merge branches
for i = 1:size(after_merge_branch_ID,2)
    ds = arrayfun(@(x) get_distance(merging_g.Nodes(x,1).ID,all_spots_D), after_merge_branch_ID{1,i});
    after_merge_branch_Ds{end+1}=ds;
end
% Loop through whole branches
for i = 1:size(whole_merge_branch_ID,2)
    ds = arrayfun(@(x) get_distance(merging_g.Nodes(x,1).ID,all_spots_D), whole_merge_branch_ID{1,i});
    whole_merge_branch_Ds{end+1}=ds;
end

% Compute contact statistics for before-merge branches
%    - Average distance per branch
%    - Contact affinity: fraction of frames with zero distance
%    - Contact durations: average run length of zero-distance sequences
%    - Contact events: number of zero-distance runs
b_average_distance = cellfun(@(x) mean(x), before_merge_branch_Ds, 'UniformOutput', false);
b_average_distance = transpose(cell2mat(b_average_distance));
b_contact_affinity = cellfun(@(x) sum(x==0)/length(x),before_merge_branch_Ds, 'UniformOutput', false);
b_contact_affinity = transpose(cell2mat(b_contact_affinity));
b_contact_duration = cellfun(@(x) contact_duration(find(x==0)),before_merge_branch_Ds, 'UniformOutput', false);
b_contact_duration_avg = cellfun(@(x) mean(x), b_contact_duration, 'UniformOutput', false);
b_contact_duration_avg = transpose(cell2mat(b_contact_duration_avg));
b_contact_event = cellfun(@(x) contact_event(find(x==0)),before_merge_branch_Ds, 'UniformOutput', false);
b_contact_frequency={};
for i=1:length(b_contact_event)
    b_contact_frequency{end+1} = b_contact_event{i}/length(before_merge_branch_Ds{i});
end

% Save before-merge summary to Excel
save('before_merge_contact_duration.mat','b_contact_duration');
save('before_merge_contact_event.mat','b_contact_event');
save('before_merge_contact_frequency.mat','b_contact_frequency');
b_contact_frequency = transpose(cell2mat(b_contact_frequency));
before_result_file = table(b_average_distance, b_contact_affinity, b_contact_duration_avg, b_contact_frequency );
writetable(before_result_file,'before_merge_result_summary.xlsx',"WriteMode","replacefile","AutoFitWidth",false);

% Identify diffusion-type merges (no contact on either paired branch)
diffusion_merge_list=[];
merged_spot=[];
whole_merging_event=0;
for i=1:length(before_merge_branch_sID)
    % Find paired branch index based on same merge spot ID
    pair_branch = find(before_merge_branch_sID==before_merge_branch_sID(i));
    if length(pair_branch)~=2
        error("Merging weired");
    end
    pair_branch=pair_branch(pair_branch~=i);
    % Count unique merging events
    if isempty(find(merged_spot==before_merge_branch_sID(i)))
        whole_merging_event=whole_merging_event+1;
        merged_spot(end+1)=before_merge_branch_sID(i);
    end
    % Check if both branches have zero affinity -> diffusion merge
    if (b_contact_affinity(i)==0) & (b_contact_affinity(pair_branch)==0)
        if isempty(find(diffusion_merge_list==pair_branch))
            diffusion_merge_list(end+1) = i;
        end
    end
end

% Save and print diffusion merge stats
save('diffusion_merge_list.mat','diffusion_merge_list');
fprintf("whole merging event: %d \n", whole_merging_event);
fprintf("diffusion merging event: %d \n", length(diffusion_merge_list));
fprintf("diffusion merging event percent: %d \n", length(diffusion_merge_list)/whole_merging_event*100);

fprintf("before_merge_average_d: %d \n", mean(b_average_distance));
fprintf("before_merge_average_affinity: %d \n", mean(b_contact_affinity)*100);
fprintf("before_merge_average_duration: %d \n", mean(b_contact_duration_avg,'omitnan'));
fprintf("before_merge_average_frequency: %d \n", mean(b_contact_frequency,'omitnan'));
fprintf("before_merge_median_affinity: %d \n", median(b_contact_affinity));
fprintf("before_merge_median_frequency: %d \n", median(b_contact_frequency,'omitnan'));

% Helper functions

% Return distance for a given spot ID
function d = get_distance(spot_id,all_spots_D)
    i=find(all_spots_D(:,1)==spot_id);
    d=all_spots_D(i,2);
end

% Compute run lengths of consecutive contact frames (zero-distance)
function d = contact_duration(ele)
    if isempty(ele)
        d=[];
        return
    end
    if length(ele)==1
        d=[1];
        return
    end
    d=[];
    temp_d = 1;
    for i =1:length(ele)
        if i == length(ele)
            d(end+1)=temp_d;
            return
        end
        temp_s=ele(i);
        if ele(i+1) == temp_s+1
            temp_d=temp_d+1;
        else
            d(end+1)=temp_d;
            temp_d=1;
        end
    end
end


% Count number of separate contact events (zero-distance runs)
function f = contact_event(ele)
    if isempty(ele)
        f=0;
        return
    end
    if length(ele)==1
        f=1;
        return
    end
    f=0;
    temp_d = 1;
    for i =1:length(ele)
        if i == length(ele)
            f=f+1;
            return
        end
        temp_s=ele(i);
        if ele(i+1) == temp_s+1
            temp_d=temp_d+1;
        else
            f=f+1;
            temp_d=1;
        end
    end
end
