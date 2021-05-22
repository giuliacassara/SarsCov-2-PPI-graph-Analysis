% GET ADJACENCY MATRIX of SARS-HUMAN-drug graph
clc
clear all

addpath(genpath(pwd))
pwd

g = importNet('./dataset/sars_human_human_drug.txt', false);
g = rmedge(g, 1:numnodes(g), 1:numnodes(g)); %remove self-loops from graph
adjMatrix = adjacency(g);

n=1;
fid = fopen('./dataset/node_names.txt','rt'); % 'rt' means "read text"
ca = cell(1, 6387);
while n < 6388
      line = fgetl(fid);
      ca{n} = line;
      n = n + 1;
      %ca{n} = n;
end
fclose(fid);

g.Nodes.Name = ca';

fileID = fopen('./output/motifs.txt','w');
%// queue - initialize to empty-keeps track of the connected nodes
current_drug = 4387;
while current_drug <= 6387
   n = neighbors(g, current_drug);
   current_sars = 1;
   while current_sars < 31
       h = neighbors(g, current_sars);
       % find the intersection between the two sets
       common_neighbors = intersect(n, h);
       idx = find(common_neighbors > 30 & common_neighbors < 4387); 
       if(numel(idx) > 1)
            %skip this sars 
            %write to file        
            %fprintf(fileID, "<--- New Motif: First row: DRUG, last row: SARS, in the middle: HUMAN---> \n");
            %fprintf(fileID,'%d\n %d\n %d\n', current_drug, common_neighbors(idx), current_sars);
            drug = g.Nodes.Name(current_drug);
            sars = g.Nodes.Name(current_sars);
            humans = common_neighbors(idx);
            humans = g.Nodes.Name(humans);
            %fprintf(fileID,'%s;%s;%s', drug{:}, humans{:} , sars{:});
            fprintf(fileID,'%s;', drug{1});
            fprintf(fileID,'[');
            for n = 1 : length(humans)
                fprintf(fileID,'%s,', humans{n});
            end
            fprintf(fileID,'];');
            fprintf(fileID,'%s\n', sars{1});


       end
       current_sars = current_sars + 1;
   end
   % aggiungo questa proteina umana alla mia coda
   current_drug = current_drug + 1;
end

fclose(fileID);
