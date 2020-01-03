function [NetFeatTable] = measureNet(BranchTable)
% BranchTable is a table the indicate measured branches corrdinates
% (X1, Y1), (X2, Y2)
net_id = unique(BranchTable.ImageName);
num_net = length(net_id);

X = [BranchTable.X1,BranchTable.X2];
Y = [BranchTable.Y1,BranchTable.Y2];

%initilise centrality variables
nb_peices_corr = zeros(num_net,1);
nb_thick_branches = zeros(num_net,1);
minspantree_length = zeros(num_net,1);
mu_closeness = zeros(num_net,1);
sd_closeness = zeros(num_net,1);
mu_weighted_closenss = zeros(num_net,1);
sd_weighted_closenss = zeros(num_net,1);
mu_betweenness = zeros(num_net,1);
sd_betweenness = zeros(num_net,1);
mu_weighted_betweenness = zeros(num_net,1);
sd_weighted_betweenness = zeros(num_net,1);
mu_eigenvector = zeros(num_net,1);
sd_eigenvector = zeros(num_net,1);

num_degree_1 = zeros(num_net,1);
num_degree_2 = zeros(num_net,1);
num_degree_3 = zeros(num_net,1);
num_degree_4 = zeros(num_net,1);
num_degree_5plus = zeros(num_net,1);

%initilise orientation variables
mu_angle_ratio1 = zeros(num_net,1);
mu_angle_ratio2 = zeros(num_net,1);
sd_angle_ratio1 = zeros(num_net,1);
sd_angle_ratio2 = zeros(num_net,1);
mu_angles = zeros(num_net,3);
sd_angles = zeros(num_net,3);


%initilise vornoi variables
mu_bthickness_mu = zeros(num_net,1);
sd_bthickness_mu = zeros(num_net,1);
mu_bthickness_sd = zeros(num_net,1);
sd_bthickness_sd = zeros(num_net,1);
mu_blength_mu = zeros(num_net,1);
sd_blength_mu = zeros(num_net,1);
mu_blength_sd = zeros(num_net,1);
sd_blength_sd = zeros(num_net,1);
mu_vornoi_area = zeros(num_net,1);
sd_vornoi_area = zeros(num_net,1);
mu_nb_vornoi_nodes = zeros(num_net,1);
sd_nb_vornoi_nodes = zeros(num_net,1);
mu_vornoi_dist_mu = zeros(num_net,1);
sd_vornoi_dist_mu = zeros(num_net,1);
mu_vornoi_dist_sd = zeros(num_net,1);
sd_vornoi_dist_sd = zeros(num_net,1);

minx = 0;
miny = 0;
maxx = 325;
maxy = 325;

for j=1:num_net %for each image
    idx = (strcmpi(BranchTable.ImageName,net_id(j)));
    num_possible_paths = (length(idx)*(length(idx)-1))/2;%n(n-1)/2
    
    
    s = [];
    t = [];
    weights = [];
    
    x = X(idx,:);
    y = Y(idx,:);
    net_weights = BranchTable.Branchlength(idx);
    net_weights2 = BranchTable.Branchthickness(idx);
    
    % Merge very close nodes
    [x, y, net_weights, net_weights2] = simplify_network(x,y,net_weights,net_weights2);
    
    %Exclude thick branches (clumps)
    thick_idx=net_weights2>9;
    nb_thick_branches(j) = sum(thick_idx);
    if sum(thick_idx)
        x(thick_idx,:) = [];
        y(thick_idx,:) = [];
        net_weights(thick_idx) = [];
        net_weights2(thick_idx) = [];
    end
    
    all_coord = [x(:,1),y(:,1);x(:,2),y(:,2)];
    u_coord = unique(all_coord,'rows');
    nodes = 1:size(u_coord,1);

    % consider networks with at least 5 branches
    nb_peices_corr(j) = size(x,1);
    if size(x,1)>=5
        % Nodes centrality
        for i=1:size(x,1)
            node1_idx = find(u_coord(:,1)==x(i,1) & u_coord(:,2)==y(i,1));
            node2_idx = find(u_coord(:,1)==x(i,2) & u_coord(:,2)==y(i,2));
            s = [s;nodes(node1_idx)];
            t = [t;nodes(node2_idx)];
            weights = [weights;net_weights(i)];
        end
        [edges ia] = unique(sort([s,t],2),'rows');
        weights = weights(ia);

        G = graph(edges(:,1),edges(:,2));
        A = adjacency(G);
        [Tree] = graphminspantree(A);
        
        minspantree_length(j) = sum(Tree(:)==1);
        degree = centrality(G,'degree');
        num_degree_1(j) = sum(degree==1);%this should be the same as the number of isolated elements+twigs
        num_degree_2(j) = sum(degree==2);
        num_degree_3(j) = sum(degree==3);
        num_degree_4(j) = sum(degree==4);
        num_degree_5plus(j) = sum(degree>=5);
        closeness = centrality(G,'closeness');
        weighted_closenss = centrality(G,'closeness','cost',weights);
        mu_closeness(j) = mean(closeness);
        sd_closeness(j) = std(closeness);
        mu_weighted_closenss(j) = mean(weighted_closenss);
        sd_weighted_closenss(j) = std(weighted_closenss);

        betweenness = centrality(G,'betweenness');
        weighted_betweenness = centrality(G,'betweenness','cost',weights);
        mu_betweenness(j) = mean(betweenness);
        sd_betweenness(j) = std(betweenness);
        mu_weighted_betweenness(j) = mean(weighted_betweenness);
        sd_weighted_betweenness(j) = std(weighted_betweenness);

        eigenvector = centrality(G,'eigenvector');
        mu_eigenvector(j) = mean(eigenvector);
        sd_eigenvector(j) = std(eigenvector);


        % Branch oreintation
        % Center at each node such that it represent (0,0) coordinates theb claculate the angle of each vector relative to x=0
        node_angles = [];
        angle_ratio1 = [];
        angle_ratio2 = [];
        for i=1:size(u_coord,1)
            vec = [];
            vec_angles = [];
            node_idx1=find(x(:,1)==u_coord(i,1) & y(:,1)==u_coord(i,2));
            vec = [vec;x(node_idx1,2),y(node_idx1,2)];
            node_idx2 = find(x(:,2)==u_coord(i,1) & y(:,2)==u_coord(i,2));
            vec = [vec;x(node_idx2,1),y(node_idx2,1)];
            vec = vec-u_coord(i,:);%center at the current node
            nvec = size(vec,1);
            if nvec==3 % As most nodes have degree =3 consider only this degree
                for k=1:nvec-1
                    for l=k+1:nvec
                        angle = abs(atan2d(vec(k,2),vec(k,1))-atan2d(vec(l,2),vec(l,1)));
                        vec_angles = [vec_angles;angle];
                    end
                end
                [s si] = sort(vec_angles);
                %replace the maximum angle to complement 360 degree with the nvec-1 angles
                vec_angles(si(nvec)) = 360-sum(vec_angles(si(1:nvec-1)));
                [s si] = sort(vec_angles);
                node_angles = [node_angles;vec_angles(si(1:2))',vec_angles(si(nvec))];
                %calculate the ratio between the two smallest angles to the largest
                %angles (symmetry of branching angle)
                angle_ratio1 = [angle_ratio1;vec_angles(si(1))./vec_angles(si(nvec))];
                angle_ratio2 = [angle_ratio2;vec_angles(si(2))./vec_angles(si(nvec))];
            end
        end
        if ~isempty(angle_ratio1)
            mu_angle_ratio1(j) = mean(angle_ratio1);
            mu_angle_ratio2(j) = mean(angle_ratio2);
            sd_angle_ratio1(j) = std(angle_ratio1);
            sd_angle_ratio2(j) = std(angle_ratio2);
            mu_angles(j,:) = mean(node_angles,1);
            sd_angles(j,:) = std(node_angles,1);
        end

        %Voronoi measurements
        mu_blength = [];
        sd_blength = [];
        mu_thickness = [];
        sd_thickness = [];
        vcell_area = [];
        num_vnodes = [];
        mu_vdist = [];
        sd_vdist = [];
        [v,c] = voronoin(u_coord);
        for i=1:size(u_coord,1)
            %for each node find half total branch lengths spining from it 
            node_idx1 = (x(:,1)==u_coord(i,1) & y(:,1)==u_coord(i,2));
            node_idx2 = (x(:,2)==u_coord(i,1) & y(:,2)==u_coord(i,2));
            blength = net_weights(node_idx1 | node_idx2);
            bthickness = net_weights2(node_idx1 | node_idx2);
            mu_blength(i) = mean(blength);
            sd_blength(i) = std(blength);
            mu_bthickness(i) = mean(bthickness);
            sd_bthickness(i) = std(bthickness);
            
            vcoord=v(c{i},:); 
            if all(c{i}~=1) & all(vcoord(:,1)>=minx) & all(vcoord(:,2)>=miny) & all(vcoord(:,1)<=maxx) & all(vcoord(:,2)<=maxy) % consider only closed cells (vertex 1 is inf)
                               
                vcell_area(i) = polyarea(v(c{i},1),v(c{i},2)); 
                %find the distance between each node and voronoi cell corner
                %points
                vdist = pdist2(u_coord(i,:),v(c{i},:));
                num_vnodes(i) = size(c{i},2);
                mu_vdist(i) = mean(vdist);
                sd_vdist(i) = std(vdist);

                %find the mean and sd of the voronoi cell polygon dimensions
                vcoord(end+1,:) = vcoord(1,:);%to meas dist between last and first nodes
                vdim = [];%added 8/12/2019
                for k=1:num_vnodes(i)
                    vdim(k) = pdist2(vcoord(k,:),vcoord(k+1,:));
                end
                mu_vdim(i) = mean(vdim(vdim>5));
                sd_vdim(i) = std(vdim(vdim>5));
                num_vnodes(i) = sum(vdim>5);%exclude nodes that are very close to each others
            end


        end
        if ~isempty(vcell_area)
            mu_bthickness_mu(j) = mean(mu_bthickness);
            sd_bthickness_mu(j) = std(mu_bthickness);
            mu_bthickness_sd(j) = mean(sd_bthickness);
            sd_bthickness_sd(j) = std(sd_bthickness);
            mu_blength_mu(j) = mean(mu_blength);
            sd_blength_mu(j) = std(mu_blength);
            mu_blength_sd(j) = mean(sd_blength);
            sd_blength_sd(j) = std(sd_blength);

            mu_vornoi_area(j) = mean(vcell_area);
            sd_vornoi_area(j) = std(vcell_area);
            mu_nb_vornoi_nodes(j) = mean(num_vnodes);
            sd_nb_vornoi_nodes(j) = std(num_vnodes);

            mu_vornoi_dist_mu(j) = mean(mu_vdist);
            sd_vornoi_dist_mu(j) = std(mu_vdist);
            mu_vornoi_dist_sd(j) = mean(sd_vdist);
            sd_vornoi_dist_sd(j) = std(sd_vdist);

            mu_vornoi_dim_mu(j) = mean(mu_vdim);
            sd_vornoi_dim_mu(j) = std(mu_vdim);
            mu_vornoi_dim_sd(j) = mean(sd_vdim);
            sd_vornoi_dim_sd(j) = std(sd_vdim);
        end
    end
end



NetFeatTable=table(net_id, nb_peices_corr,nb_thick_branches, minspantree_length, ...
mu_closeness, sd_closeness, mu_weighted_closenss, ...
sd_weighted_closenss, mu_betweenness, sd_betweenness, mu_weighted_betweenness,...
sd_weighted_betweenness, mu_eigenvector, sd_eigenvector, ...
mu_angle_ratio1, mu_angle_ratio2, sd_angle_ratio1, sd_angle_ratio2,...
mu_angles(:,1),mu_angles(:,2),mu_angles(:,3), sd_angles(:,1), sd_angles(:,2),...
sd_angles(:,3), num_degree_1, num_degree_2, num_degree_3, num_degree_4, num_degree_5plus,...
mu_bthickness_mu, sd_bthickness_mu, mu_bthickness_sd, sd_bthickness_sd, mu_blength_mu, sd_blength_mu,...
mu_blength_sd, sd_blength_sd, mu_vornoi_area, sd_vornoi_area,...
mu_nb_vornoi_nodes, sd_nb_vornoi_nodes, mu_vornoi_dist_mu,...
sd_vornoi_dist_mu, mu_vornoi_dist_sd, sd_vornoi_dist_sd,'RowNames',net_id);

NetFeatTable.Properties.VariableNames(20:25)={'mu_angle_1','mu_angle_2','mu_angle_3','sd_angle_1','sd_angle_2','sd_angle_3'};

    