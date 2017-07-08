function [RBSerialIndex,Result_C] = C_InteractionHierarchicalClustering(Result_B,InteractionDistanceDef)
%% Inputs description:
% Result_B: structure data, result produced by B_DetermineRobustClusters
% InteractionDistanceDef : string, Method to define the interaction distance between clusters
%           average = Average distances (default)
%           Pct2575 = Average between precentile 25 and 75 distances
%           median = median distance
%% Default inputs description:
PlotFigure_OrginalAndIHCSeriatedDistanceMatrix = 1; % Set to 1 to plot the colormap for distance matrix (orginal and hierarchical seriated) 
%% Outputs description:
% RBSerialIndex: 1xN integer array, the dendrogram order obtained for performing interaction-hierarchical clustering
% Result_C: structure, input for D_IdentifyClusters
% RBLinkage: K-1x3 double array, each row represents the robust cluster with column 1 index merged with
%          the robust cluster with column 2 index at distance in column 3.
%          Robust Cluster index, RCI: if RCI <= K, CI represent the index for the K robust clusters
%                                     if RCI > K, represent clusters merged at row RCI-K
%% Read Me:
% This project is published for "Cluster fusion-fission dynamics in the Singapore stock exchange", 
% by Boon Kin Teh and Siew Ann Cheong.
% Please refer to the paper for more details, and cite the paper if you are using this code to perform interaction-hierarchical clustering.
% Thank you.

%% Lastest updated date:
% 08 July 2017

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Start here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract Data
RobustClusterList = Result_B.RobustClusterList;
DistanceMatrix = Result_B.DistanceMatrix;
N = Result_B.N;
MustNotBreakCutoff = Result_B.MustNotBreakCutoff;
%% Determine the interaction distance between the robust clusters
RBDistance = zeros(size(RobustClusterList,2));
for i = 1:size(RobustClusterList,2)
    RCList1 = RobustClusterList(RobustClusterList(:,i)>0,i);
    for j = i:size(RobustClusterList,2)
        RCList2 = RobustClusterList(RobustClusterList(:,j)>0,j);
        Dist = DistanceMatrix(RCList1,RCList2);
        Dist = sort(reshape(Dist,[1,size(Dist,1)*size(Dist,2)]));
        if strcmpi(InteractionDistanceDef,'median')==1
            RBDistance(j,i) = median(Dist);
        elseif strcmpi(InteractionDistanceDef,'Pct2575')==1
            InxV = 1*(Dist>=prctile(Dist,25))+1*(Dist<=prctile(Dist,75))==2;
            RBDistance(j,i) = mean(Dist(InxV));
        else
            RBDistance(j,i) = mean(Dist);
        end
        RBDistance(i,j) = RBDistance(j,i);
    end
end
%% Interaction-Hierarchical Clustering
RBN = size(RobustClusterList,2);
RBInx = 1:RBN;
RBLinkage = zeros(RBN-1,3);
PreRBDistance = tril(ceil(2*max(max(RBDistance)))*ones(size(RBDistance)))+triu(RBDistance,1);
PreRBClusterList = zeros(N,size(RobustClusterList,2));
PreRBClusterList(1:size(RobustClusterList,1),1:size(RobustClusterList,2))=RobustClusterList;
Count = 1;
while size(PreRBClusterList,2)>1
    [Inx1,Inx2] = find(PreRBDistance==min(min(PreRBDistance)));
    RBLinkage(Count,1) = RBInx(1,Inx1(1,1));
    RBLinkage(Count,2) = RBInx(1,Inx2(1,1));
    RBLinkage(Count,3) = min(min(PreRBDistance));
    RBInx(1,Inx1(1,1)) = RBN+Count;
    RBInx(:,Inx2(1,1)) = [];
    RCList1 = unique([0;PreRBClusterList(:,Inx1(1,1));PreRBClusterList(:,Inx2(1,1))]);RCList1(1,:) = [];
    PreRBClusterList(:,Inx1(1,1)) = 0;
    PreRBClusterList(1:size(RCList1,1),Inx1(1,1)) = RCList1;
    PreRBClusterList(:,Inx2(1,1)) = [];
    PreRBDistance(:,Inx2(1,1))=[];
    PreRBDistance(Inx2(1,1),:)=[];
    for i = 1:Inx1(1,1)-1
        RCList2 = PreRBClusterList(PreRBClusterList(:,i)>0,i);
        Dist = DistanceMatrix(RCList1,RCList2);
        Dist = reshape(Dist,[1,size(Dist,1)*size(Dist,2)]);
        if strcmpi(InteractionDistanceDef,'median')==1
            PreRBDistance(i,Inx1(1,1)) = median(Dist);
        elseif strcmpi(InteractionDistanceDef,'Pct2575')==1
            InxV = 1*(Dist>=prctile(Dist,25))+1*(Dist<=prctile(Dist,75))==2;
            PreRBDistance(i,Inx1(1,1)) = mean(Dist(InxV));
        else
            PreRBDistance(i,Inx1(1,1)) = mean(Dist);
        end
    end
    for i = Inx1(1,1)+1:size(PreRBClusterList,2)
        RCList2 = PreRBClusterList(PreRBClusterList(:,i)>0,i);
        Dist = DistanceMatrix(RCList1,RCList2);
        Dist = reshape(Dist,[1,size(Dist,1)*size(Dist,2)]);
        if strcmpi(InteractionDistanceDef,'median')==1
            PreRBDistance(Inx1(1,1),i) = median(Dist);
        elseif strcmpi(InteractionDistanceDef,'Pct2575')==1
            InxV = 1*(Dist>=prctile(Dist,25))+1*(Dist<=prctile(Dist,75))==2;
            PreRBDistance(Inx1(1,1),i) = mean(Dist(InxV));
        else
            PreRBDistance(Inx1(1,1),i) = mean(Dist);
        end
    end
    Count = Count+1;
end
%% Interaction-Hierarchical Clustering Dendrogram
[~,~,PreRBSerial] = dendrogram(RBLinkage,0);
RBSerialIndex = zeros(1,N);
Count = 0;
for i = 1:size(PreRBSerial,2)
    Inx = RobustClusterList(RobustClusterList(:,PreRBSerial(1,i))>0,PreRBSerial(1,i));
    RBSerialIndex(Count+1:Count+size(Inx,1))= Inx;
    Count = Count + size(Inx,1);
end

%% Plot Original and IHC-Seriated Distance Matrix
if PlotFigure_OrginalAndIHCSeriatedDistanceMatrix == 1
    figure(2);clf;hold on;
    subplot(1,2,1);hold on;
    imagesc(-DistanceMatrix);
    colormap 'jet';
    xlim([0.5,size(DistanceMatrix,2)+0.5]);
    ylim([0.5,size(DistanceMatrix,1)+0.5]);
    axis('square');
    title('Original Distance Matrix','fontsize',16);
    
    subplot(1,2,2);hold on;
    imagesc(-DistanceMatrix(RBSerialIndex,RBSerialIndex));
    colormap 'jet';
    xlim([0.5,size(DistanceMatrix,2)+0.5]);
    ylim([0.5,size(DistanceMatrix,1)+0.5]);
    axis('square');
    title('IHC-Seriated Distance Matrix','fontsize',16);
end
%% Output Result
Result_C.RBLinkage = RBLinkage;
Result_C.RobustClusterList = RobustClusterList;
Result_C.DistanceMatrix = DistanceMatrix;
Result_C.MustNotBreakCutoff = MustNotBreakCutoff;