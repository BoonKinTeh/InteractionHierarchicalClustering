function [RobustClusterList,RobustLenght,Result_B] = B_DetermineRobustClusters(Result_A,MustBreakCutoff)
%% Inputs description:
% Result_A: structure data, result produced by A_HierarchicalClustering
% MustBreakCutoff: double, a cutoff such that no robust cluster is obtain with distance about this cutoff,
%                  define as the distance the similarity/correlation become negative
%% Default inputs description:
MustNotBreakCutoff = 0; % the distance that for elements with a cluster to be considered as
%                         strongly conected and cannot be broken up into smaller clusters
%% Outputs description:
% RobustClusterList: LxK interger array, storage for elemnt index of K robust clusters determined, the largest cluster size is L.
%                    Zeros indicate non-element
% RobustLenght: double, the robust length identified that give least number of cluster
% Result_B: structure, input for C_InteractionHierarchicalClustering
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
PDist = Result_A.PDist;
Linkage = Result_A.Linkage;
DistanceMatrix =  Result_A.DistanceMatrix;
%% Preperations
ElementDist = 0;
MaxLinkageDist = ceil(10*max(Linkage(:,3)));
Linkage(end+1,:) = [0,0,MaxLinkageDist];
%% Create Robust Length threshold array, and the robust length threshold define as the value give least number of cluster
RobustLenghtPct = linspace(0.01,49.99,50);
RobustLenghtArray = unique(prctile(PDist,50+RobustLenghtPct)-prctile(PDist,50-RobustLenghtPct));
RobustLenghtClusterSize = zeros(size(RobustLenghtArray,2),2);

%% Determine robust clusters (Pre), running through all the robust length threshold in RobustLenghtArray
N = size(DistanceMatrix,1);
for RB_i = 1:size(RobustLenghtArray,2)
    RobustLenght =RobustLenghtArray(1,RB_i);
    PreClusterInx = (1:N)';
    RobustClusterData = zeros(N,6);

    RobustClusterData(1,1) = Linkage(end-1,1);
    RobustClusterData(1,2) = Linkage(end-1,3);
    RobustClusterData(1,3) = Linkage((RobustClusterData(1,1)>N)*(RobustClusterData(1,1)-N*(RobustClusterData(1,1)>N))+N*(RobustClusterData(1,1)<=N),3);
    RobustClusterData(1,3) = (RobustClusterData(1,3)~=MaxLinkageDist)*RobustClusterData(1,3)+(RobustClusterData(1,3)==MaxLinkageDist)*(RobustClusterData(1,2)-ElementDist);
    RobustClusterData(1,4) = (RobustClusterData(1,2)>=MustNotBreakCutoff)*((RobustClusterData(1,2)-RobustClusterData(1,3))>=RobustLenght);
    RobustClusterData(2,1) = Linkage(end-1,2);
    RobustClusterData(2,2) = Linkage(end-1,3);
    RobustClusterData(2,3) = Linkage((RobustClusterData(2,1)>N)*(RobustClusterData(2,1)-N*(RobustClusterData(2,1)>N))+N*(RobustClusterData(2,1)<=N),3);
    RobustClusterData(2,3) = (RobustClusterData(2,3)~=MaxLinkageDist)*RobustClusterData(2,3)+(RobustClusterData(2,3)==MaxLinkageDist)*(RobustClusterData(2,2)-ElementDist);
    RobustClusterData(2,4) = (RobustClusterData(2,2)>=MustNotBreakCutoff)*((RobustClusterData(2,2)-RobustClusterData(2,3))>=RobustLenght);
    
    SetDiff = setdiff(PreClusterInx,unique(RobustClusterData(:,6)));
    NewClusterTF = RobustClusterData(1,4)+RobustClusterData(2,4)+(RobustClusterData(1,2)>MustBreakCutoff);
    RobustClusterData(1,5) = RobustClusterData(1,4)>0;
    RobustClusterData(1,6) = SetDiff(1,1);
    RobustClusterData(2,5) = RobustClusterData(2,4)>0;
    RobustClusterData(2,6) = SetDiff(1+(NewClusterTF>0),1);
    Index = find(RobustClusterData(:,1)>N);
    [~,Inx] = sort(RobustClusterData(Index,2),'descend');Index = Index(Inx,1);
    while size(Index,1)>0
        Inx = Index(1,1);
        DemergedCluster = RobustClusterData(Inx,:);
        RobustClusterData(Inx,:) = [];
        
        PRC = zeros(2,6);
        PRC(1,1) = Linkage(DemergedCluster(1,1)-N,1);
        PRC(1,2) = Linkage(DemergedCluster(1,1)-N,3);
        PRC(1,3) = Linkage((PRC(1,1)>N)*(PRC(1,1)-N*(PRC(1,1)>N))+N*(PRC(1,1)<=N),3);
        PRC(1,3) = (PRC(1,3)~=MaxLinkageDist)*PRC(1,3)+(PRC(1,3)==MaxLinkageDist)*(PRC(1,2)-ElementDist);
        PRC(1,4) = (PRC(1,2)>=MustNotBreakCutoff)*((PRC(1,2)-PRC(1,3))>=RobustLenght);
        PRC(2,1) = Linkage(DemergedCluster(1,1)-N,2);
        PRC(2,2) = Linkage(DemergedCluster(1,1)-N,3);
        PRC(2,3) = Linkage((PRC(2,1)>N)*(PRC(2,1)-N*(PRC(2,1)>N))+N*(PRC(2,1)<=N),3);
        PRC(2,3) = (PRC(2,3)~=MaxLinkageDist)*PRC(2,3)+(PRC(2,3)==MaxLinkageDist)*(PRC(2,2)-ElementDist);
        PRC(2,4) = (PRC(2,2)>=MustNotBreakCutoff)*((PRC(2,2)-PRC(2,3))>=RobustLenght);
        
        SetDiff = setdiff(PreClusterInx,unique(RobustClusterData(:,6)));
        NewClusterTF = PRC(1,4)+PRC(2,4)+(PRC(1,2)>MustBreakCutoff);
        PRC(1,5) = DemergedCluster(1,5)+PRC(1,4)>0;
        PRC(1,6) = (NewClusterTF==0)*DemergedCluster(1,6) + (NewClusterTF>0)*SetDiff(1,1);
        PRC(2,5) = DemergedCluster(1,5)+PRC(2,4)>0;
        PRC(2,6) = (NewClusterTF==0)*DemergedCluster(1,6) + (NewClusterTF>0)*SetDiff(2,1);
        
        RobustClusterData(1+sum(RobustClusterData(:,1)>0):2+sum(RobustClusterData(:,1)>0),:)=PRC;
        Index = find(RobustClusterData(:,1)>N);
        [~,Inx] = sort(RobustClusterData(Index,2),'descend');Index = Index(Inx,1);
    end
    RobustCluster = RobustClusterData(:,5).*RobustClusterData(:,6);
    ClusterInx = unique(RobustCluster);
    ClusterInx(ClusterInx==0,:) = [];
    SerialIndexTrail = zeros(1,N);
    RobustClusterList = zeros(N,N);
    Count = 0;Count2 = 0;
    for i = 1:size(ClusterInx,1)
        Inx = find(RobustCluster==ClusterInx(i,1));
        RobustClusterList(1:size(Inx,1),i) = RobustClusterData(Inx,1);
        SerialIndexTrail(1,Count2+1:Count2+size(Inx,1)) = RobustClusterData(Inx,1);
        Count = Count+1;
        Count2 = Count2+size(Inx,1);
    end
    Inx = find(RobustClusterData(:,5)==0);
    RobustClusterList(1,Count+1:Count+size(Inx,1)) = RobustClusterData(Inx,1);
    Inx = sum(RobustClusterList>0,1)==0;
    RobustClusterList(:,Inx) = [];
    Inx = sum(RobustClusterList>0,2)==0;
    RobustClusterList(Inx,:) = [];
    
    RobustLenghtClusterSize(RB_i,1) = size(RobustClusterList,1);
    RobustLenghtClusterSize(RB_i,2) = size(RobustClusterList,2);
end
%% Determine the robust length threshold, define as the value within RobustLenghtArray give least number of cluster
Ind = find((1*(RobustLenghtClusterSize(:,1)==max(RobustLenghtClusterSize(:,1)))+1*(RobustLenghtClusterSize(:,2)==min(RobustLenghtClusterSize(:,2))))==2);
if size(Ind,1)>0
    Ind = max(Ind);
else
    Ind = find(RobustLenghtClusterSize(:,2)==min(RobustLenghtClusterSize(:,2)),1,'last');
end
RobustLenght =RobustLenghtArray(1,Ind);
%% Determine robust clusters (final)
PreClusterInx = (1:N)';
RobustClusterData = zeros(N,6);
RobustClusterData(1,1) = Linkage(end-1,1);
RobustClusterData(1,2) = Linkage(end-1,3);
RobustClusterData(1,3) = Linkage((RobustClusterData(1,1)>N)*(RobustClusterData(1,1)-N*(RobustClusterData(1,1)>N))+N*(RobustClusterData(1,1)<=N),3);
RobustClusterData(1,3) = (RobustClusterData(1,3)~=MaxLinkageDist)*RobustClusterData(1,3)+(RobustClusterData(1,3)==MaxLinkageDist)*(RobustClusterData(1,2)-ElementDist);
RobustClusterData(1,4) = (RobustClusterData(1,2)>=MustNotBreakCutoff)*((RobustClusterData(1,2)-RobustClusterData(1,3))>=RobustLenght);
RobustClusterData(2,1) = Linkage(end-1,2);
RobustClusterData(2,2) = Linkage(end-1,3);
RobustClusterData(2,3) = Linkage((RobustClusterData(2,1)>N)*(RobustClusterData(2,1)-N*(RobustClusterData(2,1)>N))+N*(RobustClusterData(2,1)<=N),3);
RobustClusterData(2,3) = (RobustClusterData(2,3)~=MaxLinkageDist)*RobustClusterData(2,3)+(RobustClusterData(2,3)==MaxLinkageDist)*(RobustClusterData(2,2)-ElementDist);
RobustClusterData(2,4) = (RobustClusterData(2,2)>=MustNotBreakCutoff)*((RobustClusterData(2,2)-RobustClusterData(2,3))>=RobustLenght);

SetDiff = setdiff(PreClusterInx,unique(RobustClusterData(:,6)));
NewClusterTF = RobustClusterData(1,4)+RobustClusterData(2,4)+(RobustClusterData(1,2)>MustBreakCutoff);
RobustClusterData(1,5) = RobustClusterData(1,4)>0;
RobustClusterData(1,6) = SetDiff(1,1);
RobustClusterData(2,5) = RobustClusterData(2,4)>0;
RobustClusterData(2,6) = SetDiff(1+(NewClusterTF>0),1);
Index = find(RobustClusterData(:,1)>N);
[~,Inx] = sort(RobustClusterData(Index,2),'descend');Index = Index(Inx,1);
while size(Index,1)>0
    Inx = Index(1,1);
    DemergedCluster = RobustClusterData(Inx,:);
    RobustClusterData(Inx,:) = [];
    
    PRC = zeros(2,6);
    PRC(1,1) = Linkage(DemergedCluster(1,1)-N,1);
    PRC(1,2) = Linkage(DemergedCluster(1,1)-N,3);
    PRC(1,3) = Linkage((PRC(1,1)>N)*(PRC(1,1)-N*(PRC(1,1)>N))+N*(PRC(1,1)<=N),3);
    PRC(1,3) = (PRC(1,3)~=MaxLinkageDist)*PRC(1,3)+(PRC(1,3)==MaxLinkageDist)*(PRC(1,2)-ElementDist);
    PRC(1,4) = (PRC(1,2)>=MustNotBreakCutoff)*((PRC(1,2)-PRC(1,3))>=RobustLenght);
    PRC(2,1) = Linkage(DemergedCluster(1,1)-N,2);
    PRC(2,2) = Linkage(DemergedCluster(1,1)-N,3);
    PRC(2,3) = Linkage((PRC(2,1)>N)*(PRC(2,1)-N*(PRC(2,1)>N))+N*(PRC(2,1)<=N),3);
    PRC(2,3) = (PRC(2,3)~=MaxLinkageDist)*PRC(2,3)+(PRC(2,3)==MaxLinkageDist)*(PRC(2,2)-ElementDist);
    PRC(2,4) = (PRC(2,2)>=MustNotBreakCutoff)*((PRC(2,2)-PRC(2,3))>=RobustLenght);
    
    SetDiff = setdiff(PreClusterInx,unique(RobustClusterData(:,6)));
    NewClusterTF = PRC(1,4)+PRC(2,4)+(PRC(1,2)>MustBreakCutoff);
    PRC(1,5) = DemergedCluster(1,5)+PRC(1,4)>0;
    PRC(1,6) = (NewClusterTF==0)*DemergedCluster(1,6) + (NewClusterTF>0)*SetDiff(1,1);
    PRC(2,5) = DemergedCluster(1,5)+PRC(2,4)>0;
    PRC(2,6) = (NewClusterTF==0)*DemergedCluster(1,6) + (NewClusterTF>0)*SetDiff(2,1);
    
    RobustClusterData(1+sum(RobustClusterData(:,1)>0):2+sum(RobustClusterData(:,1)>0),:)=PRC;
    Index = find(RobustClusterData(:,1)>N);
    [~,Inx] = sort(RobustClusterData(Index,2),'descend');Index = Index(Inx,1);
end

RobustCluster = RobustClusterData(:,5).*RobustClusterData(:,6);
ClusterInx = unique(RobustCluster);
ClusterInx(ClusterInx==0,:) = [];
RobustClusterList = zeros(N,N);
Count = 0;Count2 = 0;
for i = 1:size(ClusterInx,1)
    Inx = find(RobustCluster==ClusterInx(i,1));
    RobustClusterList(1:size(Inx,1),i) = RobustClusterData(Inx,1);
    Count = Count+1;
    Count2 = Count2+size(Inx,1);
end
Inx = find(RobustClusterData(:,5)==0);
RobustClusterList(1,Count+1:Count+size(Inx,1)) = RobustClusterData(Inx,1);
Inx = sum(RobustClusterList>0,1)==0;
RobustClusterList(:,Inx) = [];
Inx = sum(RobustClusterList>0,2)==0;
RobustClusterList(Inx,:) = [];
%% Output Result
Result_B.RobustClusterList = RobustClusterList;
Result_B.DistanceMatrix = DistanceMatrix;
Result_B.N = N;
Result_B.MustNotBreakCutoff = MustNotBreakCutoff;