function IdentifiedClusterList = D_IdentifyClusters(Result_C)
%% Inputs description:
% Result_C: structure data, result produced by C_InteractionHierarchicalClustering
%
%% Outputs description:
% ClusterList: LxK interger array, storage for elemnt index of K identified clusters determined, the largest cluster size is L.
%                    Zeros indicate non-element
%
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
RBLinkage = Result_C.RBLinkage;
RobustClusterList = Result_C.RobustClusterList;
DistanceMatrix = Result_C.DistanceMatrix;
MustNotBreakCutoff = Result_C.MustNotBreakCutoff;
%% Identify Clusters only if number of robust cluster is more than 3
if size(RBLinkage,1)>3
Inx = min(size(RBLinkage,1)-1,sum(RBLinkage(:,3)>MustNotBreakCutoff));
IntraClusterList = zeros(Inx+1,size(DistanceMatrix,1));
SimilarityMeasure = zeros(Inx,3);
IntreClusterCoefficient= zeros(size(DistanceMatrix,1)^2,1);
IntraClusterCoefficient= reshape(DistanceMatrix,[size(DistanceMatrix,1)^2,1]);
Cluster_1 = cluster(RBLinkage,'cutoff',RBLinkage(end,3),'criterion','distance');
IntraClusterList(1,1:end) = 1:size(DistanceMatrix,1);
Count1 = 0;
Bin = linspace(min(min(DistanceMatrix)),max(max(DistanceMatrix)),50);
for Linkage_i = 1:Inx
    Cluster_2 = cluster(RBLinkage,'cutoff',RBLinkage(end-Linkage_i,3),'criterion','distance');
    Cluster_2N = zeros(size(Cluster_2,1),1);
    for Clus_i = 1:max(Cluster_1)
        Cluster_2N(Cluster_2==Cluster_2(find(Cluster_1==Clus_i,1)))=Clus_i;
    end
    Cluster_2N(Cluster_2N==0)=max(Cluster_2);
    BreakCluster = unique([Cluster_1(Cluster_1-Cluster_2N~=0);Cluster_2N(Cluster_1-Cluster_2N~=0)]);
    
    ElementShifted = zeros(1,size(DistanceMatrix,1));
    ClusterRemovedInx = find(Cluster_2N==max(Cluster_2));
    Counter = 0;
    for Rem_i = 1:size(ClusterRemovedInx,1)
        Element = RobustClusterList(RobustClusterList(:,ClusterRemovedInx(Rem_i,1))>0,ClusterRemovedInx(Rem_i,1));
        ElementShifted(Counter+1:Counter+size(Element,1))=Element;
        Counter = Counter+size(Element,1);
    end
    ElementRemains = unique([0,setdiff(IntraClusterList(BreakCluster(1,1),:),ElementShifted)]);
    ElementRemains(:,1) = [];
    ElementShifted = unique([0,ElementShifted]);
    ElementShifted(:,1) = [];
%% Compute overlapping of distribution between Shifted cluster coefficent with Intra cluster coefficient and Inter cluster coefficient
    ShiftedCoefficient = DistanceMatrix(ElementShifted(ElementShifted~=0),ElementRemains(ElementRemains~=0));
    ShiftedCoefficient = reshape(ShiftedCoefficient,[size(ShiftedCoefficient,1)*size(ShiftedCoefficient,2),1]);
    [~,D] =intersect(IntraClusterCoefficient,ShiftedCoefficient);
    IntraClusterCoefficient(D,:)=[];
    IntraClusterList(BreakCluster(1,1),:)=0;
    IntraClusterList(BreakCluster(1,1),1:size(ElementRemains,2))=ElementRemains;
    IntraClusterList(BreakCluster(2,1),1:size(ElementShifted,2))=ElementShifted;
    IntreClusterCoefficient(Count1+1:Count1+size(ShiftedCoefficient,1),1)=ShiftedCoefficient;
    Count1 = Count1+size(ShiftedCoefficient,1);
    Cluster_1 = Cluster_2N;
    
    ShiftedCoefficientFrequency = histc(ShiftedCoefficient,Bin);
    ShiftedCoefficientFrequency = ShiftedCoefficientFrequency/sum(ShiftedCoefficientFrequency);
    ShiftedCoefficientFrequency = reshape(ShiftedCoefficientFrequency,[1,size(Bin,2)]);
    IntraClusterCoefficientFrequency = histc(IntraClusterCoefficient,Bin);
    IntraClusterCoefficientFrequency = IntraClusterCoefficientFrequency/sum(IntraClusterCoefficientFrequency);
    IntraClusterCoefficientFrequency = reshape(IntraClusterCoefficientFrequency,[1,size(Bin,2)]);
    IntreClusterCoefficientFrequency = histc(IntreClusterCoefficient(1:Count1),Bin);
    IntreClusterCoefficientFrequency = IntreClusterCoefficientFrequency/sum(IntreClusterCoefficientFrequency);
    IntreClusterCoefficientFrequency = reshape(IntreClusterCoefficientFrequency,[1,size(Bin,2)]);
    SimilarityMeasure(Linkage_i,1) = sum(min(ShiftedCoefficientFrequency,IntraClusterCoefficientFrequency));
    SimilarityMeasure(Linkage_i,2) = sum(min(ShiftedCoefficientFrequency,IntreClusterCoefficientFrequency));
    SimilarityMeasure(Linkage_i,3) = sum(min(IntraClusterCoefficientFrequency,IntreClusterCoefficientFrequency));
end
%% Determine the cluster threshold, when the Shifted cluster coefficent is closer to Intra cluster coefficient than Inter cluster coefficient
ClusterCutoffThreshold = RBLinkage(find(smooth((SimilarityMeasure(:,1)-SimilarityMeasure(:,2))./SimilarityMeasure(:,3),1)>0,1),3);
Cluster = cluster(RBLinkage,'cutoff',ClusterCutoffThreshold,'criterion','distance');
IdentifiedClusterList = zeros(size(DistanceMatrix,1),max(size(unique(Cluster))));
for i = 1:size(IdentifiedClusterList,2)
    Clust = unique(RobustClusterList(:,Cluster==i));
    Clust = unique([0;reshape(Clust,[size(Clust,1)*size(Clust,2),1])]);
    Clust(1,:) = [];
    IdentifiedClusterList(1:size(Clust,1),i) = Clust;
end
IdentifiedClusterList(sum(IdentifiedClusterList,2)==0,:) = [];
else
%% If number of robust cluster is less or equal to 3, the Identify cluster list set to be equal as robust cluster list
    IdentifiedClusterList = [];
    sprintf('Cluster list set to be equal as robust cluster list')
end
