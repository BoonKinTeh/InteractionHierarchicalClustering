function [SerialIndex,Result_A] = A_HierarchicalClustering(DistanceMatrix,ClusteringMethod)
%% Inputs description:
% DistanceMatrix: NxN double array, contains distance between each pair of variables (total N variables),
%                 it is a symmetry matrix with zeros diagonal entries. 
% ClusteringMethod: string, Method to define the distance between clusters
%           complete = Maximum distance
%           average = Average distance
%           single = minimum distance
%           centroid = Euclidean distance
%% Default inputs description:
PlotFigure_OrginalAndHCSeriatedDistanceMatrix = 1; % Set to 1 to plot the colormap for distance matrix (orginal and hierarchical seriated) 
%% Outputs description:
% SerialIndex: 1xN integer array, the dendrogram order obtained for performing hierarchical clustering
% Result_A: structure, input for B_DetermineRobustClusters
% Linkage: N-1x3 double array, each row represents the cluster with column 1 index merged with
%          the cluster with column 2 index at distance in column 3.
%          Cluster index, CI: if CI <= N, CI represent the index for the N elements
%                             if CI > N, represent clusters merged at row CI-N
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
%% Hierarchical Clustering
PDist = squareform(DistanceMatrix);
Linkage = linkage(PDist,ClusteringMethod);
%% Hierarchical Clustering Dendrogram
[~,~,SerialIndex] = dendrogram(Linkage,0);

%% Plot Original and HC-Seriated Distance Matrix
if PlotFigure_OrginalAndHCSeriatedDistanceMatrix == 1
    figure(1);clf;hold on;
    subplot(1,2,1);hold on;
    imagesc(-DistanceMatrix);
    colormap 'jet';
    xlim([0.5,size(DistanceMatrix,2)+0.5]);
    ylim([0.5,size(DistanceMatrix,1)+0.5]);
    title('Original Distance Matrix','fontsize',16);
    axis('square');
    
    subplot(1,2,2);hold on;
    imagesc(-DistanceMatrix(SerialIndex,SerialIndex));
    colormap 'jet';
    xlim([0.5,size(DistanceMatrix,2)+0.5]);
    ylim([0.5,size(DistanceMatrix,1)+0.5]);
    title('HC-Seriated Distance Matrix','fontsize',16);
    axis('square');
end
%% Output Result
Result_A.PDist = PDist;
Result_A.Linkage = Linkage;
Result_A.DistanceMatrix = DistanceMatrix;