#ifndef AHC_AGGLOMERATIVE_CLUSTERING_H
#define AHC_AGGLOMERATIVE_CLUSTERING_H

#include "ahc/ahc.h"
#include <functional>

namespace ahc::impl {

using DistanceCall = std::function<Distance(const Point&, const Point&)>;
using LinkageDistanceUpdateCall = std::function<Distance(Distance, Distance, Distance, size_t, size_t, size_t)>;


class DistanceMatrix;
class LinkageTree;


class AgglomerativeClustering : public IAgglomerativeClustering {
public:
    AgglomerativeClustering(const std::vector<Point>& points, DistanceMetric metric);

    void RunClustering(Linkage linkage) final;

    std::vector<Label> ComputeLabels(size_t numClusters) const final;
    std::vector<Label> ComputeLabelsWithDistance(Distance threshold) const final;

    ~AgglomerativeClustering() override = default;

private:
    std::unique_ptr<DistanceMatrix> m_distanceMatrix;
    std::unique_ptr<LinkageTree> m_linkageTree;
};


class DistanceMatrix {
public:
    Distance GetDistance(size_t fromCluster, size_t toCluster) const;

    void UpdateDistance(size_t fromCluster, size_t toCluster, Distance value);

    size_t GetNumOfPoints() { return m_numOfPoints; }

    static std::unique_ptr<DistanceMatrix> ComputeMatrix(const std::vector<Point>& points, DistanceMetric);

    ~DistanceMatrix() = default;

private:
    explicit DistanceMatrix(std::vector<Distance> condensedData, size_t numOfPoints);
    std::vector<Distance> m_data;
    const size_t m_numOfPoints;

    static std::vector<Distance> ComputeCondensedMatrix(const std::vector<Point>& points, const DistanceCall& computeDistance);
    static void ValidatePoints(const std::vector<Point>& points);
    size_t GetCondensedIndex(size_t fromCluster, size_t toCluster) const;
};


class LinkageTree {
public:
    using Node = size_t;

    struct Children {
        Node cluster1;
        Node cluster2;
        Distance distance;
    };

    std::vector<Label> ComputeLabels(size_t numClusters) const;
    std::vector<Label> ComputeLabelsWithDistance(Distance threshold) const;

    static std::unique_ptr<LinkageTree> ComputeTree(DistanceMatrix matrix, const LinkageDistanceUpdateCall& distanceUpdate);

private:
    explicit LinkageTree(std::vector<Children> tree);

    static void CorrectLabel(std::vector<Children>& tree);
    std::vector<Label> ComputeLabels(const std::vector<Node>& nodeSet) const;
    std::vector<Node> GetAllLeaves(Node node) const;

    const std::vector<Children> m_tree;
};


}


#endif // AHC_AGGLOMERATIVE_CLUSTERING_H
