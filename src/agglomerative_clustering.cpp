#include "agglomerative_clustering.h"
#include "util.h"

#include <cmath>
#include <stack>
#include <queue>
#include <algorithm>
#include <numeric>

namespace ahc::impl {

// DistanceCall and LinkageDistanceUpdateCall

Distance UpdateDistanceAverage(Distance dxi, Distance dyi, Distance dxy, size_t size_x, size_t size_y, size_t size_i) {
    return (dxi * static_cast<Distance>(size_x) + dyi * static_cast<Distance>(size_y)) / static_cast<Distance>(size_x + size_y);
}

Distance CosineDistance(const Point& p1, const Point& p2) {
    AHC_ASSERT(p1.size() == p2.size());

    const auto dot = std::inner_product(p1.begin(), p1.end(), p2.begin(), 0.f);

    auto norm1 = std::inner_product(p1.begin(), p1.end(), p1.begin(), 0.f);
    AHC_ASSERT(norm1 > 0);
    norm1 = std::sqrt(norm1);

    auto norm2 = std::inner_product(p2.begin(), p2.end(), p2.begin(), 0.f);
    AHC_ASSERT(norm2 > 0);
    norm2 = std::sqrt(norm2);

    return 1.f - dot / (norm1 * norm2);
}


// Agglomerative Clustering


AgglomerativeClustering::AgglomerativeClustering(const std::vector<Point>& points, DistanceMetric metric) {
    m_distanceMatrix = DistanceMatrix::ComputeMatrix(points, metric);
}

void AgglomerativeClustering::RunClustering(Linkage linkage) {
    AHC_ASSERT_MSG(linkage == Linkage::AVERAGE, "Linkage not implemented");
    m_linkageTree = LinkageTree::ComputeTree(*m_distanceMatrix, UpdateDistanceAverage);
}

std::vector<Label> AgglomerativeClustering::ComputeLabels(size_t numClusters) const {
    AHC_ASSERT_MSG(m_linkageTree != nullptr, "Cannot compute levels before clustering");
    return m_linkageTree->ComputeLabels(numClusters);
}

std::vector<Label> AgglomerativeClustering::ComputeLabelsWithDistance(Distance threshold) const {
    AHC_ASSERT_MSG(m_linkageTree != nullptr, "Cannot compute levels before clustering");
    return m_linkageTree->ComputeLabelsWithDistance(threshold);
}

// Distance matrix

DistanceMatrix::DistanceMatrix(std::vector<float> data, size_t numOfPoints)
    : m_data(std::move(data)), m_numOfPoints(numOfPoints) {}



void DistanceMatrix::ValidatePoints(const std::vector<Point>& points) {
    AHC_ASSERT_MSG(!points.empty(), "Input points are empty");
    const size_t numDim = points[0].size();
    AHC_ASSERT_MSG(numDim != 0, "Input points have invalid dimensions");
    for(size_t i = 1; i < points.size(); ++i) {
        AHC_ASSERT_MSG(points[i].size() == numDim, "Input points must have the same dimensions");
    }
}

std::vector<Distance> DistanceMatrix::ComputeCondensedMatrix(const std::vector<Point>& points,
                                                          const DistanceCall& computeDistance)
{
    // The sum of an arithmetic progression
    // S = n/2 * [2a + (n-1)d] / 2; a == 1; d == 1
    // n - number of pairwise distances
    const size_t n = points.size() - 1;
    // number of non-zero distances in a upper triangle matrix
    const size_t dataSize = n * (n + 1) / 2;
    std::vector<Distance> condensedData;
    condensedData.reserve(dataSize);
    for (size_t i = 0; i < points.size(); ++i) {
        for (size_t j = i + 1; j < points.size(); ++j) {
            condensedData.push_back(computeDistance(points[i], points[j]));
        }
    }
    return condensedData;
}

std::unique_ptr<DistanceMatrix> DistanceMatrix::ComputeMatrix(const std::vector<Point>& points, DistanceMetric metric) {
    AHC_ASSERT_MSG(metric == DistanceMetric::COSINE, "Metric not implemented");
    ValidatePoints(points);
    const size_t numOfPoints = points.size();
    return std::unique_ptr<DistanceMatrix>(new DistanceMatrix(ComputeCondensedMatrix(points, CosineDistance), numOfPoints));
}

size_t DistanceMatrix::GetCondensedIndex(size_t fromCluster, size_t toCluster) const {
    AHC_ASSERT(fromCluster != toCluster);
    if (fromCluster > toCluster) std::swap(fromCluster, toCluster);
    // (fromCluster * m_dim + toCluster) is an index in a full square matrix
    // (fromCluster * (fromCluster + 1)) is a number of zeros before the row 'fromCluster' in a full upper triangle matrix
    // (fromCluster + 1) is a number of zeros in the row 'fromCluster' in a full upper triangle matrix
    const auto condensedIndex = fromCluster * m_numOfPoints + toCluster - fromCluster * (fromCluster + 1) / 2 - fromCluster - 1;
    AHC_ASSERT(condensedIndex < m_data.size());
    return condensedIndex;
}

Distance DistanceMatrix::GetDistance(size_t fromCluster, size_t toCluster) const {
    return m_data[GetCondensedIndex(fromCluster, toCluster)];
}

void DistanceMatrix::UpdateDistance(size_t fromCluster, size_t toCluster, Distance value) {
    m_data[GetCondensedIndex(fromCluster, toCluster)] = value;
}

// Linkage Tree

LinkageTree::LinkageTree(std::vector<Children> tree)
    : m_tree(std::move(tree))
{
}


std::unique_ptr<LinkageTree> LinkageTree::ComputeTree(DistanceMatrix distanceMatrix, const LinkageDistanceUpdateCall& distanceUpdate) {
    constexpr auto INF_DISTANCE = std::numeric_limits<Distance>::max();
    const auto numOfPoints = distanceMatrix.GetNumOfPoints();
    AHC_ASSERT(numOfPoints > 0);

    const auto numberOfClusterMerges = numOfPoints - 1;

    std::vector<Children> tree(numberOfClusterMerges);
    std::vector<size_t> clusterSizes(numOfPoints, 1);
    std::stack<Node> clusterChain;

    // compute the tree, finding clusters to merge (children)
    for (auto& children: tree) {

        // chose any non-empty cluster
        if (clusterChain.empty()) {
            const size_t firstNonEmptyCluster = std::find_if(clusterSizes.begin(), clusterSizes.end(),
                                                             [](auto size)  {return size > 0;}) - clusterSizes.begin();
            clusterChain.push(firstNonEmptyCluster);
        }

        // clusters to merge
        Node cluster1;
        Node cluster2;
        Distance mergeDistance;

        // build the chain until the two last clusters in it are the mutual nearest neighbours
        while (true) {

            const Node currentCluster = clusterChain.top();
            clusterChain.pop();

            // find the nearest neighbour to the currentCluster
            Node nearestNeighbour;
            Distance currentMinDistance { INF_DISTANCE };
            // start from previous cluster if it exists
            // currentCluster is the nearest cluster to the previous cluster
            if (!clusterChain.empty()) {
                nearestNeighbour = clusterChain.top();
                currentMinDistance = distanceMatrix.GetDistance(currentCluster, nearestNeighbour);
            }

            // find the nearest cluster
            for (Node cluster = 0; cluster < numOfPoints; ++cluster) {
                // skip an empty cluster and the currentCluster itself
                if (clusterSizes[cluster] == 0 || currentCluster == cluster) continue;

                const auto distance = distanceMatrix.GetDistance(currentCluster, cluster);
                if (distance < currentMinDistance) {
                    currentMinDistance = distance;
                    nearestNeighbour = cluster;
                }
            }

            // the currentCluster and the nearestNeighbour are the nearest mutual neighbors
            if (!clusterChain.empty() && nearestNeighbour == clusterChain.top()) {
                cluster1 = currentCluster;
                cluster2 = nearestNeighbour;
                mergeDistance = currentMinDistance;
                // pop the nearestNeighbour
                clusterChain.pop();
                break;
            }

            // push back the current cluster and its nearest neighbour
            // they are not the mutual nearest neighbours
            clusterChain.push(currentCluster);
            clusterChain.push(nearestNeighbour);
        }

        AHC_ASSERT(cluster1 != cluster2);
        if (cluster1 > cluster2) std::swap(cluster1, cluster2);

        // update the tree
        children.cluster1 = cluster1;
        children.cluster2 = cluster2;
        children.distance = mergeDistance;

        // merge the two clusters into the one cluster2,
        // nullify the first cluster, so it will be skipped later
        const auto cluster1Size = clusterSizes[cluster1];
        const auto cluster2Size = clusterSizes[cluster2];
        clusterSizes[cluster2] += cluster1Size;
        clusterSizes[cluster1] = 0;

        // update the distance matrix
        // update distances from all points to the new cluster
        for (size_t cluster = 0; cluster < numOfPoints; ++cluster) {
            const auto clusterSize = clusterSizes[cluster];
            // skip empty cluster of the new cluster itself
            if (clusterSize == 0 || cluster == cluster2) continue;

            const auto newDistance = distanceUpdate(
                    distanceMatrix.GetDistance(cluster, cluster1),
                    distanceMatrix.GetDistance(cluster, cluster2),
                    mergeDistance,
                    cluster1Size,
                    cluster2Size,
                    clusterSize);

            distanceMatrix.UpdateDistance(cluster, cluster2, newDistance);
        }
    }

    // sort by merge distance
    std::stable_sort(tree.begin(), tree.end(),
                     [](const auto& lhs, const auto& rhs) { return lhs.distance < rhs.distance;});

    CorrectLabel(tree);

    return std::unique_ptr<LinkageTree>(new LinkageTree(std::move(tree)));
}


void LinkageTree::CorrectLabel(std::vector<Children>& tree) {

    const auto numOfLeaves = tree.size() + 1;

    const auto NO_PARENT = std::numeric_limits<Node>::max();
    // there are 2 * numOfLeaves - 1 parent nodes in the tree;
    // leaf node is a parent for itself, as result there are numOfLeaves parents for the each leaf;
    // there are numOfLeaves - 1 merge nodes and each merge node is a parent for nodes that have been merge into it
    std::vector<Node> parent(2 * numOfLeaves - 1, NO_PARENT);

    // there are numOfLeaves labels that are leaves itself,
    // therefore, the label of the first merge node is numOfLeaves
    auto nextNode = numOfLeaves;

    // find the latest parent of the node
    auto FindParent = [&parent, NO_PARENT](Node node)-> Node {
        while (parent[node] != NO_PARENT) {
            node = parent[node];
        }
        return node;
    };

    // merge into mergeNode
    // mergeNode becomes a parent of the node1 and node2
    auto Merge = [&parent](Node node1, Node node2, Node mergeNode) {
        parent[node1] = mergeNode;
        parent[node2] = mergeNode;
    };

    for (auto& children: tree) {
        auto parent1 = FindParent(children.cluster1);
        auto parent2 = FindParent(children.cluster2);
        AHC_ASSERT(parent1 != parent2);
        if (parent1 > parent2) std::swap(parent1, parent2);
        children.cluster1 = parent1;
        children.cluster2 = parent2;
        Merge(parent1, parent2, nextNode);
        ++nextNode;
    }
}

std::vector<Label> LinkageTree::ComputeLabelsWithDistance(Distance threshold) const {
    const auto isSorted =  std::is_sorted(m_tree.begin(), m_tree.end(), [](const auto& lhs, const auto& rhs) {
       return  lhs.distance < rhs.distance;
    });
    AHC_ASSERT(isSorted);
    const auto numClusters = static_cast<size_t>(m_tree.end() - std::find_if(m_tree.begin(), m_tree.end(), [threshold](const auto& children) {
       return children.distance >= threshold;
    })) + 1;
    return ComputeLabels(numClusters);
}

std::vector<Label> LinkageTree::ComputeLabels(size_t numClusters) const {
    AHC_ASSERT_MSG(numClusters > 0, "Invalid number of clusters");
    const auto numOfPoints = m_tree.size() + 1;
    AHC_ASSERT_MSG(numClusters <= numOfPoints, "Invalid number of clusters");

    if (numClusters == 1) return std::vector<Label>(numOfPoints);

    const auto& rootChildren = m_tree.back();
    const auto treeRootNode = std::max(rootChildren.cluster1, rootChildren.cluster2) + 1;

    std::priority_queue<Node> queue;
    queue.push(treeRootNode);

    for (size_t currentNumOfClusters = 2; currentNumOfClusters <= numClusters; ++currentNumOfClusters) {
        const auto rootNode = queue.top();
        // if rootNode == numOfPoints then labels are points itself
        AHC_ASSERT(rootNode >= numOfPoints);
        const auto children = m_tree[rootNode - numOfPoints];
        AHC_ASSERT(rootNode > children.cluster1);
        AHC_ASSERT(rootNode > children.cluster2);
        AHC_ASSERT(children.cluster1 != children.cluster2);
        queue.push(children.cluster1);
        queue.push(children.cluster2);
        // remove root
        queue.pop();
    }

    std::vector<Node> nodeSet;
    while(!queue.empty()) {
        nodeSet.push_back(queue.top());
        queue.pop();
    }

    return ComputeLabels(nodeSet);
}

std::vector<Label> LinkageTree::ComputeLabels(const std::vector<Node>& nodeSet) const {
    AHC_ASSERT(nodeSet.size() > 1);
    AHC_ASSERT(!m_tree.empty());

    const auto  numOfPoints = m_tree.size() + 1;
    std::vector<Label> labels(numOfPoints);

    for (size_t i = 0; i < nodeSet.size(); ++i) {
        const auto nodeLeaves = GetAllLeaves(nodeSet[i]);
        const auto newLabel = static_cast<Label>(i);
        for (auto leaf: nodeLeaves) {
            AHC_ASSERT(leaf < labels.size());
            labels[leaf] = newLabel;
        }
    }
    return labels;
}

std::vector<LinkageTree::Node> LinkageTree::GetAllLeaves(LinkageTree::Node node) const {
    const auto  numOfPoints = m_tree.size() + 1;

    // node is a leaf itself
    if (node < numOfPoints) return {node};

    std::vector<Node> leaves;

    std::stack<Node> nodeStack;
    nodeStack.push(node);

    // descend the node
    while (!nodeStack.empty()) {
        const auto currentNode = nodeStack.top();
        nodeStack.pop();
        if (currentNode < numOfPoints) {
            // it is a leaf
            leaves.push_back(currentNode);
        } else {
            // add its children to process
            const auto& children = m_tree[currentNode - numOfPoints];
            nodeStack.push(children.cluster1);
            nodeStack.push(children.cluster2);
        }
    }
     return leaves;
}

}
