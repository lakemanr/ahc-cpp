#ifndef AHC_AHC_H
#define AHC_AHC_H

#include <vector>
#include <memory>

namespace ahc {

using Point = std::vector<float>;
using Label = int;
using Distance = float;

enum class DistanceMetric {
    COSINE
};

enum class Linkage {
    AVERAGE
};


class IAgglomerativeClustering {
public:

    static std::unique_ptr<IAgglomerativeClustering> Create(const std::vector<Point>& points, DistanceMetric metric);

    virtual void RunClustering(Linkage linkage) = 0;

    virtual std::vector<Label> ComputeLabels(size_t numOfClusters) const = 0;
    virtual std::vector<Label> ComputeLabelsWithDistance(Distance threshold) const = 0;

    virtual ~IAgglomerativeClustering() = default;
};

}
#endif // AHC_AHC_H
