#include "ahc/ahc.h"
#include "agglomerative_clustering.h"

#include <memory>

namespace ahc {

std::unique_ptr<ahc::IAgglomerativeClustering> IAgglomerativeClustering::Create(const std::vector<Point>& points, DistanceMetric metric) {
    return std::make_unique<impl::AgglomerativeClustering>(points, metric);
}

}
