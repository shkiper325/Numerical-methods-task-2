#include <vector>

#include <cmath>
#include <cstdio>

using namespace std;

#include "Point.hpp"
#include "utils.hpp"

double C_dist(vector<double> f, vector<double> g) {
    double ret = 0;

    for (int i = 0; i < f.size(); ++i) {
        ret = max(ret, fabs(f[i] - g[i]));
    }

    return ret;
}

vector<vector<double> > p2vec(vector<Point> vec) {
    auto x = vector<double>(vec.size());
    auto p = vector<double>(vec.size());

    for (int i = 0; i < vec.size(); ++i) {
        x[i] = vec[i].x;
        p[i] = vec[i].p;
    }

    return vector<vector<double> >({x, p});
}

void dump_vec(vector<double> vec, string path) {
    auto fd = fopen(path.c_str(), "w");
    auto m = vec.size();

    for (int i = 0; i < m; ++i) {
        fprintf(fd, "%30.20e\n", vec[i]);
    }

    fclose(fd);
}