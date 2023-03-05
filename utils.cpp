#include <vector>
#include <functional>

#include <cmath>
#include <cstdio>

using namespace std;

#include "Point.hpp"
#include "utils.hpp"

long double C_dist(vector<long double> f, vector<long double> g) {
    long double ret = 0;

    for (int i = 0; i < f.size(); ++i) {
        ret = max(ret, fabs(f[i] - g[i]));
    }

    return ret;
}

vector<vector<long double> > p2vec(vector<Point> vec) {
    auto x = vector<long double>(vec.size());
    auto p = vector<long double>(vec.size());

    for (int i = 0; i < vec.size(); ++i) {
        x[i] = vec[i].x;
        p[i] = vec[i].p;
    }

    return vector<vector<long double> >({x, p});
}

void dump_vec(vector<long double> vec, string path) {
    auto fd = fopen(path.c_str(), "w");
    auto m = vec.size();

    for (int i = 0; i < m; ++i) {
        fprintf(fd, "%30.20Lf\n", vec[i]);
    }

    fclose(fd);
}

long double integrate(vector<long double> t, vector<long double> x, function<long double(long double)> f) {
    long double ret = 0;

    for(int i = 0; i < t.size() - 1; ++i) {
        ret += (x[i] * f(t[i]) + x[i + 1] * f(t[i + 1])) / 2 * (t[i + 1] - t[i]);
    }

    return ret;
}

vector<long double> jacobian(function<vector<long double>(long double, long double, long double)> f, long double x, long double y, long double z, long double eps, long double delta) {
    vector<long double> ret(9);

    vector<long double> values_right = f(x + delta, y, z);
    vector<long double> values_left  = f(x - delta, y, z);

    ret[0] = (values_right[0] - values_left[0]) / (2 * eps);
    ret[3] = (values_right[1] - values_left[1]) / (2 * eps);
    ret[6] = (values_right[2] - values_left[2]) / (2 * eps);

    values_right = f(x, y + delta, z);
    values_left  = f(x, y - delta, z);

    ret[1] = (values_right[0] - values_left[0]) / (2 * eps);
    ret[4] = (values_right[1] - values_left[1]) / (2 * eps);
    ret[7] = (values_right[2] - values_left[2]) / (2 * eps);

    values_right = f(x, y, z + delta);
    values_left  = f(x, y, z - delta);

    ret[2] = (values_right[0] - values_left[0]) / (2 * eps);
    ret[5] = (values_right[1] - values_left[1]) / (2 * eps);
    ret[8] = (values_right[2] - values_left[2]) / (2 * eps);

    return ret;
}

vector<long double> inverse3x3(vector<long double> mat) {
    long double a = mat[0];
    long double b = mat[1];
    long double c = mat[2];
    long double d = mat[3];
    long double e = mat[4];
    long double f = mat[5];
    long double g = mat[6];
    long double h = mat[7];
    long double i = mat[8];

    long double det = a * e * i - a * f * h - b * d * i + b * f * g + c * d * h - c * e * g;

    vector<long double> ret(9);

    ret[0] = (e * i - f * h) / det;
    ret[1] = (c * h - b * i) / det;
    ret[2] = (b * f - c * e) / det;
    ret[3] = (f * g - d * i) / det;
    ret[4] = (a * i - c * g) / det;
    ret[5] = (c * d - a * f) / det;
    ret[6] = (d * h - e * g) / det;
    ret[7] = (b * g - a * h) / det;
    ret[8] = (a * e - b * d) / det;

    return ret;
}

 void print_mat3x3(vector<long double> mat) {
    printf("%12.6Lf %12.6Lf %12.6Lf\n", mat[0], mat[1], mat[2]);
    printf("%12.6Lf %12.6Lf %12.6Lf\n", mat[3], mat[4], mat[5]);
    printf("%12.6Lf %12.6Lf %12.6Lf\n", mat[6], mat[7], mat[8]);
}

 void test_jacobian_and_inverse3x3() {
    auto f = [](long double x, long double y, long double z) -> vector<long double> {return {x + 2 * y * z, x * x+ y * y + z, sqrt(x * y * z)};};
    
    auto mat = jacobian(f, 2, 3, 4, 1e-7, 1e-7);
    mat = inverse3x3(mat);

    print_mat3x3(mat);
}

 vector<long double> apply_mat_3x3to3 (vector<long double> mat, vector<long double> vec) {
    vector<long double> ret(3);

    ret[0] = mat[0] * vec[0] + mat[1] * vec[1] + mat[2] * vec[2];
    ret[1] = mat[3] * vec[0] + mat[4] * vec[1] + mat[5] * vec[2];
    ret[2] = mat[6] * vec[0] + mat[7] * vec[1] + mat[8] * vec[2];

    return ret;
}