#include <vector>
#include <functional>

#include <cmath>
#include <cstdio>

#include <quadmath.h>

using namespace std;

#include "Point.hpp"
#include "utils.hpp"
#include "number_type.hpp"

NT C_dist(vector<NT> f, vector<NT> g) {
    NT ret = 0;

    for (int i = 0; i < f.size(); ++i) {
        ret = maxq(ret, fabsq(f[i] - g[i]));
    }

    return ret;
}

vector<vector<NT> > p2vec(vector<Point> vec) {
    auto x = vector<NT>(vec.size());
    auto p = vector<NT>(vec.size());

    for (int i = 0; i < vec.size(); ++i) {
        x[i] = vec[i].x;
        p[i] = vec[i].p;
    }

    return vector<vector<NT> >({x, p});
}

void dump_vec(vector<NT> vec, string path) {
    auto fd = fopen(path.c_str(), "w");
    auto m = vec.size();

    for (int i = 0; i < m; ++i) {
        fprintf(fd, "%30.20lf\n", double(vec[i]));
    }

    fclose(fd);
}

NT integrate(vector<NT> t, vector<NT> x, function<NT(NT)> f) {
    NT ret = 0;

    for(int i = 0; i < t.size() - 1; ++i) {
        ret += (x[i] * f(t[i]) + x[i + 1] * f(t[i + 1])) / 2 * (t[i + 1] - t[i]);
    }

    return ret;
}

vector<NT> jacobian(function<vector<NT>(NT, NT, NT)> f, NT x, NT y, NT z, NT eps, NT delta) {
    vector<NT> ret(9);

    vector<NT> values_right = f(x + delta, y, z);
    vector<NT> values_left  = f(x - delta, y, z);

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

vector<NT> inverse3x3(vector<NT> mat) {
    NT a = mat[0];
    NT b = mat[1];
    NT c = mat[2];
    NT d = mat[3];
    NT e = mat[4];
    NT f = mat[5];
    NT g = mat[6];
    NT h = mat[7];
    NT i = mat[8];

    NT det = a * e * i - a * f * h - b * d * i + b * f * g + c * d * h - c * e * g;

    vector<NT> ret(9);

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

 void print_mat3x3(vector<NT> mat) {
    printf("%12.6lf %12.6lf %12.6lf\n", double(mat[0]), double(mat[1]), double(mat[2]));
    printf("%12.6lf %12.6lf %12.6lf\n", double(mat[3]), double(mat[4]), double(mat[5]));
    printf("%12.6lf %12.6lf %12.6lf\n", double(mat[6]), double(mat[7]), double(mat[8]));
}

 void test_jacobian_and_inverse3x3() {
    auto f = [](NT x, NT y, NT z) -> vector<NT> {return {x + 2 * y * z, x * x+ y * y + z, static_cast<NT>(sqrtq(x * y * z))};};
    
    auto mat = jacobian(f, 2, 3, 4, 1e-7, 1e-7);
    mat = inverse3x3(mat);

    print_mat3x3(mat);
}

 vector<NT> apply_mat_3x3to3 (vector<NT> mat, vector<NT> vec) {
    vector<NT> ret(3);

    ret[0] = mat[0] * vec[0] + mat[1] * vec[1] + mat[2] * vec[2];
    ret[1] = mat[3] * vec[0] + mat[4] * vec[1] + mat[5] * vec[2];
    ret[2] = mat[6] * vec[0] + mat[7] * vec[1] + mat[8] * vec[2];

    return ret;
}

NT maxq(NT a, NT b) {
    return a > b ? a : b;
}