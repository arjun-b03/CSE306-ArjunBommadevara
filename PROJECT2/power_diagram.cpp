#include <stdio.h>
#include <vector>
#include <string>
#include <cmath>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <chrono>
#include <random>
#include "stb_image.h"
#include "stb_image_write.h"
#include <cstdio>
#include "lbfgs.h"

static std::default_random_engine engine(10);
static std::uniform_real_distribution<double> uniform(0.0, 1.0);

class Vector {
public:
    explicit Vector(double x = 0, double y = 0, double z = 0) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    double norm2() const {
        return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
    }
    double norm() const {
        return sqrt(norm2());
    }
    void normalize() {
        double n = norm();
        data[0] /= n;
        data[1] /= n;
        data[2] /= n;
    }
    double operator[](int i) const { return data[i]; }
    double& operator[](int i) { return data[i]; }

private:
    double data[3];
};

Vector operator-(const Vector& a) {
    return Vector(-a[0], -a[1], -a[2]);
}

Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}

Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

Vector operator*(const double a, const Vector& b) {
    return Vector(a * b[0], a * b[1], a * b[2]);
}

Vector operator*(const Vector& a, const double b) {
    return Vector(a[0] * b, a[1] * b, a[2] * b);
}

Vector operator*(const Vector& a, const Vector& b) {
    return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}

Vector operator/(const Vector& a, const double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}

double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

class Polygon {
public:
    Polygon() = default;
    explicit Polygon(const std::vector<Vector>& vertices) : vertices(vertices) {}

    double area() const {
        if (vertices.size() < 3) return 0.0;

        double a = 0;
        int n = vertices.size();
        for (int i = 0; i < n; i++) {
            a += (vertices[i][0] * vertices[(i + 1) % n][1] - vertices[(i + 1) % n][0] * vertices[i][1]);
        }
        return 0.5 * std::abs(a);
    }

    double square_distance_integral(const Vector& point) const {
        if (vertices.size() < 3) return 0.0;

        double result = 0.0;
        int vertex_n = vertices.size();

        for (int i = 1; i < vertex_n - 1; ++i) {
            Vector triangle[3] = {vertices[0], vertices[i], vertices[i + 1]};

            double integral = 0.0;
            for (int k = 0; k < 3; ++k) {
                for (int l = k; l < 3; ++l) {
                    integral += dot(triangle[k] - point, triangle[l] - point);
                }
            }

            Vector edge1 = triangle[1] - triangle[0];
            Vector edge2 = triangle[2] - triangle[0];
            double triangleArea = 0.5 * std::abs(edge1[1] * edge2[0] - edge1[0] * edge2[1]);

            result += (integral / 6.0) * triangleArea;
        }

        return result;
    }

    void addVertex(const Vector& vertex) {
        vertices.push_back(vertex);
    }

    std::vector<Vector> vertices;
};

void save_svg(const std::vector<Polygon>& polygons, std::string filename, std::string fillcol = "none") {
    FILE* f = fopen(filename.c_str(), "w+");
    fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
    for (int i = 0; i < polygons.size(); i++) {
        fprintf(f, "<g>\n");
        fprintf(f, "<polygon points = \"");
        for (int j = 0; j < polygons[i].vertices.size(); j++) {
            fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
        }
        fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
        fprintf(f, "</g>\n");
    }
    fprintf(f, "</svg>\n");
    fclose(f);
}

Polygon clipPowerDiagram(const Polygon& subjectPolygon, const Vector& siteA, double weightA, const Vector& siteB, double weightB) {
    Polygon polygonfinal;
    int vertices_n = subjectPolygon.vertices.size();

    for (int i = 0; i < vertices_n; ++i) {
        const Vector& vertex_i = subjectPolygon.vertices[i];
        const Vector& vertex_prev = subjectPolygon.vertices[(i - 1 + vertices_n) % vertices_n];

        double powerDistanceCurrentA = (vertex_i - siteA).norm2() - weightA;
        double powerDistanceCurrentB = (vertex_i - siteB).norm2() - weightB;
        double powerDistancePreviousA = (vertex_prev - siteA).norm2() - weightA;
        double powerDistancePreviousB = (vertex_prev - siteB).norm2() - weightB;

        if (powerDistanceCurrentA <= powerDistanceCurrentB) {
            if (powerDistancePreviousA > powerDistancePreviousB) {
                Vector midPoint = (siteA + siteB) * 0.5;
                double param = dot(midPoint - vertex_prev, siteB - siteA) / dot(vertex_i - vertex_prev, siteB - siteA);
                Vector intersection = vertex_prev + param * (vertex_i - vertex_prev);
                polygonfinal.vertices.push_back(intersection);
            }
            polygonfinal.vertices.push_back(vertex_i);
        } else {
            if (powerDistancePreviousA <= powerDistancePreviousB) {
                Vector midPoint = (siteA + siteB) * 0.5;
                double param = dot(midPoint - vertex_prev, siteB - siteA) / dot(vertex_i - vertex_prev, siteB - siteA);
                Vector intersection = vertex_prev + param * (vertex_i - vertex_prev);
                polygonfinal.vertices.push_back(intersection);
            }
        }
    }

    return polygonfinal;
}

class PowerDiagram {
public:
    PowerDiagram() = default;
    PowerDiagram(const std::vector<Vector>& points, const std::vector<double>& weights) : points(points), weights(weights) {}

    void compute() {
        Polygon initialSquare;
        initialSquare.addVertex(Vector(0, 0, 0));
        initialSquare.addVertex(Vector(1, 0, 0));
        initialSquare.addVertex(Vector(1, 1, 0));
        initialSquare.addVertex(Vector(0, 1, 0));

        cells.resize(points.size());
        #pragma omp parallel for
        for (int i = 0; i < points.size(); i++) {
            Polygon cell = initialSquare;
            for (int j = 0; j < points.size(); j++) {
                if (i != j) {
                    cell = clipPowerDiagram(cell, points[i], weights[i], points[j], weights[j]);
                }
            }
            cells[i] = cell;
        }
    }

    std::vector<Vector> points;
    std::vector<double> weights;
    std::vector<Polygon> cells;
};

class SemiOptimalTransport {
public:
    std::vector<Vector> points;
    std::vector<double> lambdas;
    Polygon edges;
    std::vector<Polygon> polygons;
    std::vector<double> weights;

    SemiOptimalTransport(const std::vector<Vector>& points, const std::vector<double>& lambdas, const Polygon& edges)
        : points(points), lambdas(lambdas), edges(edges), weights(points.size(), 1.0) {}

    ~SemiOptimalTransport() {}

    void optimize(int n) {
        std::cout << "Entering optimize" << std::endl;
        double fx = 0.0;
        lbfgs(n, &weights[0], &fx, _evaluate, _progress, this, NULL);
        computePolygons();
    }

    static lbfgsfloatval_t _evaluate(
        void* instance,
        const lbfgsfloatval_t* x,
        lbfgsfloatval_t* g,
        const int n,
        const lbfgsfloatval_t step
    ) {
        return reinterpret_cast<SemiOptimalTransport*>(instance)->evaluate(x, g, n, step);
    }

    lbfgsfloatval_t evaluate(
        const lbfgsfloatval_t* x,
        lbfgsfloatval_t* g,
        const int n,
        const lbfgsfloatval_t step
    ) {
        lbfgsfloatval_t fx = 0.0;

        for (int i = 0; i < n; i++)
            weights[i] = x[i];

        computePolygons();

        double t1 = 0, t2 = 0, t3 = 0;
        for (int i = 0; i < n; i++) {
            g[i] = polygons[i].area() - lambdas[i];
            t1 += polygons[i].square_distance_integral(points[i]);
            t2 -= x[i] * polygons[i].area();
            t3 += x[i] * lambdas[i];
        }
        fx = -(t1 + t2 + t3);

        // Debugging output
        std::cout << "Evaluation: fx = " << fx << std::endl;
        std::cout << "Areas: ";
        for (const auto& poly : polygons) {
            std::cout << poly.area() << " ";
        }
        std::cout << std::endl;
        std::cout << "Gradients: ";
        for (int i = 0; i < n; i++) {
            std::cout << g[i] << " ";
        }
        std::cout << std::endl;

        return fx;
    }

    static int _progress(
        void* instance,
        const lbfgsfloatval_t* x,
        const lbfgsfloatval_t* g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
    ) {
        return reinterpret_cast<SemiOptimalTransport*>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
    }

    int progress(
        const lbfgsfloatval_t* x,
        const lbfgsfloatval_t* g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
    ) {
        std::cout << "Iteration " << k << " : fx = " << fx << ", xnorm = " << xnorm << ", gnorm = " << gnorm << ", step = " << step << std::endl;
        return 0;
    }

private:
    void computePolygons() {
        polygons.resize(points.size());
        Polygon initialPolygon = edges;

        for (size_t i = 0; i < points.size(); ++i) {
            Polygon cell = initialPolygon;
            for (size_t j = 0; j < points.size(); ++j) {
                if (i != j) {
                    cell = clipPowerDiagram(cell, points[i], weights[i], points[j], weights[j]);
                }
            }
            polygons[i] = cell;
        }
    }
};

int main() {
    auto start = std::chrono::high_resolution_clock::now();

    int n = 10;
    std::vector<double> weights(n);
    std::vector<Vector> points(n);
    std::vector<double> lambdas(n, 1.0 / n);

    for (int i = 0; i < n; i++) {
        points[i] = Vector(uniform(engine), uniform(engine), 0);
        weights[i] = 1.0;
    }

    Polygon edges({
        Vector(0., 0.), Vector(0., 1.),
        Vector(1., 1.), Vector(1., 0.)
    });

    SemiOptimalTransport ot(points, lambdas, edges);
    ot.optimize(n);

    save_svg(ot.polygons, "/Users/arjunbommadevara/Desktop/CSE306/project1/CSE306-Project-1/PROJECT2/ot_1000.svg");

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    std::cout << "Elapsed time: " << duration.count() << " milliseconds" << std::endl;

    return 0;
}
