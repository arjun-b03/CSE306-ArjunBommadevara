#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <random>
#include <cstdio>

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

Vector operator/(const Vector& a, const double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}

double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

class Polygon {
public:
    Polygon() = default;
    explicit Polygon(const std::vector<Vector>& vertices) : vertices(vertices) {}

    void addVertex(const Vector& vertex) {
        vertices.push_back(vertex);
    }

    std::vector<Vector> vertices;
};
bool inside(const Vector& p, const Vector& start_edge, const Vector& end_edge) {
    return (end_edge[0] - start_edge[0]) * (p[1] - start_edge[1]) > (end_edge[1] - start_edge[1]) * (p[0] - start_edge[0]);
}

Vector intersect(const Vector& p1, const Vector& p2, const Vector& start_edge, const Vector& end_edge) {
    Vector edge_direction = end_edge - start_edge;
    Vector line_dir = p2 - p1;
    double t = ((start_edge[0] - p1[0]) * (start_edge[1] - end_edge[1]) - (start_edge[1] - p1[1]) * (start_edge[0] - end_edge[0])) /
               (line_dir[0] * (start_edge[1] - end_edge[1]) - line_dir[1] * (start_edge[0] - end_edge[0]));
    return p1 + t * line_dir;
}

Polygon sutherlandHodgmanClip(const Polygon& subjectPolygon, const Polygon& clipPolygon) {
    Polygon polygonfinal = subjectPolygon;

    for (size_t i = 0; i < clipPolygon.vertices.size(); ++i) {
        Vector start_edge = clipPolygon.vertices[i];
        Vector end_edge = clipPolygon.vertices[(i + 1) % clipPolygon.vertices.size()];

        Polygon initial_poly = polygonfinal;
        polygonfinal.vertices.clear();

        Vector last_vertex = initial_poly.vertices.back();
        for (const Vector& vertex_i : initial_poly.vertices) {
            if (inside(vertex_i, start_edge, end_edge)) {
                if (!inside(last_vertex, start_edge, end_edge)) {
                    Vector intersection = intersect(last_vertex, vertex_i, start_edge, end_edge);
                    polygonfinal.addVertex(intersection);
                }
                polygonfinal.addVertex(vertex_i);
            } else if (inside(last_vertex, start_edge, end_edge)) {
                Vector intersection = intersect(last_vertex, vertex_i, start_edge, end_edge);
                polygonfinal.addVertex(intersection);
            }
            last_vertex = vertex_i;
        }
    }

    return polygonfinal;
}

void save_svg(const std::vector<Polygon>& polygons, const std::string& filename, const std::string& fillcol = "none") {
    FILE* f = fopen(filename.c_str(), "w+");
    fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
    for (const auto& polygon : polygons) {
        fprintf(f, "<g>\n");
        fprintf(f, "<polygon points = \"");
        for (const auto& vertex : polygon.vertices) {
            fprintf(f, "%3.3f, %3.3f ", vertex[0] * 1000, 1000 - vertex[1] * 1000);
        }
        fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
        fprintf(f, "</g>\n");
    }
    fprintf(f, "</svg>\n");
    fclose(f);
}

int main() {
    Polygon polygon1;
    polygon1.addVertex(Vector(0.1, 0.2, 0));
    polygon1.addVertex(Vector(0.4, 0.1, 0));
    polygon1.addVertex(Vector(0.5, 0.5, 0));
    polygon1.addVertex(Vector(0.3, 0.7, 0));
    polygon1.addVertex(Vector(0.1, 0.5, 0));

    Polygon polygon2;
    polygon2.addVertex(Vector(0.2, 0.3, 0));
    polygon2.addVertex(Vector(0.6, 0.2, 0));
    polygon2.addVertex(Vector(0.7, 0.6, 0));
    polygon2.addVertex(Vector(0.5, 0.8, 0));
    polygon2.addVertex(Vector(0.3, 0.6, 0));

    // Save original intertwined polygons
    std::vector<Polygon> beforeClipping = {polygon1, polygon2};
    save_svg(beforeClipping, "intertwined_polygons.svg");

    // Clip polygons
    Polygon clippedPolygon = sutherlandHodgmanClip(polygon1, polygon2);

    // Save the result after clipping
    std::vector<Polygon> afterClipping = {clipped_polygon};
    save_svg(afterClipping, "clipped_polygon.svg");

    return 0;
}