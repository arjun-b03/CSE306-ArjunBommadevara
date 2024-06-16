#include <iostream>
#include <chrono>
#include <random>
#include "stb_image.h"
#include "stb_image_write.h"
#include <cmath>

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
    
    void addVertex(const Vector& vertex) {
        vertices.push_back(vertex);
    }
    double area() const {
        if (vertices.size() < 3) return 0.0;

        double a = 0;
		int n = vertices.size();
		for (int i = 0; i < n; i++){
			a += (vertices[i][0]*vertices[(i+1)%n][1] - vertices[(i+1)%n][0]*vertices[i][1]);
		}
		return 0.5*std::abs(a);
	}

    std::vector<Vector> vertices;
};

void save_svg(const std::vector<Polygon> &polygons, std::string filename, std::string fillcol = "none") {
        FILE* f = fopen(filename.c_str(), "w+"); 
        fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
        for (int i=0; i<polygons.size(); i++) {
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
 
 
// Adds one frame of an animated svg file. frameid is the frame number (between 0 and nbframes-1).
// polygons is a list of polygons, describing the current frame.
// The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
    void save_svg_animated(const std::vector<Polygon> &polygons, std::string filename, int frameid, int nbframes) {
        FILE* f;
        if (frameid == 0) {
            f = fopen(filename.c_str(), "w+");
            fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
            fprintf(f, "<g>\n");
        } else {
            f = fopen(filename.c_str(), "a+");
        }
        fprintf(f, "<g>\n");
        for (int i = 0; i < polygons.size(); i++) {
            fprintf(f, "<polygon points = \""); 
            for (int j = 0; j < polygons[i].vertices.size(); j++) {
                fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000-polygons[i].vertices[j][1] * 1000));
            }
            fprintf(f, "\"\nfill = \"none\" stroke = \"black\"/>\n");
        }
        fprintf(f, "<animate\n");
        fprintf(f, "    id = \"frame%u\"\n", frameid);
        fprintf(f, "    attributeName = \"display\"\n");
        fprintf(f, "    values = \"");
        for (int j = 0; j < nbframes; j++) {
            if (frameid == j) {
                fprintf(f, "inline");
            } else {
                fprintf(f, "none");
            }
            fprintf(f, ";");
        }
        fprintf(f, "none\"\n    keyTimes = \"");
        for (int j = 0; j < nbframes; j++) {
            fprintf(f, "%2.3f", j / (double)(nbframes));
            fprintf(f, ";");
        }
        fprintf(f, "1\"\n   dur = \"5s\"\n");
        fprintf(f, "    begin = \"0s\"\n");
        fprintf(f, "    repeatCount = \"indefinite\"/>\n");
        fprintf(f, "</g>\n");
        if (frameid == nbframes - 1) {
            fprintf(f, "</g>\n");
            fprintf(f, "</svg>\n");
        }
        fclose(f);
    }



Polygon clipVoronoi(const Polygon& subjectPolygon, const Vector& siteA, const Vector& siteB) {
    Polygon polygonfinal;
    int vertices_n = subjectPolygon.vertices.size();

    for (int i = 0; i < vertices_n; ++i) {
        const Vector& vertex_i = subjectPolygon.vertices[i];
        const Vector& vertex_prev = subjectPolygon.vertices[(i - 1 + vertices_n) % vertices_n];

        double dist_a = (vertex_i - siteA).norm2();
        double dist_b = (vertex_i - siteB).norm2();
        double dist_a_prev = (vertex_prev - siteA).norm2();
        double dist_b_prev = (vertex_prev - siteB).norm2();

        if (dist_a <= dist_b) {
            if (dist_a_prev > dist_b_prev) {
                Vector midPoint = (siteA + siteB) * 0.5;
                double param = dot(midPoint - vertex_prev, siteB - siteA) / dot(vertex_i - vertex_prev, siteB - siteA);
                Vector intersection = vertex_prev + param * (vertex_i - vertex_prev);
                polygonfinal.vertices.push_back(intersection);
            }
            polygonfinal.vertices.push_back(vertex_i);
        } else {
            if (dist_a_prev <= dist_b_prev) {
                Vector midpoint = (siteA + siteB) * 0.5;
                double param = dot(midpoint - vertex_prev, siteB - siteA) / dot(vertex_i - vertex_prev, siteB - siteA);
                Vector intersection = vertex_prev + param * (vertex_i - vertex_prev);
                polygonfinal.vertices.push_back(intersection);
            }
        }
    }

    return polygonfinal;
}


class VoronoiDiagram {
public:
    VoronoiDiagram() = default;
    explicit VoronoiDiagram(const std::vector<Vector>& points) : points(points) {}

    void compute() {
        Polygon square1;
        square1.addVertex(Vector(0, 0, 0));
        square1.addVertex(Vector(1, 0, 0));
        square1.addVertex(Vector(1, 1, 0));
        square1.addVertex(Vector(0, 1, 0));

        cells.resize(points.size());
        #pragma omp parallel for
        for (int i = 0; i < points.size(); i++) {
            Polygon cell = square1;
            for (int j = 0; j < points.size(); j++) {
                if (i != j) {
                    cell = clipVoronoi(cell, points[i], points[j]);
                }
            }
            cells[i] = cell;
        }
    }

    std::vector<Vector> points;
    std::vector<Polygon> cells;
};

int main() {
    auto start = std::chrono::high_resolution_clock::now();

    int N = 1000;
    std::vector<Vector> points(N);
    for (int i = 0; i < N; i++) {
        points[i] = Vector(rand() / double(RAND_MAX), rand() / double(RAND_MAX), 0);
    }

    VoronoiDiagram voronoi(points);
    voronoi.compute();

    save_svg(voronoi.cells, "/Users/arjunbommadevara/Desktop/CSE306/project1/CSE306-Project-1/PROJECT2/vor.svg");

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    std::cout << "Elapsed time: " << duration.count() << " milliseconds" << std::endl;

    return 0;
}