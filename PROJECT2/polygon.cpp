#include <vector>
#include <string>
#include <cmath>
#include <cstdio>
#include "vector.cpp"
#include <cstdio>

class Polygon {
public:
    Polygon() = default;
    explicit Polygon(const std::vector<Vector>& vertices) : vertices(vertices) {}

    double computeArea() const {
        if (vertices.size() < 3) return 0.0;

        double a = 0;
		int n = vertices.size();
		for (int i = 0; i < n; i++){
			a += (vertices[i][0]*vertices[(i+1)%n][1] - vertices[(i+1)%n][0]*vertices[i][1]);
		}
		return 0.5*std::abs(a);
	}




    double SquaredDistanceIntegral(const Vector& point) const {
        if (vertices.size() < 3) return 0.0;

        double squaredDistanceIntegral = 0.0;
        int vertexCount = vertices.size();

        for (int i = 1; i < vertexCount - 1; ++i) {
            Vector triangle[3] = {vertices[0], vertices[i], vertices[i + 1]};

            double localIntegral = 0.0;
            for (int k = 0; k < 3; ++k) {
                for (int l = k; l < 3; ++l) {
                    localIntegral += dot(triangle[k] - point, triangle[l] - point);
                }
            }

            Vector edge1 = triangle[1] - triangle[0];
            Vector edge2 = triangle[2] - triangle[0];
            double triangleArea = 0.5 * std::abs(edge1[1] * edge2[0] - edge1[0] * edge2[1]);

            squaredDistanceIntegral += (localIntegral / 6.0) * triangleArea;
        }

        return squaredDistanceIntegral;
    }

    void addVertex(const Vector& vertex) {
        vertices.push_back(vertex);
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
