#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <set>
#include <unordered_map>
#include <chrono>
#include <stdio.h>
#include "stb_image_write.h"
#include "stb_image.h"

class Vector {
public:
    double data[3];
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
    double operator[](int i) const { return data[i]; };
    double& operator[](int i) { return data[i]; };
    Vector operator+(const Vector& b) const {
        return Vector(data[0] + b[0], data[1] + b[1], data[2] + b[2]);
    }
    Vector operator-(const Vector& b) const {
        return Vector(data[0] - b[0], data[1] - b[1], data[2] - b[2]);
    }
    Vector operator*(double b) const {
        return Vector(data[0] * b, data[1] * b, data[2] * b);
    }
    Vector operator/(double b) const {
        return Vector(data[0] / b, data[1] / b, data[2] / b);
    }
};

double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

class TriangleIndices {
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
    };
    int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
    int uvi, uvj, uvk;  // indices within the uv coordinates array
    int ni, nj, nk;  // indices within the normals array
    int group;       // face group
};
 
 
class TriangleMesh {
public:
  ~TriangleMesh() {}
    TriangleMesh() {};
    
    void readOBJ(const char* obj) {
 
        char matfile[255];
        char grp[255];
 
        FILE* f;
        f = fopen(obj, "r");
        int curGroup = -1;
        while (!feof(f)) {
            char line[255];
            if (!fgets(line, 255, f)) break;
 
            std::string linetrim(line);
            linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
            strcpy(line, linetrim.c_str());
 
            if (line[0] == 'u' && line[1] == 's') {
                sscanf(line, "usemtl %[^\n]\n", grp);
                curGroup++;
            }
 
            if (line[0] == 'v' && line[1] == ' ') {
                Vector vec;
 
                Vector col;
                if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
                    col[0] = std::min(1., std::max(0., col[0]));
                    col[1] = std::min(1., std::max(0., col[1]));
                    col[2] = std::min(1., std::max(0., col[2]));
 
                    vertices.push_back(vec);
                    vertexcolors.push_back(col);
 
                } else {
                    sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                    vertices.push_back(vec);
                }
            }
            if (line[0] == 'v' && line[1] == 'n') {
                Vector vec;
                sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                normals.push_back(vec);
            }
            if (line[0] == 'v' && line[1] == 't') {
                Vector vec;
                sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
                uvs.push_back(vec);
            }
            if (line[0] == 'f') {
                TriangleIndices t;
                int i0, i1, i2, i3;
                int j0, j1, j2, j3;
                int k0, k1, k2, k3;
                int nn;
                t.group = curGroup;
 
                char* consumedline = line + 1;
                int offset;
 
                nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
                if (nn == 9) {
                    if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                    if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                    if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                    if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                    if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                    if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                    if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                    if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                    if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                    indices.push_back(t);
                } else {
                    nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                    if (nn == 6) {
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                        if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                        if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                        indices.push_back(t);
                    } else {
                        nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                        if (nn == 3) {
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            indices.push_back(t);
                        } else {
                            nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                            if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                            if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                            indices.push_back(t);
                        }
                    }
                }
 
                consumedline = consumedline + offset;
 
                while (true) {
                    if (consumedline[0] == '\n') break;
                    if (consumedline[0] == '\0') break;
                    nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                    TriangleIndices t2;
                    t2.group = curGroup;
                    if (nn == 3) {
                        if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                        if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                        if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                        if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                        if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                        if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                        if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                        if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                        if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                        indices.push_back(t2);
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        k2 = k3;
                    } else {
                        nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                        if (nn == 2) {
                            if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                            if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                            if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                            if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                            if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                            if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            j2 = j3;
                            indices.push_back(t2);
                        } else {
                            nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                            if (nn == 2) {
                                if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                                if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                                if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;                             
                                consumedline = consumedline + offset;
                                i2 = i3;
                                k2 = k3;
                                indices.push_back(t2);
                            } else {
                                nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                if (nn == 1) {
                                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    indices.push_back(t2);
                                } else {
                                    consumedline = consumedline + 1;
                                }
                            }
                        }
                    }
                }
 
            }
 
        }
        fclose(f);
 
    }
 
    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
    
};

class Mesh {
public:
    std::vector<Vector> vertices;
    std::vector<std::vector<int>> adjacencyList;
    std::vector<int> boundaryIndices;
    Mesh(int n) : vertices(n), adjacencyList(n) {}
    Mesh() {}
};

void layoutBoundaryVertices(Mesh& mesh) {
    int n = mesh.boundaryIndices.size();
    double s = 0.0;
    for (int i = 0; i < n; ++i) {
        int next = (i + 1) % n;
        s += hypot(mesh.vertices[mesh.boundaryIndices[next]][0] - mesh.vertices[mesh.boundaryIndices[i]][0],
                   mesh.vertices[mesh.boundaryIndices[next]][1] - mesh.vertices[mesh.boundaryIndices[i]][1]);
    }

    double cs = 0.0;
    for (int i = 0; i < n; ++i) {
        double theta = 2 * M_PI * cs / s;
        mesh.vertices[mesh.boundaryIndices[i]] = Vector(cos(theta), sin(theta), 0);
        int next = (i + 1) % n;
        cs += hypot(mesh.vertices[mesh.boundaryIndices[next]][0] - mesh.vertices[mesh.boundaryIndices[i]][0],
                    mesh.vertices[mesh.boundaryIndices[next]][1] - mesh.vertices[mesh.boundaryIndices[i]][1]);
    }
}

void tutteEmbedding(Mesh& mesh, int iterations) {
    layoutBoundaryVertices(mesh);
    int m = mesh.vertices.size();
    std::vector<Vector> newPositions(m);

    for (int iter = 0; iter < iterations; ++iter) {
        for (int i = 0; i < m; ++i) {
            if (std::find(mesh.boundaryIndices.begin(), mesh.boundaryIndices.end(), i) == mesh.boundaryIndices.end()) {
                Vector sum(0, 0, 0);
                for (int neighbor : mesh.adjacencyList[i]) {
                    sum = sum + mesh.vertices[neighbor];
                }
                if (!mesh.adjacencyList[i].empty()) {
                    newPositions[i] = sum / mesh.adjacencyList[i].size();
                }
            } else {
                newPositions[i] = mesh.vertices[i];
            }
        }
        mesh.vertices = newPositions;
    }
}

void saveSVG(const Mesh& mesh, const std::string& filename) {
    std::ofstream file(filename);
    file << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1000\" height=\"1000\" viewBox=\"-1 -1 2 2\">\n";
    file << "<g transform=\"scale(1,-1)\">\n";
    for (const auto& pos : mesh.vertices) {
        file << "<circle cx=\"" << pos[0] << "\" cy=\"" << pos[1] << "\" r=\"0.01\" fill=\"black\" />\n";
    }
    file << "</g>\n</svg>";
    file.close();
}

Mesh convertToMesh(TriangleMesh& triMesh) {
    Mesh mesh(triMesh.vertices.size());
    mesh.vertices = triMesh.vertices;

    mesh.adjacencyList.resize(triMesh.vertices.size());
    for (const auto& face : triMesh.indices) {
        mesh.adjacencyList[face.vtxi].push_back(face.vtxj);
        mesh.adjacencyList[face.vtxi].push_back(face.vtxk);
        mesh.adjacencyList[face.vtxj].push_back(face.vtxi);
        mesh.adjacencyList[face.vtxj].push_back(face.vtxk);
        mesh.adjacencyList[face.vtxk].push_back(face.vtxi);
        mesh.adjacencyList[face.vtxk].push_back(face.vtxj);
    }

    std::unordered_map<int, int> edgeCount;
    for (const auto& face : triMesh.indices) {
        int v1 = face.vtxi;
        int v2 = face.vtxj;
        int v3 = face.vtxk;
        int edges[3][2] = {{v1, v2}, {v2, v3}, {v3, v1}};
        for (auto& edge : edges) {
            int a = std::min(edge[0], edge[1]);
            int b = std::max(edge[0], edge[1]);
            int key = a * triMesh.vertices.size() + b;
            edgeCount[key]++;
        }
    }

    for (const auto& edge : edgeCount) {
        if (edge.second == 1) {
            int a = edge.first / triMesh.vertices.size();
            int b = edge.first % triMesh.vertices.size();
            mesh.boundaryIndices.push_back(a);
            mesh.boundaryIndices.push_back(b);
        }
    }

    std::sort(mesh.boundaryIndices.begin(), mesh.boundaryIndices.end());
    mesh.boundaryIndices.erase(std::unique(mesh.boundaryIndices.begin(), mesh.boundaryIndices.end()), mesh.boundaryIndices.end());

    return mesh;
}

int main() {
    TriangleMesh triMesh;
    triMesh.readOBJ("goethe.obj");

    std::cout << "obj read done" << std::endl;


    Mesh mesh = convertToMesh(triMesh);

    std::cout << "mesh done" << std::endl;

    tutteEmbedding(mesh, 100);

    std::cout << "tutte done" << std::endl;

    saveSVG(mesh, "tutte_embedding.svg");

    return 0;
}
