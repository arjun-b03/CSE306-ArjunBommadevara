#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
 
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <iostream>
#include <chrono>
#include <random>

#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>




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
    double operator[](int i) const { return data[i]; };
    double& operator[](int i) { return data[i]; };
    double data[3];
};
 
Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector& b) {
    return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
    return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator*(const Vector &a, const Vector &b){
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


const double infinity = std::numeric_limits<double>::infinity();

class Ray {
public:
    Vector origin;
    Vector direction;

    Ray(const Vector& origin, const Vector& direction) : origin(origin), direction(direction) {}

    Vector point_at_parameter(double t) const {
        return origin + direction * t;
    }
};

class Geometry
{
public:
	Geometry(){};
	virtual bool intersect(const Ray& ray, Vector &hit, Vector &N, double &t, const Vector& vertex0, const Vector& vertex1, const Vector& vertex2, double &alpha, double &beta, double &gamma, int &meshID) = 0;
	Vector albedo;
	bool mirror,transparent;
};

class Sphere: public Geometry {
public:
    Vector center;
    double radius;

    Sphere(const Vector &center, double radius, const Vector &albedo, bool mirror = false, bool transparent = false):Geometry(){
        this->center=center;
        this->radius=radius;
        this->albedo=albedo;
        this->mirror=mirror;
        this->transparent=transparent;
    }

    bool intersect(const Ray& ray, Vector &hit, Vector &N, double &t, const Vector& vertex0, const Vector& vertex1, const Vector& vertex2, double &alpha, double &beta, double &gamma, int &meshID) {
        Vector oc = ray.origin - center;  
        double a = dot(ray.direction, ray.direction);  // Should always be 1 if direction is normalized
        double b = 2.0 * dot(oc, ray.direction);
        double c = dot(oc, oc) - radius * radius;
        double discriminant = b * b - 4 * a * c;

        if (discriminant < 0) {
            return false;  // No intersection: the ray misses the sphere
        }

        double sqrtdelta = sqrt(discriminant);
        double t1 = (-b - sqrtdelta) / (2 * a);
        double t2 = (-b + sqrtdelta) / (2 * a);

        if (t2 < 0) {
            return false;  // Both t1 and t2 are negative, ray starts after the sphere
        }


        if (t1 > 0) {
            t = t1;  // t1 is the first intersection point
        } else {
            t = t2;  // t1 is negative, use t2 (the ray starts inside the sphere)
        }

        hit = ray.origin + t * ray.direction; 
        N = (hit - center);
        N.normalize(); 

        return true;
    }
};

static std::default_random_engine engine(10) ; // random s e e d = 10
static std::uniform_real_distribution<double> uniform ( 0 , 1 );

void boxMuller ( double stdev , double &x , double &y ) {
	double r1 = uniform ( engine ) ;
	double r2 = uniform ( engine ) ;
	x = sqrt(-2 * log ( r1 ) ) *cos ( 2 * M_PI*r2 ) *stdev ;
	y = sqrt(-2 * log ( r1 ) ) *sin ( 2 * M_PI*r2 ) *stdev ;
}

Vector random_cos(const Vector &N){
    double r1,r2;
    r1 = uniform(engine);
    r2 = uniform(engine);

    double x,y,z;
    x = cos(2*M_PI* r1) * sqrt(1 - r2);
    y = sin(2*M_PI* r1) * sqrt(1 - r2);
    z = sqrt(r2);

    Vector T1,T2;
    double min_num = N[0];
    if (N[1]<min_num){
        min_num = N[1];
    }
    if (N[2]<min_num){
        min_num = N[2];
    }

    if (min_num == N[0]){
        T1 = Vector(0, N[2], -N[1]);
    }
    if (min_num == N[1]){
        T1 = Vector(N[1], 0, -N[0]);
    }
    if (min_num == N[2]){
        T1 = Vector(N[1], -N[0], 0);
    }
    T2 = cross(N, T1);

    T1.normalize();
    T2.normalize();
    return x*T1 + y*T2 + z*N;
}


bool MollerIntersection(const Ray &ray, const Vector &vertex0, const Vector &vertex1, const Vector &vertex2, double &t, double &u, double &v) {
    const double EPSILON = 1e-8;
    Vector edge1 = vertex1 - vertex0;
    Vector edge2 = vertex2 - vertex0;
    Vector h = cross(ray.direction, edge2);
    double a = dot(edge1, h);

    if (std::fabs(a) < EPSILON) {
        return false;  // This ray is parallel to this triangle.
    }

    double f = 1.0 / a;
    Vector s = ray.origin - vertex0;
    u = f * dot(s, h);
    if (u < 0.0 || u > 1.0) {
        return false;
    }

    Vector q = cross(s, edge1);
    v = f * dot(ray.direction, q);
    if (v < 0.0 || u + v > 1.0) {
        return false;
    }

    // At this stage, we can compute t to find out where the intersection point is on the line.
    t = f * dot(edge2, q);
    if (t > EPSILON) { // Ray intersection
        return true;
    } else {
        return false; // This means that there is a line intersection but not a ray intersection.
    }
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
 
 
class TriangleMesh: public Geometry {
public:
  ~TriangleMesh() {}
    TriangleMesh(const Vector albedo) {
        this->albedo = albedo;
        this->mirror = false;
        this->transparent = false;
    };

    bool intersect(const Ray& ray, Vector &hit, Vector &N, double &t, const Vector& vertex0, const Vector& vertex1, 
                    const Vector& vertex2, double &alpha, double &beta, double &gamma, int &meshID) {
        bool hasIntersection = false;
        double tMin = std::numeric_limits<double>::infinity();
        int id = -1;
        Vector n0;
        double t0, alpha0, beta0, gamma0;
        for (size_t i = 0; i < indices.size(); ++i) {
            if (MollerIntersection(ray, vertices[indices[i].vtxi], vertices[indices[i].vtxj], vertices[indices[i].vtxk], t0, beta0, gamma0)) {
                alpha0 = 1.0 - beta0 - gamma0;
                if (t0 > 0 && t0 < t) { // Only consider positive t values that are nearer than previously found intersections
                    t = t0;
                    id = i;
                    hasIntersection = true;
                    beta = beta0;
                    gamma = gamma0;
                    alpha = alpha0;
                    hit = ray.origin + ray.direction * t; // Calculate the exact intersection point
                    Vector e1 = vertices[indices[i].vtxj] - vertices[indices[i].vtxi];
                    Vector e2 = vertices[indices[i].vtxk] - vertices[indices[i].vtxi];
                    n0 = cross(e1, e2); // Compute normal at the intersection
                    n0.normalize();
                    N = n0;
                }
            }
        }
        meshID = id;
        return hasIntersection;
    }
    
    void transform(double s, const Vector& t) {
            for (int i = 0; i < vertices.size(); i++) {
                vertices[i] = vertices[i] * s + t;
            }
        }
    
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

class Scene {
public:
    Scene() {};
    void addObject(Geometry &obj){
		objects.push_back(&obj);
	}
	
	std::vector<Geometry*> objects;

    bool intersect(const Ray& ray, Vector &hit, Vector &N, double &t, int &SphereID, int &meshID) {
        bool hasIntersection = false;
        double tMin = std::numeric_limits<double>::infinity();
        for (int i =0; i<objects.size(); i++) {
            Vector hit0, n0;
            double t0;
            Vector vertex0, vertex1, vertex2;
            double alpha,beta,gamma; 
            meshID = -1;
			int meshID0;
            bool intersection = objects[i]->intersect(ray, hit0, n0, t0, vertex0, vertex1, vertex2, alpha, beta, gamma, meshID0);
            if (intersection && t0 > 0) {
                hasIntersection = true;
                if (t0 < tMin) {
                    tMin = t0;
                    hit = hit0;
                    N = n0;
                    SphereID = i;
                    meshID = meshID0;
                    
                }
            }
        }
        t = tMin;
        return hasIntersection;
    }

    Vector getColor(const Ray& ray, int ray_depth){
		Vector albedo0(0,0,0);
        Vector albedo(0,0,0);
		Vector hit,N;
		Vector Light(10,20,40);

        int SphereID, meshID;
        

		double t;
		double I = 1E10;

        if (ray_depth <= 0){
			return albedo;
		}
        
        if (intersect(ray, hit, N, t, SphereID,meshID)){

            if (objects[SphereID]->mirror){
				Vector R = ray.direction - 2*dot(ray.direction,N) *N;
				Ray reflect_ray(hit+0.001*N, R);
				return getColor(reflect_ray,ray_depth-1);
			}

            if (objects[SphereID]->transparent) {
                double n1 = 1.0; // Assume air by default
                double n2 = 1.4;
                if (dot(ray.direction, N) > 0) { // Ray inside the sphere
                    std::swap(n1, n2);
                    N = Vector(-N[0], -N[1], -N[2]); // Reverse normal if inside the sphere
                }

                double eta = n1 / n2;
                double cosi = -dot(ray.direction, N);
                double k = 1 - eta * eta * (1 - cosi * cosi);

                if (k < 0) {
                    // Total internal reflection
                    Vector R = ray.direction - 2 * dot(ray.direction, N) * N;
                    Ray reflectRay(hit + 0.001 * N, R);
                    return getColor(reflectRay, ray_depth - 1);
                } else {
                    // Refracted ray
                    Vector T = eta * ray.direction + (eta * cosi - sqrt(k)) * N;
                    Ray refractRay(hit - 0.001 * N, T); 
                    return getColor(refractRay, ray_depth - 1);
                }
            }
            Vector lightVec = Light - hit;
            double distToLight2 = lightVec.norm2();
            lightVec.normalize();

            Ray shadowRay(hit + N * 0.001, lightVec);
            Vector Plight, Nlight;
            double tlight;
            int IDlight, IDMeshLight;
            double Shadow = 1.0;

            if (intersect(shadowRay, Plight, Nlight, tlight, IDlight, IDMeshLight) && tlight * tlight < distToLight2) {
                Shadow = 0.0;
            }

            double intensity = Shadow * I / ((Light-hit).norm2() * 4 * M_PI);
            double lambertian = std::max(dot(lightVec, N), 0.0);
            albedo0 = intensity * lambertian * objects[SphereID]->albedo/M_PI;
            Ray random_ray(hit, random_cos(N)); 
            albedo = albedo0 + objects[SphereID]->albedo * getColor(random_ray, ray_depth - 1);
        }
        return albedo;
        
    }
};



int main() {
    auto start = std::chrono::high_resolution_clock::now();
    int W = 512;
    int H = 512;
    const Vector camera(0, 0, 55);  // Camera is at z = -50
    double fov = 60 * M_PI / 180;  // Field of View
    double aspect_ratio = W / static_cast<double>(H);
    double scale = tan(fov * 0.5);  
    double gamma = 2.2;

    // Sphere Sphere1(Vector(0,0,0), 10, Vector(1.0,1.0,1.0));
    // Sphere SphereMirror(Vector(-20,0,0), 10, Vector(1.0,1.0,1.0), true);
    // Sphere SphereTransparent(Vector(20,0,0), 10, Vector(1.0,1.0,1.0), false, true);

	Sphere SphereFloor(Vector(0,-1000,0), 990, Vector(0.0,0.8,0.8));
	Sphere SphereCeilling(Vector(0,1000,0), 940, Vector(0.0,1.0,0.0));
	Sphere SphereLeft(Vector(-1000,0,0), 940, Vector(0.0,0.0,1.0));
	Sphere SphereRight(Vector(1000,0,0), 940, Vector(1.0,0.0,1.0));
	Sphere SphereFront(Vector(0,0,-1000), 940, Vector(1.0,1.0,0.0));
    Sphere SphereBack(Vector(0,0,1000), 940, Vector(1.0,0.0,0.0));

	Scene scene;
	// scene.addObject(Sphere1);
    // scene.addObject(SphereMirror);
    // scene.addObject(SphereTransparent);   
	scene.addObject(SphereFloor);
	scene.addObject(SphereCeilling);
	scene.addObject(SphereLeft);
	scene.addObject(SphereRight);
	scene.addObject(SphereFront);
    scene.addObject(SphereBack);


    TriangleMesh mesh(Vector(1.0,0.4,1.0));
    mesh.readOBJ("cat.obj");
	mesh.transform(0.6, Vector(0, -10, 0));	
    scene.addObject(mesh);
    
    int nb_rays = 32;
    std::cout << nb_rays << std::endl;
    std::vector<unsigned char> image(W * H * 3, 0);
    #pragma omp parallel for schedule(dynamic, 1)
    for (int j = 0; j < H; j++) {
        std::cout << j << std::endl;
        for (int i = 0; i < W; i++) {
            Vector albedo(0.,0.,0.);
            for (int k = 0; k < nb_rays; k++){
                double x, y;
                boxMuller(1, x, y);
                x += i;
                y += j;
                Vector dir(y - W/2 + 0.5, H/2 - x - 0.5, -W / (2 * tan( fov / 2)));
                dir.normalize();
			    Ray ray(camera,dir);
			    albedo = albedo + scene.getColor(ray,5);
            }
			image[(i * W + j) * 3 + 0] = std::min(255., std::pow(albedo[0]/nb_rays, (1.0 / gamma)));
			image[(i * W + j) * 3 + 1] = std::min(255., std::pow(albedo[1]/nb_rays, (1.0 / gamma)));
			image[(i * W + j) * 3 + 2] = std::min(255., std::pow(albedo[2]/nb_rays, (1.0 / gamma)));
		}
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    std::cout << "Time taken: " << duration.count() << "ms" << std::endl;
    stbi_write_png("image.png", W, H, 3, &image[0], 0);

    return 0;
}