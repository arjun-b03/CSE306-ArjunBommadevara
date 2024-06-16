#include <stdio.h>
#include <vector>
#include <string>
#include <cmath>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <random>
#include "stb_image_write.h"
#include "lbfgs.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION

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

    double squared_distance_integral(const Vector& point) const {
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
            double triangle_area = 0.5 * std::abs(edge1[1] * edge2[0] - edge1[0] * edge2[1]);

            result += (integral / 6.0) * triangle_area;
        }

        return result;
    }

    Vector centroid() const {
        if (vertices.size() < 3) {
            std::cerr << "Error: Attempting to calculate the centroid of a polygon with fewer than 3 vertices." << std::endl;
            return Vector(0, 0, 0);
        }

        Vector c(0, 0, 0);
        double signed_area = 0.0;

        int n = vertices.size();
        for (int i = 0; i < n; ++i) {
            const Vector& vi = vertices[i];
            const Vector& vj = vertices[(i + 1) % n];
            double a = vi[0] * vj[1] - vj[0] * vi[1];
            signed_area += a;
            c = c + (vi + vj) * a;
        }

        signed_area *= 0.5;
        c = c / (6.0 * signed_area);
        return c;
    }

    void add_vertex(const Vector& vertex) {
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

int sgn(double val) {
    return (0 < val) - (val < 0);
}

void save_frame(const std::vector<Polygon>& cells, const std::vector<Vector>& points, std::string filename, int frameid = 0) {
    int W = 1000, H = 1000;
    std::vector<unsigned char> image(W * H * 3, 255);

    //Draw Voronoi cells with white edges
    for (int i = 0; i < cells.size(); i++) {
        double bminx = 1E9, bminy = 1E9, bmaxx = -1E9, bmaxy = -1E9;
        for (int j = 0; j < cells[i].vertices.size(); j++) {
            bminx = std::min(bminx, cells[i].vertices[j][0]);
            bminy = std::min(bminy, cells[i].vertices[j][1]);
            bmaxx = std::max(bmaxx, cells[i].vertices[j][0]);
            bmaxy = std::max(bmaxy, cells[i].vertices[j][1]);
        }
        bminx = std::min(W - 1., std::max(0., W * bminx));
        bminy = std::min(H - 1., std::max(0., H * bminy));
        bmaxx = std::min(W - 1., std::max(0., W * bmaxx));
        bmaxy = std::min(H - 1., std::max(0., H * bmaxy));

        for (int y = bminy; y < bmaxy; y++) {
            for (int x = bminx; x < bmaxx; x++) {
                int last_sign = 0;
                bool inside = true;
                double smallest_edge = 1E9;
                for (int j = 0; j < cells[i].vertices.size(); j++) {
                    double x0 = cells[i].vertices[j][0] * W;
                    double y0 = cells[i].vertices[j][1] * H;
                    double x1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][0] * W;
                    double y1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][1] * H;
                    double det = (x - x0) * (y1 - y0) - (y - y0) * (x1 - x0);
                    int sign = (det > 0) - (det < 0);
                    if (last_sign == 0) last_sign = sign;
                    else if (sign == 0) sign = last_sign;
                    else if (sign != last_sign) {
                        inside = false;
                        break;
                    }
                    last_sign = sign;
                    double edgeLen = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
                    double dist_edge = std::abs(det) / edgeLen;
                    double dotp = (x - x0) * (x1 - x0) + (y - y0) * (y1 - y1);
                    if (dotp < 0 || dotp > edgeLen * edgeLen) dist_edge = 1E9;
                    smallest_edge = std::min(smallest_edge, dist_edge);
                }
                if (inside) {
                    if (smallest_edge <= 2) {
                        image[((H - y - 1) * W + x) * 3] = 255;
                        image[((H - y - 1) * W + x) * 3 + 1] = 255;
                        image[((H - y - 1) * W + x) * 3 + 2] = 255;
                    }
                }
            }
        }
    }

    int dot_radius = 30;
    for (const auto& p : points) {
        int cx = static_cast<int>(p[0] * W);
        int cy = H - static_cast<int>(p[1] * H);
        for (int dy = -dot_radius; dy <= dot_radius; ++dy) {
            for (int dx = -dot_radius; dx <= dot_radius; ++dx) { 
                int x = cx + dx;
                int y = cy + dy;
                if (x >= 0 && x < W && y >= 0 && y < H && dx * dx + dy * dy <= dot_radius * dot_radius) {
                    image[(y * W + x) * 3] = 0;
                    image[(y * W + x) * 3 + 1] = 0;
                    image[(y * W + x) * 3 + 2] = 255;
                }
            }
        }
    }

    std::ostringstream os;
    os << filename << frameid << ".png";
    stbi_write_png(os.str().c_str(), W, H, 3, &image[0], 0);
}

const double GRAVITY = 9.81;
const double SPRING_CONSTANT = 0.1;
const double DELTA_TIME = 0.02;
const int NUM_FRAMES = 50;
const int NUM_PARTICLES = 100;

Polygon clip_power_diagram(const Polygon& subject_polygon, const Vector& site_a, double weight_a, const Vector& site_b, double weight_b) {
    Polygon result_polygon;
    int num_vertices = subject_polygon.vertices.size();

    for (int i = 0; i < num_vertices; ++i) {
        const Vector& current_vertex = subject_polygon.vertices[i];
        const Vector& previous_vertex = subject_polygon.vertices[(i - 1 + num_vertices) % num_vertices];

        double power_distance_current_a = (current_vertex - site_a).norm2() - weight_a;
        double power_distance_current_b = (current_vertex - site_b).norm2() - weight_b;
        double power_distance_previous_a = (previous_vertex - site_a).norm2() - weight_a;
        double power_distance_previous_b = (previous_vertex - site_b).norm2() - weight_b;

        if (power_distance_current_a <= power_distance_current_b) {
            if (power_distance_previous_a > power_distance_previous_b) {
                Vector mid_point = (site_a + site_b) * 0.5;
                double param = dot(mid_point - previous_vertex, site_b - site_a) / dot(current_vertex - previous_vertex, site_b - site_a);
                Vector intersection = previous_vertex + param * (current_vertex - previous_vertex);
                result_polygon.vertices.push_back(intersection);
            }
            result_polygon.vertices.push_back(current_vertex);
        } else {
            if (power_distance_previous_a <= power_distance_previous_b) {
                Vector mid_point = (site_a + site_b) * 0.5;
                double param = dot(mid_point - previous_vertex, site_b - site_a) / dot(current_vertex - previous_vertex, site_b - site_a);
                Vector intersection = previous_vertex + param * (current_vertex - previous_vertex);
                result_polygon.vertices.push_back(intersection);
            }
        }
    }

    return result_polygon;
}

class SemiOptimalTransport {
public:
    std::vector<Vector> point_set;
    std::vector<double> lambda_values;
    Polygon boundary;
    std::vector<Polygon> laguerre_cells;
    std::vector<double> cell_weights;
    std::vector<Vector> velocities;

    SemiOptimalTransport(const std::vector<Vector>& points, const std::vector<double>& lambdas, const Polygon& edges)
        : point_set(points), lambda_values(lambdas), boundary(edges), cell_weights(points.size(), 1.0), velocities(points.size(), Vector(0, 0, 0)) {}

    ~SemiOptimalTransport() {}

    void optimize(int num_points) {
        std::cout << "Starting optimization..." << std::endl;
        double objective_value = 0.0;
        lbfgs(num_points, &cell_weights[0], &objective_value, evaluate_wrapper, progress_wrapper, this, NULL);
        compute_laguerre_cells();
        std::cout << "Finish optimization..." << std::endl;
        
    }

    static lbfgsfloatval_t evaluate_wrapper(
        void* instance,
        const lbfgsfloatval_t* weights,
        lbfgsfloatval_t* gradient,
        const int num_variables,
        const lbfgsfloatval_t step_size
    ) {
        return reinterpret_cast<SemiOptimalTransport*>(instance)->evaluate(weights, gradient, num_variables, step_size);
    }

    lbfgsfloatval_t evaluate(
        const lbfgsfloatval_t* weights,
        lbfgsfloatval_t* gradient,
        const int num_variables,
        const lbfgsfloatval_t step_size
    ) {
        lbfgsfloatval_t objective_value = 0.0;

        for (int i = 0; i < num_variables; i++)
            cell_weights[i] = weights[i];

        compute_laguerre_cells();

        double t1 = 0, t2 = 0, t3 = 0;
        for (int i = 0; i < num_variables; i++) {
            gradient[i] = laguerre_cells[i].area() - lambda_values[i];
            t1 += laguerre_cells[i].squared_distance_integral(point_set[i]);
            t2 -= weights[i] * laguerre_cells[i].area();
            t3 += weights[i] * lambda_values[i];
        }
        objective_value = -(t1 + t2 + t3);
        return objective_value;
    }

    static int progress_wrapper(
        void* instance,
        const lbfgsfloatval_t* weights,
        const lbfgsfloatval_t* gradient,
        const lbfgsfloatval_t objective_value,
        const lbfgsfloatval_t weights_norm,
        const lbfgsfloatval_t gradient_norm,
        const lbfgsfloatval_t step_size,
        int num_variables,
        int iteration_count,
        int line_search_count
    ) {
        return reinterpret_cast<SemiOptimalTransport*>(instance)->progress(weights, gradient, objective_value, weights_norm, gradient_norm, step_size, num_variables, iteration_count, line_search_count);
    }

    int progress(
        const lbfgsfloatval_t* weights,
        const lbfgsfloatval_t* gradient,
        const lbfgsfloatval_t objective_value,
        const lbfgsfloatval_t weights_norm,
        const lbfgsfloatval_t gradient_norm,
        const lbfgsfloatval_t step_size,
        int num_variables,
        int iteration_count,
        int line_search_count
    ) {
        return 0;
    }

    void update_particles(double delta_time) {
        const double damping_factor = 0.9;
        const double bounce_factor = 0.7;
        const double min_velocity = 1e-3;

        for (size_t i = 0; i < point_set.size(); ++i) {
            Vector gravity_force = Vector(0, -GRAVITY, 0);

            Vector cell_centroid = laguerre_cells[i].centroid();
            Vector spring_force = SPRING_CONSTANT * (cell_centroid - point_set[i]);

            Vector total_force = gravity_force + spring_force;

            velocities[i] = velocities[i] + delta_time * total_force;
            velocities[i] = velocities[i] * damping_factor;

            if (velocities[i].norm() < min_velocity) {
                velocities[i].normalize();
                velocities[i] = velocities[i] * min_velocity;
            }

            point_set[i] = point_set[i] + delta_time * velocities[i];

            if (point_set[i][1] < 0) {
                point_set[i][1] = 0;
                velocities[i][1] = -velocities[i][1] * bounce_factor;
            } else if (point_set[i][1] > 1) {
                point_set[i][1] = 1;
                velocities[i][1] = -velocities[i][1] * bounce_factor;
            }

            if (point_set[i][0] < 0) {
                point_set[i][0] = 0;
                velocities[i][0] = -velocities[i][0] * bounce_factor;
            } else if (point_set[i][0] > 1) {
                point_set[i][0] = 1;
                velocities[i][0] = -velocities[i][0] * bounce_factor;
            }

            if (std::isnan(point_set[i][0]) || std::isnan(point_set[i][1]) ||
                std::isnan(velocities[i][0]) || std::isnan(velocities[i][1])) {
                std::cerr << "Error: Particle " << i << " has NaN values." << std::endl;
                continue;
            }

            //std::cout << "Particle " << i << ": Position(" << point_set[i][0] << ", " << point_set[i][1] << "), " << "Velocity(" << velocities[i][0] << ", " << velocities[i][1] << ")" << std::endl;
        }
    }

    void validate_particles() {
        for (size_t i = 0; i < point_set.size(); ++i) {
            if (std::isnan(point_set[i][0]) || std::isnan(point_set[i][1]) || 
                std::isnan(velocities[i][0]) || std::isnan(velocities[i][1])) {
                std::cerr << "Error: Particle " << i << " has invalid position or velocity." << std::endl;
                exit(EXIT_FAILURE);
            }
        }
    }

private:
    void compute_laguerre_cells() {
        laguerre_cells.resize(point_set.size());
        Polygon initial_polygon = boundary;

        for (size_t i = 0; i < point_set.size(); ++i) {
            Polygon cell = initial_polygon;
            for (size_t j = 0; j < point_set.size(); ++j) {
                if (i != j) {
                    cell = clip_power_diagram(cell, point_set[i], cell_weights[i], point_set[j], cell_weights[j]);
                }
            }
            laguerre_cells[i] = cell;
            if (cell.vertices.size() < 3) {
                //std::cerr << "Warning: Computed Laguerre cell " << i << " has fewer than 3 vertices." << std::endl;
            }
        }
    }
};

int main() {
    auto start = std::chrono::high_resolution_clock::now();

    Polygon bounds({
        Vector(0., 0.), Vector(0., 1.),
        Vector(1., 1.), Vector(1., 0.)
    });
    std::vector<Vector> points(NUM_PARTICLES);
    std::vector<double> lambdas(NUM_PARTICLES, 1.0 / NUM_PARTICLES);

    for (int i = 0; i < NUM_PARTICLES; i++) {
        points[i] = Vector(uniform(engine), uniform(engine));
    }

    SemiOptimalTransport sd_ot(points, lambdas, bounds);

    for (int frame = 0; frame < NUM_FRAMES; frame++) {
        std::cout << "Frame " << frame << std::endl;
        sd_ot.optimize(NUM_PARTICLES);
        std::cout << "optimization done" << std::endl;
        sd_ot.update_particles(DELTA_TIME);
        std::cout << "update particle done" << std::endl;
        sd_ot.validate_particles();
        save_frame(sd_ot.laguerre_cells, sd_ot.point_set, "frames/frame_", frame);
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    std::cout << "Elapsed time: " << duration.count() << " milliseconds" << std::endl;

    return 0;
}
