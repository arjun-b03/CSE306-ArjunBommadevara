#define _CRT_SECURE_NO_WARNINGS 1
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define ITERATIONS 100

class Vector {
public:
    double x, y, z;
    Vector() : x(0), y(0), z(0) {}
    Vector(double x, double y, double z) : x(x), y(y), z(z) {}

    double dot(const Vector& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    Vector operator+(const Vector& other) const {
        return Vector(x + other.x, y + other.y, z + other.z);
    }

    Vector operator-(const Vector& other) const {
        return Vector(x - other.x, y - other.y, z - other.z);
    }

    Vector operator*(double scalar) const {
        return Vector(x * scalar, y * scalar, z * scalar);
    }

    Vector& operator+=(const Vector& other) {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }

    static Vector random_direction() {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_real_distribution<> dis(0, 1);

        double r1 = dis(gen);
        double r2 = dis(gen);
        double theta = 2 * M_PI * r1;
        double phi = acos(1 - 2 * r2);
        double x = sin(phi) * cos(theta);
        double y = sin(phi) * sin(theta);
        double z = cos(phi);
        return Vector(x, y, z);
    }
};

int clamp(int value, int min, int max) {
    if (value < min) return min;
    if (value > max) return max;
    return value;
}

void match_color(int width, int height, int channels, unsigned char* sourceImage, unsigned char* modelImage) {
    size_t pixels = width * height;
    std::vector<std::pair<double, int>> source(pixels);
    std::vector<std::pair<double, int>> model(pixels);

    for (int iter = 0; iter < ITERATIONS; ++iter) {
        Vector random_dir = Vector::random_direction();

        for (size_t i = 0; i < pixels; ++i) {
            Vector source_pixel(sourceImage[i * channels], sourceImage[i * channels + 1], sourceImage[i * channels + 2]);
            Vector model_pixel(modelImage[i * channels], modelImage[i * channels + 1], modelImage[i * channels + 2]);

            source[i] = { source_pixel.dot(random_dir), static_cast<int>(i) };
            model[i] = { model_pixel.dot(random_dir), static_cast<int>(i) };
        }

        std::sort(source.begin(), source.end());
        std::sort(model.begin(), model.end());

        for (size_t i = 0; i < pixels; ++i) {
            int sourceIndex = source[i].second;
            int modelIndex = model[i].second;

            Vector source_pixel(sourceImage[sourceIndex * channels], sourceImage[sourceIndex * channels + 1], sourceImage[sourceIndex * channels + 2]);
            Vector model_pixel(modelImage[modelIndex * channels], modelImage[modelIndex * channels + 1], modelImage[modelIndex * channels + 2]);

            Vector adjustment = random_dir * (model[i].first - source[i].first);

            source_pixel += adjustment;

            sourceImage[sourceIndex * channels] = clamp(static_cast<int>(source_pixel.x), 0, 255);
            sourceImage[sourceIndex * channels + 1] = clamp(static_cast<int>(source_pixel.y), 0, 255);
            sourceImage[sourceIndex * channels + 2] = clamp(static_cast<int>(source_pixel.z), 0, 255);
        }
    }
}

int main(int argc, char** argv) {
    int width, height, channels;
    int modelWidth, modelHeight, modelChannels;

    unsigned char* sourceImage = stbi_load("imgA.jpg", &width, &height, &channels, 0);
    unsigned char* modelImage = stbi_load("redim.jpg", &modelWidth, &modelHeight, &modelChannels, 0);


    match_color(width, height, channels, sourceImage, modelImage);

    stbi_image_free(sourceImage);
    stbi_image_free(modelImage);

    return 0;
}
