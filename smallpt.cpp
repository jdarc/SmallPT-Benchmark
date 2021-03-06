#include <cmath>
#include <cstdlib>
#include <cstdio>

struct Vec {
    double x, y, z;

    explicit Vec(double x_ = 0, double y_ = 0, double z_ = 0) {
        x = x_;
        y = y_;
        z = z_;
    }

    Vec operator+(const Vec &b) const { return Vec(x + b.x, y + b.y, z + b.z); }

    Vec operator-(const Vec &b) const { return Vec(x - b.x, y - b.y, z - b.z); }

    Vec operator*(double b) const { return Vec(x * b, y * b, z * b); }

    Vec mul(const Vec &b) const { return Vec(x * b.x, y * b.y, z * b.z); }

    Vec &norm() { return *this = *this * (1.0 / sqrt(x * x + y * y + z * z)); }

    double dot(const Vec &b) const { return x * b.x + y * b.y + z * b.z; }

    Vec operator%(Vec &b) const { return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }
};

struct Ray {
    Vec o, d;

    Ray(Vec o_, Vec d_) : o(o_), d(d_) {}
};

enum Refl_t {
    DIFF, SPEC, REFR
};

struct Sphere {
    double rad;
    Vec p, e, c;
    Refl_t refl;

    Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_) : rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {}

    double intersect(const Ray &r) const {
        Vec op = p - r.o;
        double t, eps = 1e-4, b = op.dot(r.d), det = b * b - op.dot(op) + rad * rad;
        if (det < 0) return 0; else det = sqrt(det);
        return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
    }
};

Sphere spheres[] = {
        Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(), Vec(.75, .25, .25), DIFF),
        Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(), Vec(.25, .25, .75), DIFF),
        Sphere(1e5, Vec(50, 40.8, 1e5), Vec(), Vec(.75, .75, .75), DIFF),
        Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Vec(), Vec(), DIFF),
        Sphere(1e5, Vec(50, 1e5, 81.6), Vec(), Vec(.75, .75, .75), DIFF),
        Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(), Vec(.75, .75, .75), DIFF),
        Sphere(16.5, Vec(27, 16.5, 47), Vec(), Vec(1, 1, 1) * .999, SPEC),
        Sphere(16.5, Vec(73, 16.5, 78), Vec(), Vec(0.700, 0.999, 0.900), REFR),
        Sphere(600, Vec(50, 681.6 - .27, 81.6), Vec(12, 12, 12), Vec(), DIFF)
};

inline double clamp(double x) { return x < 0 ? 0 : x > 1 ? 1 : x; }

inline int toInt(double x) { return (int)(pow(clamp(x), 1.0 / 2.2) * 255.0 + 0.5); }

inline bool intersect(const Ray &r, double &t, int &id) {
    double n = sizeof(spheres) / sizeof(Sphere);
    double d;
    double inf = t = 1e20;
    for (int i = int(n); i--;)
        if ((d = spheres[i].intersect(r)) > 0.0 && d < t) {
            t = d;
            id = i;
        }
    return t < inf;
}

Vec radiance(const Ray &r, int depth, unsigned short *Xi) {
    double t;
    int id = 0;
    if (!intersect(r, t, id)) return Vec();

    const Sphere &obj = spheres[id];
    Vec x = r.o + r.d * t;
    Vec n = (x - obj.p).norm();
    Vec nl = n.dot(r.d) < 0.0 ? n : n * -1.0;

    Vec f = obj.c;
    if (++depth > 5) {
        double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z;
        if (depth < 256 && erand48(Xi) < p) f = f * (1 / p); else return obj.e;
    }

    if (obj.refl == DIFF) {
        double r2 = erand48(Xi);
        double r2s = sqrt(r2);
        double r1 = 2.0 * M_PI * erand48(Xi);
        Vec w = nl;
        Vec u = ((fabs(w.x) > 0.1 ? Vec(0, 1) : Vec(1)) % w).norm();
        Vec v = w % u;
        Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();
        return obj.e + f.mul(radiance(Ray(x, d), depth, Xi));
    }

    if (obj.refl == SPEC) {
        return obj.e + f.mul(radiance(Ray(x, r.d - n * 2 * n.dot(r.d)), depth, Xi));
    }

    Ray reflRay(x, r.d - n * 2.0 * n.dot(r.d));
    bool into = n.dot(nl) > 0.0;
    double nc = 1.0;
    double nt = 1.5;
    double nnt = into ? nc / nt : nt / nc;
    double ddn = r.d.dot(nl);
    double cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);
    if (cos2t < 0.0) {
        return obj.e + f.mul(radiance(reflRay, depth, Xi));
    }

    Vec tdir = (r.d * nnt - n * ((into ? 1.0 : -1.0) * (ddn * nnt + sqrt(cos2t)))).norm();
    double a = nt - nc;
    double b = nt + nc;
    double R0 = a * a / (b * b);
    double c = 1.0 - (into ? -ddn : tdir.dot(n));
    double Re = R0 + (1.0 - R0) * pow(c, 5);
    double Tr = 1.0 - Re;
    double P = 0.25 + 0.5 * Re;
    double RP = Re / P;
    double TP = Tr / (1.0 - P);
    return obj.e + f.mul(depth > 2 ? (erand48(Xi) < P ?
        radiance(reflRay, depth, Xi) * RP : radiance(Ray(x, tdir), depth, Xi) * TP) :
        radiance(reflRay, depth, Xi) * Re + radiance(Ray(x, tdir), depth, Xi) * Tr);
}

int main(int argc, char *argv[]) {
    int samples = argc == 2 ? atoi(argv[1]) / 4 : 1;
    unsigned short w = 1024;
    unsigned short h = 768;
    Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm());
    Vec cx = Vec(w * .5135 / h), cy = (cx % cam.d).norm() * .5135, r, *c = new Vec[w * h];

#pragma omp parallel for schedule(dynamic, 1) private(r)
    for (int y = 0; y < h; y++) {
        fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samples * 4, 100. * y / (h - 1));
        for (unsigned short x = 0, Xi[3] = {0, 0, (unsigned short) (y * y * y)}; x < w; x++) {
            int i = (h - y - 1) * w + x;
            for (int sy = 0; sy < 2; sy++) {
                for (int sx = 0; sx < 2; sx++) {
                    Vec r = Vec();
                    for (int s = 0; s < samples; s++) {
                        double r1 = 2 * erand48(Xi), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                        double r2 = 2 * erand48(Xi), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                        Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) + cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;
                        r = r + radiance(Ray(cam.o + d * 140, d.norm()), 0, Xi) * (1. / samples);
                    }
                    c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * 0.25;
                }
            }
        }
    }

    fprintf(stderr, "\rRendering (%d spp) 100.00%%\n\n", samples * 4);

    FILE *f = fopen("image-cpp.ppm", "w");
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    for (int i = 0; i < w * h; i++) {
        fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
    }
}
