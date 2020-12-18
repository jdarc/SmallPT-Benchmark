import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.Arrays;
import java.util.concurrent.Executors;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.ThreadLocalRandom;
import java.util.stream.Collectors;

import static java.lang.Math.abs;
import static java.lang.Math.ceil;
import static java.lang.Math.cos;
import static java.lang.Math.max;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static java.util.stream.IntStream.range;

final class SmallPT {
    private static final int WIDTH = 1024;
    private static final int HEIGHT = 768;

    private final Sphere[] spheres = {
        new Sphere(100000.0, new Vec3(100001.0, 40.8, 81.6), Vec3.ZERO, new Vec3(0.75, 0.25, 0.25), ReflectanceType.DIFFUSE),
        new Sphere(100000.0, new Vec3(-99901.0, 40.8, 81.6), Vec3.ZERO, new Vec3(0.25, 0.25, 0.75), ReflectanceType.DIFFUSE),
        new Sphere(100000.0, new Vec3(50.0, 40.8, 100000.0), Vec3.ZERO, new Vec3(0.75, 0.75, 0.75), ReflectanceType.DIFFUSE),
        new Sphere(100000.0, new Vec3(50.0, 40.8, -99830.0), Vec3.ZERO, Vec3.ZERO, ReflectanceType.DIFFUSE),
        new Sphere(100000.0, new Vec3(50.0, 100000.0, 81.6), Vec3.ZERO, new Vec3(0.75, 0.75, 0.75), ReflectanceType.DIFFUSE),
        new Sphere(100000.0, new Vec3(50.0, -99918.4, 81.6), Vec3.ZERO, new Vec3(0.75, 0.75, 0.75), ReflectanceType.DIFFUSE),
        new Sphere(16.5, new Vec3(27.0, 16.5, 47.0), Vec3.ZERO, new Vec3(0.999, 0.999, 0.999), ReflectanceType.MIRROR),
        new Sphere(16.5, new Vec3(73.0, 16.5, 78.0), Vec3.ZERO, new Vec3(0.700, 0.999, 0.900), ReflectanceType.GLASS),
        new Sphere(600.0, new Vec3(50.0, 681.33, 81.6), new Vec3(12.0, 12.0, 12.0), Vec3.ZERO, ReflectanceType.DIFFUSE)
    };

    private Vec3 radiance(final Ray ray, int depth) {
        Sphere obj = null;
        var nearest = Double.POSITIVE_INFINITY;
        for (var sphere : spheres) {
            var dist = sphere.intersect(ray);
            if (dist > 0.0 && dist < nearest) {
                obj = sphere;
                nearest = dist;
            }
        }

        if (obj == null) {
            return Vec3.ZERO;
        }

        var col = obj.color;
        var x = ray.origin.plus(ray.direction.times(nearest));
        var n = Vec3.normalize(x.minus(obj.position));
        var nl = Vec3.dot(n, ray.direction) < 0.0 ? n : n.times(-1.0);
        if (++depth > 5) {
            var p = max(col.x, max(col.y, col.z));
            if (rnd() < p) {
                col = col.times(1.0 / p);
            } else {
                return obj.emission;
            }
        }

        switch (obj.type) {
            case DIFFUSE:
                var r = rnd();
                var phi = 2.0 * Math.PI * rnd();
                var u = Vec3.normalize(abs(nl.x) > 0.1 ? new Vec3(-nl.z, 0.0, nl.x) : new Vec3(0.0, -nl.z, nl.y));
                var d = u.times(cos(phi)).plus(Vec3.cross(nl, u).times(sin(phi))).times(sqrt(r)).plus(nl.times(sqrt(1.0 - r)));
                return obj.emission.plus(col.times(radiance(new Ray(x, d), depth)));
            case MIRROR:
                return obj.emission.plus(col.times(radiance(new Ray(x, Vec3.reflect(ray.direction, n)), depth)));
            default:
                var reflected = new Ray(x, ray.direction.minus(n.times(2.0 * Vec3.dot(n, ray.direction))));
                var into = Vec3.dot(n, nl) > 0.0;
                var nc = 1.0;
                var nt = 1.5;
                var nnt = into ? nc / nt : nt / nc;
                var ddn = Vec3.dot(ray.direction, nl);
                var cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);
                if (cos2t < 0.0) {
                    return obj.emission.plus(col.times(radiance(reflected, depth)));
                } else {
                    var tdir = Vec3.normalize(ray.direction.times(nnt).minus(n.times((into ? 1.0 : -1.0) * (ddn * nnt + Math.sqrt(cos2t)))));
                    var R0 = Math.pow(nt - nc, 2) / Math.pow(nt + nc, 2);
                    var c = 1.0 - (into ? -ddn : Vec3.dot(tdir, n));
                    var Re = R0 + (1.0 - R0) * c * c * c * c * c;
                    var Tr = 1.0 - Re;
                    var P = 0.25 + 0.5 * Re;
                    if (depth > 2) {
                        Vec3 v;
                        if (Math.random() < P) {
                            v = radiance(reflected, depth).times(Re / P);
                        } else {
                            v = radiance(new Ray(x, tdir), depth).times(Tr / (1.0 - P));
                        }
                        return obj.emission.plus(col.times(v));
                    }
                    var rad = col.times(radiance(reflected, depth).times(Re).plus(radiance(new Ray(x, tdir), depth).times(Tr)));
                    return obj.emission.plus(rad);
                }
        }
    }

    private void run(int samples) throws IOException {
        var cam = new Ray(new Vec3(50.0, 52.0, 295.6), Vec3.normalize(new Vec3(0.0, -0.042612, -1.0)));
        var cx = new Vec3(WIDTH * 0.5135 / HEIGHT, 0.0, 0.0);
        var cy = Vec3.normalize(Vec3.cross(cx, cam.direction)).times(0.5135);

        var pixels = new Vec3[WIDTH * HEIGHT];
        Arrays.fill(pixels, Vec3.ZERO);
        ForkJoinPool.commonPool().invokeAll(range(0, HEIGHT).mapToObj(y -> Executors.callable(() -> {
            System.out.printf("\rRendering (%d spp) %5.2f%%", samples * 4, 100.0 * y / (HEIGHT - 1));
            for (var x = 0; x < WIDTH; ++x) {
                var acc = Vec3.ZERO;
                for (var sy = 0; sy < 2; ++sy) {
                    for (var sx = 0; sx < 2; ++sx) {
                        var rad = Vec3.ZERO;
                        for (var s = 0; s < samples; ++s) {
                            var r1 = 2.0 * rnd();
                            var r2 = 2.0 * rnd();
                            var dx = r1 < 1.0 ? sqrt(r1) - 1.0 : 1.0 - sqrt(2.0 - r1);
                            var dy = r2 < 1.0 ? sqrt(r2) - 1.0 : 1.0 - sqrt(2.0 - r2);
                            var tx = cx.times(((sx + 0.5 + dx) / 2.0 + x) / WIDTH - 0.5);
                            var ty = cy.times(((sy + 0.5 + dy) / 2.0 + y) / HEIGHT - 0.5);
                            var d = cam.direction.plus(tx.plus(ty));
                            var ray = new Ray(cam.origin.plus(d.times(137)), Vec3.normalize(d));
                            rad = rad.plus(radiance(ray, 0));
                        }
                        rad = rad.times(1.0 / samples);
                        acc = acc.plus(new Vec3(clamp(rad.x), clamp(rad.y), clamp(rad.z)).times(0.25));
                    }
                }
                pixels[(HEIGHT - 1 - y) * WIDTH + x] = acc;
            }
        })).collect(Collectors.toList()));
        System.out.println();

        writeImage(pixels);
    }

    enum ReflectanceType {
        DIFFUSE,
        MIRROR,
        GLASS
    }

    static class Vec3 {
        final static Vec3 ZERO = new Vec3(0.0, 0.0, 0.0);

        final double x;
        final double y;
        final double z;

        Vec3(double x, double y, double z) {
            this.x = x;
            this.y = y;
            this.z = z;
        }

        Vec3 plus(Vec3 v) {
            return new Vec3(x + v.x, y + v.y, z + v.z);
        }

        Vec3 minus(Vec3 v) {
            return new Vec3(x - v.x, y - v.y, z - v.z);
        }

        Vec3 times(double s) {
            return new Vec3(x * s, y * s, z * s);
        }

        Vec3 times(Vec3 v) {
            return new Vec3(x * v.x, y * v.y, z * v.z);
        }

        static Vec3 normalize(Vec3 v) {
            return v.times(1.0 / sqrt(dot(v, v)));
        }

        static double dot(Vec3 a, Vec3 b) {
            return a.x * b.x + a.y * b.y + a.z * b.z;
        }

        static Vec3 cross(Vec3 a, Vec3 b) {
            return new Vec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
        }

        static Vec3 reflect(Vec3 v, Vec3 n) {
            return v.minus(n.times(2.0 * dot(v, n)));
        }
    }

    static class Ray {
        final Vec3 origin;
        final Vec3 direction;

        Ray(Vec3 origin, Vec3 direction) {
            this.origin = origin;
            this.direction = direction;
        }
    }

    static class Sphere {
        private static final double EPSILON = 0.0001;

        final double radius;
        final Vec3 position;
        final Vec3 emission;
        final Vec3 color;
        final ReflectanceType type;

        Sphere(double radius, Vec3 position, Vec3 emission, Vec3 color, ReflectanceType type) {
            this.radius = radius;
            this.position = position;
            this.emission = emission;
            this.color = color;
            this.type = type;
        }

        double intersect(Ray r) {
            var op = position.minus(r.origin);
            var b = Vec3.dot(op, r.direction);
            var det = b * b - Vec3.dot(op, op) + radius * radius;
            if (det < 0.0) {
                return 0.0;
            }
            det = sqrt(det);
            if (b - det > EPSILON) {
                return b - det;
            }
            if (b + det > EPSILON) {
                return b + det;
            }
            return 0.0;
        }
    }

    private static double rnd() {
        return ThreadLocalRandom.current().nextDouble();
    }

    private static double clamp(double n) {
        return max(0.0, Math.min(1.0, n));
    }

    private static int toInt(double n) {
        return (int)(255.0 * Math.pow(clamp(n), 1.0 / 2.2) + 0.5);
    }

    private void writeImage(Vec3[] pixels) throws IOException {
        try (var writer = new OutputStreamWriter(new FileOutputStream("image-java.ppm"))) {
            var sb = new StringBuilder();
            sb.append("P3\n");
            sb.append(WIDTH);
            sb.append(" ");
            sb.append(HEIGHT);
            sb.append("\n255\n");
            for (var col : pixels) {
                sb.append(toInt(col.x));
                sb.append(" ");
                sb.append(toInt(col.y));
                sb.append(" ");
                sb.append(toInt(col.z));
                sb.append(" ");
            }
            writer.write(sb.toString());
        }
    }

    public static void main(String[] args) throws IOException {
        var samples = args.length > 0 ? (int)ceil(Integer.parseInt(args[0]) / 4.0) : 1;
        new SmallPT().run(samples);
    }
}

