using System;
using System.IO;
using System.Runtime.CompilerServices;
using System.Threading;
using System.Threading.Tasks;
using static System.Math;

namespace SmallPT
{
    internal enum ReflectanceType
    {
        Diffuse,
        Mirror,
        Glass
    }

    internal readonly struct Vec3
    {
        public static readonly Vec3 Zero = new(0.0, 0.0, 0.0);

        public readonly double X;
        public readonly double Y;
        public readonly double Z;

        public Vec3(in double x, in double y, in double z)
        {
            X = x;
            Y = y;
            Z = z;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vec3 operator +(in Vec3 a, in Vec3 b) => new(a.X + b.X, a.Y + b.Y, a.Z + b.Z);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vec3 operator -(in Vec3 v) => new(-v.X, -v.Y, -v.Z);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vec3 operator -(in Vec3 a, in Vec3 b) => new(a.X - b.X, a.Y - b.Y, a.Z - b.Z);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vec3 operator *(in Vec3 a, in double s) => new(a.X * s, a.Y * s, a.Z * s);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vec3 operator *(in Vec3 a, in Vec3 b) => new(a.X * b.X, a.Y * b.Y, a.Z * b.Z);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vec3 operator /(in Vec3 a, in double s) => new(a.X / s, a.Y / s, a.Z / s);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Dot(in Vec3 a, in Vec3 b) => a.X * b.X + a.Y * b.Y + a.Z * b.Z;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vec3 Cross(in Vec3 a, in Vec3 b) => new(a.Y * b.Z - a.Z * b.Y, a.Z * b.X - a.X * b.Z, a.X * b.Y - a.Y * b.X);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vec3 Normalize(in Vec3 v) => v / Sqrt(Dot(v, v));

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vec3 Reflect(in Vec3 l, in Vec3 n) => l - n * 2.0 * Dot(l, n);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vec3 Refract(Vec3 i, Vec3 n, double eta)
        {
            var k = 1.0 - eta * eta * (1.0 - Dot(n, i) * Dot(n, i));
            return k < 0.0 ? Zero : i * eta - (n * (eta * Dot(n, i) + Sqrt(k)));
        }
    }

    internal readonly struct Ray
    {
        public readonly Vec3 Origin;
        public readonly Vec3 Direction;

        public Ray(in Vec3 origin, in Vec3 direction)
        {
            Origin = origin;
            Direction = direction;
        }
    }

    internal sealed class Sphere
    {
        private const double Epsilon = 0.0001;

        public readonly ReflectanceType Type;
        public readonly double Radius;
        public Vec3 Position;
        public Vec3 Emission;
        public Vec3 Color;

        public Sphere(in double radius, in Vec3 position, in Vec3 emission, in Vec3 color, in ReflectanceType type)
        {
            Radius = radius;
            Position = position;
            Emission = emission;
            Color = color;
            Type = type;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public double Intersect(in Ray ray)
        {
            var op = Position - ray.Origin;
            var b = Vec3.Dot(op, ray.Direction);
            var det = b * b - Vec3.Dot(op, op) + Radius * Radius;
            if (det < 0.0) return 0.0;
            det = Sqrt(det);
            if (b - det > Epsilon) return b - det;
            if (b + det > Epsilon) return b + det;
            return 0.0;
        }
    }

    internal sealed class Program
    {
        private readonly Sphere[] _spheres =
        {
            new(100000.0, new Vec3(100001.0, 40.8, 81.6), Vec3.Zero, new Vec3(0.75, 0.25, 0.25), ReflectanceType.Diffuse),
            new(100000.0, new Vec3(-99901.0, 40.8, 81.6), Vec3.Zero, new Vec3(0.25, 0.25, 0.75), ReflectanceType.Diffuse),
            new(100000.0, new Vec3(50.0, 40.8, 100000.0), Vec3.Zero, new Vec3(0.75, 0.75, 0.75), ReflectanceType.Diffuse),
            new(100000.0, new Vec3(50.0, 40.8, -99830.0), Vec3.Zero, Vec3.Zero, ReflectanceType.Diffuse),
            new(100000.0, new Vec3(50.0, 100000.0, 81.6), Vec3.Zero, new Vec3(0.75, 0.75, 0.75), ReflectanceType.Diffuse),
            new(100000.0, new Vec3(50.0, -99918.4, 81.6), Vec3.Zero, new Vec3(0.75, 0.75, 0.75), ReflectanceType.Diffuse),
            new(16.5, new Vec3(27.0, 16.5, 47.0), Vec3.Zero, new Vec3(0.999, 0.999, 0.999), ReflectanceType.Mirror),
            new(16.5, new Vec3(73.0, 16.5, 78.0), Vec3.Zero, new Vec3(0.700, 0.999, 0.900), ReflectanceType.Glass),
            new(600.0, new Vec3(50.0, 681.33, 81.6), new Vec3(12.0, 12.0, 12.0), Vec3.Zero, ReflectanceType.Diffuse)
        };

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private double Intersect(in Ray ray, out Sphere hit)
        {
            hit = null;
            var nearest = double.PositiveInfinity;
            foreach (var sphere in _spheres)
            {
                var dist = sphere.Intersect(ray);
                if (dist > 0.0 && dist < nearest)
                {
                    hit = sphere;
                    nearest = dist;
                }
            }

            return nearest;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining | MethodImplOptions.AggressiveOptimization)]
        private Vec3 Radiance(in Ray ray, int depth)
        {
            var nearest = Intersect(ray, out var obj);
            if (obj == null) return Vec3.Zero;

            var x = ray.Origin + ray.Direction * nearest;
            var n = Vec3.Normalize(x - obj.Position);
            var nl = Vec3.Dot(n, ray.Direction) < 0.0 ? n : -n;
            
            var color = obj.Color;
            if (++depth > 5)
            {
                var maxComp = Max(color.X, Max(color.Y, color.Z));
                if (depth < 256 && Rnd() < maxComp) color /= maxComp; else return obj.Emission;
            }

            if (obj.Type == ReflectanceType.Diffuse)
            {
                var r = Rnd();
                var phi = 2.0 * PI * Rnd();
                var u = Vec3.Normalize(Abs(nl.X) > 0.1 ? new Vec3(-nl.Z, 0.0, nl.X) : new Vec3(0.0, -nl.Z, nl.Y));
                var nr = new Ray(x, (u * Cos(phi) + Vec3.Cross(nl, u) * Sin(phi)) * Sqrt(r) + nl * Sqrt(1.0 - r));
                return obj.Emission + color * Radiance(nr, depth);
            }

            if (obj.Type == ReflectanceType.Mirror)
            {
                return obj.Emission + color * Radiance(new Ray(x, Vec3.Reflect(ray.Direction, n)), depth);
            }

            var into = Vec3.Dot(n, nl) > 0.0;
            const double nc = 1.0;
            const double nt = 1.5;
            var nnt = into ? nc / nt : nt / nc;
            var ddn = Vec3.Dot(ray.Direction, nl);
            var cos2T = 1.0 - nnt * nnt * (1.0 - ddn * ddn);
            if (cos2T < 0.0)
            {
                return obj.Emission + color * Radiance(new Ray(x, Vec3.Reflect(ray.Direction, n)), depth);
            }

            var tDir = Vec3.Refract(ray.Direction, nl, nnt);
            var r0 = Pow((nt - nc) / (nt + nc), 2);
            var re = r0 + (1.0 - r0) * Pow(1.0 - (into ? -ddn : Vec3.Dot(tDir, n)), 5);
            var tr = 1.0 - re;
            if (depth > 2)
            {
                var p = 0.25 + 0.5 * re;
                return Rnd() < p
                    ? obj.Emission + color * Radiance(new Ray(x, Vec3.Reflect(ray.Direction, n)), depth) * (re / p)
                    : obj.Emission + color * Radiance(new Ray(x, tDir), depth) * (tr / (1.0 - p));
            }

            return obj.Emission + color * (Radiance(new Ray(x, Vec3.Reflect(ray.Direction, n)), depth) * re +
                                           Radiance(new Ray(x, tDir), depth) * tr);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining | MethodImplOptions.AggressiveOptimization)]
        private void Run(int samples)
        {
            const int width = 1024;
            const int height = 768;

            var cam = new Ray(new Vec3(50.0, 52.0, 295.6), Vec3.Normalize(new Vec3(0.0, -0.042612, -1.0)));
            var cx = new Vec3(width * 0.5135 / height, 0.0, 0.0);
            var cy = Vec3.Normalize(Vec3.Cross(cx, cam.Direction)) * 0.5135;

            var pixels = new Vec3[width * height];

            var progress = 0;
            Parallel.For(0, height, y =>
            {
                Console.Write($"\rRendering ({samples * 4} spp) {100.0 * Interlocked.Increment(ref progress) / height:0.0#}%");
                var offset = (height - y - 1) * width;
                for (var x = 0; x < width; ++x)
                {
                    for (var sy = 0; sy < 2; ++sy)
                    {
                        for (var sx = 0; sx < 2; ++sx)
                        {
                            var rad = Vec3.Zero;
                            for (var s = 0; s < samples; ++s)
                            {
                                var r1 = 2.0 * Rnd();
                                var r2 = 2.0 * Rnd();
                                var dx = r1 < 1.0 ? Sqrt(r1) - 1.0 : 1.0 - Sqrt(2.0 - r1);
                                var dy = r2 < 1.0 ? Sqrt(r2) - 1.0 : 1.0 - Sqrt(2.0 - r2);
                                var tx = cx * (((sx + 0.5 + dx) / 2.0 + x) / width - 0.5);
                                var ty = cy * (((sy + 0.5 + dy) / 2.0 + y) / height - 0.5);
                                var d = cam.Direction + tx + ty;
                                rad += Radiance(new Ray(cam.Origin + d * 137.0, Vec3.Normalize(d)), 0);
                            }

                            rad /= samples;
                            pixels[offset + x] += new Vec3(Clamp(rad.X), Clamp(rad.Y), Clamp(rad.Z)) * 0.25;
                        }
                    }
                }
            });
            Console.WriteLine();

            var ppm = new StreamWriter("image-cs.ppm", false);
            ppm.Write($"P3\n{width} {height}\n255\n");
            foreach (var xyz in pixels) ppm.Write($"{ToInt(xyz.X)} {ToInt(xyz.Y)} {ToInt(xyz.Z)} ");
            ppm.Close();
        }

        private static readonly ThreadLocal<Random> Random = new(() => new Random(Environment.TickCount));

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double Rnd() => Random.Value.NextDouble();

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double Clamp(double n) => Math.Clamp(n, 0.0, 1.0);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static int ToInt(double n) => (int) (Pow(Clamp(n), 1.0 / 2.2) * 255.0 + 0.5);

        public static void Main(string[] args) => new Program().Run(args.Length > 0 ? (int) Ceiling(int.Parse(args[0]) / 4.0) : 1);
    }
}
