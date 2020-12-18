import java.io.BufferedWriter
import java.io.FileWriter
import java.util.concurrent.Executors
import java.util.concurrent.ThreadLocalRandom
import kotlin.math.*

class SmallPT {

    private val spheres = arrayOf(
        Sphere(100000.0, Vec3(100001.0, 40.8, 81.6), Vec3.ZERO, Vec3(0.75, 0.25, 0.25), ReflectanceType.DIFFUSE),
        Sphere(100000.0, Vec3(-99901.0, 40.8, 81.6), Vec3.ZERO, Vec3(0.25, 0.25, 0.75), ReflectanceType.DIFFUSE),
        Sphere(100000.0, Vec3(50.0, 40.8, 100000.0), Vec3.ZERO, Vec3(0.75, 0.75, 0.75), ReflectanceType.DIFFUSE),
        Sphere(100000.0, Vec3(50.0, 40.8, -99830.0), Vec3.ZERO, Vec3.ZERO, ReflectanceType.DIFFUSE),
        Sphere(100000.0, Vec3(50.0, 100000.0, 81.6), Vec3.ZERO, Vec3(0.75, 0.75, 0.75), ReflectanceType.DIFFUSE),
        Sphere(100000.0, Vec3(50.0, -99918.4, 81.6), Vec3.ZERO, Vec3(0.75, 0.75, 0.75), ReflectanceType.DIFFUSE),
        Sphere(16.5, Vec3(27.0, 16.5, 47.0), Vec3.ZERO, Vec3(0.999, 0.999, 0.999), ReflectanceType.MIRROR),
        Sphere(16.5, Vec3(73.0, 16.5, 78.0), Vec3.ZERO, Vec3(0.700, 0.999, 0.900), ReflectanceType.GLASS),
        Sphere(600.0, Vec3(50.0, 681.33, 81.6), Vec3(12.0, 12.0, 12.0), Vec3.ZERO, ReflectanceType.DIFFUSE)
    )

    private fun radiance(ray: Ray, dep: Int): Vec3 {
        var obj = Sphere.NO_HIT
        var nearest = Double.POSITIVE_INFINITY
        for (sphere in spheres) {
            val dist = sphere.intersect(ray)
            if (dist > 0.0 && dist < nearest) {
                obj = sphere
                nearest = dist
            }
        }
        if (obj == Sphere.NO_HIT) return Vec3.ZERO

        var col = obj.color
        val x = ray.origin + ray.direction * nearest
        val n = Vec3.normalize(x - obj.position)
        val a = Vec3.dot(n, ray.direction)
        val nl = if (a < 0.0) n else -n
        val depth = dep + 1
        if (depth > 5) {
            val p = max(col.x, max(col.y, col.z))
            if (rnd() < p) col /= p else return obj.emission
        }

        return obj.emission + col * when (obj.type) {
            ReflectanceType.DIFFUSE -> {
                val r = rnd()
                val phi = 2.0 * Math.PI * rnd()
                val u = Vec3.normalize(if (abs(nl.x) > 0.1) Vec3(-nl.z, 0.0, nl.x) else Vec3(0.0, -nl.z, nl.y))
                radiance(Ray(x, (u * cos(phi) + Vec3.cross(nl, u) * sin(phi)) * sqrt(r) + nl * sqrt(1.0 - r)), depth)
            }
            ReflectanceType.MIRROR -> radiance(Ray(x, Vec3.reflect(ray.direction, n)), depth)
            else -> {
                val reflected = Ray(x, ray.direction - n * (2.0 * Vec3.dot(n, ray.direction)))
                val into = Vec3.dot(n, nl) > 0.0
                val nc = 1.0
                val nt = 1.5
                val nnt = if (into) nc / nt else nt / nc
                val ddn = Vec3.dot(ray.direction, nl)
                val cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn)
                return if (cos2t < 0.0) {
                    obj.emission + col * radiance(reflected, depth)
                } else {
                    val tdir = Vec3.normalize(ray.direction * nnt - n * ((if (into) 1.0 else -1.0) * (ddn * nnt + sqrt(cos2t))))
                    val r0 = (nt - nc).pow(2.0) / (nt + nc).pow(2.0)
                    val c = 1.0 - if (into) -ddn else Vec3.dot(tdir, n)
                    val re = r0 + (1.0 - r0) * c * c * c * c * c
                    val tr = 1.0 - re
                    val p = 0.25 + 0.5 * re
                    if (depth > 2) {
                        return obj.emission + col * if (Math.random() < p) {
                            radiance(reflected, depth) * (re / p)
                        } else {
                            radiance(Ray(x, tdir), depth) * (tr / (1.0 - p))
                        }
                    }
                    obj.emission + col * (radiance(reflected, depth) * re + radiance(Ray(x, tdir), depth) * tr)
                }
            }
        }
    }

    fun run(samples: Int) {
        val cam = Ray(Vec3(50.0, 52.0, 295.6), Vec3.normalize(Vec3(0.0, -0.042612, -1.0)))
        val cx = Vec3(WIDTH * 0.5135 / HEIGHT, 0.0, 0.0)
        val cy = Vec3.normalize(Vec3.cross(cx, cam.direction)).times(0.5135)

        val pixels = Array(WIDTH * HEIGHT) { Vec3.ZERO }
        val tasks = (0 until HEIGHT).map {
            Executors.callable {
                print("\rRendering (${samples * 4} spp) ${"%5.2f".format(100.0 * it / (HEIGHT - 1))}%")
                val y = (HEIGHT - 1) - it
                for (x in 0 until WIDTH) {
                    for (sy in 0..1) {
                        for (sx in 0..1) {
                            var rad = Vec3.ZERO
                            for (s in 0 until samples) {
                                val r1 = 2.0 * rnd()
                                val r2 = 2.0 * rnd()
                                val dx = if (r1 < 1.0) sqrt(r1) - 1.0 else 1.0 - sqrt(2.0 - r1)
                                val dy = if (r2 < 1.0) sqrt(r2) - 1.0 else 1.0 - sqrt(2.0 - r2)
                                val tx = cx * (((sx + 0.5 + dx) / 2.0 + x) / WIDTH - 0.5)
                                val ty = cy * (((sy + 0.5 + dy) / 2.0 + y) / HEIGHT - 0.5)
                                val d = cam.direction + tx + ty
                                rad += radiance(Ray(cam.origin + d * 137.0, Vec3.normalize(d)), 0)
                            }
                            rad /= samples.toDouble()
                            pixels[it * WIDTH + x] += Vec3(clamp(rad.x), clamp(rad.y), clamp(rad.z)) * 0.25
                        }
                    }
                }
            }
        }

        Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors()).apply {
            invokeAll(tasks)
            shutdown()
        }
        println()

        BufferedWriter(FileWriter("image-kotlin.ppm")).use {
            it.write("P3\n$WIDTH $HEIGHT\n255\n")
            for (col in pixels) it.write("${toInt(col.x)} ${toInt(col.y)} ${toInt(col.z)} ")
        }
    }

    companion object {
        const val WIDTH = 1024
        const val HEIGHT = 768

        @JvmStatic
        fun main(args: Array<String>) {
            val samples = if (args.isEmpty()) 1 else ceil(args[0].toInt() / 4.0).toInt()
            SmallPT().run(samples)
        }

        fun rnd() = ThreadLocalRandom.current().nextDouble()
        fun clamp(x: Double) = max(0.0, min(1.0, x))
        fun toInt(x: Double) = (clamp(x).pow(1.0 / 2.2) * 255.0 + 0.5).toInt()

        enum class ReflectanceType { DIFFUSE, MIRROR, GLASS }

        data class Vec3(val x: Double, val y: Double, val z: Double) {

            operator fun unaryMinus() = Vec3(-x, -y, -z)
            operator fun plus(v: Vec3) = Vec3(x + v.x, y + v.y, z + v.z)
            operator fun minus(v: Vec3) = Vec3(x - v.x, y - v.y, z - v.z)
            operator fun times(s: Double) = Vec3(x * s, y * s, z * s)
            operator fun times(v: Vec3) = Vec3(x * v.x, y * v.y, z * v.z)
            operator fun div(s: Double) = Vec3(x / s, y / s, z / s)

            companion object {
                val ZERO = Vec3(0.0, 0.0, 0.0)

                fun normalize(v: Vec3) = v / sqrt(dot(v, v))
                fun dot(a: Vec3, b: Vec3) = a.x * b.x + a.y * b.y + a.z * b.z
                fun cross(a: Vec3, b: Vec3) = Vec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x)
                fun reflect(v: Vec3, n: Vec3) = v - n * (2.0 * dot(v, n))
            }
        }

        data class Ray(val origin: Vec3, val direction: Vec3)

        class Sphere(val radius: Double, val position: Vec3, val emission: Vec3, val color: Vec3, val type: ReflectanceType) {

            fun intersect(ray: Ray): Double {
                val op = position - ray.origin
                val b = Vec3.dot(op, ray.direction)
                var det = b * b - Vec3.dot(op, op) + radius * radius
                if (det < 0.0) return 0.0
                det = sqrt(det)
                if (b - det > EPSILON) return b - det
                if (b + det > EPSILON) return b + det
                return 0.0
            }

            companion object {
                private const val EPSILON = 0.0001

                val NO_HIT = Sphere(0.0, Vec3.ZERO, Vec3.ZERO, Vec3.ZERO, ReflectanceType.DIFFUSE)
            }
        }
    }
}
