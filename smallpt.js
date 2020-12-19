const fs = require('fs');
const os = require('os');
const {Worker, isMainThread, workerData, parentPort} = require("worker_threads");

const ReflectanceType = {
    DIFFUSE: 0,
    MIRROR: 1,
    GLASS: 2
}

class ObjectPool {

    constructor(size, generator) {
        this.index = 0;
        this.pool = new Array(size);
        this.generator = generator;
        for (let i = 0; i < this.pool.length; i++) {
            this.pool[i] = generator();
        }
    }

    next = () => {
        this.#ensureCapacity();
        return this.pool[this.index++];
    };

    reset = () => this.index = 0;

    #ensureCapacity = () => {
        if (this.index >= this.pool.length) {
            const more = new Array(this.pool.length);
            for (let i = 0, length = more.length; i < length; ++i) more[i] = this.generator();
            this.pool = [...this.pool, ...more];
        }
    };
}

class Vec3 {
    static POOL = new ObjectPool(16384, () => new Vec3(0.0, 0.0, 0.0));

    static UP = new Vec3(0.0, 1.0, 0.0);
    static ZERO = new Vec3(0.0, 0.0, 0.0);
    static RIGHT = new Vec3(1.0, 0.0, 0.0);

    constructor(x = 0.0, y = 0.0, z = 0.0) {
        this.set(x, y, z);
    }

    set = (x, y, z) => {
        this.x = x;
        this.y = y;
        this.z = z;
        return this;
    };

    length = () => Math.sqrt(Vec3.dot(this, this));

    plus = (v, dst = Vec3.POOL.next()) => dst.set(this.x + v.x, this.y + v.y, this.z + v.z);
    minus = (v, dst = Vec3.POOL.next()) => dst.set(this.x - v.x, this.y - v.y, this.z - v.z);
    times = (v, dst = Vec3.POOL.next()) => dst.set(this.x * v.x, this.y * v.y, this.z * v.z);
    scale = (s, dst = Vec3.POOL.next()) => dst.set(this.x * s, this.y * s, this.z * s);
    clamp = (dst = Vec3.POOL.next()) => dst.set(Math.min(1.0, this.x), Math.min(1.0, this.y), Math.min(1.0, this.z));

    static dot = (a, b) => a.x * b.x + a.y * b.y + a.z * b.z
    static cross = (a, b, dst = Vec3.POOL.next()) => dst.set(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x)
    static normalize = (v, dst = Vec3.POOL.next()) => v.scale(1.0 / v.length(), dst);
    static reflect = (v, n, dst = Vec3.POOL.next()) => v.minus(n.scale(2.0 * Vec3.dot(v, n)), dst);
    static refract = (i, n, eta, dst = Vec3.POOL.next()) => {
        const k = 1.0 - eta * eta * (1.0 - Vec3.dot(n, i) * Vec3.dot(n, i));
        return k < 0.0 ? Vec3.ZERO : i.scale(eta).minus(n.scale(eta * Vec3.dot(n, i) + Math.sqrt(k)), dst);
    }
}

class Ray {
    static POOL = new ObjectPool(2048, () => new Ray(new Vec3(), new Vec3()));

    constructor(origin, direction) {
        this.set(origin, direction);
    }

    set(origin, direction) {
        this.origin = origin;
        this.direction = direction;
        return this;
    }

    static build = (origin, direction) => Ray.POOL.next().set(origin, direction);
}

class Sphere {
    static EPSILON = 0.0001;

    constructor(radius, position, emission, color, type) {
        this.radius = radius;
        this.position = position;
        this.emission = emission;
        this.color = color;
        this.type = type;
    }

    intersect = ray => {
        const x = this.position.x - ray.origin.x;
        const y = this.position.y - ray.origin.y;
        const z = this.position.z - ray.origin.z;
        const b = x * ray.direction.x + y * ray.direction.y + z * ray.direction.z;
        const lenSqr = x * x + y * y + z * z;
        let det = b * b - lenSqr + this.radius * this.radius;
        if (det < 0.0) {
            return 0.0;
        }
        det = Math.sqrt(det);
        let t;
        return (t = b - det) > Sphere.EPSILON ? t : (t = b + det) > Sphere.EPSILON ? t : 0.0;
    };
}

const mainProgram = () => {
    (async () => {
        const width = 1024;
        const height = 768;
        const args = process.argv.slice(2);
        const samples = args.length > 0 ? ~~Math.ceil(parseInt(args[0]) / 4.0) : 1;
        const cpuCount = os.cpus().length;

        let progress = 0;
        const tasks = [];
        for (let y = 0; y < height; y += cpuCount) {
            const workerData = {width, height, samples, y1: y, y2: y + cpuCount};
            tasks.push(new Promise((resolve, reject) => {
                const worker = new Worker('./smallpt.js', {workerData});
                worker.on('message', evt => {
                    if (evt.msg === "processing") {
                        process.stdout.write(`\rRendering (${samples * 4} spp) ${Math.min(100.0, ++progress * 100.0 / (height - 1)).toFixed(2)}%`);
                    } else {
                        resolve(evt.data);
                    }
                });
                worker.on('error', reject);
                worker.on('exit', code => {
                    if (code !== 0) {
                        reject(new Error(`Worker stopped with exit code ${code}`));
                    }
                });
            }));
        }

        const scanLines = await Promise.all(tasks);

        process.stdout.write(`\rRendering (${samples * 4} spp) 100.00%\n\n`)

        const toInt = n => ~~(Math.pow(Math.max(0.0, Math.min(1.0, n)), 1.0 / 2.2) * 255.0 + 0.5);
        let sb = `P3\n${width} ${height}\n255\n`;
        for (let y = 0; y < scanLines.length; y++) {
            const scanline = scanLines[y];
            for (let x = 0; x < scanline.length; x++) {
                sb += `${toInt(scanline[x])} `;
            }
        }
        fs.writeFile("image-js.ppm", sb, err => {
            if (err) console.log(err);
        });

    })().catch(err => console.log(err));
};

const workerProgram = () => {
    const radiance = (ray, spheres, depth) => {
        let hit = null;
        let nearest = Number.POSITIVE_INFINITY;
        for (const sphere of spheres) {
            const dist = sphere.intersect(ray);
            if (dist > 0.0 && dist < nearest) {
                hit = sphere;
                nearest = dist;
            }
        }

        if (hit == null) {
            return Vec3.ZERO;
        }

        const x = ray.origin.plus(ray.direction.scale(nearest));
        const n = Vec3.normalize(x.minus(hit.position));
        const nl = Vec3.dot(n, ray.direction) < 0.0 ? n : n.scale(-1.0);
        let f = hit.color;
        const p = Math.max(f.x, Math.max(f.y, f.z));
        if (++depth > 5) {
            if (depth < 256 && Math.random() < p) {
                f = f.scale(1.0 / p);
            } else {
                return hit.emission;
            }
        }

        if (hit.type === ReflectanceType.DIFFUSE) {
            const r1 = 2.0 * Math.PI * Math.random();
            const r2 = Math.random();
            const r2s = Math.sqrt(r2);
            const vector3 = Math.abs(nl.x) > 0.1 ? Vec3.UP : Vec3.RIGHT;
            const u = Vec3.normalize(Vec3.cross(vector3, nl));
            const d = Vec3.normalize(u.scale(Math.cos(r1) * r2s).plus(Vec3.cross(nl, u).scale(Math.sin(r1) * r2s)).plus(nl.scale(Math.sqrt(1.0 - r2))));
            return hit.emission.plus(f.times(radiance(Ray.build(x, d), spheres, depth)));
        }

        if (hit.type === ReflectanceType.MIRROR) {
            return hit.emission.plus(f.times(radiance(Ray.build(x, Vec3.reflect(ray.direction, n)), spheres, depth)));
        }

        const into = Vec3.dot(n, nl) > 0.0;
        const nc = 1.0;
        const nt = 1.5;
        const nnt = into ? nc / nt : nt / nc;
        const ddn = Vec3.dot(ray.direction, nl);
        const cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);
        if (cos2t < 0.0) {
            return hit.emission.plus(f.times(radiance(Ray.build(x, Vec3.reflect(ray.direction, n)), spheres, depth)));
        }

        const tdir = Vec3.refract(ray.direction, nl, nnt);
        const R0 = Math.pow((nt - nc) / (nt + nc), 2);
        const Re = R0 + (1.0 - R0) * Math.pow(1.0 - (into ? -ddn : Vec3.dot(tdir, n)), 5);
        const Tr = 1.0 - Re;
        if (depth > 2) {
            const P = 0.25 + 0.5 * Re;
            const v = Math.random() < P ?
                radiance(Ray.build(x, Vec3.reflect(ray.direction, n)), spheres, depth).scale(Re / P) :
                radiance(Ray.build(x, tdir), spheres, depth).scale(Tr / (1.0 - P));
            return hit.emission.plus(f.times(v));
        }
        const reflRay = Ray.build(x, Vec3.reflect(ray.direction, n));
        const rad = radiance(reflRay, spheres, depth).scale(Re).plus(radiance(Ray.build(x, tdir), spheres, depth).scale(Tr));
        return hit.emission.plus(f.times(rad));
    }

    const spheres = [
        new Sphere(100000.0, new Vec3(100001.0, 40.8, 81.6), Vec3.ZERO, new Vec3(0.75, 0.25, 0.25), ReflectanceType.DIFFUSE),
        new Sphere(100000.0, new Vec3(-99901.0, 40.8, 81.6), Vec3.ZERO, new Vec3(0.25, 0.25, 0.75), ReflectanceType.DIFFUSE),
        new Sphere(100000.0, new Vec3(50.0, 40.8, 100000.0), Vec3.ZERO, new Vec3(0.75, 0.75, 0.75), ReflectanceType.DIFFUSE),
        new Sphere(100000.0, new Vec3(50.0, 40.8, -99830.0), Vec3.ZERO, Vec3.ZERO, ReflectanceType.DIFFUSE),
        new Sphere(100000.0, new Vec3(50.0, 100000.0, 81.6), Vec3.ZERO, new Vec3(0.75, 0.75, 0.75), ReflectanceType.DIFFUSE),
        new Sphere(100000.0, new Vec3(50.0, -99918.4, 81.6), Vec3.ZERO, new Vec3(0.75, 0.75, 0.75), ReflectanceType.DIFFUSE),
        new Sphere(16.5, new Vec3(27.0, 16.5, 47.0), Vec3.ZERO, new Vec3(0.999, 0.999, 0.999), ReflectanceType.MIRROR),
        new Sphere(16.5, new Vec3(73.0, 16.5, 78.0), Vec3.ZERO, new Vec3(0.700, 0.999, 0.900), ReflectanceType.GLASS),
        new Sphere(600.0, new Vec3(50.0, 681.33, 81.6), new Vec3(12.0, 12.0, 12.0), Vec3.ZERO, ReflectanceType.DIFFUSE)
    ];

    const width = workerData.width;
    const height = workerData.height;
    const samples = workerData.samples;
    const y1 = workerData.y1;
    const y2 = workerData.y2;
    const data = new Float64Array(new SharedArrayBuffer(width * (y2 - y1) * 3 * Float64Array.BYTES_PER_ELEMENT));

    const cam = new Ray(new Vec3(50.0, 52.0, 295.6), Vec3.normalize(new Vec3(0.0, -0.042612, -1.0), new Vec3()));
    const cx = new Vec3(width * 0.5135 / height, 0.0, 0.0);
    const cy = Vec3.normalize(Vec3.cross(cx, cam.direction)).scale(0.5135, new Vec3());

    let r = new Vec3();
    let acc = new Vec3();
    let mem = 0;
    for (let i = y1; i < y2; i++) {
        parentPort.postMessage({msg: "processing", data: i});
        const y = height - 1.0 - i;
        for (let x = 0; x < width; x++) {
            acc.set(0.0, 0.0, 0.0);
            for (let sy = 0; sy < 2; sy++) {
                for (let sx = 0; sx < 2; sx++) {
                    r.set(0.0, 0.0, 0.0);
                    for (let s = 0; s < samples; s++) {
                        Vec3.POOL.reset();
                        Ray.POOL.reset();
                        const r1 = 2.0 * Math.random();
                        const r2 = 2.0 * Math.random();
                        const dx = r1 < 1.0 ? Math.sqrt(r1) - 1.0 : 1.0 - Math.sqrt(2.0 - r1);
                        const dy = r2 < 1.0 ? Math.sqrt(r2) - 1.0 : 1.0 - Math.sqrt(2.0 - r2);
                        const tx = cx.scale(((sx + 0.5 + dx) / 2.0 + x) / width - 0.5);
                        const ty = cy.scale(((sy + 0.5 + dy) / 2.0 + y) / height - 0.5);
                        const d = cam.direction.plus(tx.plus(ty));
                        r.plus(radiance(Ray.build(cam.origin.plus(d.scale(137.0)), Vec3.normalize(d)), spheres, 0), r);
                    }
                    acc.plus(r.scale(1.0 / samples, r).clamp(r).scale(0.25, r), acc);
                }
            }
            data[mem++] = acc.x;
            data[mem++] = acc.y;
            data[mem++] = acc.z;
        }
    }

    parentPort.postMessage({msg: "done", data});
};

if (isMainThread) mainProgram(); else workerProgram();
