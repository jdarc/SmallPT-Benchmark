extern crate rand;

use std::{env, f64};
use std::f64::consts::PI;
use std::io::Write;
use std::ops::{Add, Div, Mul, Neg, Sub};
use std::sync::atomic::{AtomicUsize, Ordering};

use rayon::prelude::*;

const EPSILON: f64 = 0.0001;
const WIDTH: usize = 1024;
const HEIGHT: usize = 768;

#[inline]
fn clamp(x: f64) -> f64 { x.max(0.0).min(1.0) }

#[inline]
fn to_int(x: f64) -> u32 { (clamp(x).powf(1.0 / 2.2) * 255.0 + 0.5) as u32 }

enum ReflectanceType {
    Diffuse,
    Mirror,
    Glass,
}

#[derive(Copy, Clone)]
struct Vector3 {
    x: f64,
    y: f64,
    z: f64,
}

impl Vector3 {
    fn zero() -> Vector3 { Vector3::new(0.0, 0.0, 0.0) }

    fn new(x: f64, y: f64, z: f64) -> Vector3 { Vector3 { x, y, z } }

    #[inline]
    fn length(self) -> f64 {
        Vector3::dot(self, self).sqrt()
    }

    #[inline]
    fn dot(a: Vector3, b: Vector3) -> f64 {
        a.x * b.x + a.y * b.y + a.z * b.z
    }

    #[inline]
    fn cross(a: Vector3, b: Vector3) -> Vector3 {
        Vector3::new(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x)
    }

    #[inline]
    fn normalize(v: Vector3) -> Vector3 { v / v.length() }

    #[inline]
    fn reflect(l: Vector3, n: Vector3) -> Vector3 { l - n * 2.0 * Vector3::dot(l, n) }

    #[inline]
    fn refract(i: Vector3, n: Vector3, eta: f64) -> Vector3 {
        let k = 1.0 - eta * eta * (1.0 - Vector3::dot(n, i) * Vector3::dot(n, i));
        return if k < 0.0 {
            Vector3::zero()
        } else {
            i * eta - (n * (eta * Vector3::dot(n, i) + k.sqrt()))
        };
    }
}

impl Neg for Vector3 {
    type Output = Vector3;
    fn neg(self) -> Vector3 { Vector3::new(-self.x, -self.y, -self.z) }
}

impl Add for Vector3 {
    type Output = Vector3;
    fn add(self, v: Vector3) -> Vector3 { Vector3::new(self.x + v.x, self.y + v.y, self.z + v.z) }
}

impl Sub for Vector3 {
    type Output = Vector3;
    fn sub(self, v: Vector3) -> Vector3 { Vector3::new(self.x - v.x, self.y - v.y, self.z - v.z) }
}

impl Mul for Vector3 {
    type Output = Vector3;
    fn mul(self, v: Vector3) -> Vector3 { Vector3::new(self.x * v.x, self.y * v.y, self.z * v.z) }
}

impl Mul<f64> for Vector3 {
    type Output = Vector3;
    fn mul(self, v: f64) -> Vector3 { Vector3::new(self.x * v, self.y * v, self.z * v) }
}

impl Div<f64> for Vector3 {
    type Output = Vector3;
    fn div(self, v: f64) -> Vector3 { Vector3::new(self.x / v, self.y / v, self.z / v) }
}

struct Ray {
    origin: Vector3,
    direction: Vector3,
}

impl Ray {
    fn new(origin: Vector3, direction: Vector3) -> Ray { Ray { origin, direction } }
}

struct Sphere {
    radius: f64,
    position: Vector3,
    emission: Vector3,
    color: Vector3,
    refl: ReflectanceType,
}

impl Sphere {
    fn new(rad: f64, pos: Vector3, emit: Vector3, col: Vector3, refl: ReflectanceType) -> Sphere {
        Sphere { radius: rad, position: pos, emission: emit, color: col, refl }
    }

    #[inline]
    fn intersect(&self, ray: &Ray) -> f64 {
        let op = self.position - ray.origin;
        let b = Vector3::dot(op, ray.direction);
        let det = b * b - Vector3::dot(op, op) + self.radius * self.radius;
        if det < 0.0 { return 0.0; }
        let det = det.sqrt();
        if b - det > EPSILON { return b - det; }
        if b + det > EPSILON { return b + det; }
        return 0.0;
    }
}

#[inline]
fn intersect(ray: &Ray, spheres: &[Sphere]) -> (usize, f64) {
    let mut index = 0;
    let mut nearest = f64::INFINITY;
    for i in 0..spheres.len() {
        let distance = spheres[i].intersect(ray);
        if distance > 0.0 && distance < nearest {
            nearest = distance;
            index = i;
        }
    }
    (index, nearest)
}

#[inline]
fn diffuse(point: Vector3, normal: Vector3) -> Ray {
    let r1 = rand::random::<f64>();
    let r2 = 2.0 * PI * rand::random::<f64>();
    let u = Vector3::normalize(if normal.x.abs() > 0.1 {
        Vector3::new(-normal.z, 0.0, normal.x)
    } else {
        Vector3::new(0.0, -normal.z, normal.y)
    });
    Ray::new(point, (u * r2.cos() + Vector3::cross(normal, u) * r2.sin()) * r1.sqrt() + normal * (1.0 - r1).sqrt())
}

#[inline]
fn mirror(point: Vector3, normal: Vector3, direction: Vector3) -> Ray {
    Ray::new(point, Vector3::reflect(direction, normal))
}

fn radiance(ray: &Ray, spheres: &[Sphere], depth: i32) -> Vector3 {
    let (index, nearest) = intersect(ray, spheres);
    return if nearest < f64::INFINITY {
        let obj = &spheres[index];
        let x = ray.origin + ray.direction * nearest;
        let n = Vector3::normalize(x - obj.position);
        let nl = if Vector3::dot(n, ray.direction) < 0.0 { n } else { -n };

        let mut f = obj.color;
        if depth > 4 {
            let p = f.x.max(f.y.max(f.z));
            if depth < 256 && rand::random::<f64>() < p { f = f / p; } else { return obj.emission; }
        }

        obj.emission + f * match obj.refl {
            ReflectanceType::Diffuse => radiance(&diffuse(x, nl), spheres, depth + 1),
            ReflectanceType::Mirror => radiance(&mirror(x, n, ray.direction), spheres, depth + 1),
            _ => {
                let into = Vector3::dot(n, nl) > 0.0;
                let nc = 1.0;
                let nt = 1.5;
                let nnt = if into { nc / nt } else { nt / nc };
                let ddn = Vector3::dot(ray.direction, nl);
                let cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);
                if cos2t < 0.0 {
                    let reflected = Ray::new(x, Vector3::reflect(ray.direction, n));
                    radiance(&reflected, spheres, depth + 1)
                } else {
                    let t_dir = Vector3::refract(ray.direction, nl , nnt);
                    let r0 = (nt - nc).powi(2) / (nt + nc).powi(2);
                    let re = r0 + (1.0 - r0) * (1.0 - (if into { -ddn } else { Vector3::dot(t_dir, n) })).powi(5);
                    let tr = 1.0 - re;
                    if depth > 1 {
                        let p = 0.25 + 0.5 * re;
                        if rand::random::<f64>() < p {
                            let reflected = Ray::new(x, Vector3::reflect(ray.direction, n));
                            radiance(&reflected, spheres, depth + 1) * re / p
                        } else {
                            radiance(&Ray::new(x, t_dir), spheres, depth + 1) * tr / (1.0 - p)
                        }
                    } else {
                        let reflected = Ray::new(x, Vector3::reflect(ray.direction, n));
                        radiance(&reflected, spheres, depth + 1) * re +
                        radiance(&Ray::new(x, t_dir), spheres, depth + 1) * tr
                    }
                }
            }
        }
    } else {
        Vector3::zero()
    };
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let samples = if args.len() == 2 { (args[1].parse::<u32>().unwrap() as f64 / 4.0).ceil() as u32 } else { 1 };

    let spheres = [
        Sphere::new(100000.0, Vector3::new(100001.0, 40.8, 81.6), Vector3::zero(), Vector3::new(0.75, 0.25, 0.25), ReflectanceType::Diffuse),
        Sphere::new(100000.0, Vector3::new(-99901.0, 40.8, 81.6), Vector3::zero(), Vector3::new(0.25, 0.25, 0.75), ReflectanceType::Diffuse),
        Sphere::new(100000.0, Vector3::new(50.0, 40.8, 100000.0), Vector3::zero(), Vector3::new(0.75, 0.75, 0.75), ReflectanceType::Diffuse),
        Sphere::new(100000.0, Vector3::new(50.0, 40.8, -99830.0), Vector3::zero(), Vector3::zero(), ReflectanceType::Diffuse),
        Sphere::new(100000.0, Vector3::new(50.0, 100000.0, 81.6), Vector3::zero(), Vector3::new(0.75, 0.75, 0.75), ReflectanceType::Diffuse),
        Sphere::new(100000.0, Vector3::new(50.0, -99918.4, 81.6), Vector3::zero(), Vector3::new(0.75, 0.75, 0.75), ReflectanceType::Diffuse),
        Sphere::new(16.5, Vector3::new(27.0, 16.5, 47.0), Vector3::zero(), Vector3::new(0.999, 0.999, 0.999), ReflectanceType::Mirror),
        Sphere::new(16.5, Vector3::new(73.0, 16.5, 78.0), Vector3::zero(), Vector3::new(0.700, 0.999, 0.900), ReflectanceType::Glass),
        Sphere::new(600.0, Vector3::new(50.0, 681.33, 81.6), Vector3::new(12.0, 12.0, 12.0), Vector3::zero(), ReflectanceType::Diffuse)
    ];

    let cam = Ray::new(Vector3::new(50.0, 52.0, 295.6), Vector3::normalize(Vector3::new(0.0, -0.042612, -1.0)));
    let cx = Vector3 { x: WIDTH as f64 * 0.5135 / HEIGHT as f64, y: 0.0, z: 0.0 };
    let cy = Vector3::normalize(Vector3::cross(cx, cam.direction)) * 0.5135;

    let atomic_counter: AtomicUsize = AtomicUsize::new(0);

    let mut pixels = std::vec::from_elem(Vector3::zero(), WIDTH * HEIGHT);
    let bands: Vec<(usize, &mut [Vector3])> = pixels.chunks_mut(WIDTH).enumerate().collect();
    bands.into_par_iter().for_each(|(i, band)| {
        let y = HEIGHT - i;
        let progress = atomic_counter.fetch_add(1, Ordering::SeqCst);
        eprint!("\rRendering ({} spp) {:.2}% ", samples * 4, 100.0 * progress as f64 / (HEIGHT - 1) as f64);
        for x in 0..WIDTH {
            for sy in 0..2 {
                for sx in 0..2 {
                    let mut r = Vector3::zero();
                    for _ in 0..samples {
                        let r1 = 2.0 * rand::random::<f64>();
                        let r2 = 2.0 * rand::random::<f64>();
                        let dx = if r1 < 1.0 { r1.sqrt() - 1.0 } else { 1.0 - (2.0 - r1).sqrt() };
                        let dy = if r2 < 1.0 { r2.sqrt() - 1.0 } else { 1.0 - (2.0 - r2).sqrt() };
                        let tx = cx * (((sx as f64 + 0.5 + dx) / 2.0 + x as f64) / WIDTH as f64 - 0.5);
                        let ty = cy * (((sy as f64 + 0.5 + dy) / 2.0 + y as f64) / HEIGHT as f64 - 0.5);
                        let direction = cam.direction + tx + ty;
                        let origin = cam.origin + direction * 137.0;
                        r = r + radiance(&Ray::new(origin, Vector3::normalize(direction)), &spheres, 0);
                    }
                    r = r / samples as f64;
                    band[x] = band[x] + Vector3::new(clamp(r.x), clamp(r.y), clamp(r.z)) * 0.25;
                }
            }
        }
    });
    eprintln!();

    let ppm: String = pixels.iter().map(|pix| {
        format!("{} {} {} ", to_int(pix.x), to_int(pix.y), to_int(pix.z))
    }).collect();
    let data = (format!("P3\n{} {}\n{}\n", WIDTH, HEIGHT, 255) + &ppm);

    std::fs::File::create(&std::path::Path::new("image-rust.ppm"))
        .expect("image file creation failed")
        .write_all(data.as_bytes())
        .expect("writing image data failed");
}
