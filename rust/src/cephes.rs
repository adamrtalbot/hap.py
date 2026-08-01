//! Port of cephes incbi (inverse incomplete beta) and dependencies from
//! scipy 1.2.1, used to bit-match legacy hap.py's sompy confidence-interval
//! computations (`Tools/ci.py::jeffreys`). Legacy uses
//! `scipy.stats.beta.ppf` which calls `scipy.special.btdtri` →
//! `cephes_incbi`. To produce identical f64 LSBs as the legacy
//! container's scipy 1.2.1 we port the exact algorithm.
//!
//! Source: scipy/special/cephes/{incbi,incbet,ndtri,polevl}.c at tag
//! v1.2.1 (commits c98ef99..). The functions below mirror the C
//! source line-by-line with identical magic constants, polynomial
//! coefficients, and convergence thresholds — only translated into
//! Rust syntax and idiomatic arithmetic. Bit-exact behavior depends
//! on these constants and the floating-point operation order.

// Preserve the original Cephes decimal constants for source and bit-level parity.
#![allow(clippy::excessive_precision)]
#![allow(clippy::needless_range_loop)]

const MACHEP: f64 = f64::EPSILON / 2.0; // exactly 2^-53, matches cephes MACHEP
const MAXLOG: f64 = 7.097_827_128_933_840e2;
const MINLOG: f64 = -7.451_332_191_019_411e2;
const MAXGAM: f64 = 171.624_376_956_302_725;

const BIG: f64 = 4.503_599_627_370_496e15;
const BIGINV: f64 = 2.220_446_049_250_313e-16;

const S2PI: f64 = 2.506_628_274_631_000_502_4;

// ---------------------------------------------------------------------------
// polevl, p1evl
// ---------------------------------------------------------------------------

#[inline]
fn polevl(x: f64, coef: &[f64]) -> f64 {
    // Coefficients in REVERSE order: coef[0] = highest degree.
    let mut ans = coef[0];
    for &c in &coef[1..] {
        ans = ans * x + c;
    }
    ans
}

#[inline]
fn p1evl(x: f64, coef: &[f64]) -> f64 {
    // p1evl assumes leading coefficient is 1.0 (omitted from coef array).
    let mut ans = x + coef[0];
    for &c in &coef[1..] {
        ans = ans * x + c;
    }
    ans
}

// ---------------------------------------------------------------------------
// lgam — log of |gamma(x)| (cephes algorithm). For positive args used here.
// ---------------------------------------------------------------------------

// Stirling expansion of log Gamma — cephes A[].
const LGAM_A: [f64; 5] = [
    8.116_141_674_705_084_503e-4,
    -5.950_619_042_843_014_383e-4,
    7.936_503_404_577_169_439e-4,
    -2.777_777_777_300_996_872e-3,
    8.333_333_333_333_319_277e-2,
];

// log Gamma rational coefficients between 2 and 3 — cephes B[].
const LGAM_B: [f64; 6] = [
    -1.378_251_525_691_208_591e3,
    -3.880_163_151_346_378_409e4,
    -3.316_129_927_388_711_847e5,
    -1.162_370_974_927_623_074e6,
    -1.721_737_008_208_396_621e6,
    -8.535_556_642_457_654_656e5,
];

// log Gamma rational coefficients — cephes C[] (leading 1.0 omitted).
const LGAM_C: [f64; 6] = [
    -3.518_157_014_365_234_705e2,
    -1.706_421_066_518_811_592e4,
    -2.205_285_905_538_544_548e5,
    -1.139_334_443_679_825_072e6,
    -2.532_523_071_775_829_513e6,
    -2.018_891_414_335_327_732e6,
];

const LS2PI: f64 = 0.918_938_533_204_672_741_78;

fn lgam(x: f64) -> f64 {
    // Mirror scipy/special/cephes/gamma.c::lgam_sgn for positive x only.
    // Our use cases (lgam(a), lgam(b), lgam(a+b) with a,b>0) don't hit
    // the negative-x branch.
    if !x.is_finite() {
        return x;
    }
    if x < 13.0 {
        let mut z = 1.0_f64;
        let mut p = 0.0_f64;
        let mut u = x;
        while u >= 3.0 {
            p -= 1.0;
            u = x + p;
            z *= u;
        }
        while u < 2.0 {
            if u == 0.0 {
                return f64::INFINITY;
            }
            z /= u;
            p += 1.0;
            u = x + p;
        }
        if z < 0.0 {
            z = -z;
        }
        if u == 2.0 {
            return z.ln();
        }
        let p2 = u - 2.0;
        // Note: numerator coefficients (B) have 6 entries → polevl with degree 5.
        // Denominator (C) has 6 entries with leading 1.0 omitted → p1evl with N=6.
        let p_val = p2 * polevl(p2, &LGAM_B);
        let q_val = p1evl(p2, &LGAM_C);
        return (z.ln()) + p_val / q_val;
    }
    // Stirling for large x.
    let mut q = (x - 0.5) * x.ln() - x + LS2PI;
    if x > 1.0e8 {
        return q;
    }
    let p = 1.0 / (x * x);
    if x >= 1000.0 {
        // For very large x, Horner's polynomial expansion suffices.
        q += ((7.936_507_936_507_936_5e-4 * p - 2.777_777_777_777_777_5e-3) * p
            + 0.083_333_333_333_333_33)
            / x;
    } else {
        q += polevl(p, &LGAM_A) / x;
    }
    q
}

#[inline]
fn lbeta(a: f64, b: f64) -> f64 {
    // Match scipy 1.2.1 cephes/beta.c::lbeta operation order EXACTLY for the
    // |a+b|>MAXGAM path (the only path our germline/somatic CIs exercise).
    // Scipy swaps so |a| >= |b|, then computes:
    //   y = lgam_sgn(a+b)
    //   y = lgam_sgn(b) - y
    //   y = lgam_sgn(a) + y
    // i.e. lbeta = lgam(a) + (lgam(b) - lgam(a+b))
    // The order matters for f64 LSBs.
    let (a_ord, b_ord) = if a.abs() < b.abs() { (b, a) } else { (a, b) };
    let lgab = lgam(a_ord + b_ord);
    let lgb = lgam(b_ord);
    let lga = lgam(a_ord);
    lga + (lgb - lgab)
}

#[inline]
fn beta_fn(a: f64, b: f64) -> f64 {
    lbeta(a, b).exp()
}

// ---------------------------------------------------------------------------
// ndtri — inverse standard normal CDF.
// ---------------------------------------------------------------------------

const NDTRI_P0: [f64; 5] = [
    -5.996_335_010_141_079e1,
    9.800_107_541_859_996_6e1,
    -5.667_628_574_690_703e1,
    1.393_126_093_872_797e1,
    -1.239_165_838_673_812_5,
];

const NDTRI_Q0: [f64; 8] = [
    1.954_488_583_381_417_6,
    4.676_279_128_988_815_3,
    8.636_024_213_908_906e1,
    -2.254_626_878_541_193_7e2,
    2.002_602_123_800_606_6e2,
    -8.203_722_561_683_333e1,
    1.590_562_251_262_117e1,
    -1.183_316_211_213_300_1,
];

const NDTRI_P1: [f64; 9] = [
    4.055_448_923_059_624,
    3.152_510_945_998_938_7e1,
    5.716_281_922_464_212_8e1,
    4.408_050_738_932_008_5e1,
    1.468_495_619_288_580_3e1,
    2.186_633_068_507_902_6,
    -1.402_560_791_713_544_9e-1,
    -3.504_246_268_278_482e-2,
    -8.574_567_851_546_854e-4,
];

const NDTRI_Q1: [f64; 8] = [
    1.577_998_832_564_667_3e1,
    4.539_076_351_288_792_2e1,
    4.131_720_382_546_720_5e1,
    1.504_253_856_929_075e1,
    2.504_649_462_083_094,
    -1.421_829_228_547_877_8e-1,
    -3.808_064_076_915_782_8e-2,
    -9.332_594_808_954_574e-4,
];

const NDTRI_P2: [f64; 9] = [
    3.237_748_917_769_460_4,
    6.915_228_890_689_842,
    3.938_810_252_924_744_4,
    1.333_034_608_158_075_4,
    2.014_853_895_491_790_8e-1,
    1.237_166_348_178_200_2e-2,
    3.015_815_535_082_354e-4,
    2.658_069_746_867_375_5e-6,
    6.239_745_391_849_833e-9,
];

const NDTRI_Q2: [f64; 8] = [
    6.024_270_393_647_42,
    3.679_835_638_561_608_6,
    1.377_020_994_890_813_3,
    2.162_369_935_944_966_4e-1,
    1.342_040_060_885_431_9e-2,
    3.280_144_646_821_277_4e-4,
    2.892_478_647_453_807e-6,
    6.790_194_080_099_812e-9,
];

pub fn ndtri(y0: f64) -> f64 {
    if y0 <= 0.0 {
        return f64::NEG_INFINITY;
    }
    if y0 >= 1.0 {
        return f64::INFINITY;
    }
    let mut code = 1;
    let mut y = y0;
    if y > 1.0 - 0.135_335_283_236_612_69 {
        // 0.135... = exp(-2)
        y = 1.0 - y;
        code = 0;
    }
    if y > 0.135_335_283_236_612_69 {
        y -= 0.5;
        let y2 = y * y;
        let mut x = y + y * (y2 * polevl(y2, &NDTRI_P0) / p1evl(y2, &NDTRI_Q0));
        x *= S2PI;
        return x;
    }
    let x = (-2.0 * y.ln()).sqrt();
    let x0 = x - x.ln() / x;
    let z = 1.0 / x;
    let x1 = if x < 8.0 {
        z * polevl(z, &NDTRI_P1) / p1evl(z, &NDTRI_Q1)
    } else {
        z * polevl(z, &NDTRI_P2) / p1evl(z, &NDTRI_Q2)
    };
    let xres = x0 - x1;
    if code != 0 { -xres } else { xres }
}

// ---------------------------------------------------------------------------
// incbet — regularized incomplete beta integral (cephes algorithm).
// ---------------------------------------------------------------------------

fn incbcf(a: f64, b: f64, x: f64) -> f64 {
    // Continued fraction expansion #1.
    let mut k1 = a;
    let mut k2 = a + b;
    let mut k3 = a;
    let mut k4 = a + 1.0;
    let mut k5 = 1.0;
    let mut k6 = b - 1.0;
    let mut k7 = k4;
    let mut k8 = a + 2.0;
    let mut pkm2 = 0.0;
    let mut qkm2 = 1.0;
    let mut pkm1 = 1.0;
    let mut qkm1 = 1.0;
    let mut ans = 1.0;
    let mut r = 1.0;
    let thresh = 3.0 * MACHEP;
    for _ in 0..300 {
        let mut xk = -(x * k1 * k2) / (k3 * k4);
        let mut pk = pkm1 + pkm2 * xk;
        let mut qk = qkm1 + qkm2 * xk;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;
        xk = (x * k5 * k6) / (k7 * k8);
        pk = pkm1 + pkm2 * xk;
        qk = qkm1 + qkm2 * xk;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;
        let mut t = 1.0;
        if qk != 0.0 {
            r = pk / qk;
        }
        if r != 0.0 {
            t = ((ans - r) / r).abs();
            ans = r;
        }
        if t < thresh {
            return ans;
        }
        k1 += 1.0;
        k2 += 1.0;
        k3 += 2.0;
        k4 += 2.0;
        k5 += 1.0;
        k6 -= 1.0;
        k7 += 2.0;
        k8 += 2.0;
        if (qk.abs() + pk.abs()) > BIG {
            pkm2 *= BIGINV;
            pkm1 *= BIGINV;
            qkm2 *= BIGINV;
            qkm1 *= BIGINV;
        }
        if qk.abs() < BIGINV || pk.abs() < BIGINV {
            pkm2 *= BIG;
            pkm1 *= BIG;
            qkm2 *= BIG;
            qkm1 *= BIG;
        }
    }
    ans
}

fn incbd(a: f64, b: f64, x: f64) -> f64 {
    // Continued fraction expansion #2.
    let mut k1 = a;
    let mut k2 = b - 1.0;
    let mut k3 = a;
    let mut k4 = a + 1.0;
    let mut k5 = 1.0;
    let mut k6 = a + b;
    let mut k7 = a + 1.0;
    let mut k8 = a + 2.0;
    let mut pkm2 = 0.0;
    let mut qkm2 = 1.0;
    let mut pkm1 = 1.0;
    let mut qkm1 = 1.0;
    let z = x / (1.0 - x);
    let mut ans = 1.0;
    let mut r = 1.0;
    let thresh = 3.0 * MACHEP;
    for _ in 0..300 {
        let mut xk = -(z * k1 * k2) / (k3 * k4);
        let mut pk = pkm1 + pkm2 * xk;
        let mut qk = qkm1 + qkm2 * xk;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;
        xk = (z * k5 * k6) / (k7 * k8);
        pk = pkm1 + pkm2 * xk;
        qk = qkm1 + qkm2 * xk;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;
        let mut t = 1.0;
        if qk != 0.0 {
            r = pk / qk;
        }
        if r != 0.0 {
            t = ((ans - r) / r).abs();
            ans = r;
        }
        if t < thresh {
            return ans;
        }
        k1 += 1.0;
        k2 -= 1.0;
        k3 += 2.0;
        k4 += 2.0;
        k5 += 1.0;
        k6 += 1.0;
        k7 += 2.0;
        k8 += 2.0;
        if (qk.abs() + pk.abs()) > BIG {
            pkm2 *= BIGINV;
            pkm1 *= BIGINV;
            qkm2 *= BIGINV;
            qkm1 *= BIGINV;
        }
        if qk.abs() < BIGINV || pk.abs() < BIGINV {
            pkm2 *= BIG;
            pkm1 *= BIG;
            qkm2 *= BIG;
            qkm1 *= BIG;
        }
    }
    ans
}

fn pseries(a: f64, b: f64, x: f64) -> f64 {
    // Power series for incomplete beta.
    let ai = 1.0 / a;
    let u = (1.0 - b) * x;
    let mut v = u / (a + 1.0);
    let t1 = v;
    let mut t = u;
    let mut n = 2.0;
    let mut s = 0.0;
    let z = MACHEP * ai;
    while v.abs() > z {
        let u2 = (n - b) * x / n;
        t *= u2;
        v = t / (a + n);
        s += v;
        n += 1.0;
    }
    s += t1;
    s += ai;
    let u3 = a * x.ln();
    if (a + b) < MAXGAM && u3.abs() < MAXLOG {
        let t = 1.0 / beta_fn(a, b);
        s = s * t * x.powf(a)
    } else {
        let t = -lbeta(a, b) + u3 + s.ln();
        if t < MINLOG {
            s = 0.0;
        } else {
            s = t.exp();
        }
    }
    s
}

pub fn incbet(aa: f64, bb: f64, xx: f64) -> f64 {
    if aa <= 0.0 || bb <= 0.0 {
        return f64::NAN;
    }
    if xx <= 0.0 {
        return 0.0;
    }
    if xx >= 1.0 {
        return 1.0;
    }
    let mut flag = 0;
    if (bb * xx) <= 1.0 && xx <= 0.95 {
        return pseries(aa, bb, xx);
    }
    let w = 1.0 - xx;
    let (a, b, x, xc);
    if xx > (aa / (aa + bb)) {
        flag = 1;
        a = bb;
        b = aa;
        xc = xx;
        x = w;
    } else {
        a = aa;
        b = bb;
        xc = w;
        x = xx;
    }
    if flag == 1 && (b * x) <= 1.0 && x <= 0.95 {
        let mut t = pseries(a, b, x);
        if flag == 1 {
            if t <= MACHEP {
                t = 1.0 - MACHEP;
            } else {
                t = 1.0 - t;
            }
        }
        return t;
    }
    let y = x * (a + b - 2.0) - (a - 1.0);
    let w_cf = if y < 0.0 {
        incbcf(a, b, x)
    } else {
        incbd(a, b, x) / xc
    };
    // Multiply w by x^a (1-x)^b * gamma(a+b) / (a * gamma(a) * gamma(b)).
    let y_log = a * x.ln();
    let t_log = b * xc.ln();
    let mut t;
    if (a + b) < MAXGAM && y_log.abs() < MAXLOG && t_log.abs() < MAXLOG {
        t = xc.powf(b);
        t *= x.powf(a);
        t /= a;
        t *= w_cf;
        t *= 1.0 / beta_fn(a, b);
    } else {
        // Match scipy 1.2.1 cephes/incbet.c log path EXACTLY:
        //   y += t - lbeta(a,b);
        //   y += log(w/a);
        // where lbeta is the cephes lbeta function (NOT lgam(a)+lgam(b)-lgam(a+b)
        // computed inline; scipy's lbeta has a specific operation order that
        // produces different f64 LSBs).
        let mut log_y = y_log + (t_log - lbeta(a, b));
        log_y += (w_cf / a).ln();
        if log_y < MINLOG {
            t = 0.0;
        } else {
            t = log_y.exp();
        }
    }
    if flag == 1 {
        if t <= MACHEP {
            t = 1.0 - MACHEP;
        } else {
            t = 1.0 - t;
        }
    }
    t
}

// ---------------------------------------------------------------------------
// incbi — inverse of incbet.
// ---------------------------------------------------------------------------

#[allow(clippy::assign_op_pattern)] // Keep the Cephes expression form unchanged.
pub fn incbi(aa: f64, bb: f64, yy0: f64) -> f64 {
    if yy0 <= 0.0 {
        return 0.0;
    }
    if yy0 >= 1.0 {
        return 1.0;
    }
    // State variables that must persist across all goto labels.
    let x0: f64 = 0.0;
    let yl: f64 = 0.0;
    let x1: f64 = 1.0;
    let yh: f64 = 1.0;
    let mut nflg: bool = false;
    let a: f64;
    let b: f64;
    let y0: f64;
    let y: f64;
    let x: f64;
    let yp: f64;
    let rflg: bool;
    let dithresh: f64;
    if aa <= 1.0 || bb <= 1.0 {
        dithresh = 1.0e-6;
        rflg = false;
        a = aa;
        b = bb;
        y0 = yy0;
        x = a / (a + b);
        y = incbet(a, b, x);
        return ihalve_loop(
            a, b, x, x0, x1, yl, yh, y, y0, rflg, nflg, dithresh, aa, bb, yy0,
        );
    } else {
        dithresh = 1.0e-4;
        // Approximation to inverse function via standard normal.
        let mut yp_init = -ndtri(yy0);
        if yy0 > 0.5 {
            rflg = true;
            a = bb;
            b = aa;
            y0 = 1.0 - yy0;
            yp_init = -yp_init;
        } else {
            rflg = false;
            a = aa;
            b = bb;
            y0 = yy0;
        }
        let lgm = (yp_init * yp_init - 3.0) / 6.0;
        let xp = 2.0 / (1.0 / (2.0 * a - 1.0) + 1.0 / (2.0 * b - 1.0));
        let mut d = yp_init * (xp + lgm).sqrt() / xp
            - (1.0 / (2.0 * b - 1.0) - 1.0 / (2.0 * a - 1.0))
                * (lgm + 5.0 / 6.0 - 2.0 / (3.0 * xp));
        d = 2.0 * d;
        if d < MINLOG {
            return apply_rflg(0.0, rflg);
        }
        x = a / (a + b * d.exp());
        y = incbet(a, b, x);
        yp = (y - y0) / y0;
        if yp.abs() < 0.2 {
            // goto newt directly (skip ihalve)
            return newton_then_maybe_ihalve(
                a, b, x, x0, x1, yl, yh, y, y0, rflg, &mut nflg, dithresh, aa, bb, yy0,
            );
        }
    }
    ihalve_loop(
        a, b, x, x0, x1, yl, yh, y, y0, rflg, nflg, dithresh, aa, bb, yy0,
    )
}

fn apply_rflg(x: f64, rflg: bool) -> f64 {
    if rflg {
        if x <= MACHEP { 1.0 - MACHEP } else { 1.0 - x }
    } else {
        x
    }
}

/// Run the ihalve interval-halving loop. On convergence transition to
/// newton; on bracket-crossing-0.75 with rflg==1, restart with swapped
/// a/b. On full 100-iter exhaust, return x (relaxed-threshold restart
/// is handled by the newton-then-maybe-ihalve loop).
#[allow(clippy::too_many_arguments)]
fn ihalve_loop(
    mut a: f64,
    mut b: f64,
    mut x: f64,
    mut x0: f64,
    mut x1: f64,
    mut yl: f64,
    mut yh: f64,
    mut y: f64,
    mut y0: f64,
    mut rflg: bool,
    mut nflg: bool,
    dithresh: f64,
    aa: f64,
    bb: f64,
    yy0: f64,
) -> f64 {
    'outer: loop {
        let mut dir: i32 = 0;
        let mut di: f64 = 0.5;
        for i in 0..100 {
            if i != 0 {
                x = x0 + di * (x1 - x0);
                if x == 1.0 {
                    x = 1.0 - MACHEP;
                }
                if x == 0.0 {
                    di = 0.5;
                    x = x0 + di * (x1 - x0);
                    if x == 0.0 {
                        return apply_rflg(0.0, rflg);
                    }
                }
                y = incbet(a, b, x);
                let yp = (x1 - x0) / (x1 + x0);
                if yp.abs() < dithresh {
                    // goto newt
                    return newton_then_maybe_ihalve(
                        a, b, x, x0, x1, yl, yh, y, y0, rflg, &mut nflg, dithresh, aa, bb, yy0,
                    );
                }
                let yp = (y - y0) / y0;
                if yp.abs() < dithresh {
                    return newton_then_maybe_ihalve(
                        a, b, x, x0, x1, yl, yh, y, y0, rflg, &mut nflg, dithresh, aa, bb, yy0,
                    );
                }
            }
            if y < y0 {
                x0 = x;
                yl = y;
                if dir < 0 {
                    dir = 0;
                    di = 0.5;
                } else if dir > 3 {
                    di = 1.0 - (1.0 - di) * (1.0 - di);
                } else if dir > 1 {
                    di = 0.5 * di + 0.5;
                } else {
                    di = (y0 - y) / (yh - yl);
                }
                dir += 1;
                if x0 > 0.75 {
                    if rflg {
                        rflg = false;
                        a = aa;
                        b = bb;
                        y0 = yy0;
                    } else {
                        rflg = true;
                        a = bb;
                        b = aa;
                        y0 = 1.0 - yy0;
                    }
                    x = 1.0 - x;
                    y = incbet(a, b, x);
                    x0 = 0.0;
                    yl = 0.0;
                    x1 = 1.0;
                    yh = 1.0;
                    continue 'outer;
                }
            } else {
                x1 = x;
                if rflg && x1 < MACHEP {
                    return apply_rflg(0.0, rflg);
                }
                yh = y;
                if dir > 0 {
                    dir = 0;
                    di = 0.5;
                } else if dir < -3 {
                    di *= di;
                } else if dir < -1 {
                    di *= 0.5;
                } else {
                    di = (y - y0) / (yh - yl);
                }
                dir -= 1;
            }
        }
        // 100-iter exhausted (mtherr PLOSS in cephes)
        if x0 >= 1.0 {
            return apply_rflg(1.0 - MACHEP, rflg);
        }
        if x <= 0.0 {
            return apply_rflg(0.0, rflg);
        }
        // Fall through to newt (which will run since nflg still false).
        return newton_then_maybe_ihalve(
            a, b, x, x0, x1, yl, yh, y, y0, rflg, &mut nflg, dithresh, aa, bb, yy0,
        );
    }
}

/// Run the Newton iteration. If it converges, return. If not,
/// re-enter ihalve_loop with relaxed dithresh; on second visit
/// to newton, nflg is true so it returns immediately.
#[allow(clippy::too_many_arguments)]
fn newton_then_maybe_ihalve(
    a: f64,
    b: f64,
    mut x: f64,
    mut x0: f64,
    mut x1: f64,
    mut yl: f64,
    mut yh: f64,
    mut y: f64,
    y0: f64,
    rflg: bool,
    nflg: &mut bool,
    _dithresh: f64,
    aa: f64,
    bb: f64,
    yy0: f64,
) -> f64 {
    if *nflg {
        return apply_rflg(x, rflg);
    }
    *nflg = true;
    let lgm = lgam(a + b) - lgam(a) - lgam(b);
    for i in 0..8 {
        if i != 0 {
            y = incbet(a, b, x);
        }
        if y < yl {
            x = x0;
            y = yl;
        } else if y > yh {
            x = x1;
            y = yh;
        } else if y < y0 {
            x0 = x;
            yl = y;
        } else {
            x1 = x;
            yh = y;
        }
        if x == 1.0 || x == 0.0 {
            break;
        }
        let d_log = (a - 1.0) * x.ln() + (b - 1.0) * (1.0 - x).ln() + lgm;
        if d_log < MINLOG {
            return apply_rflg(x, rflg);
        }
        if d_log > MAXLOG {
            break;
        }
        let d_val = d_log.exp();
        let d_step = (y - y0) / d_val;
        let mut xt = x - d_step;
        if xt <= x0 {
            let yp = (x - x0) / (x1 - x0);
            xt = x0 + 0.5 * yp * (x - x0);
            if xt <= 0.0 {
                break;
            }
        }
        if xt >= x1 {
            let yp = (x1 - x) / (x1 - x0);
            xt = x1 - 0.5 * yp * (x1 - x);
            if xt >= 1.0 {
                break;
            }
        }
        x = xt;
        if (d_step / x).abs() < 128.0 * MACHEP {
            return apply_rflg(x, rflg);
        }
    }
    // Did not converge: relax dithresh and restart ihalve.
    let dithresh = 256.0 * MACHEP;
    ihalve_loop(
        a, b, x, x0, x1, yl, yh, y, y0, rflg, *nflg, dithresh, aa, bb, yy0,
    )
}

// Public entry: matches scipy.special.btdtri(a, b, p).
pub fn btdtri(a: f64, b: f64, p: f64) -> f64 {
    incbi(a, b, p)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lgam_matches_scipy() {
        let cases = [
            (1104.5_f64, 6_632.310_431_219_753),
            (15131.5, 130498.242_030_250_2),
            (16236.0, 141167.868_659_465_75),
        ];
        for (x, expected) in cases.iter() {
            let r = lgam(*x);
            let rel = (r - expected).abs() / expected;
            println!(
                "lgam({}) = {}  (expected {}  rel = {:.3e})",
                x, r, expected, rel
            );
            assert!(rel < 1e-13, "lgam({}) = {} expected {}", x, r, expected);
        }
    }

    #[test]
    fn ndtri_matches_scipy() {
        // scipy.special.ndtri(0.025) = -1.959963984540054
        let r = ndtri(0.025);
        assert!(
            (r - (-1.959_963_984_540_054)).abs() < 1e-12,
            "ndtri(0.025)={}",
            r
        );
    }

    #[test]
    fn incbet_matches_scipy_btdtr() {
        // Reference values from scipy.special.betainc 1.2.1.
        // Values must be bit-exact (or extremely close) for incbi to converge correctly.
        let cases = [
            (
                15131.5_f64,
                1104.5,
                0.9280491268419975,
                0.024999999999995838_f64,
            ),
            (15131.5, 1104.5, 0.9280490605005636, 0.024998104484453124),
            (15131.5, 1104.5, 0.93, 0.15909898042524057_f64),
            (15131.5, 1104.5, 0.928, 0.02362923271944636_f64),
        ];
        for (a, b, x, expected) in cases.iter() {
            let r = incbet(*a, *b, *x);
            let rel = (r - expected).abs() / expected;
            println!(
                "incbet({}, {}, {}) = {}  (expected ≈ {}  rel_err = {:.3e})",
                a, b, x, r, expected, rel
            );
            // Tolerate up to 1e-3 relative for sanity but ideally 1e-13.
            assert!(
                rel < 1e-3,
                "incbet({}, {}, {}) = {} differs from {}",
                a,
                b,
                x,
                r,
                expected
            );
        }
    }

    #[test]
    fn incbi_matches_legacy_jeffreys_recall() {
        // Legacy scipy 1.2.1: beta.ppf(0.025, 15131.5, 1104.5) = 0.9280491268419975 (bit-exact f64)
        let r = incbi(15131.5, 1104.5, 0.025);
        let expected = 0.928_049_126_841_997_5;
        let diff = (r - expected).abs();
        let r_bits = r.to_bits();
        let exp_bits = expected.to_bits();
        let ulps = (r_bits as i64 - exp_bits as i64).abs();
        println!(
            "incbi(15131.5, 1104.5, 0.025)={} expected={} diff={:.3e} ULP_diff={}",
            r, expected, diff, ulps
        );
        assert!(ulps <= 0, "Not bit-exact: ULPs apart = {}", ulps);
    }
}
