//! Shared, filesystem-free statistical calculations.

use super::cephes;

/// Modified Jeffreys interval used by legacy Tools/ci.py.
pub(crate) fn jeffreys_interval(x: usize, n: usize, alpha: f64) -> (f64, f64) {
    if n == 0 {
        return (0.0, 1.0);
    }
    let lower = if x == n {
        cephes::legacy_pow(alpha / 2.0, 1.0 / n as f64)
    } else if x <= 1 {
        0.0
    } else {
        cephes::incbi(x as f64 + 0.5, (n - x) as f64 + 0.5, alpha / 2.0)
    };
    let upper = if x == 0 {
        1.0 - cephes::legacy_pow(alpha / 2.0, 1.0 / n as f64)
    } else if x >= n - 1 {
        1.0
    } else {
        cephes::incbi(x as f64 + 0.5, (n - x) as f64 + 0.5, 1.0 - alpha / 2.0)
    };
    (lower.max(0.0), upper.min(1.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_observation_has_full_interval() {
        assert_eq!(jeffreys_interval(0, 0, 0.05), (0.0, 1.0));
    }
}
