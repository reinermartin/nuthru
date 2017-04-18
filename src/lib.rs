//! The crate `nathru` (NUmber THeory in RUst) implements a few simple number theory functions. 


use std::cmp;

#[cfg(test)]
mod tests {
	use super::gcd;

    #[test]
    fn test1_gcd() {
    	assert_eq!(gcd(945, 165), 15);
    	assert_eq!(gcd(945, -165), -15);
    }
}



/// This function complements the % operator, as it behaves like the normal mathematical mod operator even for negative inputs.
pub fn modulo(n: i64, m: i64) -> i64 {
    if n >= 0 { n % m } else { - (-n) % m }
}

/// Computes powers modulo m
pub fn power_mod(b: i64, e: i64, m: i64) -> i64 {
    if e < 0 { unimplemented!() }
    let mut bm = b % m;
    let mut res = 1;
    let mut e = e;
    while e > 0 {
        if e & 1 != 0 { res = (res * bm) % m }
        bm = (bm * bm) % m;
        e >>= 1
    }
    modulo(res, m)
}

/// The greatest common divisor of a and b
pub fn gcd(a: i64, b: i64) -> i64
{
    if b == 0 { a } else { gcd(b, a % b) }
}

/// The integer square root
pub fn int_sqrt(n: i64) -> Option<i64> {
	let t = (n as f64).sqrt() as i64;
	if t*t == n { return Some(t) }
	if (t+1)*(t+1) == n { return Some(t+1) }
    None
}


/// Legendre symbol, returns 1, 0, or -1 mod p
pub fn legendre_symbol(a: i64, p: i64) -> i64 {
    power_mod(a, (p-1)/2, p)
}


/// The Tonelli–Shanks algorithm finds solutions to x^2 = n mod p, where p is an odd prime.
// We are following the notation in https://en.wikipedia.org/wiki/Tonelli–Shanks_algorithm (WP)
pub fn ts(n: i64, p: i64) -> (i64, i64, bool) {

    if legendre_symbol(n, p) != 1 { return (0, 0, false) }

    // WP step 1, factor out powers two.
    // variables Q, S named as at WP.
    let mut q = p - 1;
    let mut s = 0;
    while q & 1 == 0 {
        s += 1;
        q >>= 1
    }

    // WP step 1, direct solution
    if s == 1 {
        let r1 = power_mod(n, (p+1)/4, p);
        return (r1, p - r1, true)
    }

    // WP step 2, select z, assign c
    let mut z = 2;
    while legendre_symbol(z, p) != p-1 { z += 1 }
    let mut c = power_mod(z, q, p);

    // WP step 3, assign R, t, M
    let mut r = power_mod(n, (q+1)/2, p);
    let mut t = power_mod(n, q, p);
    let mut m = s;

    // WP step 4, loop
    loop {
        // WP step 4.1, termination condition
        if t == 1 { return (r, p - r, true) }

        // WP step 4.2, find lowest i...
        let mut i = 0;
        let mut z = t;
        while z != 1 && i < m-1 {
            z = z * z % p;
            i += 1
        }

        // WP step 4.3, using a variable b, assign new values of R, t, c, M
        let mut b = c;
        let mut e = m - i - 1;
        while e > 0 {
            b = b * b % p;
            e -= 1
        }
        r = r * b % p;
        c = b * b % p;
        t = t * c % p;
        m = i;
    }
}

/// Finds integer solution to x^2 + y^2 = p
pub fn cornacchia(p: i64) -> (i64, i64) {
    if p == 1 { return (1, 0)}
    if p == 2 { return (1, 1)}
    if p % 4 != 1 { panic!(""); }

    let res = ts(p-1, p);
    let mut a = p;
    let mut b = cmp::max(res.0, res.1);
    let l = (p as f64).sqrt() as i64;
    while b > l {
        let r = a % b;
        a = b;
        b = r;
    }
    let c = p - b*b;
    (b, int_sqrt(c).unwrap())
}
