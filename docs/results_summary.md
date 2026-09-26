# Compact results summary

## Positive coefficient peak

At q-power exponent `j=100` (coefficient of `q^100`):

- `Q_peak=28`
- `max coefficient = 498517883639602510421134799165663534780829757383208`

Fit over q-power exponents `j=50..100`:

- `Q_peak ~ 1.291882 j^(2/3)`
- `log c_max ~ 6.112350 j^(2/3) - 2.135663 log j - 5.118475`

## Direct fixed-Q crossings

- Q=4/3: `2J_cross=22.003694`; quoted Appendix-A prediction `21.9575`.
- Q=3: `2J_cross=103.251119`; nearest image `2J_ctr=103.776695`.
- Q=6: `2J_cross=328.862367`; nearest image `2J_ctr=328.302031`.

## Q<=40 ratio test

Selected integer charges:

| Q | 2J_cross | 2J_ctr closest | ratio |
|---:|---:|---:|---:|
| 20 | 2123.982740 | 2104.288279 | 1.009359 |
| 25 | 2971.684802 | 2953.996144 | 1.005988 |
| 30 | 3904.989538 | 3899.973638 | 1.001286 |
| 35 | 4916.131688 | 4934.094271 | 0.996359 |
| 40 | 5998.898966 | 5984.125145 | 1.002469 |

The ratio is the preferred macroscopic convergence diagnostic.

## Vertical mismatch

Directly measured through Q=18:

`|slope| ~ 0.449 Q_ctr^(-0.495)`.

Examples:

| Q | Delta_log | slope | measured delta(2J) | -Delta/slope |
|---:|---:|---:|---:|---:|
| 3 | -0.049650 | -0.094483 | -0.525576 | -0.525490 |
| 8 | 0.780400 | -0.058749 | 13.282020 | 13.283736 |
| 18 | 0.701643 | -0.039820 | 17.620703 | 17.620400 |

This demonstrates that an O(1) vertical correction can produce an O(sqrt(Q)) horizontal displacement because the crossing function becomes flat.
