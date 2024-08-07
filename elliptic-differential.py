from math import factorial, ceil
from sympy import symbols, expand, collect, O, Poly, Matrix, diff

z, w = symbols('z w') # variables
[a1, a2, a3, a4, a6] = [0, 0, 0, 1, 1] # coefficients
n = 15 # highest power in trunacted series

# Return truncated Laurent series for w in terms of z
def w_trunc():
    f = [0]*n
    f[0] = z**3 + a1*z*w + a2*z**2*w + a3*w**2 + a4*z*w**2 + a6*w**3
    for i in range(1, n):
        new_f = expand(f[i-1].subs(w, f[0]))
        f[i] = (new_f + O(w**ceil((n+3)/3))).removeO()
    return (collect(expand(f[n-1].subs(w, 0)), z) + O(z**(n+3))).removeO()

# Compute the inverse of a truncated power series with nonzero constant term
def invert_trunc(series):
    coeffs = Poly(series).all_coeffs()
    coeffs.reverse()
    coeffs += [0]*(n-len(coeffs))

    a0 = coeffs[0]
    nm_coeffs = [m/a0 for m in coeffs]

    inv_coeffs = [0]*n
    inv_coeffs[0] = 1

    for i in range(1,n):
        inv_coeffs[i] = -sum([inv_coeffs[j]*nm_coeffs[i-j] for j in range(i)])
    
    inv = 0

    for i in range(n):
        inv += (1/a0)*inv_coeffs[i]*z**i
        
    return inv

# Return truncated Laurent series for x in terms of z
def x_trunc(wt):
    return expand(z**(-2) * invert_trunc(wt/z**3))

# Return truncated Laurent series for y in terms of z
def y_trunc(wt):
    return expand(-z**(-3) * invert_trunc(wt/z**3))

# Return truncated Laurent series for the invariant differential in terms of z
def diff_trunc(xt, yt):
    top = diff(xt, z)
    bottom = 2*yt + a1*xt + a3
    bottom_inv = z**3 * invert_trunc(expand(bottom*z**3))
    return (expand(top*bottom_inv) + O(z**n)).removeO()

def get_series():
    wt = w_trunc()
    xt = x_trunc(wt)
    yt = y_trunc(wt)
    return expand(z*diff_trunc(xt, yt))

def main():
    print(get_series())

if __name__ == '__main__':
    main()
