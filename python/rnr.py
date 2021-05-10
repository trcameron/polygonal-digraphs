# Restricted Numerical Range Functions
from numpy import argsort, array, conj, cos, diag, dot, exp, imag, pi, real, sin, sum, transpose, zeros
from numpy.linalg import eigh, eigvals, eigvalsh, norm
from nauty_directg_reader import digraph6
from networkx import DiGraph, draw_shell
from sys import stdin
from matplotlib import pyplot as plt

# Tolerance Parameters
EPS = 2**(-52)
TOL = 2**(-42)
TOL2 = 2**(-12)

###############################################
###             Cmplx Convex Hull           ###
###############################################
def cmplxConvHull(points):
    # Sort the points lexicographically (tuples are compared lexicographically).
    # Remove duplicates to detect the case we have just one unique point.
    points = sorted(set(points))
    # Boring case: no points or a single point, possibly repeated multiple times.
    if(len(points) <= 1):
        return points
    # 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
    # Returns a positive value, if OAB makes a counter-clockwise turn,
    # negative for clockwise turn, and zero if the points are collinear.
    def _cross(o, a, b):
        return (a.real - o.real)*(b.imag - o.imag) - (a.imag - o.imag)*(b.real-o.real)
    # Build lower hull
    lower = []
    for p in points:
        while((len(lower)>=2) and (_cross(lower[-2],lower[-1],p)<EPS)):
            lower.pop()
        lower.append(p)
    # Build upper hull
    upper = []
    for p in reversed(points):
        while((len(upper)>=2) and (_cross(upper[-2],upper[-1],p)<EPS)):
            upper.pop()
        upper.append(p)
    # Concatenation of the lower and upper hulls gives the convex hull.
    # Last point of each list is omitted because it is repeated at the beginning of the other list.
    return lower[:-1] + upper[:-1]
###############################################
###             Polygon Area                ###
###############################################
def poly_area(v):
    s = 0
    for k in range(len(v)-1):
        s += conj(v[k])*v[k+1]
    return 0.5*s.imag
###############################################
###             Numerical Range Main        ###
###############################################
def nr_main(a):
    # set size and norm
    n,_ = a.shape
    nrm = norm(a,ord='fro')
    # Hermitian test
    err = norm(a - conj(transpose(a)),ord='fro')
    if(err < max(TOL*nrm,TOL)):
        eig = eigvalsh(a)
        ind = argsort(eig)
        eig = eig[ind]
        return [eig[0],eig[-1]], eig
    # Normal test
    err = norm(dot(conj(transpose(a)),a) - dot(a,conj(transpose(a))),ord='fro')
    if(err < max(TOL*nrm,TOL)):
        eig = eigvals(a)
        convHull = cmplxConvHull(eig)
        convHull.append(convHull[0])
        return convHull, eig
    # eigenvalues and inscribed polygon vertices (for non-normal matrices)
    r = 3
    e = [0 for k in range(r+1)]
    p = [0 for k in range(r+1)]
    for k in range(r):
        # multiply by complex exponential
        a1 = exp(2*pi*1j*k/r)*a
        # compute hermitian part
        a2 = (a1 + transpose(conj(a1)))/2
        # compute eigenvalues and eigenvectors of hermitian part
        eigval, eigvec = eigh(a2)
        # sort eigenvalues and store maximum eigenvalue and associated inscribed polygonal vertex
        ind = argsort(eigval)
        eigval = eigval[ind]
        e[r-k] = eigval[-1]
        eigvec = eigvec[:,ind]
        p[r-k] = dot(conj(eigvec[:,-1]),dot(a,eigvec[:,-1]))
    # complete cycle
    e[0] = e[r]
    p[0] = p[r]
    # compute circumscribed polygon vertices
    q = [0 for k in range(r+1)]
    for k in range(r):
        q[r-k] = exp(-2*pi*1j*k/r)*(e[k]+1j*(e[k]*cos(2*pi/r)-e[k+1])/sin(2*pi/r))
    # complete cycle
    q[0] = q[r]
    # compute area
    s1 = poly_area(p)
    s2 = poly_area(q)
    # refine approximation
    return _nr_sub(a,e,p,q,n,r,s1,s2)
###############################################
###             Numerical Range Sub         ###
###############################################
def _nr_sub(a,e,p,q,n,r,s1,s2):
    while((s2-s1)>TOL2*s2):
        # update r
        r = 2*r
        e = [0 for k in range(r+1)]
        p = [0 for k in range(r+1)]
        for k in range(r):
            # multiply by complex exponential
            a1 = exp(2*pi*1j*k/r)*a
            # compute hermitian part
            a2 = (a1 + transpose(conj(a1)))/2
            # compute eigenvalues and eigenvectors of hermitian part
            eigval, eigvec = eigh(a2)
            # sort eigenvalues and store maximum eigenvalue and associated inscribed polygonal vertex
            ind = argsort(eigval)
            eigval = eigval[ind]
            e[r-k] = eigval[-1]
            eigvec = eigvec[:,ind]
            p[r-k] = dot(conj(eigvec[:,-1]),dot(a,eigvec[:,-1]))
        # complete cycle
        e[0] = e[r]
        p[0] = p[r]
        # compute circumscribed polygon vertices
        q = [0 for k in range(r+1)]
        for k in range(r):
            q[r-k] = exp(-2*pi*1j*k/r)*(e[k]+1j*(e[k]*cos(2*pi/r)-e[k+1])/sin(2*pi/r))
        # complete cycle
        q[0] = q[r]
        # compute area
        s1 = poly_area(p)
        s2 = poly_area(q)
    # return
    return q, eigvals(a)
###############################################
###         Restricted Laplacian            ###
###############################################
def ql(l):
    n, _ = l.shape
    q = zeros((n,n-1),dtype=float)
    for j in range(n-1):
        q[0:j+1,j] = 1
        q[j+1,j] = -(j+1)
        q[:,j] = q[:,j]/norm(q[:,j])
    a = dot(transpose(q),dot(l,q))
    return a
###############################################
###         Restricted Numerical Range      ###
###############################################
def qnr(l):
    return nr_main(ql(l))
###############################################
###             Normality Tests             ###
###############################################
def normality(l):
    # check if l is normal
    ln = False
    nrm = norm(l,ord='fro')
    err = norm(dot(conj(transpose(l)),l) - dot(l,conj(transpose(l))),ord='fro')
    if(err < max(TOL,TOL*nrm)):
        ln = True
    # build a = q^{T}lq
    a = ql(l)
    # check if a is normal
    an = False
    nrm = norm(a,ord='fro')
    err = norm(dot(conj(transpose(a)),a) - dot(a,conj(transpose(a))),ord='fro')
    if(err < max(TOL*nrm,TOL)):
        an = True
    # return
    return ln, an
###############################################
###             main                        ###
###############################################
def main():
    try:
        for line in stdin:
            a = digraph6(bytes(line.rstrip(),'utf-8'))
            f,e = qnr(diag(sum(a,axis=1))-a)
            g = DiGraph(a)
            
            fig, axes = plt.subplots(nrows=2, ncols=1)
            fig.tight_layout()
            axes[1].plot(real(f),imag(f),color='#000000',fillstyle='right')
            axes[1].fill(real(f),imag(f),color='#C0C0C0')
            axes[1].plot(real(e),imag(e),color='#606060',linestyle='none',marker='*',markersize=10.0)
            draw_shell(g,node_color='#606060',ax=axes[0])
            plt.show()
    except Exception as e:
        print(e)
if __name__ == '__main__':
    main()
